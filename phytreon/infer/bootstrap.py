"""Bootstrap branch support.

Resample alignment columns with replacement, rebuild the tree for each
replicate, and score every clade of the reference tree by how often it
recurs.  Works on any builder ``Alignment -> Tree``; the default is
p-distance + neighbor-joining.
"""
from __future__ import annotations

import random
from typing import Callable, Dict, List, Optional, Tuple

from ..core.tree import Tree
from .align import Alignment


def distance_matrix_model(aln: Alignment, model: str = "jc69"
                          ) -> Tuple[List[str], List[List[float]]]:
    """Pairwise evolutionary distances with a substitution-model correction.

    ``model``:
      * ``"raw"``     -- uncorrected p-distance,
      * ``"jc69"``    -- Jukes-Cantor: d = -3/4 ln(1 - 4p/3),
      * ``"k2p"``     -- Kimura 2-parameter (transition/transversion separated),
      * ``"poisson"`` -- Poisson correction (protein analogue of JC69):
        d = -ln(1 - p).

    ``"jc69"``/``"k2p"`` need nucleotide data; on non-nucleotide alignments
    they fall back to raw p-distance (unchanged).  ``"poisson"`` is an
    explicit opt-in for protein data and is applied as requested.  Gaps,
    missing data and ambiguity codes (IUPAC ``N``/``R``/``Y``... in DNA,
    ``X``/``B``/``Z`` in protein) are ignored pairwise rather than counted as
    differences.  Saturated pairs, where the correction diverges, get the
    model's largest finite distance, so a more divergent pair never comes out
    closer than a less divergent one.
    """
    import math

    import numpy as np

    from .align import guess_type
    n = aln.nseq
    seqs = [s.upper() for s in aln.seqs]
    # judged on the whole sample, not "every character is ACGTN": one IUPAC
    # ambiguity code (R, Y, ...) used to flip a DNA alignment to protein and
    # silently drop the JC69/K2P correction for every pair
    nucleotide = guess_type(seqs) == "nucleotide"
    if nucleotide:
        seqs = [s.replace("U", "T") for s in seqs]
    if not nucleotide and model not in ("raw", "poisson"):
        model = "raw"
    codes = np.array([np.frombuffer(s.encode("latin-1", "replace"), np.uint8)
                      for s in seqs]) if n else np.zeros((0, 0), np.uint8)
    if nucleotide:
        # anything but ACGT is a gap, missing data or an ambiguity code
        usable = np.isin(codes, np.frombuffer(b"ACGT", np.uint8))
    else:
        # protein, or a character matrix (0/1 traits, lineage alleles): any
        # state counts except gaps, missing data and protein ambiguity codes
        usable = ~np.isin(codes, np.frombuffer(b"-.?*XBZJ", np.uint8))
    purine = np.isin(codes, np.frombuffer(b"AG", np.uint8))
    # the log arguments are floored where the correction diverges, so a
    # saturated pair gets the model's largest finite distance (~3.24 for
    # JC69/K2P, ~9.2 for Poisson) and distances stay monotone in divergence
    floor_jc = 1 - 4 * 0.74 / 3
    floor_poisson = 1e-4
    D = [[0.0] * n for _ in range(n)]
    for a in range(n - 1):
        both = usable[a] & usable[a + 1:]               # (n-a-1, ncol)
        differ = both & (codes[a] != codes[a + 1:])
        comps = both.sum(axis=1)
        diffs = differ.sum(axis=1)
        # a transition keeps purine (A<->G) or pyrimidine (C<->T) status
        tss = (differ & (purine[a] == purine[a + 1:])).sum(axis=1)
        for k, b in enumerate(range(a + 1, n)):
            comp, diff = int(comps[k]), int(diffs[k])
            ts = int(tss[k])
            tv = diff - ts
            if comp == 0:
                d = 1.0
            elif model == "k2p":
                P, Q = ts / comp, tv / comp
                d = (-0.5 * math.log(max(1 - 2 * P - Q, floor_jc))
                     - 0.25 * math.log(max(1 - 2 * Q, floor_jc)))
            elif model == "jc69":
                p = diff / comp
                d = -0.75 * math.log(max(1 - 4 * p / 3, floor_jc))
            elif model == "poisson":
                p = diff / comp
                d = -math.log(max(1 - p, floor_poisson))
            else:
                d = diff / comp
            D[a][b] = D[b][a] = d
    return list(aln.names), D


def p_distance_matrix(aln: Alignment) -> Tuple[List[str], List[List[float]]]:
    """Uncorrected pairwise p-distance (kept for callers that want it)."""
    return distance_matrix_model(aln, model="raw")


def nj_builder(aln: Alignment, model: str = "jc69") -> Tree:
    from .distance import neighbor_joining
    names, D = distance_matrix_model(aln, model)
    return neighbor_joining(names, D)


def upgma_builder(aln: Alignment, model: str = "jc69") -> Tree:
    from .distance import upgma
    names, D = distance_matrix_model(aln, model)
    return upgma(names, D)


def _bipartitions(tree: Tree, all_leaves: frozenset, anchor: str) -> List[frozenset]:
    """Non-trivial bipartitions (splits) of the tree, rooting-independent.

    Each internal edge splits the tips in two; we represent the split by the
    side that does *not* contain ``anchor`` so it is canonical regardless of
    where the tree is rooted.
    """
    n = len(all_leaves)
    out = []
    for node in tree.traverse():
        if node.is_leaf or node.is_root:
            continue
        side = frozenset(node.leaf_names())
        if not (2 <= len(side) <= n - 2):
            continue
        out.append(side if anchor not in side else (all_leaves - side))
    return out


def bootstrap_support(aln: Alignment,
                      builder: Optional[Callable[[Alignment], Tree]] = None,
                      n: int = 100, as_percent: bool = True,
                      seed: Optional[int] = None,
                      reference: Optional[Tree] = None
                      ) -> Tuple[Tree, Dict[frozenset, float]]:
    """Return ``(reference_tree_with_support, {split: support})``.

    Support is computed on **bipartitions**, so it is correct even if the
    reference tree is (midpoint-)rooted differently from the replicates.
    Each internal node gets ``node.support`` set to its split's frequency
    (0-100 if ``as_percent``, else 0-1).
    """
    builder = builder or nj_builder
    ref = reference or builder(aln)
    all_leaves = frozenset(ref.leaf_names())
    anchor = min(all_leaves)
    nt = len(all_leaves)

    ref_map: Dict[frozenset, List] = {}
    for node in ref.traverse():
        if node.is_leaf or node.is_root:
            continue
        side = frozenset(node.leaf_names())
        if not (2 <= len(side) <= nt - 2):
            continue
        canon = side if anchor not in side else (all_leaves - side)
        ref_map.setdefault(canon, []).append(node)

    counts: Dict[frozenset, int] = {c: 0 for c in ref_map}
    rng = random.Random(seed)
    ncol = aln.ncol
    for _ in range(n):
        idx = [rng.randrange(ncol) for _ in range(ncol)]
        rep = builder(aln.select_columns(idx))
        for split in set(_bipartitions(rep, all_leaves, anchor)):
            if split in counts:
                counts[split] += 1

    scale = 100.0 if as_percent else 1.0
    support = {c: counts[c] / n * scale for c in counts}
    for canon, nodes in ref_map.items():
        for node in nodes:
            node.support = round(support[canon], 1)
    return ref, support
