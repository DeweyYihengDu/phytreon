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


_TRANSITION = {frozenset("AG"), frozenset("CT"), frozenset("AU"), frozenset("CU")}


def distance_matrix_model(aln: Alignment, model: str = "jc69"
                          ) -> Tuple[List[str], List[List[float]]]:
    """Pairwise evolutionary distances with a substitution-model correction.

    ``model``:
      * ``"raw"``  -- uncorrected p-distance,
      * ``"jc69"`` -- Jukes-Cantor: d = -3/4 ln(1 - 4p/3),
      * ``"k2p"``  -- Kimura 2-parameter (transition/transversion separated).

    Corrections need nucleotide data; non-nucleotide alignments fall back to
    raw p-distance.  Gapped positions are ignored pairwise.
    """
    import math
    n = aln.nseq
    seqs = [s.upper() for s in aln.seqs]
    nucleotide = set("".join(seqs[:2])[:500]) <= set("ACGTUN-.")
    if not nucleotide and model != "raw":
        model = "raw"
    D = [[0.0] * n for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            comp = ts = tv = diff = 0
            for ca, cb in zip(seqs[a], seqs[b]):
                if ca in "-.N" or cb in "-.N":
                    continue
                comp += 1
                if ca != cb:
                    diff += 1
                    if frozenset((ca, cb)) in _TRANSITION:
                        ts += 1
                    else:
                        tv += 1
            if comp == 0:
                d = 1.0
            elif model == "k2p":
                P, Q = ts / comp, tv / comp
                try:
                    d = -0.5 * math.log((1 - 2 * P - Q) * math.sqrt(1 - 2 * Q))
                except ValueError:
                    d = 3.0                      # saturated
            elif model == "jc69":
                p = diff / comp
                d = -0.75 * math.log(1 - 4 * p / 3) if p < 0.74 else 3.0
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
