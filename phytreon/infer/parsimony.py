"""Maximum-parsimony phylogenetics (pure Python).

Fitch parsimony scoring, vectorised over site patterns using integer state
masks (numpy bitwise ops -- fast), plus NNI hill-climbing to find the
minimum-changes topology.  Reuses the NNI machinery from :mod:`ml_native`.

States are derived independently for each site from whatever characters
actually appear there, so this works for nucleotide sequences, amino-acid
sequences, and arbitrary discrete character/trait matrices alike (see
:func:`phytreon.read_character_matrix`) -- not just A/C/G/T. ``-``, ``.``,
``?`` and ``N`` are always treated as missing/ambiguous: they match
whichever state a site's neighbours settle on, so missing data never
forces a spurious change.
"""
from __future__ import annotations

from typing import Dict, Optional

from ..core.tree import Tree
from .align import Alignment

_MISSING = frozenset({"-", ".", "?", "N"})
_MAX_STATES = 32                        # bits available in the uint32 mask


def _normalize(ch: str) -> str:
    return "T" if ch == "U" else ch     # RNA U == DNA T for scoring purposes


def _encode(aln: Alignment):
    """Compressed patterns: (ntip, npat) uint32 state masks + pattern weights.

    Each pattern's state space is derived independently from the symbols
    observed at that site (excluding ``_MISSING``); missing symbols get the
    bitwise OR of every state seen at that site, so they are compatible with
    any resolution.
    """
    import numpy as np
    from collections import Counter
    ncol = aln.ncol
    counts: Counter = Counter()
    order = []
    for j in range(ncol):
        col = tuple(_normalize(s[j].upper()) for s in aln.seqs)
        if col not in counts:
            order.append(col)
        counts[col] += 1
    masks = np.zeros((aln.nseq, len(order)), dtype=np.uint32)
    for p, col in enumerate(order):
        states = sorted({ch for ch in col if ch not in _MISSING})
        if not states:
            masks[:, p] = 1              # no information at this site -> free
            continue
        if len(states) > _MAX_STATES:
            raise ValueError(
                f"site has {len(states)} distinct states {states!r}; "
                f"parsimony supports at most {_MAX_STATES} per site")
        bit = {ch: np.uint32(1) << i for i, ch in enumerate(states)}
        all_bits = np.uint32((1 << len(states)) - 1)
        for i, ch in enumerate(col):
            masks[i, p] = bit.get(ch, all_bits)
    weights = np.array([counts[c] for c in order], dtype=float)
    return masks, weights


def parsimony_score(tree: Tree, aln: Alignment, data=None) -> float:
    """Total Fitch parsimony score (number of changes) of ``tree``."""
    import numpy as np
    masks, weights = data if data is not None else _encode(aln)
    idx = {n: i for i, n in enumerate(aln.names)}
    npat = masks.shape[1]
    changes = np.zeros(npat, dtype=float)
    cache: Dict[int, "np.ndarray"] = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf:
            cache[id(node)] = masks[idx[node.name]].astype(np.uint32)
        else:
            acc = cache[id(node.children[0])].copy()
            for c in node.children[1:]:
                cm = cache[id(c)]
                inter = acc & cm
                empty = inter == 0
                changes += empty * weights          # a change where sets disjoint
                acc = np.where(empty, acc | cm, inter).astype(np.uint32)
            cache[id(node)] = acc
    return float(changes.sum())


def parsimony_tree(aln: Alignment, start: Optional[Tree] = None,
                   search: bool = True, max_sweeps: int = 50) -> Tree:
    """Minimum-parsimony tree by NNI hill-climbing (from an NJ start).

    Result carries ``tree.data['parsimony_score']`` plus the consistency index
    (``ci``) and retention index (``ri``).  Note: NNI returns one locally
    optimal tree; equally-parsimonious alternatives are not enumerated
    (``tree.data['mp_search'] = 'NNI-local'``).
    """
    from .bootstrap import nj_builder
    from ..treeops import midpoint_root
    from ._search import internal_edges, nni_neighbors

    data = _encode(aln)
    tree = start or midpoint_root(nj_builder(aln))
    best = parsimony_score(tree, aln, data)

    if search:
        for _ in range(max_sweeps):
            improved = False
            for node in internal_edges(tree):
                for swap in list(nni_neighbors(node)):
                    swap()
                    sc = parsimony_score(tree, aln, data)
                    if sc < best - 1e-9:
                        best = sc
                        improved = True
                        break
                    swap()                              # involution -> undo
            if not improved:
                break

    m, g = _parsimony_bounds(aln)                       # min & max possible steps
    tree.data["parsimony_score"] = best
    tree.data["ci"] = (m / best) if best > 0 else None              # consistency
    tree.data["ri"] = ((g - best) / (g - m)) if g > m else None     # retention
    tree.data["mp_search"] = "NNI-local"
    return tree


def _parsimony_bounds(aln):
    """(min, max) possible parsimony steps summed over characters."""
    from collections import Counter
    m = g = 0
    for j in range(aln.ncol):
        col = [_normalize(c) for c in (s[j].upper() for s in aln.seqs)
               if c not in _MISSING]
        if not col:
            continue
        cnt = Counter(col)
        m += len(cnt) - 1               # min: each extra state = >=1 change
        g += len(col) - max(cnt.values())   # max: star-tree changes
    return m, g
