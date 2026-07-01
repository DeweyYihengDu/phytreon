"""Maximum-parsimony phylogenetics (pure Python).

Fitch parsimony scoring, vectorised over site patterns using 4-bit state
masks (numpy bitwise ops -- fast), plus NNI hill-climbing to find the
minimum-changes topology.  Reuses the NNI machinery from :mod:`ml_native`.
"""
from __future__ import annotations

from typing import Dict, Optional

from ..core.tree import Tree
from .align import Alignment

# nucleotide -> 4-bit mask (A,C,G,T); ambiguous/gap -> all states
_MASK = {"A": 1, "C": 2, "G": 4, "T": 8, "U": 8}
_ALL = 15


def _encode(aln: Alignment):
    """Compressed patterns: (ntip, npat) uint8 masks + pattern weights."""
    import numpy as np
    from collections import Counter
    ncol = aln.ncol
    counts: Counter = Counter()
    order = []
    for j in range(ncol):
        col = tuple(s[j].upper() for s in aln.seqs)
        if col not in counts:
            order.append(col)
        counts[col] += 1
    masks = np.zeros((aln.nseq, len(order)), dtype=np.uint8)
    for p, col in enumerate(order):
        for i, ch in enumerate(col):
            masks[i, p] = _MASK.get(ch, _ALL)
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
            cache[id(node)] = masks[idx[node.name]].astype(np.uint8)
        else:
            acc = cache[id(node.children[0])].copy()
            for c in node.children[1:]:
                cm = cache[id(c)]
                inter = acc & cm
                empty = inter == 0
                changes += empty * weights          # a change where sets disjoint
                acc = np.where(empty, acc | cm, inter).astype(np.uint8)
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
        col = [("T" if c == "U" else c) for c in
               (s[j].upper() for s in aln.seqs) if c in "ACGTU"]
        if not col:
            continue
        cnt = Counter(col)
        m += len(cnt) - 1               # min: each extra state = >=1 change
        g += len(col) - max(cnt.values())   # max: star-tree changes
    return m, g
