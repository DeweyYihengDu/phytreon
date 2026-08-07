"""Distance-based tree inference (neighbour-joining / UPGMA).

Thin wrappers over Biopython's ``DistanceTreeConstructor`` that accept a
plain labels + matrix pair (or an alignment) and return a phytreon
:class:`~phytreon.core.tree.Tree`.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

from ..core.tree import Node, Tree
from ..core.io import from_biopython


def _to_bio_dm(names: Sequence[str], matrix):
    """Build a Biopython lower-triangular DistanceMatrix."""
    from Bio.Phylo.TreeConstruction import DistanceMatrix
    n = len(names)
    lower = [[float(matrix[i][j]) for j in range(i + 1)] for i in range(n)]
    return DistanceMatrix(list(names), lower)


def _clamp_negative(tree: Tree) -> Tree:
    """Set negative branch lengths to 0 (a standard fix for the NJ artefact)."""
    for node in tree.traverse():
        if node.length is not None and node.length < 0:
            node.length = 0.0
    return tree


def neighbor_joining(names: Sequence[str], matrix, nonneg: bool = True) -> Tree:
    """Neighbor-joining tree from a square distance matrix.

    ``nonneg=True`` (default) clamps NJ's negative branch lengths to 0.
    """
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
    dm = _to_bio_dm(names, matrix)
    tree = DistanceTreeConstructor().nj(dm)
    _strip_inner_names(tree)
    out = from_biopython(tree)
    return _clamp_negative(out) if nonneg else out


def upgma(names: Sequence[str], matrix, nonneg: bool = True) -> Tree:
    """UPGMA (ultrametric) tree from a square distance matrix."""
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
    dm = _to_bio_dm(names, matrix)
    tree = DistanceTreeConstructor().upgma(dm)
    _strip_inner_names(tree)
    out = from_biopython(tree)
    return _clamp_negative(out) if nonneg else out


def _mean_between(matrix, idx_a: List[int], idx_b: List[int]) -> float:
    return sum(matrix[i][j] for i in idx_a for j in idx_b) / (len(idx_a) * len(idx_b))


def constrained_nj(names: Sequence[str], matrix, groups: Dict[str, object],
                   nonneg: bool = True) -> Tree:
    """Neighbour-joining with every named group forced monophyletic.

    ``groups`` maps a tip name to a label (its genus, say); a tip missing
    from ``groups``, or mapped to ``None``, is left free -- on its own,
    joined to whatever NJ finds closest, exactly as if there were no
    constraint for it at all.

    Two passes of ordinary NJ, not one search under a constraint. NJ first
    runs *inside* each group, on only the distances between that group's own
    tips, giving each group its own internally-resolved subtree. A second NJ
    pass then places the groups relative to each other -- and any free tip,
    still its own trivial group of one -- on a reduced matrix, one row per
    group, each entry the mean of the real distances between the two groups'
    tips. Every group's subtree, midpoint-rooted for a defensible attachment
    point, then replaces its placeholder leaf on that backbone.

    This *forces* the monophyly a constrained ML search
    (:func:`~phytreon.infer.constraint.constraint_tree`, for IQ-TREE's ``-g``
    /RAxML-NG's ``--tree-constraint``) only asks for: there is no way for the
    result to show a group as anything but monophyletic, because tips from
    two different groups never end up on opposite sides of an NJ split in
    either pass -- not even where the sequence data itself disagrees. Reach
    for this when the taxonomy should win outright, such as a display tree
    organised by genus; reach for a constrained ML search when the sequence
    data should still have the final word on whether a genus really is
    monophyletic, and the constraint only needs to settle the ties.
    """
    index = {n: i for i, n in enumerate(names)}
    by_group: Dict[object, List[str]] = {}
    for n in names:
        label = groups.get(n)
        by_group.setdefault(n if label is None else label, []).append(n)

    if len(by_group) < 2:
        return neighbor_joining(names, matrix, nonneg=nonneg)   # nothing to graft

    subtrees: Dict[object, Tree] = {}
    for label, tips in by_group.items():
        if len(tips) == 1:
            continue                                             # the leaf itself
        if len(tips) == 2:
            a, b = tips
            d = matrix[index[a]][index[b]]
            root = Node()
            root.add_child(Node(name=a, length=d / 2))
            root.add_child(Node(name=b, length=d / 2))
            subtrees[label] = Tree(root=root)
        else:
            idx = [index[t] for t in tips]
            sub_d = [[matrix[i][j] for j in idx] for i in idx]
            from ..treeops import midpoint_root
            subtrees[label] = midpoint_root(neighbor_joining(tips, sub_d, nonneg=nonneg))

    labels = list(by_group)
    if len(labels) == 2:
        a, b = labels                       # NJ needs >= 3 taxa; place directly
        d = _mean_between(matrix, [index[t] for t in by_group[a]],
                          [index[t] for t in by_group[b]])
        backbone_root = Node()
        backbone_root.add_child(Node(name=str(a), length=d / 2))
        backbone_root.add_child(Node(name=str(b), length=d / 2))
        backbone = Tree(root=backbone_root)
    else:
        reduced = [[0.0 if gi == gj else
                   _mean_between(matrix, [index[t] for t in by_group[gi]],
                                [index[t] for t in by_group[gj]])
                   for gj in labels] for gi in labels]
        backbone = neighbor_joining([str(g) for g in labels], reduced, nonneg=nonneg)

    by_str = {str(label): sub for label, sub in subtrees.items()}
    for leaf in backbone.leaves():
        sub = by_str.get(leaf.name)
        if sub is None:
            continue                                             # free tip / singleton
        sub.root.length = leaf.length
        sub.root.parent = leaf.parent
        leaf.parent.children[leaf.parent.children.index(leaf)] = sub.root
    return backbone


def distance_matrix(alignment, model: str = "identity"):
    """Compute (names, matrix) from a Biopython ``MultipleSeqAlignment``.

    ``model`` is any name accepted by Biopython's ``DistanceCalculator``
    (e.g. ``"identity"``, ``"blastn"``, ``"blosum62"``).
    """
    from Bio.Phylo.TreeConstruction import DistanceCalculator
    dm = DistanceCalculator(model).get_distance(alignment)
    names = list(dm.names)
    n = len(names)
    mat = [[dm[i, j] for j in range(n)] for i in range(n)]
    return names, mat


def tree_from_alignment(alignment, method: str = "nj", model: str = "identity") -> Tree:
    """One-shot: alignment -> distances -> NJ/UPGMA tree."""
    names, mat = distance_matrix(alignment, model)
    if method == "nj":
        return neighbor_joining(names, mat)
    if method == "upgma":
        return upgma(names, mat)
    raise ValueError(f"unknown method {method!r}; use 'nj' or 'upgma'")


def _strip_inner_names(bp_tree) -> None:
    """Biopython names internal nodes 'Inner1'...; drop those for clean plots."""
    for clade in bp_tree.find_clades():
        if clade.name and clade.name.startswith("Inner"):
            clade.name = None
