"""Distance-based tree inference (neighbour-joining / UPGMA).

Thin wrappers over Biopython's ``DistanceTreeConstructor`` that accept a
plain labels + matrix pair (or an alignment) and return a phytreon
:class:`~phytreon.core.tree.Tree`.
"""
from __future__ import annotations

from typing import List, Sequence

from ..core.tree import Tree
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
