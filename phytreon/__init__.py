"""phytreon -- phylogenetic trees and publication figures in Python.

Quick start
-----------
>>> import phytreon as pt
>>> tr = pt.datasets.primates()
>>> (pt.TreeFigure(tr)
...     .tip_labels()
...     .support_labels()).save("tree.pdf")

The package is layered:

* ``phytreon.core``        -- ``Tree`` / ``Node`` data model and I/O
* ``phytreon.layout``      -- topology -> display coordinates
* ``phytreon.infer``       -- distance-based inference (NJ / UPGMA), ML, parsimony
  (from sequences or a discrete character/trait matrix)
* ``phytreon.comparative`` -- ancestral state reconstruction
* ``phytreon.plot``        -- the ``TreeFigure`` builder + matplotlib/plotly backends
"""
from __future__ import annotations

from . import core, layout, infer, comparative, plot, datasets, treeops
from .core import Tree, Node
from .plot import TreeFigure
from .infer import (
    neighbor_joining, upgma, tree_from_alignment, distance_matrix,
    Alignment, align, read_fasta, read_character_matrix, trim, bootstrap_support,
    build_tree, infer_ml, ml_tree, log_likelihood, model_finder,
    parsimony_tree, parsimony_score,
    read_allele_table, read_mutation_matrix, sankoff_score, camin_sokal_score, lineage_tree,
    reconstruct_ancestral_mutations,
    expression_distance_matrix, expression_dendrogram,
)
from .comparative import ace_parsimony, ace_ml, ace_continuous, stochastic_map
from .treeops import (
    rotate, flip, swap_children, ladderize, collapse_low_support,
    scale_clade, cut_tree, midpoint_root, group_clade, group_otu,
    robinson_foulds, prune_to_taxa,
)

__version__ = "0.1.1"

__all__ = [
    "core", "layout", "infer", "comparative", "plot", "datasets", "treeops",
    "Tree", "Node",
    "TreeFigure",
    "neighbor_joining", "upgma", "tree_from_alignment", "distance_matrix",
    "Alignment", "align", "read_fasta", "read_character_matrix", "trim", "bootstrap_support",
    "build_tree", "infer_ml", "ml_tree", "log_likelihood", "model_finder",
    "parsimony_tree", "parsimony_score", "robinson_foulds",
    "read_allele_table", "read_mutation_matrix", "sankoff_score", "camin_sokal_score", "lineage_tree",
    "reconstruct_ancestral_mutations",
    "expression_distance_matrix", "expression_dendrogram",
    "ace_parsimony", "ace_ml", "ace_continuous", "stochastic_map",
    "rotate", "flip", "swap_children", "ladderize", "collapse_low_support",
    "scale_clade", "cut_tree", "midpoint_root", "group_clade", "group_otu",
    "prune_to_taxa",
]
