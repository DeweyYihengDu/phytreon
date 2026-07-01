"""Tree inference: alignment, trimming, distance/ML methods, bootstrap."""
from .align import Alignment, align, align_external, read_fasta
from .trim import trim, trim_terminal_gaps, column_gap_fraction, column_conservation
from .distance import (
    neighbor_joining,
    upgma,
    distance_matrix,
    tree_from_alignment,
)
from .bootstrap import (bootstrap_support, p_distance_matrix, distance_matrix_model,
                        nj_builder, upgma_builder)
from .ml import infer_ml, infer_iqtree, infer_fasttree
from .ml_native import ml_tree, log_likelihood, model_finder
from .parsimony import parsimony_tree, parsimony_score
from .pipeline import build_tree

__all__ = [
    "Alignment", "align", "align_external", "read_fasta",
    "trim", "trim_terminal_gaps", "column_gap_fraction", "column_conservation",
    "neighbor_joining", "upgma", "distance_matrix", "tree_from_alignment",
    "bootstrap_support", "p_distance_matrix", "nj_builder", "upgma_builder",
    "infer_ml", "infer_iqtree", "infer_fasttree",
    "ml_tree", "log_likelihood", "model_finder",
    "parsimony_tree", "parsimony_score",
    "distance_matrix_model",
    "build_tree",
]
