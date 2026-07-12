"""Expression-similarity dendrograms -- explicitly not phylogenetic trees.

See phytreon/infer/expression.py's module docstring for why: this groups
cells by transcriptional similarity (cell state), which is a fundamentally
different question from lineage_tree()'s genotype-based common ancestry
(tests/test_lineage.py). Purely additive; reuses the existing
neighbor_joining()/upgma() with a new distance metric for continuous data.
"""
from __future__ import annotations

import pandas as pd
import pytest

from phytreon.infer.expression import expression_distance_matrix, expression_dendrogram


def test_single_gene_euclidean_distance_is_exact():
    # hand-computable: |1.0-1.2|=0.2, |8.0-8.5|=0.5, cross-cluster ~7
    df = pd.DataFrame({"CD3D": [1.0, 1.2, 8.0, 8.5]}, index=["a", "b", "c", "d"])
    names, mat = expression_distance_matrix(df, genes=["CD3D"])
    by_name = dict(zip(names, mat))
    assert by_name["a"][names.index("b")] == pytest.approx(0.2)
    assert by_name["c"][names.index("d")] == pytest.approx(0.5)


def test_multi_gene_correlation_distance_matches_hand_calculation():
    # cellB's expression is a scaled copy of cellA's (r=1, distance=0);
    # cellC's is cellA's exact reversal (r=-1, distance=2) -- both
    # hand-verifiable Pearson correlation cases (pdist distance = 1 - r).
    # row cellA=[1,2,3], cellB=[2,4,6] (scaled copy), cellC=[3,2,1] (reversal)
    df = pd.DataFrame({
        "geneP": [1., 2., 3.],
        "geneQ": [2., 4., 2.],
        "geneR": [3., 6., 1.],
    }, index=["cellA", "cellB", "cellC"])
    names, mat = expression_distance_matrix(df, metric="correlation")
    by_name = dict(zip(names, mat))
    assert by_name["cellA"][names.index("cellB")] == pytest.approx(0.0, abs=1e-9)
    assert by_name["cellA"][names.index("cellC")] == pytest.approx(2.0)


def test_single_gene_correlation_auto_falls_back_to_euclidean():
    df = pd.DataFrame({"CD3D": [1.0, 1.2, 8.0, 8.5]}, index=["a", "b", "c", "d"])
    tree = expression_dendrogram(df, genes=["CD3D"], metric="correlation")
    assert tree.data["metric"] == "euclidean"     # not the requested "correlation"


def test_expression_dendrogram_recovers_clusters_and_is_flagged():
    df = pd.DataFrame({"CD3D": [1.0, 1.2, 8.0, 8.5]}, index=["a", "b", "c", "d"])
    tree = expression_dendrogram(df, genes=["CD3D"])
    assert set(tree.leaf_names()) == {"a", "b", "c", "d"}
    leaf_sets = [frozenset(n.leaf_names()) for n in tree.traverse() if not n.is_leaf]
    assert frozenset({"a", "b"}) in leaf_sets
    assert frozenset({"c", "d"}) in leaf_sets
    assert tree.data["tree_type"] == "expression_similarity_dendrogram"
    assert tree.data["genes"] == ["CD3D"]


def test_expression_dendrogram_nj_method():
    df = pd.DataFrame({"g1": [1.0, 1.1, 9.0, 9.2], "g2": [2.0, 2.1, 8.0, 8.3]},
                      index=["a", "b", "c", "d"])
    tree = expression_dendrogram(df, method="nj")
    assert set(tree.leaf_names()) == {"a", "b", "c", "d"}


def test_expression_dendrogram_rejects_unknown_method():
    df = pd.DataFrame({"g1": [1.0, 2.0, 3.0]}, index=["a", "b", "c"])
    with pytest.raises(ValueError, match="method"):
        expression_dendrogram(df, method="ml")


def test_expression_distance_matrix_rejects_no_genes_selected():
    df = pd.DataFrame({"g1": [1.0, 2.0]}, index=["a", "b"])
    with pytest.raises(ValueError, match="genes"):
        expression_distance_matrix(df, genes=[])
