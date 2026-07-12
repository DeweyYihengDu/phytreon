"""Transcriptional similarity dendrograms from single-cell expression data.

**This is explicitly not phylogenetics.** Two cells with similar expression
of a gene (or a small combination of genes) are similar in cell *state*, not
necessarily in ancestry -- two unrelated cells of the same type/state look
identical here, while two truly sibling cells that have diverged in
expression can end up far apart. For real single-cell lineage
reconstruction (actual evolutionary/genealogical relationships between
cells) use :mod:`phytreon.infer.lineage` instead -- CRISPR scars
(:func:`~phytreon.infer.lineage.read_allele_table`) or somatic mutations
(:func:`~phytreon.infer.lineage.read_mutation_matrix`) feeding
:func:`~phytreon.infer.lineage.camin_sokal_score`/
:func:`~phytreon.infer.lineage.lineage_tree`.

Grouping cells by expression similarity (a marker gene, or a small
combination) is still a common and legitimate thing to visualize -- it's
just a hierarchical-clustering *dendrogram* of cell state, not a tree of
common descent, and is named/documented accordingly throughout: the
function names avoid "tree"/"phylogen-", matching how scipy/seaborn name the
exact same distinction, and the resulting :class:`~phytreon.core.tree.Tree`
carries ``tree.data["tree_type"] = "expression_similarity_dendrogram"`` as a
machine-readable flag.

Purely additive: reuses the existing, alphabet-agnostic
:func:`phytreon.infer.distance.neighbor_joining`/
:func:`phytreon.infer.distance.upgma` (they already accept any raw square
distance matrix) -- no new tree-building algorithm, just a distance metric
for continuous expression data that didn't exist in the package before.
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple, Union

from ..core.tree import Tree


def _effective_metric(ncols: int, metric: str) -> str:
    """"correlation" is mathematically undefined for a single gene (there's
    nothing to correlate *across*) -- a lone gene automatically falls back
    to "euclidean". Shared by expression_distance_matrix/expression_dendrogram
    so the metric actually used can never drift out of sync with what gets
    reported."""
    return "euclidean" if ncols == 1 and metric == "correlation" else metric


def expression_distance_matrix(expr: Union[str, "object"], *,
                               genes: Optional[Sequence[str]] = None,
                               metric: str = "correlation"
                               ) -> Tuple[List[str], List[List[float]]]:
    """Pairwise transcriptional dissimilarity between cells/samples, from a
    numeric expression matrix (cells as rows, genes as columns; CSV path or
    an existing :class:`pandas.DataFrame`).

    ``genes`` optionally restricts the comparison to one gene or a small
    combination (rather than the whole transcriptome) -- e.g.
    ``genes=["CD3D"]`` groups cells purely by that one marker's expression.
    ``metric`` is any value accepted by
    :func:`scipy.spatial.distance.pdist` (``"correlation"``, ``"euclidean"``,
    ``"cosine"``, ...); see :func:`_effective_metric` for the single-gene
    fallback.
    """
    import pandas as pd
    from scipy.spatial.distance import pdist, squareform

    df = expr.copy() if isinstance(expr, pd.DataFrame) else pd.read_csv(expr, index_col=0)
    if genes is not None:
        df = df[list(genes)]
    if df.shape[1] == 0:
        raise ValueError("no genes selected")
    metric = _effective_metric(df.shape[1], metric)

    names = [str(n) for n in df.index]
    if len(names) < 2:
        return names, [[0.0] * len(names) for _ in names]
    mat = squareform(pdist(df.to_numpy(dtype=float), metric=metric))
    return names, mat.tolist()


def expression_dendrogram(expr: Union[str, "object"], *,
                          genes: Optional[Sequence[str]] = None,
                          metric: str = "correlation",
                          method: str = "upgma") -> Tree:
    """Hierarchical-clustering dendrogram of transcriptional similarity for
    one gene or a small combination of genes.

    **This is not a phylogenetic tree.** It groups cells by how alike their
    expression is (cell state/type), not by shared ancestry. For real
    single-cell lineage reconstruction, use
    :func:`phytreon.infer.lineage.lineage_tree` on CRISPR scar data
    (:func:`~phytreon.infer.lineage.read_allele_table`) or somatic-mutation
    data (:func:`~phytreon.infer.lineage.read_mutation_matrix`) instead.

    ``method`` is ``"upgma"`` (default -- ultrametric, the conventional
    choice for expression dendrograms) or ``"nj"``. Result carries
    ``tree.data["tree_type"] = "expression_similarity_dendrogram"``.
    """
    import pandas as pd

    from .distance import neighbor_joining, upgma

    builder = {"upgma": upgma, "nj": neighbor_joining}.get(method)
    if builder is None:
        raise ValueError(f"unknown method {method!r}; use 'upgma' or 'nj'")

    df = expr if isinstance(expr, pd.DataFrame) else pd.read_csv(expr, index_col=0)
    ncols = len(genes) if genes is not None else df.shape[1]
    effective_metric = _effective_metric(ncols, metric)

    names, mat = expression_distance_matrix(df, genes=genes, metric=metric)
    tree = builder(names, mat)
    tree.data["tree_type"] = "expression_similarity_dendrogram"
    tree.data["metric"] = effective_metric
    if genes is not None:
        tree.data["genes"] = list(genes)
    return tree
