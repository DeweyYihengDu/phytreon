# Tutorial: expression-similarity dendrograms

!!! warning "This is not a phylogenetic tree"
    Grouping cells by how similar their expression of a gene (or genes) is
    reflects **cell state/type**, not shared ancestry. Two unrelated cells in
    the same state can end up as "sisters" here; two truly related cells
    that have diverged in expression can end up far apart. For real
    single-cell lineage reconstruction (actual evolutionary relationships
    between cells), see the
    [lineage tracing tutorial](lineage_tracing.md) instead -- CRISPR scars
    or somatic mutations, not expression.

That caveat aside, a hierarchical-clustering dendrogram of transcriptional
similarity is still a common and legitimate thing to visualize. phytreon
builds one from a numeric expression matrix (cells as rows, genes as
columns), reusing the same `neighbor_joining`/`upgma` builders used
everywhere else in the package -- there's no new tree-building algorithm
here, just a distance metric for continuous data.

## From an expression matrix to a dendrogram

```python
import pandas as pd
import phytreon as pt

expr = pd.read_csv("expression.csv", index_col=0)   # cells x genes
tree = pt.expression_dendrogram(expr, genes=["CD3D"])   # one marker gene

print(tree.data["tree_type"])   # "expression_similarity_dendrogram"
(pt.TreeFigure(tree).tip_labels()).save("cd3d_similarity.pdf")
```

`genes` restricts the comparison to one gene or a small combination, rather
than the whole transcriptome -- matching the common case of grouping cells
by a marker gene or a small panel. Omit it to use every column in `expr`.

## Choosing a distance metric

`metric` is any value accepted by `scipy.spatial.distance.pdist`
(`"correlation"` by default, or `"euclidean"`, `"cosine"`, ...). Pearson
correlation is mathematically undefined for a single gene (there's nothing
to correlate *across*), so a lone `genes=[...]` of length 1 automatically
falls back to `"euclidean"` -- check `tree.data["metric"]` to see which one
actually ran:

```python
tree = pt.expression_dendrogram(expr, genes=["CD3D"])
print(tree.data["metric"])   # "euclidean", not the "correlation" default
```

## Lower-level access

`expression_distance_matrix(expr, genes=..., metric=...)` returns the
`(names, matrix)` pair directly, if you want the distances without a tree
attached -- e.g. to feed a different clustering method.
