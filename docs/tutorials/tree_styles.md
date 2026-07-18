# Tutorial: drawing styles beyond the plain tree

Four figure types that show up constantly in the literature but need more than
a bare phylogram. All four are in `examples/tree_styles_demo.py`.

## Collapsed clades

Compress a big tree down to the clades you actually want to discuss. Each
collapsed clade becomes a triangle whose two sides reach its **nearest** and
**farthest** hidden leaf, so the wedge shows how deep and how ragged the hidden
group is — the convention iTOL uses.

```python
node = tree.get_mrca(cyanobacteria)
pt.collapse_clade(tree, node, name="Cyanobacteria (37)")

(pt.TreeFigure(tree)
    .collapsed_clades()
    .tip_labels()
).save("collapsed.pdf")
```

`collapse_clade()` **modifies the tree in place** — the clade's children are
dropped and a summary (`n`, `near`, `far`, `leaves`) is written to
`node.data["_collapsed"]`. Work on a copy to keep the original:

```python
display = pt.Tree.from_newick(tree.write())
```

!!! warning "Check monophyly first"
    `get_mrca(taxa)` returns the *smallest clade containing* those taxa, which
    may contain others too. Collapsing it then silently swallows taxa that do
    not belong to the group you named. Verify before collapsing:

    ```python
    node = tree.get_mrca(taxa)
    if set(node.leaf_names()) == set(taxa):
        pt.collapse_clade(tree, node, name="…")
    else:
        print("not monophyletic; also contains",
              sorted(set(node.leaf_names()) - set(taxa)))
    ```

    The demo does exactly this, and on the bundled 16S tree it correctly
    refuses to collapse Euryarchaeota (*Saccharolobus*, a Thermoproteota, falls
    inside it).

Options: `color=` (literal or a data column), `scale_height=True` to let the
triangle's width grow with the number of hidden tips, `edgecolor=`. Works on
circular layouts too.

## Node interval bars (95% HPD)

The standard annotation on a dated Bayesian tree: a bar across each node
spanning its age credible interval, as FigTree's "node bars" and ggtree's
`geom_range` draw it.

```python
(pt.TreeFigure(dated_tree)
    .node_bars(lower="height_95_lower", upper="height_95_upper")
    .time_axis(geo=True)
    .tip_labels()
).save("dated.pdf")
```

Values are read as **ages** on the same scale as `time_axis()` — distance back
from `present`, increasing toward the root — so a node's bar lines up with the
time axis underneath it. Pass `as_age=False` if your two keys already hold plot
x coordinates. Rectangular layouts only (the bar runs along the time axis).

## Connections between tips

Curved links between arbitrary pairs — horizontal gene transfer, gene sharing,
co-occurrence, host-symbiont pairings. This is iTOL's `DATASET_CONNECTION`.

```python
pairs = [("Escherichia_coli", "Salmonella_enterica", 0.82), ...]

(pt.TreeFigure(tree, layout="circular")
    .connections(pairs, color="value")     # or a literal colour
    .tip_labels()
).save("hgt.pdf")
```

`pairs` is an iterable of `(name1, name2)` or `(name1, name2, value)`, or a
DataFrame with those columns. On a **circular** layout the curves bend toward
the centre (iTOL's `CENTER_CURVES`), which is what keeps a dense set readable;
on a **rectangular** one they bow out past the tips so they clear the tree.
`curvature=0` gives straight lines. Unknown names raise rather than being
silently dropped.

## DensiTree

A whole set of trees drawn on top of each other, so topological uncertainty is
visible instead of hidden behind one summary tree — a posterior sample, a
bootstrap set, or a collection of gene trees.

```python
pt.DensiTreeFigure(trees).titled("500 posterior trees").save("cloud.pdf")
```

Shared structure accumulates into dark, confident edges; conflicting
placements stay faint and fan out. Opacity defaults to a value that scales with
the size of the set.

The trees only line up if their tips share an order, so each is first rotated
to match a reference (`consensus=`, default the first tree) via
[`untangle`](tanglegram.md) — rotation reorders children without touching
topology or branch lengths, so aligning the cloud never changes what any tree
says. Pass `align=False` to skip it.

## Scale bar

A compact branch-length scale, ggtree's `geom_treescale`:

```python
pt.TreeFigure(tree).tip_labels().scale_bar()          # auto round length
pt.TreeFigure(tree).scale_bar(length=0.05, label="0.05 subs/site")
```

Unlike `time_axis()` this assumes nothing about branch lengths being time and
works on any layout, which is what a plain substitutions/site phylogram needs.
