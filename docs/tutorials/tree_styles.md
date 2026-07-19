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

Read a BEAST/MrBayes summary tree with `fmt="beast"` and it works directly —
that reader keeps the per-node estimates, and its `{lower, upper}` intervals
land on the keys `node_bars()` reads by default:

```python
tree = pt.Tree.read("beast_summary.tre", fmt="beast")
tree.root.data["posterior"]          # 0.97
tree.root.data["height_95_lower"]    # 7.0

(pt.TreeFigure(tree)
    .node_bars()                     # picks up height_95_lower/_upper
    .time_axis(geo=True)
    .tip_labels()
).save("dated.pdf")
```

A plain `fmt="nexus"` read keeps only the topology — the annotations are
dropped — so the bars would have nothing to draw. Other keys work too:

```python
.node_bars(lower="length_95_lower", upper="length_95_upper")
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

## Grey the default state, colour the exceptions

When one level of a category covers most of the tree — `GTDB_reference`,
`wild_type`, `present` — giving it a saturated colour spends most of the
figure's ink on its least informative part and buries the handful of tips that
actually carry the finding. `baseline=` renders those levels neutral grey, and
they stop consuming palette slots so the remaining levels keep the strongest
hues:

```python
context = ["GTDB reference", "Near-known target", "Intermediate",
           "Deep or mixed", "Unresolved"]

(pt.TreeFigure(tree, layout="circular")
    .branches(color="#707070", size=0.4)     # neutral baseline for the tree
    .ring(meta, columns=["Phylogenetic context"],
          width=0.07,                        # a narrow band, not a heavy hoop
          baseline="GTDB reference",         # the expected state recedes
          order=context,                     # known -> unresolved, not A-Z
          leaders=True)
).save("context.pdf")
```

`order=` matters more than it looks: categorical levels are otherwise sorted
alphabetically, which almost never matches a meaningful progression.

Legend keys follow the mark they stand for — filled layers (rings, heatmaps)
get square swatches, marker layers get dots.

### Prefer clade annotation over a per-tip ring

If a variable is essentially clade-structured (a phylum that forms one
monophyletic block), a per-tip ring repeats the same colour across hundreds of
sectors to convey one fact. A shaded sector plus an arc label says it once:

```python
for name, taxa in phyla.items():
    node = tree.get_mrca(taxa)
    if set(node.leaf_names()) == set(taxa):        # verify monophyly first
        fig.highlight(node=node, fill=colour[name], alpha=0.10)
        fig.clade_label(name, node=node)
```

That frees a whole ring, and the tree usually gets much cleaner. Where the
group is *not* monophyletic a ring is the honest choice — a shaded sector
would imply a clade that the tree does not support.

## Scale bar

A compact branch-length scale, ggtree's `geom_treescale`:

```python
pt.TreeFigure(tree).tip_labels().scale_bar()          # auto round length
pt.TreeFigure(tree).scale_bar(length=0.05, label="0.05 subs/site")
```

Unlike `time_axis()` this assumes nothing about branch lengths being time and
works on any layout, which is what a plain substitutions/site phylogram needs.
