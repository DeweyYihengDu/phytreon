# Tutorial: annotated trees

phytreon uses a fluent builder: start with `TreeFigure(tree)` and chain methods
to add elements.

## Rectangular tree with metadata

```python
import phytreon as pt
tr = pt.datasets.primates()
tr.join_data(pt.datasets.primates_metadata().reset_index(), on="name")
apes = tr.get_mrca(["Human", "Gibbon"])

(pt.TreeFigure(tr)
    .highlight(node=apes, fill="#cfe8f3")           # shade a clade
    .tip_points(color="habitat", shape="habitat")   # colour + shape mapping
    .tip_labels()
    .support_labels()
    .clade_label("Apes", node=apes)
).save("tree.pdf")   # .html -> interactive plotly
```

## Right-side tracks (rectangular)

```python
(pt.TreeFigure(tr).tip_labels()
    .heatmap(meta[["phylum"]])          # categorical tile track (own scale)
    .heatmap(meta[["group"]])           # another track, stacks rightward
    .bar_track(meta, "length")          # horizontal bar track
    .alignment(alignment))              # residue matrix (raster)
```

## Circular tree with metadata rings

```python
(pt.TreeFigure(tr, layout="circular", extent=320)
    .branches(color="lineage")          # branches coloured by clade
    .ring(meta, columns=["serotype"], geom="tile")
    .ring(meta, columns=["AMR_score"], geom="bar")
    .tip_labels(max_labels=60)          # thin labels on big trees
).save("rings.png")
```

Continuous columns get a colorbar; categorical columns get a legend.
Tracks, labels and legends are placed so nothing overlaps.
