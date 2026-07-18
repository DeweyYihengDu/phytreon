# Tutorial: comparing two trees (tanglegrams)

A **tanglegram** faces two trees of the same taxa at each other and joins the
shared tips, so disagreements show up as crossing links. It is the standard
picture for questions like *"does the transcriptome tree agree with the
genome tree?"*, *"do these two inference methods give the same answer?"* or
*"do the parasites track their hosts?"*.

```python
import phytreon as pt

genome = pt.build_tree("genome.fasta", method="nj")
transcriptome = pt.build_tree("transcriptome.fasta", method="nj")

fig = pt.TangleFigure(genome, transcriptome,
                      titles=("genome", "transcriptome"))
fig.untangle()                              # line the tips up first
fig.connect(highlight_discordant=True)      # red = this taxon disagrees
fig.save("tanglegram.pdf")
```

## Untangle before you interpret

Rotating a node swaps the order of its children. That slides a clade up or
down the page **without changing the topology or any branch length** — the
tree says exactly the same thing, it just reads differently. So the raw number
of crossings is partly an artefact of how each tree happened to be drawn.

`untangle()` hill-climbs over single rotations to remove the crossings that
are mere drawing accidents; whatever still crosses afterwards is real
disagreement.

```python
fig = pt.TangleFigure(t1, t2)
print(fig.crossings())          # e.g. 111  -- mostly drawing artefact
fig.untangle()
print(fig.crossings())          # e.g. 17   -- real conflict
```

| `fix` | effect |
|---|---|
| `"left"` (default) | keep the left tree's tip order, rotate only the right |
| `"right"` | the reverse |
| `None` | alternate between the two; usually untangles further, but neither tree keeps its original order |

!!! warning "Zero crossings does not mean the trees are identical"
    Crossings measure *tip order*. Two trees can resolve splits differently
    and still admit one common tip order, giving zero crossings. Always read
    the Robinson-Foulds distance alongside it:

    ```python
    pt.robinson_foulds(t1, t2, normalized=True)   # rotation-independent
    pt.crossing_number(t1, t2)                    # what the plot shows
    ```

    In `examples/tanglegram_demo.py`, neighbour joining and parsimony untangle
    to **0 crossings** while RF is still 0.067 — they agree on the order, not
    on every split.

## Styling each side

`fig.left` and `fig.right` are ordinary
[`TreeFigure`](../api.md) objects, so every element works on either tree:

```python
fig.left.tip_points(color="phylum", size=6)
fig.right.support_labels()
fig.left.highlight(node=mrca, fill="#eef3f8")
```

Pass ready-made figures when you want full control of both sides:

```python
left = pt.TreeFigure(t1).tip_labels(italic=True).tip_points(color="host")
right = pt.TreeFigure(t2).tip_points(color="host")
fig = pt.TangleFigure(left, right)
```

## Links

```python
fig.connect(color="#b0b0b0", width=0.7, dash="dot")   # plain styling
fig.connect(color="phylum")                           # colour by a data column
fig.connect(highlight_discordant=True)                # red = crossing link
```

`color=` accepts a column on the **left** tree's tips and emits a legend, which
answers a different question from `highlight_discordant`: not *which* taxa move,
but whether the movement is confined to one group or cuts across the tree.

## Layout knobs

Tip labels and links share the middle band. Text is measured in points and the
tree in data units, and the two are only related once the figure exists, so the
band is **sized automatically** from the longest tip label — long species names
widen it rather than overrunning the links. Three knobs override that:

- `tip_labels=` — `"both"` (default), `"left"`, `"right"` or `False`. Both
  sides are named by default so each tree can be read on its own.
- `gap=` — set the band yourself, as a fraction of the two trees' combined
  depth, instead of the automatic estimate.
- `connect(inset=(left, right))` — fraction of the band to keep clear at each
  end. `(0, 0)` runs the links tip to tip.

The automatic estimate assumes the default figure width; if you pass a very
different `figsize` to `save()`, set `gap=` too.

Trees with only partly overlapping taxa are fine: unmatched tips are still
drawn, just left unlinked.

```python
fig = pt.TangleFigure(t1, t2, tip_labels="left", gap=0.8)
```
