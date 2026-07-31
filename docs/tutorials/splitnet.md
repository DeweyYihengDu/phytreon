# Tutorial: split networks

A tree can only say that one grouping is right. Real data often supports two
groupings at once — recombination, hybridisation and incomplete lineage sorting
all leave that signature — and a tree resolves the conflict silently, by
picking a winner.

A **split network** refuses to. Every conflicting pair of splits is drawn as a
pair of parallel edges, so the conflict appears as a **box**. Data with a clean
tree signal draws as a tree; reticulate data draws a lattice, and the size of
the boxes is the size of the disagreement.

```python
net = pt.SplitNetwork.from_trees(bootstrap_trees)
net.color_by(domain, title="domain")
net.titled("bootstrap splits").save("network.pdf")
```

## Where the splits come from

```python
# from a bootstrap or posterior sample: weight = fraction of trees containing it
pt.SplitNetwork.from_trees(trees)

# from a distance matrix: NeighborNet — fit every circular split to the distances
pt.SplitNetwork.from_distances(names, matrix)

# from an alignment
pt.SplitNetwork.from_alignment(aln)
```

`from_trees` turns a bootstrap set straight into a picture of *which groupings
the replicates disagree about* — a split found in every tree draws long, one
found in half of them draws short and boxed against its rival.

`from_distances` does something that looks similar and is not. It builds a
neighbour-joining tree only to get a circular ordering out of it, then throws
the tree's splits away and fits **every** split that ordering can draw — all
`n(n-1)/2` of them — to the distances:

$$\min_w \lVert Aw - d Vert \quad	ext{subject to}\quad w \ge 0$$

where `A[pair, split]` says whether that split separates that pair of taxa. The
non-negativity is what does the work: a split the data does not support is
driven to exactly zero rather than to a small negative number, so the result is
sparse and every surviving split has earned its place.

!!! warning "Why the fit is not optional"
    Splits read off a single tree are compatible with one another **by
    construction**, so taking them and stopping produces a drawing with no
    boxes in it however conflicted the data is. Measured on the 18-taxon 16S
    distance matrix:

    | | splits | conflicts | boxes |
    |---|---|---|---|
    | splits from the NJ tree | 33 | 0 | **0** |
    | all circular splits fitted | 40 | 20 | **20** |

    Same matrix, same ordering. The conflict was in the distances the whole
    time and only the fit can report it. The 40 fitted splits reproduce the
    distance matrix to a relative residual of 4.6%.

The fit is square in the number of taxon *pairs*, so it grows as the fourth
power of the taxon count: 0.02 s at 30 taxa, 0.8 s at 60, 4.5 s at 80. Past
`FIT_MAX_TAXA` (80) `from_distances` reads the tree instead and **warns**,
because a picture with no boxes is a claim about the data and should not be
made on the quiet. At that size use `from_trees` with a bootstrap sample
instead. Which route ran is recorded in `net.estimated`.

## Reading one

- a **box** is two groupings the data supports at once; its size is how much
  support the losing one has
- a **tree-like** region is agreement
- **terminal edges** carry each taxon out to the rim, so the names have room;
  they are in every tree by definition and say nothing about topology
- taxa that no retained split separates share a node and are fanned around it
  — being together at one node *is* the statement

`net.conflicts()` lists the conflicting split pairs, so a box you think you see
can be confirmed rather than assumed. `net.dropped` lists any split the
circular ordering could not lay out as one arc, so nothing goes missing
silently.

## `max_splits` is the readability knob

It caps the *informative* splits. Terminal splits are exempt — they can never
conflict, so they can never be the point of the picture. Each conflict opens a
box, and past a couple of dozen boxes the drawing is a mesh. Measured on a
60-replicate 16S bootstrap set (47 distinct splits, 7 of which no circular
ordering can draw):

| `max_splits` | conflicts | boxes | nodes | time |
|---|---|---|---|---|
| 12 | 0 | 0 | 31 | 0.1 s |
| 16 | 2 | 2 | 37 | 0.2 s |
| **20** (default) | 9 | 9 | 48 | 0.3 s |
| 28 (all) | 13 | 13 | 54 | 0.6 s |

Boxes and conflicts match exactly at every setting. That is the check worth
knowing about: the drawing neither invents a box nor swallows a conflict.

!!! note "Why the selection is not simply 'the strongest N'"
    Splits present in more than half the trees are exactly the majority-rule
    consensus — and a majority consensus is compatible **by construction**. A
    plain top-N cut therefore draws a tree no matter how reticulate the data
    is. So most of the budget goes to the strongest splits and the rest is
    reserved for the strongest splits that actually conflict with something
    already kept. If the data really is tree-like, nothing conflicts, the
    reserve goes unused, and the drawing is a tree because the data says so
    rather than because the selection said so.

## Why it does not cross itself

Three steps, and the middle one is what makes the picture readable.

1. **Splits**, from neighbour joining or from your tree set.
2. **A circular ordering** of the taxa, so that as many splits as possible cut
   the circle as a single arc — which makes each of them a *chord*. This is the
   idea NeighborNet turns on. Without it the split directions are arbitrary and
   the drawing crosses itself: 48 edges and 87 crossings on the 16S set above.
3. **The chord arrangement.** The chords cut the disc into cells; the network
   is the dual — one node per cell, one edge per shared chord segment. Two
   chords that cross leave a cell on each of their four sides, and those four
   cells are the box. Each edge is drawn perpendicular to its own chord, and
   the result is planar.

Step 3 is the one that has to be got right. Taking the median closure of the
taxon signatures instead — the Buneman graph — overshoots: for three mutually
conflicting splits it returns the whole 3-cube, eight nodes, which in two
dimensions can only be drawn as a wireframe cube with its hidden edges crossing
the visible ones. The chord arrangement returns seven cells, which is the
hexagon of three rhombi SplitsTree draws.

## When a split simply cannot be drawn

Four taxa admit three ways of splitting two against two, and a circle can show
only two of them — the third pair would have to sit opposite each other. So a
dataset supporting all three resolutions of some quartet is **not circular
under any ordering**, and no amount of searching will fix it.

That is not hypothetical. On the 60-replicate 16S bootstrap set, 392 of the
3060 quartets carry all three resolutions. Seven splits therefore have to go;
they hold 0.8% of the total weight, and in every one of those quartets the
losing resolution appears in 1–3 replicates out of 60 against 45 for the
winner. A search from 31 independent starting orderings — swap, move and
reverse moves run to convergence — returned the same ordering every time, so
40 of 47 is the ceiling rather than a search failure.

Those splits are listed in `net.dropped` rather than drawn wrong.

## What this is and is not

The circular ordering, the weight fit, the planarity and the boxes are all the
same construction SplitsTree uses. What differs is the ordering search:
NeighborNet agglomerates on the distance matrix directly, while this nests the
compatible splits and then reorders their children. If you would rather bring
your own split system — from SplitsTree, or any other source — hand it over
directly:

```python
pt.SplitNetwork(names, [(frozenset(side), weight), ...])
```

## Labels

Past a dozen taxa the names move out to a ring with hairline leaders back to
their nodes, because a split network bunches taxa wherever the splits between
them are short and a long species name then lies across its neighbour. Force it
either way with `label_ring=True` / `False`.
