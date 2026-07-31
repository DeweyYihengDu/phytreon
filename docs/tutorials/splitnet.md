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

# from a distance matrix (splits via neighbour joining, weight = branch length)
pt.SplitNetwork.from_distances(names, matrix)

# from an alignment
pt.SplitNetwork.from_alignment(aln)
```

`from_trees` is the one that turns a bootstrap set straight into a picture of
*which groupings the replicates disagree about* — a split found in every tree
draws long, one found in half of them draws short and boxed against its rival.

## Reading one

- a **box** is two groupings the data supports at once; its size is how much
  support the losing one has
- a **tree-like** region is agreement
- **pendant edges** fan taxa that no retained split separates around their
  shared vertex — being together at one vertex *is* the statement

`net.conflicts()` lists the conflicting split pairs, so a box you think you see
can be confirmed rather than assumed.

## `max_splits` is the readability knob

This matters more than it looks. A split network stays legible only while the
conflict is modest; past that it degenerates into a mesh, and the median
closure that draws it grows steeply. Measured on a 60-replicate 16S bootstrap
set (31 distinct splits, 82 conflicting pairs in total):

| `max_splits` | conflicts | boxes | vertices | time |
|---|---|---|---|---|
| 16 | 1 | 1 | 18 | instant |
| **20** (default) | 11 | 15 | 34 | 0.1 s |
| 24 | 31 | 82 | 80 | 1.1 s |
| 31 (all) | 82 | 474 | 284 | 80 s |

Raise it to chase weaker conflicting signal and expect both the picture and the
wait to degrade; lower it for a cleaner figure showing only the conflicts among
the strongest splits.

!!! note "Why the selection is not simply 'the strongest N'"
    Splits present in more than half the trees are exactly the majority-rule
    consensus — and a majority consensus is compatible **by construction**. A
    plain top-N cut therefore draws a tree no matter how reticulate the data
    is; on the 16S set above it discarded all 82 conflicting pairs. So most of
    the budget goes to the strongest splits and the rest is reserved for the
    strongest splits that actually conflict with something already kept. If the
    data really is tree-like, nothing conflicts, the reserve goes unused, and
    the drawing is a tree because the data says so rather than because the
    selection said so.

## What this is and is not

Splits are extracted by neighbour joining (or from your tree set) and drawn by
the split-decomposition convention: each split becomes a displacement shared by
every taxon on one side, so conflicting splits open into boxes. The vertex set
is closed under coordinate-wise medians, which supplies the box corners.

This is **not** a reimplementation of SplitsTree's NeighborNet. It recovers the
boxes for the conflicts an NJ tree plus a tree set expose; NeighborNet's
circular ordering finds splits this will not. For a publication figure where
the split system itself is the claim, compute it in SplitsTree and bring the
splits here to draw.
