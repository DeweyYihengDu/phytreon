# Tutorial: sequence-similarity networks

A phylogenetic tree assumes you can align the sequences well enough that
branch order means something. For a protein family whose members are only
detectable by profile searches, that assumption fails: the alignment carries
real error, and the deep branching order a tree reports is then an artefact of
that error rather than a record of descent.

The standard alternative is to drop the tree and draw the **sequence space**
itself — one node per sequence, an edge wherever a pairwise search finds
significant similarity, laid out by a force-directed algorithm so that groups
of mutually similar sequences fall into visible clusters. This is what CLANS
produces, and it is a common sight in comparative-genomics papers where a
family's internal structure is real but its deep branching is not recoverable.

```python
net = pt.SequenceNetwork.from_alignment(aln, cutoff=0.85)
net.color_by(phylum, title="phylum", baseline="other")
net.save("clusters.pdf")
```

## Reading one

- a **tight ball** of nodes is a group whose members all detect each other
- a **thin bridge** between two balls means the groups are related, but only
  distantly
- an **isolated node** found nothing above the cutoff

Distance on the page is **not** an evolutionary distance. It is the layout's
compromise between many pairwise attractions, so read clusters and
connections — never a ruler.

## Getting the edges in

Three constructors, depending on what you already have:

```python
# from an alignment: all-against-all pairwise identity
pt.SequenceNetwork.from_alignment(aln, cutoff=0.85)

# from a distance matrix you computed elsewhere (similarity = 1 - distance)
pt.SequenceNetwork.from_distances(names, matrix, cutoff=0.5)

# from a real search — e.g. rows parsed out of a BLAST tabular report
pt.SequenceNetwork.from_pairs([("q1", "s1", 0.82), ("q1", "s2", 0.44), ...])
```

!!! note "`from_alignment` needs an alignment to exist"
    It measures identity down alignment columns, so it only works where the
    sequences *were* globally alignable. For the case this method exists for —
    a family too divergent to align — run a real all-against-all search
    (BLASTP, DIAMOND, `jackhmmer`) and feed the hits to `from_pairs`.

## The cutoff is the figure

The cutoff decides how much of the sequence space you see, exactly like the
E-value slider in CLANS. Too high and everything falls apart into singletons;
too low and the whole family fuses into one ball. Neither is wrong — but a
single panel can always be tuned to tell a tidy story, and the reader cannot
see that from one picture. Report the cutoff, and show more than one if the
clustering is your claim:

```python
for cutoff in (0.78, 0.85):
    net = pt.SequenceNetwork.from_alignment(aln, cutoff=cutoff)
    net.titled(f"identity > {cutoff}").save(f"net_{cutoff}.pdf")
```

`examples/network_demo.py` does exactly this on the bundled 106-sequence 16S
set: at 0.78 there is one 90-sequence hairball, at 0.85 the phyla resolve into
separate clean clusters.

## Checking what you are looking at

Clusters that look obvious to the eye should be confirmed, not assumed:

```python
for comp in net.components():        # connected components, largest first
    print(len(comp), comp[:5])
```

## Styling

```python
net.color_by(groups, title="family",
             baseline="uncharacterised",   # grey out the background category
             order=["known", "novel", "uncharacterised"])
net.label_clusters({"Photoglobin": photoglobin_ids,
                    "Phycocyanin": phycocyanin_ids})
net.label_nodes(["the_one_you_care_about"])
```

Cluster labels are placed outside their own cluster and anchored on the side
facing it, so long names grow outward instead of back over the nodes.

## Layout knobs

`gravity` pulls everything gently toward the centre. Without it a graph with
disconnected pieces blows apart — nothing but repulsion acts between the
pieces, isolated nodes drift off, and since they then set the scale of the
picture the connected core collapses into an unreadable dot. Sequence networks
are disconnected by construction, so this matters more here than in the
textbook algorithm.

The default (0.8) was chosen by measuring the trade-off on a real 16S graph and
a synthetic three-family one: going from 0.25 to 0.8 grows the main cluster
from 16% to 27% of the frame while cluster separation only falls from 4.0x to
3.6x. Raise it if isolated sequences are pushing your clusters into a corner;
lower it if distinct families are being squeezed together.

```python
pt.SequenceNetwork.from_pairs(pairs, seed=3, iterations=500)
```

The layout is randomised but **deterministic for a given `seed`** — the same
seed always gives the same picture, and a different seed gives a different
(equally valid) arrangement of the same graph.
