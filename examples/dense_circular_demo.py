"""A dense circular tree in the layered style journals use for large families.

The recipe behind most published trees of a big protein or gene family is the
same one every time, and it is a *layout* decision rather than a phylogenetics
one: far too many tips to label, so the branches carry only the topology and
everything else moves out to concentric rings, each with its own key.

  * thin branches, no tip labels -- at a few thousand tips a name is narrower
    than the sector it would sit in
  * a handful of named reference sequences, and nothing else labelled, so the
    reader has somewhere to start (``tip_labels(only=...)``)
  * a tile ring for what each tip *is* (taxonomy, group)
  * a composition ring for a mixture that no single value can carry -- what
    fraction of the sequences at this tip came from each domain
    (``ring(geom="stack")``)
  * a bar ring for a number, with its own axis
  * one clade called out by name, shaded behind the branches
  * a scale bar, because branch lengths in substitutions/site read as nothing
    without one

Run:  python examples/dense_circular_demo.py
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import phytreon as pt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

meta = pd.read_csv(os.path.join(DATA, "big16S_metadata.csv"))
tree = pt.build_tree(os.path.join(DATA, "big16S_aligned.fasta"),
                     aligner="none", trim_kw=dict(max_gap=0.5),
                     method="nj", dist_model="k2p", root="midpoint")
tree.join_data(meta, on="name")

# the annotation table: one row per tip, one column per ring
ann = pd.DataFrame({"name": tree.leaf_names()})
by_name = meta.set_index("name")
for col in ("phylum", "domain", "length"):
    ann[col] = [by_name.loc[n, col] for n in ann["name"]]

# a composition ring needs a mixture per tip. Real ones come from counting the
# sequences behind each collapsed clade; here they are made up, so that the
# figure shows the layer rather than claims the numbers.
rng = np.random.default_rng(7)
for part in ("Archaea", "Bacteria", "Eukaryota"):
    ann[part] = rng.random(len(ann))

# two references named out of 106 tips
refs = ["Escherichia_coli", "Methanocaldococcus_jannaschii"]
cyanos = [n for n in tree.leaf_names()
          if by_name.loc[n, "phylum"] == "Cyanobacteriota"]

fig = (pt.TreeFigure(tree, layout="circular", extent=340)
       .branches(size=0.7)
       .highlight(taxa=cyanos, fill="#f2c14e", reach=1.0)
       .clade_label("Cyanobacteriota", taxa=cyanos, size=9)
       .ring(ann, columns=["phylum"], width=0.10, gap=0.015, colnames=False)
       .ring(ann, columns=["domain"], width=0.06, gap=0.015, colnames=False)
       .ring(ann, columns=["Archaea", "Bacteria", "Eukaryota"], geom="stack",
             title="Domain of origin", width=0.16, gap=0.02)
       .ring(ann, columns=["length"], geom="bar", width=0.20, gap=0.02,
             colnames=False)
       .tip_labels(only=refs, size=9, italic="taxa")
       .scale_bar()
       .titled("16S rRNA, 106 taxa: four annotation rings"))
fig.save(os.path.join(OUT, "dense_circular.png"))
print("[ok] dense_circular.png")

# The same tree with every tip named, for comparison: this is the figure the
# rings replace, and at 106 tips it is already spending most of the page on
# names. Past a few hundred it stops being possible at all.
(pt.TreeFigure(tree, layout="circular", extent=340)
    .branches(size=0.7)
    .tip_labels(size=6, italic="taxa")
    .ring(ann, columns=["phylum"], width=0.10, colnames=False)
    .scale_bar()
    .titled("The same tree with every tip named")
    ).save(os.path.join(OUT, "dense_circular_named.png"))
print("[ok] dense_circular_named.png")
