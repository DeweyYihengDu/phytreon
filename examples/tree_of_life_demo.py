"""Real-data example: a microbial "tree of life" from 16S rRNA (NCBI).

Builds a tree from the bundled 16S sequences (common model organisms across
the bacterial phyla plus four archaea) and renders it as a rectangular tree
and as a circular tree with metadata rings.

Data prep (run once):  python examples/data/fetch_example_data.py
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib
matplotlib.use("Agg")
import pandas as pd
import phytreon as pt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

meta = pd.read_csv(os.path.join(DATA, "tol_metadata.csv")).set_index("name")

# build from the pre-aligned fasta (fast); use the raw fasta to run MSA too
tree = pt.build_tree(
    os.path.join(DATA, "tol_16S_aligned.fasta"),
    aligner="none",                       # already aligned
    trim_kw=dict(max_gap=0.5),            # cut very gappy columns
    method="nj",
    root="midpoint",
    bootstrap=200,
    seed=1,
)
for tip in tree.leaves():
    for col in ("domain", "phylum"):
        tip.data[col] = meta.loc[tip.name, col]

print("tree:", tree.write()[:90], "...")
archaea = tree.get_mrca([n for n in tree.leaf_names()
                         if meta.loc[n, "domain"] == "Archaea"])
mono = set(archaea.leaf_names()) == set(meta[meta.domain == "Archaea"].index)
print("Archaea MRCA spans", len(archaea.leaf_names()), "tips "
      f"({'monophyletic' if mono else 'paraphyletic here'})")

# 1. rectangular -----------------------------------------------------------
(pt.TreeFigure(tree.ladderize())
    .tip_points(color="domain", size=8)
    .tip_labels(italic=True)
    .support_labels(size=7)).save(os.path.join(OUT, "tol_rect.png"))
print("[ok] tol_rect.png")

# 2. circular with metadata rings -----------------------------------------
# Use compact "G. species" tip labels on the circular figure so the long
# binomials never reach the legend column (rename tips + metadata index
# together so the rings still match each tip).
abbrev = {}
for tip in tree.leaves():
    parts = tip.name.split("_")
    abbrev[tip.name] = f"{parts[0][0]}. {parts[1]}" if len(parts) > 1 else tip.name
for tip in tree.leaves():
    tip.name = abbrev[tip.name]
meta_c = meta.rename(index=abbrev)

(pt.TreeFigure(tree, layout="circular", extent=345)
    .tip_points(size=4)
    .ring(meta_c, columns=["domain", "phylum", "length"], width=0.11, gap=0.022,
          colnames=False)                       # legend already names each ring
    .tip_labels(size=9, italic=True)).save(os.path.join(OUT, "tol_rings.png"),
                                           figsize=(12.5, 9))
print("[ok] tol_rings.png")
