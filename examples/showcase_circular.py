"""Showcase: a richly annotated circular tree.

Demonstrates the elements needed for pathogen-phylogeny figures: lineage-
coloured branches, shaped tip points, a categorical ring, an outer bar ring,
tip labels, and several non-overlapping legends.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib
matplotlib.use("Agg")
import random
import pandas as pd
import phytreon as pt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)

tr = pt.datasets.random_tree(60, seed=3).ladderize()

# --- define ~6 lineages by expanding the shallowest splits, colour branches ---
roots = [tr.root]
while len(roots) < 6:
    internal = [x for x in roots if not x.is_leaf]
    node = min(internal, key=lambda x: x.depth(use_lengths=False))
    roots.remove(node)
    roots.extend(node.children)
pt.group_clade(tr, {r: f"Lineage {i+1}" for i, r in enumerate(roots)},
               key="lineage", default="other")

# --- per-tip metadata ---
rng = random.Random(7)
tips = tr.leaf_names()
meta = pd.DataFrame({
    "name": tips,
    "city": [rng.choice(["HCMC", "Hue", "KH"]) for _ in tips],
    "serotype": [rng.choice(["Ogawa", "Inaba", "Hikojima"]) for _ in tips],
    "AMR_score": [round(rng.uniform(0, 1), 2) for _ in tips],
}).set_index("name")
tr.join_data(meta.reset_index(), on="name")     # so tip points can map city/serotype

p = (pt.TreeFigure(tr, layout="circular_slanted", extent=350)
     .branches(color="lineage", size=1.1)
     .tip_points(shape="city", size=5, color="#333333")
     .ring(meta, columns=["serotype"], geom="tile", width=0.07, offset=0.03)
     .ring(meta, columns=["AMR_score"], geom="bar",
           width=0.28, offset=0.03)
     .tip_labels(size=4)
     .titled("Annotated circular tree (lineage / city / serotype / AMR)"))
p.save(os.path.join(OUT, "showcase_circular.png"))
print("[ok] showcase_circular.png")
