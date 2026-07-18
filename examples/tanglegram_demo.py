"""Tanglegram: where do two trees of the same taxa disagree?

Builds trees from the *same* 16S alignment by two different methods and faces
them off, so every taxon they place differently shows up as a crossing link.

The headline figure uses the large bundled set -- 106 taxa across 25
prokaryotic phyla -- because the comparison only really becomes worth drawing
once a tree has enough taxa to disagree about:

* neighbour joining vs UPGMA -- UPGMA assumes a molecular clock and NJ does
  not, so they genuinely conflict; the crossings that survive untangling are
  the taxa the clock assumption moves.
* neighbour joining vs parsimony (on the small 18-taxon set) -- these agree on
  the tip *order* while still differing in two splits, the cautionary case:
  zero crossings does not mean two trees are identical.

The same recipe compares trees from two different *datasets* rather than two
methods: build one tree from genomic data and one from transcriptomic data
over the same taxa, then hand both to ``TangleFigure``.

Run from the repo root:  python examples/tanglegram_demo.py
Outputs land in examples/out/.

Data prep (run once):  python examples/data/fetch_example_data.py
                       python examples/data/fetch_large_16S.py
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

COMMON = dict(aligner="none", trim_kw=dict(max_gap=0.5), root="midpoint")


def load(aln, meta_csv, methods):
    meta = pd.read_csv(os.path.join(DATA, meta_csv))
    trees = []
    for method, kw in methods:
        tree = pt.build_tree(os.path.join(DATA, aln), method=method,
                             **dict(COMMON, **kw))
        tree.join_data(meta, on="name")
        trees.append(tree)
    return trees


# -- 1. the big one: 106 taxa, NJ vs UPGMA -------------------------------
nj, upgma = load("big16S_aligned.fasta", "big16S_metadata.csv",
                 [("nj", dict(dist_model="k2p")),
                  ("upgma", dict(dist_model="k2p"))])

# RF counts conflicting splits and does not care how the trees are drawn;
# crossings count links that tangle, which *does* depend on clade rotation --
# so untangle first, then read the number.
rf = pt.robinson_foulds(nj, upgma, normalized=True)
before = pt.crossing_number(nj, upgma)

fig = pt.TangleFigure(nj, upgma, titles=("neighbour joining", "UPGMA (clock)"))
after = fig.untangle().crossings()
fig.connect(highlight_discordant=True)
# colour by domain, not phylum: 25 phyla is past what any categorical palette
# can keep distinguishable, and the two domains are the readable signal here
fig.left.tip_points(color="domain", size=5)
fig.right.tip_points(color="domain", size=5)
fig.titled(f"16S, {nj.n_leaves} taxa: neighbour joining vs UPGMA (RF={rf:.2f})")
for ext in ("pdf", "png"):
    fig.save(os.path.join(OUT, f"tanglegram.{ext}"))

print(f"NJ vs UPGMA  ({nj.n_leaves} taxa)")
print(f"  normalised RF                  : {rf:.3f}")
print(f"  crossings before/after untangle: {before} -> {after}")
print("  the links still crossing are the taxa the clock assumption moves")

# links coloured by phylum instead, to ask a different question: is the
# conflict confined to one group, or does it cut across the whole tree?
(pt.TangleFigure(nj, upgma, titles=("neighbour joining", "UPGMA (clock)"))
    .untangle()
    .connect(color="phylum")
    .titled("16S: method comparison, links coloured by phylum")
    .save(os.path.join(OUT, "tanglegram_by_phylum.png")))

# -- 2. the cautionary case: agree on order, still not the same tree -----
nj_s, mp_s = load("tol_16S_aligned.fasta", "tol_metadata.csv",
                  [("nj", dict(dist_model="k2p")), ("parsimony", {})])
rf_mp = pt.robinson_foulds(nj_s, mp_s, normalized=True)
mp_fig = pt.TangleFigure(nj_s, mp_s,
                         titles=("neighbour joining", "parsimony"))
after_mp = mp_fig.untangle().crossings()
mp_fig.connect(highlight_discordant=True)
mp_fig.titled(f"16S: neighbour joining vs parsimony (RF={rf_mp:.2f})")
mp_fig.save(os.path.join(OUT, "tanglegram_parsimony.png"))

print(f"\nNJ vs parsimony  ({nj_s.n_leaves} taxa)")
print(f"  normalised RF                  : {rf_mp:.3f}")
print(f"  crossings after untangle       : {after_mp}")
print("  no crossings, yet RF > 0: the two trees admit the same tip order")
print("  while still resolving some splits differently -- read RF too.")

print("\nwrote tanglegram.pdf/.png, tanglegram_by_phylum.png and "
      "tanglegram_parsimony.png to examples/out/")
