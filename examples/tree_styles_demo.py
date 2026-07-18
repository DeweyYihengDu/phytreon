"""Drawing styles beyond the plain tree.

Four figure types that show up constantly in the literature but need more than
a bare phylogram:

1. collapsed clades  -- compress a big tree to the clades you want to discuss,
   each drawn as a triangle whose sides reach its nearest and farthest hidden
   leaf (iTOL's convention)
2. node interval bars -- 95% HPD divergence-time uncertainty, the standard
   annotation on a dated Bayesian tree
3. connections        -- curved links between arbitrary tips: horizontal gene
   transfer, gene sharing, co-occurrence, host-symbiont pairs
4. DensiTree          -- a whole set of trees overlaid, so topological
   uncertainty is visible instead of hidden behind one summary tree

Run from the repo root:  python examples/tree_styles_demo.py
Outputs land in examples/out/.

Data prep (run once):  python examples/data/fetch_example_data.py
"""
import os
import random
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

ALN = os.path.join(DATA, "tol_16S_aligned.fasta")
meta = pd.read_csv(os.path.join(DATA, "tol_metadata.csv"))
COMMON = dict(aligner="none", trim_kw=dict(max_gap=0.5), root="midpoint")


def fresh_tree():
    tree = pt.build_tree(ALN, method="nj", dist_model="k2p", **COMMON)
    tree.join_data(meta, on="name")
    return tree


# -- 1. collapse the phylum-level clades ---------------------------------
tree = fresh_tree()
for phylum in ("Euryarchaeota", "Pseudomonadota", "Bacillota"):
    taxa = [r["name"] for _, r in meta.iterrows() if r["phylum"] == phylum]
    if len(taxa) < 2:
        continue
    node = tree.get_mrca(taxa)
    # only collapse where the phylum really is a clade, otherwise collapsing
    # would silently swallow unrelated taxa sitting inside it
    if set(node.leaf_names()) == set(taxa):
        pt.collapse_clade(tree, node, name=f"{phylum} ({len(taxa)})")
    else:
        extra = sorted(set(node.leaf_names()) - set(taxa))
        print(f"  ! {phylum} is not monophyletic here (also contains {extra}) "
              f"-- left expanded")

(pt.TreeFigure(tree)
    .collapsed_clades(color="#8494a8")
    .tip_labels()
    .scale_bar()
    .titled("16S: phylum-level clades collapsed to triangles")
    ).save(os.path.join(OUT, "style_collapsed.png"), figsize=(9, 6))
print(f"collapsed tree: {tree.n_leaves} rows")

# -- 2. dated tree with 95% HPD node bars --------------------------------
# The bundled 16S tree is not dated, so scale it to a plausible depth and
# attach *illustrative* intervals -- this panel demonstrates the drawing, it
# is not a divergence-time estimate. Real use reads the height_95%_HPD
# annotations off a BEAST/TreeAnnotator summary tree.
dated = fresh_tree()
pt.scale_clade(dated, dated.root,
               3500.0 / max(n.depth() for n in dated.leaves()))
maxd = max(n.depth(use_lengths=True) for n in dated.leaves())
for node in dated.traverse():
    if node.is_leaf:
        continue
    age = maxd - node.depth(use_lengths=True)
    spread = 0.15 * age + 40.0
    node.data["height_95_lower"] = max(0.0, age - spread)
    node.data["height_95_upper"] = age + spread

(pt.TreeFigure(dated)
    .node_bars()
    .time_axis()
    .tip_labels()
    .titled("Illustrative node-age intervals (not a dating analysis)")
    ).save(os.path.join(OUT, "style_nodebars.png"), figsize=(10, 6))

# -- 3. connections between tips -----------------------------------------
# Stand-in pairs; in practice these come from an HGT caller, a gene-sharing
# matrix or a co-occurrence analysis.
conn_tree = fresh_tree()
rng = random.Random(7)
names = conn_tree.leaf_names()
pairs = []
while len(pairs) < 10:
    a, b = rng.choice(names), rng.choice(names)
    if a != b:
        pairs.append((a, b, round(rng.uniform(0.1, 1.0), 2)))

(pt.TreeFigure(conn_tree, layout="circular")
    .connections(pairs, color="value")
    .tip_points(color="domain", size=6)
    .tip_labels()
    .titled("Links between tips (illustrative pairs)")
    ).save(os.path.join(OUT, "style_connections.png"), figsize=(9, 9))

# -- 4. DensiTree over bootstrap replicates ------------------------------
# A genuine tree set: NJ trees from bootstrap resamples of the alignment, so
# the cloud shows real support -- dark where the replicates agree.
aln = pt.read_fasta(ALN)
seq_names = [n for n, _ in aln]
seqs = [s for _, s in aln]
ncol = len(seqs[0])
boot = []
for _ in range(60):
    cols = [rng.randrange(ncol) for _ in range(ncol)]
    resampled = [(nm, "".join(s[c] for c in cols)) for nm, s in zip(seq_names, seqs)]
    boot.append(pt.build_tree(resampled, aligner="none", trim_kw=None,
                              method="nj", dist_model="k2p", root="midpoint"))

(pt.DensiTreeFigure(boot)
    .titled(f"{len(boot)} bootstrap NJ trees overlaid")
    ).save(os.path.join(OUT, "style_densitree.png"), figsize=(8, 6))

print("wrote style_collapsed.png, style_nodebars.png, style_connections.png "
      "and style_densitree.png to examples/out/")
