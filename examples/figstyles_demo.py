"""Figure styles borrowed from the comparative-genomics literature.

Five layouts that recur in papers where a tree alone would not carry the
argument:

1. ribbon tanglegram  -- two facing trees joined by filled bands, so a group
   that the two trees place differently shows as a twist rather than as a
   bundle of crossing lines
2. multi-panel grid   -- many small trees under one colour key, to compare
   several gene families at a glance
3. domain track       -- each protein's architecture beside its tip
4. stacked support    -- several support values per branch, as printed
5. split network      -- conflicting splits drawn as boxes

Run from the repo root:  python examples/figstyles_demo.py
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
rng = random.Random(11)


def build(method="nj", **kw):
    tree = pt.build_tree(ALN, method=method, **dict(COMMON, **kw))
    tree.join_data(meta, on="name")
    return tree


# -- 1. ribbon tanglegram -------------------------------------------------
# NJ vs UPGMA again, but with the taxa grouped by domain: the question is not
# "which tip moved" but "did the clock assumption move a whole group".
nj = build("nj", dist_model="k2p")
upgma = build("upgma", dist_model="k2p")
fig = pt.TangleFigure(nj, upgma, titles=("neighbour joining", "UPGMA (clock)"))
fig.untangle()
fig.ribbons("phylum", title="phylum")
fig.titled("Ribbons: which groups the clock assumption moves")
fig.save(os.path.join(OUT, "style_ribbons.png"))
print("ribbons: crossings after untangling =", fig.crossings())

# -- 2. multi-panel grid --------------------------------------------------
# Stand-in "gene families": bootstrap replicates of the same alignment, which
# is honest about what they are -- resampled versions of one dataset, not
# independent genes -- while still showing the layout.
seqs = pt.read_fasta(ALN)
names = [n for n, _ in seqs]
cols = [s for _, s in seqs]
ncol = len(cols[0])
panels = []
for i in range(8):
    idx = [rng.randrange(ncol) for _ in range(ncol)]
    tree = pt.build_tree([(nm, "".join(s[c] for c in idx))
                          for nm, s in zip(names, cols)],
                         aligner="none", trim_kw=None, method="nj",
                         dist_model="k2p", root="midpoint")
    tree.join_data(meta, on="name")
    panels.append(pt.TreeFigure(tree, layout="unrooted")
                  .tip_points(color="domain", size=5)
                  .titled(f"replicate {i + 1}"))
(pt.panels(panels, ncols=4, share_legend=True, label_panels=True)
    .titled("Eight bootstrap replicates under one key")
    .save(os.path.join(OUT, "style_panels.png")))

# -- 3. domain architecture beside the tips -------------------------------
# Illustrative architectures, clearly labelled as such: the point of the panel
# is the drawing, and inventing plausible-looking real annotations would be
# worse than obviously schematic ones.
small = pt.Tree.from_newick(
    "(((NsnA_1:.1,NsnA_2:.1):.1,NsnA_3:.15):.1,(NsnA_4:.2,NsnA_5:.2):.1);")
arch = {
    "NsnA_1": [("wHTH", 60), ("ParB", 180), ("DUF262", 210), ("ParBDB", 90)],
    "NsnA_2": [("wHTH", 55), ("ParB", 175), ("DUF262", 205), ("ParBDB", 88)],
    "NsnA_3": [("B3-like", 80), ("ParB", 180), ("DUF262", 210), ("ParBDB", 90)],
    "NsnA_4": [("ParB", 180), ("DUF262", 210), ("ParBDB", 90)],
    "NsnA_5": [("PUA-like", 95), ("ParB", 180), ("DUF262", 210),
               ("ParBDB", 90), ("TRD", 70)],
}
(pt.TreeFigure(small).tip_labels().domains(arch, width=0.9, labels=True)
    .titled("Schematic domain architectures (illustrative)")
    .save(os.path.join(OUT, "style_domains.png"), figsize=(11, 4)))

operons = {
    "NsnA_1": [("nsnA", 300), ("nsnB", 250), ("nsnC", -280)],
    "NsnA_2": [("nsnA", 300), ("nsnB", 250), ("nsnC", -280)],
    "NsnA_3": [("tnp", -200), ("nsnA", 300), ("nsnB", 250), ("nsnC", -280)],
    "NsnA_4": [("nsnA", 300), ("nsnC", -280)],
    "NsnA_5": [("nsnA", 300), ("nsnB", 250), ("nsnC", -280), ("hyp", 150)],
}
(pt.TreeFigure(small).tip_labels()
    .domains(operons, width=0.9, arrows=True, labels=True)
    .titled("Schematic gene neighbourhoods (arrows show strand)")
    .save(os.path.join(OUT, "style_operons.png"), figsize=(11, 4)))

# -- 4. several support values per branch ---------------------------------
boot = pt.build_tree(ALN, method="nj", dist_model="k2p", bootstrap=100,
                     **COMMON)
for node in boot.traverse():
    if node.is_leaf or node.support is None:
        continue
    # a second and third "method" would come from FastTree / BEAST runs; here
    # they are derived from the bootstrap so the panel stays self-contained
    node.data["b"] = round(node.support)
    node.data["s"] = round(min(100, node.support * 1.02))
    node.data["p"] = round(min(1.0, node.support / 100 * 1.01), 2)
(pt.TreeFigure(boot).tip_labels()
    .support_labels(attr=["b", "s", "p"], stack=True,
                    prefixes=["b", "s", "p"], size=5.5)
    .titled("Three support values per branch")
    .save(os.path.join(OUT, "style_support.png"), figsize=(9, 6)))

# -- 5. split network -----------------------------------------------------
replicates = []
for _ in range(60):
    idx = [rng.randrange(ncol) for _ in range(ncol)]
    replicates.append(pt.build_tree(
        [(nm, "".join(s[c] for c in idx)) for nm, s in zip(names, cols)],
        aligner="none", trim_kw=None, method="nj", dist_model="k2p",
        root="midpoint"))
net = pt.SplitNetwork.from_trees(replicates, label_size=7)
net.color_by({r["name"]: r["domain"] for _, r in meta.iterrows()},
             title="domain")
net.titled("60 bootstrap replicates as a split network")
net.save(os.path.join(OUT, "style_splitnet.png"), figsize=(11, 8))
verts, edges = net._network()
print(f"split network: {len(net.splits)} splits, "
      f"{len(net.conflicts())} conflicting pairs, "
      f"{len(edges) - len(verts) + 1} boxes")

# -- 6. the same data as a NeighborNet ------------------------------------
# Straight from the distance matrix, with no tree set at all. The point of the
# comparison is what estimate= does: splits read off the NJ tree are compatible
# with one another by construction, so that route cannot draw a box however
# conflicted the data is, while fitting every circular split to the distances
# finds the conflict that was in the matrix the whole time.
aln = pt.Alignment(names, cols)
taxa, dmat = pt.infer.distance_matrix_model(aln, "k2p")
for label, estimate in (("from the NJ tree", False), ("all splits fitted", True)):
    nn = pt.SplitNetwork.from_distances(taxa, dmat, estimate=estimate,
                                        label_size=7)
    v, e = nn._network()
    print(f"neighbornet ({label}): {len(nn.splits)} splits, "
          f"{len(nn.conflicts())} conflicting pairs, "
          f"{len(e) - len(v) + 1} boxes")
nn.color_by({r["name"]: r["domain"] for _, r in meta.iterrows()}, title="domain")
nn.titled("16S distances: every circular split fitted (NeighborNet)")
nn.save(os.path.join(OUT, "style_neighbornet.png"), figsize=(11, 8))

print("wrote style_ribbons / style_panels / style_domains / style_operons / "
      "style_support / style_splitnet / style_neighbornet to examples/out/")
