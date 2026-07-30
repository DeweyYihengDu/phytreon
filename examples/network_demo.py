"""Sequence-similarity network -- when a tree would not be honest.

A phylogenetic tree assumes you can align the sequences well enough that
branch order means something. For a protein family whose members are only
detectable by profile searches, that assumption fails: the alignment carries
real error, and the deep branching order a tree reports is then an artefact of
that error rather than a record of descent.

The standard alternative is to draw the sequence space itself -- one node per
sequence, an edge wherever a pairwise search finds significant similarity,
laid out by a force-directed algorithm. Groups of mutually similar sequences
fall into visible clusters, and the *absence* of a connection is as
informative as its presence. This is what CLANS produces.

Run from the repo root:  python examples/network_demo.py
Outputs land in examples/out/.

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

# Real 16S sequences: alignable, so a tree is perfectly defensible here. The
# point of the panel is the method, and using real data keeps the clusters
# honest -- they come out of the sequences, not out of a random generator.
aln_path = os.path.join(DATA, "big16S_aligned.fasta")
if not os.path.exists(aln_path):
    aln_path = os.path.join(DATA, "tol_16S_aligned.fasta")
    meta_path = os.path.join(DATA, "tol_metadata.csv")
else:
    meta_path = os.path.join(DATA, "big16S_metadata.csv")

pairs = pt.read_fasta(aln_path)
aln = pt.Alignment(names=[n for n, _ in pairs], seqs=[s for _, s in pairs])
meta = pd.read_csv(meta_path).set_index("name")

# The cutoff is the knob that decides how much of the sequence space you see:
# too high and everything falls apart into singletons, too low and the whole
# family fuses into one ball. Reporting it is part of reporting the figure.
CUTOFF = 0.78
net = pt.SequenceNetwork.from_alignment(aln, cutoff=CUTOFF, seed=4)

phylum = {nm: meta.loc[nm, "phylum"] for nm in net.names if nm in meta.index}
counts = phylum.value_counts() if hasattr(phylum, "value_counts") else None

# grey out everything except the few phyla worth following, so the eye lands
# on those rather than on 25 competing hues
focus = ["Pseudomonadota", "Bacillota", "Actinomycetota", "Euryarchaeota"]
shown = {nm: (p if p in focus else "other") for nm, p in phylum.items()}

net.color_by(shown, title="phylum", baseline="other",
             order=[*focus, "other"])
net.titled(f"16S sequence space, identity > {CUTOFF} "
           f"({len(net.names)} sequences)")
net.save(os.path.join(OUT, "network_16S.png"), figsize=(9, 7))

comps = net.components()
print(f"sequences        : {len(net.names)}")
print(f"edges (>{CUTOFF}) : {len(net.edges)}")
print(f"components       : {len(comps)}  sizes {[len(c) for c in comps][:8]}")
print("largest component is dominated by:",
      pd.Series([phylum.get(n, "?") for n in comps[0]]).value_counts().head(3).to_dict())

# The same data at a stricter cutoff: fewer edges, the family fragments. Worth
# showing side by side, because a single cutoff can always be made to tell a
# tidy story and the reader cannot see that from one panel.
strict = pt.SequenceNetwork.from_alignment(aln, cutoff=0.85, seed=4)
strict.color_by(shown, title="phylum", baseline="other", order=[*focus, "other"])
strict.titled("Same sequences, identity > 0.85 (stricter cutoff)")
strict.save(os.path.join(OUT, "network_16S_strict.png"), figsize=(9, 7))
print(f"at cutoff 0.85   : {len(strict.edges)} edges, "
      f"{len(strict.components())} components")

print("wrote network_16S.png and network_16S_strict.png to examples/out/")
