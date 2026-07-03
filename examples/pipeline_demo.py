"""Raw sequences -> tree, the whole configurable pipeline."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib
matplotlib.use("Agg")
import phytreon as pt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)

# --- a bunch of UNALIGNED sequences (two groups, varying lengths/indels) ---
seqs = [
    ("A1", "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"),
    ("A2", "ATGGCCATTGTTATGGGCCGCTGAAAGGGTGCCCGATAG"),
    ("A3", "ATGGCCATTGTAATGGGCCGCTGTAAGGGTGCCGATAG"),     # small deletion
    ("B1", "ATGTCGATTCTAATGAACCGCTGAAAGCGTGACCTTTAG"),
    ("B2", "ATGTCGATTCTAATGAACCGCTGTAAGCGTGACCTTTAG"),
    ("B3", "ATGTCGATTCTAATGAACCGGCTGAAAGCGTGACCTTTAG"),    # small insertion
]

# one call, every stage configured: builtin MSA -> trim -> NJ -> 200 boots
tree, aln = pt.build_tree(
    seqs,
    aligner="builtin",
    align_kw=dict(seqtype="nucleotide", match=2, mismatch=-1, gap=-3),
    trim_kw=dict(max_gap=0.4, min_occupancy=0.5),
    method="nj",
    root="midpoint",          # sensible root for the unrooted NJ tree
    bootstrap=200,
    seed=1,
    return_alignment=True,
)
print("alignment:", aln, "->", aln.ncol, "cols after trim")
print("tree:", tree.write())
print("supports:", [round(n.support, 0) for n in tree.traverse()
                    if not n.is_leaf and n.support is not None])

# --- tree manipulation: ladderize, rotate a clade, cut into clusters ---
pt.ladderize(tree)
clusters = pt.cut_tree(tree, k=2)
print("2 clusters:", clusters)

for tip in tree.leaves():
    tip.data["cluster"] = f"clade{clusters[tip.name]}"

(pt.TreeFigure(tree)
    .tip_points(color="cluster", size=8)
    .tip_labels()
    .support_labels()).save(os.path.join(OUT, "pipeline_tree.png"))
print("[ok] wrote pipeline_tree.png")

# --- show rotate/flip actually move branches ---
order_before = tree.leaf_names()
pt.rotate(tree, tree.root)
print("leaf order before rotate:", order_before)
print("leaf order after  rotate:", tree.leaf_names())
