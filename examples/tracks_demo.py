"""Rectangular multi-column tracks + alignment track, on the 16S data."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib
matplotlib.use("Agg")
import pandas as pd
import phytreon as pt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA, OUT = os.path.join(HERE, "data"), os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

meta = pd.read_csv(os.path.join(DATA, "tol_metadata.csv")).set_index("name")
tree = pt.build_tree(os.path.join(DATA, "tol_16S_aligned.fasta"),
                     aligner="none", method="nj", root="midpoint").ladderize()
tree.join_data(meta.reset_index(), on="name")

# 1. multi-column tracks: two categorical tiles + a bar ----------------------
(pt.TreeFigure(tree)
    .tip_points(size=6)
    .tip_labels(italic=True)
    .heatmap(meta[["domain"]], width=0.05)      # each track its own scale
    .heatmap(meta[["phylum"]], width=0.05)
    .bar_track(meta, "length", width=0.30, fill="#4c78a8")
 ).save(os.path.join(OUT, "tol_tracks.png"))
print("[ok] tol_tracks.png  (tile + tile + bar tracks)")

# 2. alignment track -------------------------------------------------------
aln = pt.Alignment.from_fasta(os.path.join(DATA, "tol_16S_aligned.fasta"))
(pt.TreeFigure(tree)
    .tip_labels(italic=True)
    .alignment(aln, width=1.4)                      # full alignment as a raster
 ).save(os.path.join(OUT, "tol_msa.png"))
print(f"[ok] tol_msa.png  (alignment {aln.ncol} cols x {aln.nseq} seqs)")
