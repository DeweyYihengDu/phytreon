"""Native (pure-Python) maximum-likelihood tree from the 16S data.

This evaluates the model and optimises branch lengths on the NJ starting tree
(``ml_search=False``) for a quick demo. Set ``ml_search=True`` to also run the
NNI topology search (much slower on this many taxa). For rate heterogeneity add
``ml_gamma=4``; to pick a model by AIC use ``pt.model_finder(...)``; for large
data use ``ml_engine="iqtree"``.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib; matplotlib.use("Agg")
import pandas as pd
import phytreon as pt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA, OUT = os.path.join(HERE, "data"), os.path.join(HERE, "out")

meta = pd.read_csv(os.path.join(DATA, "tol_metadata.csv")).set_index("name")

# pure-Python ML: HKY85 on the NJ tree (set ml_search=True for NNI search)
tree = pt.build_tree(os.path.join(DATA, "tol_16S_aligned.fasta"),
                     aligner="none", method="ml", ml_model="HKY85",
                     ml_search=False)
print("ML logL = %.1f   model = %s   kappa = %.3f"
      % (tree.data["logL"], tree.data["model"], tree.data["model_params"][0]))

tree = pt.midpoint_root(tree)
tree.join_data(meta.reset_index(), on="name")
(pt.TreeFigure(tree.ladderize())
    .tip_points(color="domain", size=8)
    .tip_labels(italic=True)).save(os.path.join(OUT, "tol_ml.png"))
print("[ok] tol_ml.png")
