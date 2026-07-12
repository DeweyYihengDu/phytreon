"""Single-cell CRISPR lineage-tracing tree reconstruction, on real data.

Real mouse tumor clonal sample (226 cells, 10 CRISPR integration barcodes)
from the Cassiopeia package (see examples/data/SOURCES.md) -- reconstructs
the within-clone cell-division tree from shared indel "scars" using the
irreversible (Camin-Sokal) parsimony model, and validates the result
against Cassiopeia's own published reconstruction of the same sample via
Robinson-Foulds distance (both trees cover the same 226 cells, confirmed
below, so no leaf-subsetting is needed).

This runs the fast path by default (``search=False`` -- branch/model fit on
the NJ starting tree, no NNI topology search: ~15s total). Set
``LINEAGE_DEMO_SEARCH=1`` for the NNI hill-climbing search too; on this
226-taxon dataset that took ~153s (2.5 minutes) and improved the
reconstruction from camin_sokal_score=901 (NJ start only) to 720, with the
Robinson-Foulds distance to Cassiopeia's published tree correspondingly
improving from 0.673 to 0.632 -- both numbers last measured on the machine
this demo was developed on; see the printed wall-clock time for your own run.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import time
import matplotlib
matplotlib.use("Agg")
import phytreon as pt
from phytreon.core.tree import Tree

HERE = os.path.dirname(os.path.abspath(__file__))
DATA, OUT = os.path.join(HERE, "data"), os.path.join(HERE, "out")

SEARCH = os.environ.get("LINEAGE_DEMO_SEARCH", "0") == "1"

# --- load real allele table -> character-matrix Alignment ------------------
t0 = time.time()
aln = pt.read_allele_table(os.path.join(DATA, "lineage_alleletable.txt"))
print(f"[ok] parsed {aln.nseq} cells x {aln.ncol} CRISPR sites "
      f"({time.time() - t0:.1f}s)")

# --- reconstruct the lineage tree -------------------------------------------
t0 = time.time()
tree = pt.lineage_tree(aln, search=SEARCH)
elapsed = time.time() - t0
print(f"[ok] lineage_tree(search={SEARCH}) in {elapsed:.1f}s")
print(f"     camin_sokal_score={tree.data['camin_sokal_score']:.0f}  "
      f"min_possible_score={tree.data['min_possible_score']:.0f}  "
      f"excess_origins={tree.data['excess_origins']:.0f}")

# --- validate against Cassiopeia's own published reconstruction ------------
ref = Tree.from_newick(open(os.path.join(DATA, "lineage_reference_tree.nwk")).read())
assert set(tree.leaf_names()) == set(ref.leaf_names()), \
    "reconstructed tree and reference tree should cover the same 226 cells"
rf = pt.robinson_foulds(tree, ref, normalized=True)
print(f"[ok] Robinson-Foulds distance to Cassiopeia's published tree: "
      f"{rf:.3f} (normalized, 0=identical, 1=maximally different)")
print("     Cassiopeia's own solver (greedy/ILP) is a different, specialized "
      "algorithm -- exact agreement isn't the bar; this number is reported, "
      "not gated.")

# --- visualize the reconstructed lineage tree -------------------------------
(pt.TreeFigure(tree.ladderize())
    .branches()
    .tip_labels(max_labels=40)).save(os.path.join(OUT, "lineage_tree.png"))
print("[ok] lineage_tree.png")
