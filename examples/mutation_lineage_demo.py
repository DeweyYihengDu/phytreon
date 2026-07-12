"""Single-cell somatic-mutation lineage tree reconstruction, on real data.

Real clonal-evolution cancer dataset (Hou et al. 2012, 18 genes, 58 cells --
see examples/data/SOURCES.md) -- reconstructs the clonal cell-division tree
from shared somatic mutations using the same irreversible (Camin-Sokal)
parsimony model as the CRISPR lineage-tracing demo (examples/lineage_demo.py),
generalized beyond CRISPR via read_mutation_matrix(). Also demonstrates
reconstruct_ancestral_mutations(): which gene mutation arose on which branch.

Unlike the CRISPR demo, this does not compare against an independent
published tree -- I did not find a conveniently-parseable reconstruction
for this exact dataset, so results are reported honestly on their own
terms (reconstruction cost, timing) rather than against an unearned
external benchmark. Expect a large ``excess_origins``: this 2012-era
single-cell exome data has ~45% missing calls and substantial allele-dropout
noise (exactly what SCITE's own error-rate model was built to address), so
a plain parsimony reconstruction sees real homoplasy from sequencing error,
not just true convergent mutation.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import time
import matplotlib
matplotlib.use("Agg")
import pandas as pd
import phytreon as pt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA, OUT = os.path.join(HERE, "data"), os.path.join(HERE, "out")

# --- load real genotype matrix -> character-matrix Alignment ---------------
t0 = time.time()
genotypes = pd.read_csv(os.path.join(DATA, "mutation_genotypes.csv"), index_col=0)
aln = pt.read_mutation_matrix(genotypes)
print(f"[ok] parsed {aln.nseq} cells x {aln.ncol} genes "
      f"({time.time() - t0:.2f}s)")

# --- reconstruct the clonal lineage tree ------------------------------------
t0 = time.time()
tree = pt.lineage_tree(aln, search=True)
elapsed = time.time() - t0
print(f"[ok] lineage_tree(search=True) in {elapsed:.1f}s")
print(f"     camin_sokal_score={tree.data['camin_sokal_score']:.0f}  "
      f"min_possible_score={tree.data['min_possible_score']:.0f}  "
      f"excess_origins={tree.data['excess_origins']:.0f}")

# --- trace back which mutation arose on which branch ------------------------
pt.reconstruct_ancestral_mutations(tree, aln, site_names=list(genotypes.columns))
acquisitions = [(n.name or "(clade)", n.data["mutations_acquired"])
               for n in tree.traverse() if n.data["mutations_acquired"]]
print(f"[ok] {len(acquisitions)} branches acquired a new mutation; first 10:")
for name, genes in acquisitions[:10]:
    print(f"     {name}: {', '.join(genes)}")

# --- visualize the reconstructed lineage tree -------------------------------
(pt.TreeFigure(tree.ladderize())
    .branches()
    .tip_labels(max_labels=40)).save(os.path.join(OUT, "mutation_lineage_tree.png"))
print("[ok] mutation_lineage_tree.png")
