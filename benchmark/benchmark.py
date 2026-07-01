"""Timing + scaling benchmark for the phytreon validated core.

Times the core operations on the bundled 16S data and on random trees of
increasing size, and writes benchmark/benchmark_report.md.  The "validated
core" (NJ, ML+G, rectangular/circular plotting, heatmap) is what we recommend
for production; pure-Python ML/MSA do not scale to hundreds of taxa.
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib
matplotlib.use("Agg")
import phytreon as pt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "examples", "data")
REPORT = os.path.join(HERE, "benchmark_report.md")


def timed(fn):
    t = time.perf_counter()
    out = fn()
    return out, time.perf_counter() - t


def main():
    rows = []
    aln = pt.Alignment.from_fasta(os.path.join(DATA, "tol_16S_aligned.fasta"))
    n, ncol = aln.nseq, aln.ncol

    _, t = timed(lambda: pt.build_tree(aln, aligner="none", method="nj"))
    rows.append(("NJ (JC69)", f"{n} taxa x {ncol} bp", t))
    _, t = timed(lambda: pt.parsimony_tree(aln, search=True))
    rows.append(("Parsimony (Fitch+NNI)", f"{n} x {ncol}", t))
    _, t = timed(lambda: pt.ml_tree(aln, model="JC69", search=False))
    rows.append(("ML JC69 (branch opt)", f"{n} x {ncol}", t))
    _, t = timed(lambda: pt.ml_tree(aln, model="HKY85", gamma=4, search=False))
    rows.append(("ML HKY85+G (branch opt)", f"{n} x {ncol}", t))
    _, t = timed(lambda: pt.ml_tree(aln, model="HKY85", search=True))
    rows.append(("ML HKY85 (full NNI)", f"{n} x {ncol}", t))

    # plotting scaling on random trees
    for nt in (50, 200, 500):
        tr = pt.datasets.random_tree(nt, seed=1)
        _, t = timed(lambda tr=tr: (pt.TreeFigure(tr, layout="circular")
                                    .tip_points(size=2)
                                    ).save(os.path.join(HERE, "_tmp.png")))
        rows.append(("circular plot", f"{nt} tips", t))

    # builtin aligner scaling (small n, growing length)
    import random
    rng = random.Random(0)
    base = "".join(rng.choice("ACGT") for _ in range(300))
    seqs = [(f"s{i}", "".join(c if rng.random() > 0.1 else rng.choice("ACGT")
                              for c in base)) for i in range(8)]
    _, t = timed(lambda: pt.align(seqs, seqtype="nucleotide"))
    rows.append(("builtin MSA", "8 seqs x 300 bp", t))

    if os.path.exists(os.path.join(HERE, "_tmp.png")):
        os.remove(os.path.join(HERE, "_tmp.png"))

    lines = ["# phytreon benchmark (timings)", "",
             "| operation | size | time (s) |", "|---|---|---|"]
    for name, size, t in rows:
        lines.append(f"| {name} | {size} | {t:.2f} |")
    lines += [
        "", "## Guidance",
        "- **Validated core (recommended for production):** NJ/UPGMA, native ML",
        "  +G on small/medium alignments, rectangular/circular plotting, heatmap,",
        "  rings. Correctness checks in `validation/` (pure Python).",
        "- Pure-Python **ML/MSA scale to ~tens of taxa**; for hundreds+ use",
        "  `aligner='mafft'` and `ml_engine='iqtree'`.",
        "- Plotting handles hundreds of tips; use `tip_labels(max_labels=...)`",
        "  to thin labels on large trees.",
    ]
    with open(REPORT, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
