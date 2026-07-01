# phytreon benchmark (timings)

| operation | size | time (s) |
|---|---|---|
| NJ (JC69) | 12 taxa x 1626 bp | 0.25 |
| Parsimony (Fitch+NNI) | 12 x 1626 | 0.05 |
| ML JC69 (branch opt) | 12 x 1626 | 21.39 |
| ML HKY85+G (branch opt) | 12 x 1626 | 179.56 |
| ML HKY85 (full NNI) | 12 x 1626 | 221.24 |
| circular plot | 50 tips | 2.33 |
| circular plot | 200 tips | 1.96 |
| circular plot | 500 tips | 4.55 |
| builtin MSA | 8 seqs x 300 bp | 0.63 |

## Guidance
- **Validated core (recommended for production):** NJ/UPGMA, native ML
  +G on small/medium alignments, rectangular/circular plotting, heatmap,
  rings. Correctness checks in `validation/` (pure Python).
- Pure-Python **ML/MSA scale to ~tens of taxa**; for hundreds+ use
  `aligner='mafft'` and `ml_engine='iqtree'`.
- Plotting handles hundreds of tips; use `tip_labels(max_labels=...)`
  to thin labels on large trees.
