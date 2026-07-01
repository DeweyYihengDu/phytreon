# Changelog

All notable changes to phytreon are documented here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/).

## [0.1.0] — 2026-06-24

First working version: a pure-Python phylogenetics + tree-visualisation
library covering inference, comparative methods, and publication figures.

### Core
- `Tree`/`Node` data model; Newick/Nexus/PhyloXML I/O (Biopython bridge +
  native Newick parser); metadata join (`Tree.join_data`).

### Layouts
- rectangular, slanted, dendrogram, circular, fan, radial, inward-circular,
  unrooted (equal-angle), equal-daylight. Backend-agnostic scene graph.

### Inference
- Distance: NJ / UPGMA with JC69/K2P-corrected distances (negative branches
  clamped); built-in progressive MSA aligner; configurable alignment trimming.
- Maximum likelihood (native, pure Python): JC69/K80/HKY85/GTR, discrete-Γ
  rate heterogeneity, NNI search, AIC/BIC + `model_finder`; external IQ-TREE/
  FastTree adapters.
- Maximum parsimony (Fitch + NNI). Bipartition bootstrap for NJ/UPGMA/ML/MP.

### Comparative methods
- Ancestral states: Fitch parsimony, Mk marginal ML, Brownian (continuous).
- Stochastic character mapping (`stochastic_map`) with painted branches.
- Tree ops: rotate/flip/ladderize/collapse/scale_clade/cut_tree/midpoint_root/
  group_clade/group_otu; Robinson-Foulds distance.

### Plotting
- `TreeFigure` fluent builder with matplotlib (static) and plotly (interactive)
  backends. Elements: branches, tip_labels, node_labels, support_labels, points
  (colour/size/shape mapping), highlight, clade_label, heatmap, ring (tile/bar),
  bar_track, alignment, painted_branches, time_axis (with geological scale),
  node_pies.
- HCL hue-wheel palette; continuous colorbars; non-overlapping label/track/
  legend placement.

### Validation
- `validation/validate.py` (pure Python): likelihood engine matches an
  independent naive implementation to machine precision; NJ recovers a tree
  from its own additive distances (RF = 0). 37 pytest tests.
