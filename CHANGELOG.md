# Changelog

All notable changes to phytreon are documented here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Fixed
- `build_tree()` silently fell back to UPGMA for any unrecognised `method`
  value instead of raising; now validates against an explicit whitelist.
  Removed the unused, dead `model` kwarg.
- Maximum parsimony scoring assumed a fixed A/C/G/T/U alphabet, so any other
  alphabet (protein sequences, or a discrete character/trait matrix such as
  a 0/1 gene presence/absence table) silently scored every tree 0.0. States
  are now derived per-site from whatever characters actually appear, so
  parsimony works correctly for nucleotide, amino-acid, and arbitrary
  discrete character matrices alike.
- `heatmap()` only matched rows by DataFrame index, despite its docstring
  promising a `name` column would also work (matching `ring()`/`bar_track()`);
  now uses the same name-column lookup.
- The Plotly backend did not shift aligned `Path` primitives (e.g.
  `clade_label()`'s bracket bar) past the tip labels the way aligned
  Polygon/Label/Raster already were, causing interactive HTML output to
  diverge from the matplotlib rendering.
- CI lint (`pyflakes ... || true`) could never fail the build. Switched to
  `ruff` (added to the `dev` extra) with lint now enforced.

### Added
- `read_character_matrix()`: build an `Alignment` directly from a discrete
  character/trait matrix (CSV/TSV file or DataFrame; taxa as rows, one
  column per character), ready for `parsimony_tree()` /
  `build_tree(..., method="parsimony")`.

## [0.1.1] — 2026-07-01

Renamed the project **phytree → phytreon**; first release published to PyPI.
Added GitHub Actions workflows for automated PyPI publishing (trusted
publishing / OIDC) and docs deployment to GitHub Pages.

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
