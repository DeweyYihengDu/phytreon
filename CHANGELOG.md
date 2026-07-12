# Changelog

All notable changes to phytreon are documented here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Fixed
- The native Newick writer/parser (`to_newick()`/`parse_newick()`, used by
  `Tree.write()`/`Tree.from_newick()` whenever no file path is given) never
  quoted or unquoted taxon names containing reserved Newick punctuation
  (`()[]{}/\,;:=*'` or whitespace). A name like `"weird(name),here"` wrote
  out unquoted and silently split into three unrelated leaves on
  read-back -- no error, just a corrupted tree. Now quotes such names on
  write (doubling any embedded `'`) and correctly parses quoted labels,
  including embedded reserved characters, back out again.
- `Tree.ladderize()` recomputed each node's subtree size from scratch inside
  its own sort comparator, so every level of nesting re-triggered a full
  recursive re-descent through everything beneath it -- exponential blowup
  on deep/unbalanced trees (harmless at the small scale previously
  exercised; non-terminating after 50+ minutes on a real 226-taxon,
  depth-37 tree). Now computes every node's size once and sorts from that.
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
- `Tree.get_mrca()` silently computed the MRCA of whichever requested taxa
  *were* found, so a typo'd or missing name (e.g. `get_mrca(["Human",
  "NotExisting"])`) returned a misleadingly small clade instead of an error.
  Now raises `ValueError` listing the missing taxa by default (`strict=True`);
  pass `strict=False` for the old lenient behaviour.
- K2P distance treated `A<->U` as a transition; it is a purine<->pyrimidine
  transversion (only `A<->G` and `C<->T`/`C<->U` are transitions), so K2P
  distances on RNA data were biased.
- The Plotly backend rendered every point as a circle regardless of a
  `shape=` mapping -- the marker→symbol table existed but was only wired up
  for legend swatches, not the actual marker trace.
- `ace_parsimony()` computed the Fitch state set for internal nodes as the
  intersection of *all* children at once, undercounting steps at polytomies
  (3+ children, as produced by `collapse_low_support()`): 3 children with 3
  disjoint states scored 1 instead of the correct 2. Now combines children
  sequentially.
- `robinson_foulds(normalized=True)` divided by `2*(n-3)`, which is `-2` (not
  0) for `n=2`, so `or 1` didn't guard it and small trees could return a
  negative distance. Now returns `0.0` for `n<4`.
- `Alignment(names, seqs)` accepted mismatched-length sequences, a
  names/seqs count mismatch, or duplicate names without error; each now
  raises `ValueError` in `__post_init__`.
- `stochastic_map()` recorded only `params[0]` as `"rate"`, discarding the
  other fitted rates for `SYM`/`ARD` models (which have more than one).  Now
  records `"model"` and the full `"rates"` list.
- `_resolve_size()` used `isinstance(v, (int, float))`, which is `False` for
  `numpy.int64`/`float64` (the dtype pandas normally produces), so a `size=`
  mapping from a DataFrame column silently fell back to a constant size.
  Now shares the `numbers.Number`-based check already used for colour scales.

### Added
- `read_character_matrix()`: build an `Alignment` directly from a discrete
  character/trait matrix (CSV/TSV file or DataFrame; taxa as rows, one
  column per character), ready for `parsimony_tree()` /
  `build_tree(..., method="parsimony")`.
- Protein (amino acid) support for the native ML engine, purely additive
  alongside the existing nucleotide-only code path: pass `ml_model="JTT"`
  / `"WAG"` / `"LG"` to `ml_tree()`/`build_tree(..., method="ml")` for
  empirical protein substitution models (each with its own published
  equilibrium frequencies), and `model_finder()` now ranks JTT/WAG/LG
  instead of the nucleotide model set when it detects protein data. A new
  alphabet-mismatch guard raises `ValueError` if a nucleotide model is
  used on protein data or vice versa, rather than silently miscoding
  amino acid letters that happen to coincide with nucleotide codes.
  `distance_matrix_model()` gains an explicit opt-in `dist_model="poisson"`
  (the protein analogue of the Jukes-Cantor correction). None of this
  changes any existing nucleotide default: `ml_model` still defaults to
  `"HKY85"` and `dist_model` still defaults to `"jc69"` (which still falls
  back to raw p-distance on non-nucleotide data unless you opt in to
  `"poisson"`).
- Single-cell CRISPR lineage-tracing tree reconstruction
  (`phytreon/infer/lineage.py`), purely additive alongside the existing
  reversible Fitch parsimony: `read_allele_table()` turns a Cassiopeia-style
  allele table into an `Alignment` (reusing `read_character_matrix()`,
  handling allele dropout and near-saturated sites correctly);
  `sankoff_score()`/`camin_sokal_score()` add a general Sankoff parsimony
  engine plus the irreversible preset appropriate for CRISPR scars (a
  derived state can arise independently more than once, but never reverts
  or converts directly to a different derived state); `lineage_tree()` is
  the NNI hill-climbing search, also reachable via
  `build_tree(..., method="parsimony", parsimony_model="camin_sokal")`.
  Validated against a real published dataset in `examples/lineage_demo.py`
  (Robinson-Foulds distance to Cassiopeia's own reconstruction of the same
  226-cell sample, reported honestly rather than gated against a threshold).
- `prune_to_taxa()` (`phytreon/treeops.py`): restrict a tree to a leaf
  subset, collapsing now-redundant unary nodes.
- `read_mutation_matrix()` (`phytreon/infer/lineage.py`): generalizes
  lineage-tracing reconstruction beyond CRISPR allele tables to any
  single-gene or multi-gene somatic-mutation/genotype matrix -- the same
  irreversible-mutation model applies (a mutated gene doesn't spontaneously
  revert to wild-type), so it feeds directly into the existing
  `camin_sokal_score()`/`lineage_tree()` with no changes needed there. The
  "phantom ancestral row" correctness fix from `read_allele_table()`
  (guarantees wild-type always codes as `"0"`, even at a gene mutated in
  100% of profiled cells) is now a shared helper both readers use.
- `expression_dendrogram()`/`expression_distance_matrix()`
  (`phytreon/infer/expression.py`): hierarchical-clustering dendrograms of
  transcriptional similarity for one gene or a small combination of genes,
  reusing the existing alphabet-agnostic `neighbor_joining()`/`upgma()`
  with a new distance metric for continuous expression data
  (`scipy.spatial.distance.pdist`, no new dependency). **Explicitly not
  phylogenetic** -- expression similarity reflects cell state, not
  ancestry -- so it's named/documented distinctly from `lineage_tree()`,
  and the result carries `tree.data["tree_type"] =
  "expression_similarity_dendrogram"` as a machine-readable flag.

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
