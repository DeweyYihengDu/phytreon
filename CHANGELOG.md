# Changelog

All notable changes to phytreon are documented here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Added
- **Community phylogenetics: `patristic_distances`, `mpd`, `mntd`, `ses_mpd`,
  `ses_mntd`, `beta_mntd`, `beta_nti`, plus `permanova` and `mantel`.** The
  existing diversity functions answer "how much evolutionary history is in this
  sample" (Faith's PD) and "how different are two samples' branches" (UniFrac).
  These answer the question after that one: *given* that samples differ, are the
  taxa that co-occur more closely related than a random draw from the same pool
  would be, and is the turnover between samples more phylogenetic than chance can
  explain -- which is what separates habitat filtering and selection from drift.

  `patristic_distances` was the missing primitive underneath all of it and is
  worth having on its own: there was no way to get tip-to-tip distances *through
  the tree* at all, only `distance_matrix`, which compares sequences without
  reference to any tree. It is independent of rooting, checked against a
  hand-computed four-tip matrix and, separately, against `phylo_vcv`'s identity
  `d_ij = V_ii + V_jj - 2 V_ij` on a random tree, which is a different code path
  to the same numbers.

  `permanova` (Anderson 2001) and `mantel` are the two tests that actually get
  run on a `unifrac_matrix` once samples carry labels or environmental data --
  "do these groups differ in composition" and "does dissimilarity track
  environment" -- and both get their p-values by permutation, since distances
  over shared samples are not independent observations. Both were held to the
  same Type-I error standard as `pgls` was: on a homogeneous cloud with random
  group labels, `permanova` rejects 3.2-5.5% of true nulls at a nominal 5% across
  four configurations of sample size and group count (600 replicates each), i.e.
  correctly calibrated and slightly conservative. An initial reading of 8.0%
  looked like inflation and was noise from 200 replicates -- worth measuring
  properly rather than either ignoring it or "fixing" a problem that was not
  there.

  **Sign conventions get their own tests, because inverting one flips the
  biological conclusion while leaving every magnitude plausible.** NRI and NTI
  carry a factor of -1 (Webb 2000; they are the negations of the standardised
  effect sizes, Kembel 2009), so positive is clustered. betaNTI does not (Stegen
  et al. 2012, *ISME J* 6:1653-1664), so positive is more turnover than expected.
  Both directions were confirmed against the primary sources before implementing
  rather than from memory, and are pinned by tests on communities whose answer is
  known by construction: a real clade of the tree comes out at NRI +5.2 / NTI
  +3.3, two halves of one clade give a negative betaNTI and two different clades
  +5.8.

  The null model is a tip relabelling (Stegen's randomisation, picante's
  `taxa.labels`), which keeps community richness and tree shape exactly as
  observed. Implemented as a permutation of *positions* rather than of names,
  which is the same operation for every metric here and made one trap explicit:
  a pair's two samples must share one relabelling, or taxa present in both stop
  being present in both and the null inflates. Verified directly -- with the
  shared permutation, two identical samples give a betaMNTD of exactly 0 in all
  200 null draws, so betaNTI is `NaN` and that is the correct answer; permuting
  independently instead gives the null a spread (sd 0.070 about a mean of 0.302)
  and would return a confident-looking -4.3.

  Beyond the closed-form checks, the standardised indices are validated by the
  property any correctly specified effect size must have: on random samples they
  centre on zero **with unit standard deviation** (measured NRI -0.05, NTI +0.08,
  sd 1.01 and 1.06), which tests the null's spread and not merely its centre. A
  first attempt to demonstrate NTI used a tree regular enough that patristic
  distances took only four distinct values, inflating the null's sd and capping
  the index at +1.8; the metric was right and the test design was wrong.
- **`tests/test_reference_crosscheck.py`: the comparative methods checked against
  something other than themselves.** Every other test of these functions checks
  a statistic against the property it is *defined* to have -- K averaging to 1
  under Brownian motion, PGLS not inflating Type I error -- which is strong but
  circular in one respect: the same implementation produces both the number and
  the behaviour being used to judge it. Two outside sources now, both pure
  Python, since a Python package that needed R to verify itself would be a
  contradiction:
  - **`statsmodels.GLS`.** PGLS at a fixed `lambda` *is* generalised least
    squares with a known error covariance, so a mature independent
    implementation of that exact computation can be held against ours term by
    term. Agreement is 1.4e-14 relative at worst, on coefficients, standard
    errors, t-values, p-values and residual degrees of freedom, across
    `lambda` of 0.0/0.3/0.5/1.0 and one to three predictors.
  - **`sympy` exact rationals**, on a four-tip tree small enough that the
    answers can also be read off it by hand. `phylo_vcv` reproduces the matrix
    written out by inspection; Blomberg's K comes to exactly 70014/93775, which
    the implementation matches to 1.5e-16; Faith's PD and UniFrac match hand
    arithmetic (PD of every taxon = the tree's 10 total branch length, UniFrac
    of two disjoint halves = 1, `{A}` against `{A,B}` = 3/5).

  Both are `dev`-extra only and the tests `importorskip`, so a bare install
  skips them rather than failing -- verified by running the file with both
  blocked, where the six reference tests skip and the two hand-arithmetic ones
  still run.

### Changed
- **`unifrac_matrix` is 6-15x faster.** Both paths now build each sample's
  per-branch data into a row of a `samples x edges` array once, so every pair is
  a numpy expression instead of a Python loop over a dict; the weighted
  denominator also separates (`sum(L * (F_i + F_j))` is just `r_i + r_j`), so it
  costs one matrix-vector product for all pairs rather than a pass per pair.
  A 967-sample table on a 500-tip tree went from about 5 minutes to under 10
  seconds. Verified numerically identical to the untouched pairwise functions
  across all three modes on random trees and tables (worst disagreement
  4.4e-16, i.e. floating-point rounding), which is now a test rather than a
  one-off check, since a rewrite for speed is exactly the kind that can change
  answers quietly. Explicitly **not** a change of complexity: the cost is still
  proportional to `samples^2 x edges` -- measured at 10-13x for a 10x bigger
  tree, before and after -- because comparing two samples branch by branch has
  to look at every branch. A first version of the accompanying test asserted the
  scaling had flattened; it had not, and the test was wrong rather than the code.

### Fixed
- **`unifrac_matrix(weighted=True)` silently returned wrong distances when the
  table had a column the tree does not have as a tip.** The realistic case, not
  an exotic one: tree building drops sequences (too short, failed alignment,
  chimeras filtered after the table was counted), so an ASV table normally has
  *more* ASVs than the tree has tips. The extra taxa's abundance still went into
  each sample's total, so every real taxon's subtree fraction came out too small
  -- 0.846 where the answer was 0.5 on a four-tip test case, with the fractions
  summing to 0.18 instead of 1.0 -- and those distances then feed an ordination
  with nothing having gone visibly wrong. The same mistake produced three
  different behaviours across the module: this silent one, a bare `KeyError`
  from `unifrac_matrix(weighted=False)`, and a clear `ValueError` from every
  single-pair function. Now one `ValueError` from all of them, naming the
  offending columns and the one-liner that subsets the table. `faiths_pd_table`
  checked all columns up front too, since it had been validating only the
  nonzero entries of each row -- so whether the same table raised depended on
  the data in it.
- **`ace_parsimony`, `ace_ml` and `stochastic_map` silently ignored a trait name
  that is not a tip of the tree.** A misspelling, or a tree and a metadata table
  that disagree on how names are formatted -- and not harmless when ignored: the
  tip the entry was meant for is left unlabelled, so it gets reconstructed as
  fully ambiguous and resolved from its parent instead. One misspelling among the
  seven bundled primates moved the reconstructed state of seven nodes with
  nothing reported. And since all three read their state space off
  `trait.values()`, a misspelled key can add a state that no tip in the tree
  carries, which `ace_ml` will then estimate rates over. Now a `ValueError`
  naming the offending keys, from these three and from `ace_continuous`, which
  did fail loudly but blamed the tip left without a value rather than the
  misspelled key that stranded it. The *other* direction is deliberately still
  allowed: a tip with no entry in `trait` is documented behaviour for the
  discrete methods (treated as fully ambiguous -- partial data is a normal thing
  to reconstruct from), and only `ace_continuous` requires every tip, a weighted
  average having nothing to average over an unknown.
- **`blomberg_k`, `pagels_lambda` and `pgls` failed three different ways on a
  tree whose covariance matrix is singular.** Also not exotic: any two tips at
  zero patristic distance have identical rows in `phylo_vcv`'s matrix, and a
  zero-length terminal branch does exactly that -- IQ-TREE and RAxML emit them
  routinely. `blomberg_k` raised a bare `LinAlgError: Singular matrix` with
  nothing to act on; `pgls` continued, having landed on a `lambda` where the
  matrix happened to invert; and `pagels_lambda` was the bad one -- its
  optimiser escaped to `lambda ~ 0`, where the off-diagonals vanish and the
  matrix inverts again, and reported `p = 1.0`, which reads as "no phylogenetic
  signal" for what is really "not computable on this tree". All three now raise
  one `ValueError` naming the tied tips and what to do about them. The guard is
  deliberately only on the paths that estimate `lambda`: a `lambda_` handed over
  as a number is the caller asserting the error structure, and `lambda_=0.0`
  ignores the tree entirely, so it is well defined on exactly the trees this
  check rejects and still works.

### Changed
- **`pgls` estimates `lambda` by REML by default, and charges it a degree of
  freedom.** A Type-I error sweep over tree shape (ultrametric and not), tree
  size (10-80 taxa) and the true lambda that generated the traits -- two traits
  simulated independently on one tree, so every rejection is a false positive
  against a nominal 5% -- found the old `lambda_="ML"` default rejecting 7.9%
  of true nulls at 10-20 taxa (6400 replicates; the same quantity reads 8.1% on
  a separate 4800-replicate set, and every figure below is quoted from whichever
  run its comparison arm was paired against on identical data). Two causes,
  both now addressed. Plain ML pulls a
  variance component towards zero when the mean structure is estimated from the
  same data, and here that means a systematically too-small `lambda` (measured
  at 10 tips on data whose real lambda was 1.0: ML averaged 0.60, REML 0.78) --
  a too-small lambda understates how dependent the tips are and so understates
  the standard errors. And `lambda` is a parameter read out of the data like
  any other, so the t-test now spends a degree of freedom on it rather than
  treating it as given. Together: 7.9% -> 6.8% pooled over the small-tree
  cells, and 7.3-9.5% -> 5.2-7.2% in the cells where the true lambda really was
  1.0. The rejections removed are strictly the right ones -- paired against ML
  on identical data, "REML called it significant and ML did not" happened 0
  times in 1200 replicates, while the reverse happened 33 times. `lambda_="ML"`
  is still available for reproducing software that does it that way
  (`caper::pgls`). Two smaller things fell out of the same sweep: `pgls` now
  reports `lambda_method` so a result says how its lambda was obtained, and it
  refuses a fit with no residual degrees of freedom left instead of returning
  meaningless p-values from one.
- **Two README examples were not valid Python.** A bare `...` inside a dict
  literal (`{"Human": 1.4, "Chimp": 1.35, ...}`), so copy-pasting either gave a
  `SyntaxError`. The ellipsis moved into a comment, where it still says "one per
  tip" without breaking the code. Found by a check that walks every `pt.*` call
  in the README and verifies the attribute exists and the keyword arguments are
  real parameters -- worth noting because the two unparseable blocks were being
  *skipped* by that check, and one of them was the block documenting every
  function added this release. With them parsing, coverage went from 22 calls to
  39, and all 39 resolve. That check is now `tests/test_readme_api.py` rather
  than something run once: the README is the first thing most people read and
  nothing else verified it, so a renamed parameter could stay documented
  indefinitely with the reader finding out via `TypeError`. It cannot check that
  the examples *run* -- they use placeholder filenames and variables on purpose
  -- but "this function accepts that argument" is the part that rots silently,
  and one of its four tests exists only to fail if the block extraction itself
  ever stops matching, so the others cannot pass by checking nothing.
- **`pagels_lambda`'s `logLik` is now an actual log-likelihood.** It was the
  concentrated criterion with the additive constants dropped for optimisation:
  correct to a constant, which cancels inside the reported likelihood ratio but
  not in anything else built on it, so an AIC computed from it was wrong.
  Checked against the multivariate normal density at the fitted parameters,
  evaluated independently by scipy.

### Added
- **`pgls(..., n_boot=)`, a parametric bootstrap p-value.** For the residual
  the two changes above do not reach: the t-test treats the estimated `lambda`
  as if it were known exactly, and at 10-20 taxa it is being read out of very
  few points, so when it comes out too low the standard errors come out too
  small. The bootstrap simulates from the fitted reduced model (that predictor
  dropped, the others kept and refitted) and re-estimates `lambda` on every
  replicate, which prices that uncertainty in rather than conditioning on one
  point estimate of it. Paired against the t-test on identical data over 3200
  replicates: 7.3% -> 5.5%, which is the difference between six standard errors
  above the nominal 5% and one and a half -- inflated versus not distinguishable
  from correct. Worth the `n_boot` extra fits per predictor below roughly 20
  taxa and not otherwise.
  Also measured and **rejected** along the way: a likelihood-ratio test on the
  predictor instead of the Wald t-test, which sounds like the more principled
  option and is markedly worse calibrated here. The first comparison of the two
  was confounded -- Wald with a REML lambda against an LRT with an ML one, since
  an LRT on fixed effects has to use ML (restricted likelihoods for models with
  different fixed effects are not comparable) -- so it mixed the test statistic
  and the lambda estimator into one 3.4-point gap. Re-run with all three arms on
  the same 6400 replicates and the lambda estimator held fixed, it decomposes,
  and the test statistic turns out to be the larger of the two effects:
  Wald+ML 7.9%, LRT+ML 10.1% (+2.3 points from the statistic alone),
  Wald+REML 6.8% (-1.1 from the estimator). The LRT is worst at the smallest
  trees (11.3-11.4% at 10 taxa against 8.4-10.2% at 20), which is the mechanism
  showing through: it leans on an asymptotic chi2_1 null, where the Wald t has
  a finite-sample reference distribution to lean on instead. What survives all of this is a real
  limit, not a bug to keep chasing: below ~20 taxa a PGLS p-value is mildly
  anti-conservative whatever test you use, and one of 0.04 from 10 species
  should not be leaned on. Two other p-values in the module were audited on the
  same terms and left alone: `pagels_lambda`'s likelihood-ratio test rejects
  0.5-2.0% of true nulls, conservative because lambda=0 sits on the boundary of
  [0,1] and the LR's true null there is a 50:50 chi2_0/chi2_1 mixture rather
  than the plain chi2_1 it is tested against -- the safe direction, and what
  phytools' `phylosig` reports too; `blomberg_k`'s permutation test came out at
  3.0-5.6%, i.e. calibrated, as a permutation test should be by construction.
- **Phylogenetic diversity: `faiths_pd`/`faiths_pd_table`,
  `unweighted_unifrac`/`weighted_unifrac`/`unifrac_matrix`.** Not "what did
  the ancestors look like" (the rest of `comparative`), but "how much
  evolutionary history does this sample's -- or these two samples' -- taxa
  collectively represent": the standard alpha/beta diversity metrics for a
  community sitting on a tree (a 16S ASV table, say), and the natural next
  step once that tree is built rather than an edge case of a general
  phylogenetics toolbox. Faith's PD is a single tree walk (the total branch
  length of the smallest subtree connecting a sample's taxa to the root);
  UniFrac's weighted form needs one post-order pass per sample for the
  per-branch abundance fractions, reused across every pair `unifrac_matrix`
  computes rather than redone per pair -- the difference between
  `O(samples)` and `O(samples^2)` tree walks once a table has hundreds of
  rows. Verified against the properties each metric is defined to have,
  not just "runs without raising": PD of every taxon in the tree equals the
  tree's own total branch length exactly; UniFrac between a sample and
  itself is 0 and between a clean bipartition of all taxa is exactly 1,
  weighted or not.
- **Phylogenetic signal and PGLS: `phylo_vcv`, `blomberg_k`, `pagels_lambda`,
  `pgls`.** All three read a continuous trait's covariance across tips
  against the Brownian-motion expectation implied by the tree itself
  (`phylo_vcv`, ape's `vcv.phylo`) -- Blomberg et al. (2003)'s K asks how
  the trait's *actual* fit compares to that expectation (K=1 matches BM);
  Pagel (1999)'s lambda finds, by maximum likelihood, how much of the
  tree's shared history the trait's covariance actually supports (0 = none,
  1 = all of it), more robust than K to polytomies and uncertain branch
  lengths (Freckleton, Harvey & Pagel 2002); PGLS regresses one trait on
  another using that same covariance as the error structure, rather than
  the plain least-squares assumption that every tip is an independent data
  point, which is false by construction for any two related tips and
  inflates false positives the more of the tree they share (Felsenstein
  1985). None of the three formulas has a simple closed-form answer to
  assert against for an arbitrary tree, so each is verified by the
  statistical property that *does* have a known answer: K averages to 1.0
  over many Brownian-motion simulations on the same tree it is being asked
  about (measured: 0.996-1.034 across two tree sizes and 300+ simulations
  each); lambda's maximum-likelihood estimate recovers the true value used
  to simulate the data on average (mean 0.98 for lambda=1 data, 0.03 for
  lambda=0 data, at 60 tips -- lambda is a genuinely bimodal estimator at
  very small sample sizes, a documented property of the statistic itself,
  not a bug, and not the regime these numbers were measured in); and PGLS
  reproduces the Felsenstein (1985) demonstration that founded this whole
  field of methods -- two independently-evolved traits sharing one **60-tip**
  tree, 400 replicates, gave a 36-37% false-positive rate under plain OLS
  against 5.0-6.5% under PGLS on the same data, at a nominal 5%. An earlier
  draft of this entry said "41% and 4%" and gave no tree size; both numbers
  came from a 150-replicate run and were within a standard error of these,
  but the omission was the real problem. OLS's inflation comes from shared
  ancestry, so it scales with the tree (16% at 10 tips, 43% at 80) rather
  than being one constant of the method, and an unqualified 41% reads as
  exactly that. See the `Changed` entry above for the full sweep.
- **Bootstrap now respects `constraint=`.** Every replicate is rebuilt under
  the same constraint as the main tree (`constrained_nj` for `method="nj"`,
  the same `-g`/`--tree-constraint` for an external `ml_engine`) instead of
  without it, which used to score how well a *forced* split does against data
  that was never made to honour it -- not a fair reading of the tree actually
  reported. A low number now means what it should: the resampled data
  disagree with the grouping even once the constraint is applied, not that an
  unconstrained rebuild would have found a different grouping on its own
  (which the constraint had already decided regardless of the data).
  Uncovered while wiring this up: `build_tree(method="ml", ml_engine="iqtree",
  bootstrap=N)` built the tree successfully and silently produced **no**
  support values at all -- `bootstrap` only ever reached the generic
  resampling loop, which is skipped for every external ML engine (one
  subprocess search per replicate is a non-starter), and was never forwarded
  to IQ-TREE's own `-bb` either. Fixed by forwarding it at build time instead,
  the same as the constraint file.
- **`infer_raxmlng()`, RAxML-NG as a second constrained-ML engine**
  (`build_tree(method="ml", ml_engine="raxml-ng", ...)`) alongside IQ-TREE --
  same `constraint=` support (`--tree-constraint`, needs at least 4 taxa, a
  RAxML-NG restriction not relaxed here) and the same own-bootstrap-not-a-
  rebuild-loop treatment (`--all`/`--bs-trees` in one run, reading
  `<prefix>.raxml.support` instead of `.raxml.bestTree` when it ran). Unlike
  IQ-TREE's `-m MFP`, RAxML-NG has no built-in "find the model for me", so
  `model=` needs a real one (`"GTR+G"` default for nucleotides).
- **`infer_rapidnj()`, fast approximate NJ for alignments too large for the
  builtin engine** (`build_tree(method="nj", nj_engine="rapidnj")`).
  `neighbor_joining()` wraps Biopython's textbook O(n^3) implementation,
  impractical somewhere in the low thousands of tips -- exactly the range a
  large 16S ASV table reaches. RapidNJ computes its own distance matrix from
  the alignment (a smaller model menu than `dist_model`'s `k2p`/`poisson`/...,
  the tradeoff for not needing phytreon's own matrix built in Python first)
  and is incompatible with `constraint=` for the same reason. Bootstrap goes
  through phytreon's own rebuild loop here rather than RapidNJ's own `-b` --
  unlike an external ML search, a fast approximate NJ rebuild is cheap enough
  that the loop is not the bottleneck, and RapidNJ's own flag does not
  document how it reports support well enough to trust guessing at it.
- **`highlight(..., label=...)` prints a clade's name centred inside its own
  band**, instead of (or alongside) a separate `clade_label()` outside the
  tree -- the way a dense published tree names a filled sector directly, with
  no room to spare on naming every tip. A single `node=`/`taxa=` highlight
  takes the text itself; `by=` takes `label=True` and uses each group's own
  value, one label per band including every run of a scattered group. Ink is
  chosen for contrast against the band as actually rendered -- its pale wash
  (alpha 0.3), not the saturated fill colour underneath it -- which
  `readable_on()` could not do before this: it read the fill's raw hex, so a
  colour dark enough to want white ink at full strength but pale enough to
  want black once washed out the way a highlight always is would have picked
  wrong. Only draws for `shape="wedge"` on a circular layout -- a ring is thin
  and curved, and a label centred on it would need to run tangentially, which
  nothing here typesets, so asking for one there warns and leaves it out.
- **`clade_label(..., leader=...)` parks the name further out and draws a
  thin dotted line back**, for a clade whose own bracket is too narrow --
  angularly, or too few rows -- to hold its name beside it without crowding a
  neighbour's, the way a dense published tree calls a small clade out with a
  name on a leader rather than jammed against its bracket. A number, the same
  convention as `reach`/`gap`/`offset` elsewhere: a fraction of the tree's own
  depth. Getting the line to land correctly needed `_push_radius` (the
  mechanism that moves a circular tree's clade bracket/label out past
  wherever the tip names turned out to end, since only the renderer knows
  that) generalised from *the path's outermost point, moved as a whole* to
  *every point moved by the same amount*: a bracket's own points already
  share one radius so the two gave the same answer there, but a leader line's
  two ends do not, and the old version would have collapsed them onto a
  single radius -- erasing the very gap the line exists to show.
- **`outgroup_root(tree, taxa)`**, and `build_tree(..., root=taxa)` to reach it
  in one call. `midpoint_root` assumes every lineage drifts at about the same
  rate, which is often not even approximately true; an outgroup only assumes
  it split off before everything else did, usually a claim already settled
  independently of the gene tree being rooted. Finding the right branch is not
  simply "the MRCA of `taxa`": an unrooted method's own display root is an
  arbitrary implementation choice -- NJ roots wherever its last join happened
  to land -- and can seat some of `taxa` as direct children of that root while
  the rest sit nested deep in what looks like an unrelated subtree, which
  hides a real split from a plain ancestor search even though the two are
  still one clade in the *unrooted* tree. Checking every node's own leaf set
  against both `taxa` and its complement finds the branch regardless of how
  the tree happened to arrive rooted (a real case, not a hypothetical: the
  bundled 6-taxon NJ demo tree does exactly this to one of its two groups).
  Raises rather than silently rooting somewhere plausible-looking if `taxa`
  does not correspond to any single branch at all.
- **`build_tree(..., constraint={tip: label})`** forces a taxonomy grouping --
  a genus column, say -- to come out monophyletic in the tree itself, at two
  different strengths depending on `method`, rather than only rearranging a
  tree already built (`sort_by`, above, is display-only and cannot do this: it
  can bring two clades a tree already keeps apart no closer, it can only
  rotate what is already a fork). A tip left out of the mapping, or mapped to
  `None`, is placed exactly as if there were no constraint for it.
  - `method="ml"` with `ml_engine` set to an external engine turns the mapping
    into a **constrained search**: `constraint_tree()` builds one polytomy per
    label -- unresolved internally, so nothing about *how* those tips relate
    is dictated -- and IQ-TREE's `-g` (RAxML-NG's `--tree-constraint` is the
    same idea, not wrapped here) still picks the topology within and between
    groups by likelihood; only "no other tip may land inside this group" is
    fixed. Pass a path to an existing constraint file, or a `Tree`, instead of
    a mapping to skip that conversion.
  - `method="nj"` instead runs **`constrained_nj`**, two plain NJ passes with
    no search or external engine involved: once inside each group on only
    that group's own tips, then again on a reduced matrix (each pair of
    groups averaged over the real distances between their tips) to place the
    groups -- and any free tip -- against each other, grafting each group's
    own midpoint-rooted subtree onto its slot on that backbone. This *forces*
    monophyly outright: there is no way for the result to disagree, even
    where NJ run without a constraint would put two of a group's tips on
    opposite sides of a split (verified against a real case in the bundled
    16S data: *Bacillus subtilis* and *B. cereus* do not come out sister
    under plain NJ -- four other genera's tips fall between them -- and do
    under `constrained_nj`, every other leaf unmoved).
  - Reach for the constrained search when the sequence data should still have
    the final word on whether a genus really is monophyletic and the
    constraint only needs to settle the ties it leaves open; reach for
    `constrained_nj` when the taxonomy should win outright, such as a display
    tree organised by genus. Native ML and the distance/parsimony methods
    other than NJ have no constrained search to hook into and raise rather
    than silently ignore `constraint`.
- **`sort_by(tree, key)`** reorders sibling branches so tips sharing a category
  -- a genus column on a big 16S ASV tree, say -- sit together, without moving
  a single branch or adding/dropping a split: it only ever changes which side
  of a fork a subtree is drawn on, the same move `ladderize` already makes
  with subtree size as the key. That bounds what it can do honestly: rotating
  a fork brings two clades that are already siblings closer together, it
  cannot join clades the tree keeps apart. A genus that a single 16S
  hypervariable region does not resolve as one clade still comes out in more
  than one run after this -- correctly, since forcing it into one run would
  require moving a branch, which would tell the reader the opposite of what
  the gene supports (`highlight(by=key)` shows those remaining runs as
  separate bands rather than hiding the split). Searches for the arrangement
  minimizing the number of label changes tip-to-tip, the way `untangle`
  already searches rotations to minimize crossings against a reference tree,
  and keeps a rotation only when it measurably helps -- so the result is never
  worse than what the tree started with. The first cut at this sorted every
  node once by its own subtree's majority category and skipped the search;
  measured on a real 106-taxon 16S tree it made the phylum grouping *worse*
  (33 runs against a plain ladderize's 31), because a summary taken in
  isolation at one node knows nothing about what its neighbours need.
- **`highlight(by="column")`** shades every group in a joined column behind the
  branches — each in its own colour, with a key — instead of one
  `highlight(node=...)` call per clade and no way to tell the shades apart.
  Works on rectangular and circular layouts — a band behind the branches, and
  on a circular tree a **ring** around it. Not sectors filled from the centre:
  a sector's area grows with the square of its radius, so filling from the
  middle swamps the drawing and leaves the tree the smallest thing in it, which
  is why iTOL draws its coloured ranges as a band too. `shape="wedge"` asks for
  the filled sector anyway; `width`/`offset` size the ring, and it claims a slot
  from the same running cursor the `ring()` tracks use, so the two stack outward
  instead of colliding. The ring is drawn solid (0.85) rather than at the band's
  pale 0.3 — nothing is printed over it, and at 0.3 it read as a faded copy of
  its own legend swatch.
  A group whose taxa are **not** monophyletic comes out as several bands, one
  per run of adjacent tips, rather than one band over their common ancestor:
  that ancestor reaches down over other groups' taxa, so a single band there
  would colour in tips that do not belong to the group — tidier, and false.
  The shades stay pale enough to read a name on (measured: the darkest leaves
  black text at 12.5:1).
  - `gap` is the white space between neighbouring bands, in tip rows. The
    default `0` butts them together into one continuous field of colour, which
    reads as a single stratified panel rather than a set of stripes; raise it if
    two adjacent groups drew colours close enough that the boundary stops being
    obvious.
  - `reach` decides how wide the colour block is. A **number** is a multiple of
    the tree's own depth — `0.5` reaches half way in from the root, `1.0` to
    the tips, `1.3` a third past them. The reference is the branch-length axis
    the tree is drawn on, so the same number means the same width in any figure
    at any font size. `"labels"` runs out past the longest tip label instead, so
    no species name hangs off the end of the block it belongs to; that end is
    found while rendering, since how much room a name takes depends on the font
    and the figure rather than on the tree. Left unset it picks per layout:
    `"labels"` on a rows layout, `1.0` on a round one — out past the names a
    wedge's area grows with the square of its radius, so carrying the colour
    there floods the whole disc and leaves the tree a smudge in the middle.
    Rendered side by side, the version stopping at the tips is the readable one,
    with the names outside on white.
  - `anchor` says which end holds still as the width changes. `"root"` (the
    default) pins the inner edge where the clade starts and moves the outer
    one, so the block grows outward from the tree. `"tips"` pins the outer edge
    out past the names and retreats the inner one, so the species names stay on
    the colour at every width and what shortens is the end nearest the root.
  - `span` decides where the bands start. `"aligned"` (the default) lines every
    left edge up at the shallowest clade's own start: flush, and no band
    reaches further rootward than the drawing already did somewhere, so the
    deep backbone stays outside the colour. `"clade"` is the old behaviour —
    each band hugs its own ancestor, which makes every left edge a fact about
    that clade at the cost of a ragged column. `"full"` runs from the root,
    also flush, and reads as *these rows* rather than *this clade* (iTOL's
    coloured ranges work this way); the trunk is then under colour too.
    Aligning needs to know every group at once, so a lone `node=`/`taxa=`
    highlight still hugs its clade whatever the default.

- **`tip_labels(only=...)`** names just the tips you ask for and leaves the rest
  unlabelled. `max_labels` thins *evenly*, which is right when the point is to
  sample a dense tree and wrong when the point is to call out particular
  sequences — a 5,000-tip figure that names the two references the paper is
  about cannot get there by thinning, because the tip it needs is not on the
  every-*n*th grid. Takes one name or a list; it overrides the thinning rather
  than fighting it. `max_labels=0` now means *name none of them*, which is what
  a dense circular tree wants once the ring tracks carry the story (it was
  silently read as "no limit").
- **`tip_points(only=...)`** marks just the tips a measurement exists for — a
  structure prediction, a validated hit — instead of every tip, which says the
  opposite. The colour scale is still built from the ones drawn, so the key
  matches what is on the figure.
- **A bar ring has a scale** (`ring(geom="bar")`, off with `axis=False`):
  circles at round values, drawn *over* the bars so the reference is there both
  where a bar falls short and where it runs past, with the baseline and the
  outermost circle numbered in the fan opening. Only those two are numbered, and
  they read *across* the spoke rather than along it — a ring is a fifth of the
  tree's radius wide and a 7 pt number a tenth of it tall in the same figure, so
  numbering all three circles along the spoke stacks them into one blot. The
  circles in between are what a grid is for: they divide the interval the two
  numbers give. Same argument as the rectangular `bar_track(axis=True)`: without
  it the bars are decoration, and it matters most exactly where they look alike.
- **`ring(geom="stack")`** draws one ring in which every column is a *segment*
  of each tip's bar, scaled so the segments fill the ring's width — a
  composition ring: *of the sequences at this tip, what fraction came from each
  domain*. Neither of the other geoms can say that: a `"bar"` carries one
  number and a `"tile"` only which category won. The parts are normalised per
  tip, so the ring reads as a composition however the raw counts are scaled, and
  a tip whose parts are all zero is left blank rather than normalised to
  something. `title=` names the key, since the columns are its entries and none
  of them can double as its heading.

- **`tip_labels(italic=...)` decides per name, not per figure.** Genus and
  species are italicised by convention and a catch-all label like `others` or
  `unclassified` is not, so one flag for the whole figure is wrong wherever both
  appear in it. `True`/`False` still set every label; `"taxa"` italicises the
  ones that read as an organism; a function receives the name and decides.
  `pt.looks_like_a_taxon` is the test `"taxa"` uses, exported so it can be
  wrapped rather than reimplemented — it is deliberately shallow (a capitalised
  first letter, some letters, and not one of a short list of catch-all words),
  because nothing short of a taxonomy lookup does better and a rule that can be
  read off the name is easier to predict than a clever one.

### Fixed
- **A failed external engine run now says why.** Every `subprocess.run(...,
  check=True)` across `infer_iqtree`/`infer_fasttree`/`infer_raxmlng`/
  `infer_rapidnj`, and MAFFT/MUSCLE in `align_external`, raised
  `CalledProcessError` on a non-zero exit, whose default message is "Command
  [...] returned non-zero exit status 1" -- the actual reason (a model name
  the engine did not recognise, a constraint with too few taxa, an alignment
  it could not parse) was captured in `stderr` and then silently discarded
  every single time a run failed, contradicting `ml.py`'s own claim to be
  "graceful... a clear, actionable error rather than failing obscurely."
  Both now surface it.
- **`build_tree(..., seed=)` reaches an external ML engine's own seed.**
  `infer_raxmlng` was hardcoding `--seed 1` regardless of what `seed=` the
  caller passed -- reproducible, since it was always the same fixed value,
  but not *user*-reproducible: the parameter existed and looked like it
  should apply, the same shape of bug as the dead `ml_tool` parameter above.
  `infer_iqtree` gained a `seed=`/`-seed` it never had at all. Left unset,
  RAxML-NG's default is still the fixed `1` rather than its own random one --
  a tree that changes between identical-looking calls is its own kind of
  surprise -- while IQ-TREE keeps picking its own, matching prior behaviour.
- **`ml_model` reaches an external ML engine's `model=` too, when set.** It
  was silently native-engine-only: `build_tree(ml_engine="iqtree",
  ml_model="GTR")` built successfully with IQ-TREE's own default (`"MFP"`)
  regardless, no error, nothing to notice -- the exact docstring claim
  ("protein sequences work the same way -- pass `ml_model=`... for
  `method='ml'`") was simply false for every engine but the native one.
  `ml_model`'s default changed from the literal `"HKY85"` to `None`, so
  `build_tree` can tell "unset" from "the caller chose HKY85" -- forwarded
  verbatim only when actually set, so an unset `ml_model` still leaves
  IQ-TREE's `"MFP"`/RAxML-NG's `"GTR+G"` alone rather than overwriting either
  engine's own, better default with phytreon's native-engine one. IQ-TREE and
  RAxML-NG do not share one model-string dialect, so a string valid for one
  is not guaranteed valid for the other -- now the caller's problem when they
  opt into both an external engine and an explicit model at once, not
  something this can validate for them. `ml_gamma` has no equivalent for an
  external engine (it encodes rate variation in the model string itself,
  e.g. `"GTR+G4"`) and stays native-only; asking for it or `ml_model` on
  FastTree, which does not accept a model string of its own, now raises
  rather than a bare `TypeError` from a stray keyword argument.
- **`save()` no longer leaks the figure it draws.** It builds a fresh
  matplotlib `Figure` via `draw()`, saves it, and returns only the output
  path -- never the figure itself -- so nothing else could ever close it
  either. Harmless for one figure; most of `examples/` (and any real script
  saving one figure per sample) calls `save()` in a loop, and every past
  call stayed registered with pyplot, artists and all, for the rest of the
  process. Found from the same `RuntimeWarning: More than 20 figures have
  been opened` that this whole test suite has been quietly printing on
  every run. `PanelFigure.save()` had the identical leak independently and
  gets the identical fix.
- **A dead `ml_tool` parameter is gone from `build_tree()`.** It was declared
  in the signature, defaulted to `"iqtree"`, and read nowhere in the
  function body -- `ml_engine` is what actually selects native vs. an
  external tool, so passing `ml_tool="fasttree"` silently ran the native
  engine instead, no error, nothing to notice. Grepped the repo first:
  nothing outside its own declaration ever referenced it.
- **A clade bracket on a circular tree goes outside the rings and the names.**
  It was drawn at the tip circle, so on any figure with ring tracks the bracket
  ran *under* them and its name across them — illegible, and it hid the data it
  was put there to point at. It now sits outside whatever the rings claimed and
  outside the tip labels **in its own sector**, and reads outward so the clades
  on the left half are not upside down. Per-sector rather than per-figure on
  purpose: measured against the whole drawing, one long name on the far side
  pushes every bracket out to that radius and leaves the rest floating in white
  space. Where the names end depends on the font and the figure size, so the
  renderer resolves it after drawing them (`push_out`/`push_span` on a scene
  item — the radial counterpart of the `align` shift the right-side tracks use).
- **A closed fan no longer seats its last leaf on its first.** At `extent=360`
  the two ends of the arc are the same angle, and spacing the leaf rows over
  `n-1` gaps put leaf *n* exactly where leaf 1 already was — one branch, one
  label and one stack of ring tiles drawn over another, at the seam, in every
  full circle. A closed fan now spaces over `n` gaps, so the last leaf lands one
  slot short of coming round; an open fan keeps `n-1`, which is what puts a leaf
  on each end of the arc it was given.
- **The scale bar on a circular tree sits in the fan opening.** It was placed
  from the whole drawing's bounding box, which on a round figure is several
  times the tree's radius: the tick marks and the label offset came out
  bounding-box sized, so the number landed on the rule it labels, and the bar
  itself — correct in branch-length units — read as a typo next to a figure that
  wide. It now measures against the tree's radius and goes in the one wedge a
  fan leaves empty of tips, branches and ring segments, just outside the tip
  circle where that wedge is widest (near the centre it is narrower than the
  number, and the two ends of the number spill over the leaves flanking the
  opening). A closed circle has no such wedge, so there it goes below the
  outermost ring. The default length is now a tenth of what the tree *spans*
  rather than of its radius — on a circle the radius is half the span, which
  made the default bar a stub shorter than its own label.
- **Keys stack by measurement instead of by guess.** Each legend block's height
  was estimated from its entry count, which left a ragged column and, on the
  tree-of-life figure, pushed the colour bar a quarter of the figure's height
  below the key it belonged with. Each block is now measured and the next sits
  directly beneath it.
- **A bar track has a scale.** Without one the bars were decoration — you could
  compare two of them and read nothing off either. It matters most exactly
  where the bars look alike: 16S lengths run 1238 to 1584, and against the zero
  baseline a bar chart owes them every bar is between 78% and 100% of full
  width. The axis is what makes "nearly equal" read as the finding it is rather
  than as a drawing that failed. Turn it off with `bar_track(..., axis=False)`.
- **A packed network centres on its bounding box**, not on the origin the
  packer started from, so a graph with one big piece and a few stragglers is
  not sized for the far side and margined on the near one: the content now
  fills 92–95% of the square, against 62% before.
- **Figures no longer print names on top of each other.** Three defects, all
  found by counting real glyph collisions rather than by eye — measured off
  glyph ink and off *oriented* boxes, since a circular tree's labels are
  rotated and matplotlib's window extent for rotated text is the axis-aligned
  envelope, nearly twice the ink in each direction.
  - The default canvas grows with the tip count for round layouts, as it
    always has for rectangular ones. It was a flat 8×8 whatever the tree, so a
    circular tree ran its names together past about thirty tips: **215
    collisions at 106 taxa, now none**. The rule — `(0.105 × labels + 4.85) ×
    size / 10` inches, floored at the old 8 — comes from the smallest square
    that measured clean at 30, 50, 75 and 106 tips, at two label sizes.
    Thinning with `max_labels` lowers the count and shrinks the figure with it.
  - Unrooted trees point each name away from the middle of the drawing instead
    of along its own branch. A terminal branch can point anywhere, including
    back across the figure, so neighbouring names ran parallel and straight
    through one another: **164 collisions at 106 taxa, now 84** at the same
    canvas. The rest cannot be fixed by placement — the layout puts taxa a gene
    cannot separate at the same point, so their names land there too, at any
    figure size — so it now says so and points at `max_labels` or
    `layout="circular"`, which stays clean at every size tested.
  - Ring names run along the spoke rather than across it. Across the spoke a
    name takes up its own width in angle, which is more than the fan opening
    leaves, so it printed itself over the ring it was naming.
- **A dendrogram now grows sideways.** Its leaves run along x, so the *width*
  has to follow the tip count; it was being given the tall, narrow default a
  rectangular tree wants, which crushed 106 names into eight inches — **151
  collisions, now none**.
- **A ring's name is only drawn when there is room for it.** It goes in the fan
  opening, but the first and last sectors each hang half a sector into that
  opening, and on a tree with few tips the sectors swallow it whole. Rather
  than print the name across the ring it is naming, it is now left to the
  legend, which already carries the column name. Measured over tip counts 10 to
  106 and sizes 6 to 14 pt: no name lands on data anywhere.
- **A name inside a coloured block picks its own ink.** It was white whatever
  the block, which measures 2.2:1 against the paler half of the palette, below
  the 3:1 floor for large text. Choosing black or white by the block's
  luminance guarantees at least 4.58:1.
- **Nothing is drawn thinner than a press can hold.** `MIN_STROKE_PT` (0.3 pt)
  is now a package-wide floor, applied in both backends, so no element and no
  combination of options can emit a rule that shows on screen and disappears
  in the journal — verified with deliberately abusive settings
  (`branches(size=0.05)`, `edge_width=0.02`, `leader_width=0.02`): the
  thinnest stroke rendered anywhere is exactly the floor.
  - The elements that scale width by data now map *into* a range beginning at
    the floor rather than being clipped onto it. Clipping was the first fix
    and it was the wrong one: it gave every edge below the floor the same
    width, so the weakest hits — the ones a cutoff exists to include — all
    came out identical. Mapping keeps them distinguishable and printable at
    once, with opacity carrying the rest of the fade.
- **A faint connection line is still on the page.** `connections(color="value")`
  ran the default sequential ramp, whose pale end measures **1.23:1 against
  white paper** — and no opacity rescues that, since fully opaque the line is
  still that colour. A filled mark may legitimately fade to nearly white (the
  cell has area, its neighbours give it an edge, pale reads as "low"); a line
  has no area, so it simply is not there. Links now start the ramp where it
  first reaches the 3:1 floor, 48% along, and the colour bar is built from the
  same mapper so the key still matches the lines.
- **A disconnected sequence network is no longer squashed into a corner.**
  Each connected piece is laid out on its own and the pieces are then packed,
  sized by node count so density stays comparable. Laid out in one pass the
  stragglers set the scale: on the 106-sequence 16S graph the component
  holding 90 of those sequences occupied **33% of the frame's width, now
  79%**, and at the stricter cutoff 7% → 24%. A graph in one piece is
  unaffected.
- **Taxa a gene cannot separate no longer hide each other** in a split network:
  markers drawn closer than a marker's width are nudged just far enough apart
  to both be seen.
- New `tests/test_figure_quality.py` holds all of this as regression tests: a
  106-taxon tree drawn rectangular, circular, circular with a ring, as a
  dendrogram and as a tanglegram must print every name without a single
  collision; a round layout must not be scaled into an oval; no default figure
  may set text under 5 pt; a block label must reach 4.5:1 contrast; network
  edges must stay printable.

### Added
- **Sequence-similarity networks** (`SequenceNetwork`) -- the CLANS-style
  cluster map. A tree assumes the sequences align well enough for branch order
  to mean something; for a family only detectable by profile searches that
  assumption fails, and the deep branching a tree reports is then an artefact
  of alignment error rather than a record of descent. This draws the sequence
  space instead: one node per sequence, an edge above a similarity cutoff, laid
  out by Fruchterman-Reingold so mutually-similar groups fall into visible
  clusters.
  - Build from an alignment (`from_alignment`), a distance matrix
    (`from_distances`), or `(name1, name2, similarity)` triples straight out of
    a BLAST report (`from_pairs`).
  - `color_by` (with the same `baseline` / `order` controls as the tree
    figures), `label_clusters`, `label_nodes`, and `components()` to confirm
    that the clusters the eye picks out really are connected components.
  - The layout adds a gravity term to the textbook algorithm: sequence networks
    are disconnected by construction, and without it the isolated nodes drift
    off, set the scale of the picture, and squash the connected core into an
    unreadable dot. The default was chosen by measuring the trade-off against
    cluster separation on a real 16S graph.
  - New `docs/tutorials/network.md` and `examples/network_demo.py`.
- `support_labels(attr=[...])` combines several support values into one label
  (`"88/95/0.98"`). A topology cross-checked by ML bootstrap, SH-aLRT and a
  Bayesian posterior carries three of them, and that is how such trees are
  annotated in print; calling `support_labels()` three times stacked all three
  on exactly the same point. `min_value=` hides weakly supported nodes.
- Reading a BEAST/MrBayes tree now also copies the posterior probability onto
  `node.support` when that is empty -- it *is* the clade support of a Bayesian
  tree, and without it the default `support_labels()` drew nothing on a BEAST
  tree while working fine on a bootstrapped one. An existing support value is
  never overwritten.
- **Ribbon tanglegrams.** `TangleFigure.ribbons(groups)` joins the two trees
  with one filled band per group instead of a line per tip. A line answers
  "where did this taxon go"; a band answers "where did this *group* go", which
  is usually the question -- a flat band is a group both trees agree on, a
  twisted one is a group they place differently, and that does not emerge from
  a bundle of individual lines.
- **Multi-panel figures.** `pt.panels([...], ncols=4)` lays any figures out as
  one grid -- trees, tanglegrams, DensiTree clouds, sequence networks -- with
  `share_legend=True` drawing each distinct colour key once beside the grid
  instead of repeating it in every cell, and optional a/b/c panel labels.
- **Domain architectures and gene neighbourhoods.** `TreeFigure.domains(data)`
  draws each tip's architecture beside it, so a domain gained, lost or swapped
  along a clade can be read off the figure. Plain names give evenly-spaced
  blocks, `(name, length)` pairs draw to scale, `arrows=True` gives block
  arrows for an operon, and a negative length flips the arrow for a gene on the
  other strand.
- **Split networks.** `SplitNetwork` draws conflicting splits as boxes rather
  than resolving them away: build from a bootstrap/posterior tree set
  (`from_trees`, weight = fraction of trees containing the split), a distance
  matrix, or an alignment. `conflicts()` lists the conflicting pairs so a box
  can be confirmed rather than assumed. New `docs/tutorials/splitnet.md`.
  - The drawing is **planar**, by the same construction SplitsTree uses: the
    taxa are put in a circular ordering so that each split cuts the circle as
    a single arc — one chord — and the network is the dual of the resulting
    chord arrangement, each edge drawn perpendicular to its own chord. A
    regression test asserts zero crossings over random circular split systems.
  - Boxes and conflicting split pairs now match exactly at every `max_splits`
    setting, on random systems and on real bootstrap sets: the drawing neither
    invents a box nor swallows a conflict. `dropped` lists any split no
    circular ordering can lay out as an arc, rather than drawing it wrong.
  - Terminal splits are drawn and are exempt from `max_splits`. They can never
    conflict, so they never cost a box; without them every taxon sits on an
    internal node and the names land on each other.
  - Past a dozen taxa the labels move out to a ring with hairline leaders
    (`label_ring`), because a split network bunches taxa wherever the splits
    between them are short.
- **Neighbor-Net.** `pt.neighbor_net(names, matrix)` — a split network
  straight from a distance matrix, both halves of the Bryant–Moulton method.
  - **The ordering** (`neighbornet_ordering`) is built by agglomeration on the
    distances. Neighbour joining picks the two nodes to *merge*, and merging is
    what costs it: from then on the pair is one subtree and nothing can come
    between them. This picks the two to stand *next to each other* and merges
    nothing, so clusters are chains growing at both ends, and a chain can seat
    two taxa together that no tree groups. Measured over 40 distance matrices
    built from known circular split systems: every generating split came back
    drawable in **40 of 40**, against **3 of 40** for a neighbour-joining
    tree's leaf order, which left a fifth of the split weight undrawable. With
    10% noise on the distances it still manages 34 of 40 and keeps 99.7% of
    the weight. Select with `ordering="neighbornet"` (default) or `"tree"`.
  - **The weight estimation** fits every split the ordering can draw — all
    `n(n-1)/2` of them — to the distances by non-negative least squares. This
    is not a refinement, it is the difference between a network and a tree:
    splits read off one tree are compatible with each other by construction,
    so the previous behaviour could not draw a box however conflicted the data
    was. On the 18-taxon 16S matrix — same matrix all three times — 33 splits
    and **0 boxes** reading the tree, 38 splits and 11 boxes fitting against
    the tree's ordering (4.6% residual), 38 splits and **24 boxes** fitting
    against the agglomerative ordering (**2.7% residual**).
  - The fit grows as the fourth power of the taxon count (0.02 s at 30, 0.8 s
    at 60, 4.5 s at 80); past `FIT_MAX_TAXA` it reads the tree and warns,
    because a drawing with no boxes is a claim about the data and should not
    be made on the quiet. `net.estimated` records which route ran.
  - A split and its complement are now recognised as one split. They were
    drawn as two chords lying on top of each other, which lost the cells that
    belong between them — and a rooted tree hands over both sides of its root
    edge, so this was the ordinary case rather than a rare one.
- `support_labels(stack=True, prefixes=[...])` writes several support values on
  separate lines (`p:1.00` / `b:100` / `n:0.98`) as well as the joined form,
  which stays readable with four values where a slash-joined string does not.
- New `examples/figstyles_demo.py` and four more gallery figures.
- README gallery now shows the 0.3.0 drawing styles (collapsed clades, node
  interval bars, connections, DensiTree) that had no visual example before.

### Changed
- The ruff rule set is now stated explicitly in `pyproject.toml`
  (`select = ["E4", "E7", "E9", "F"]`). `ruff` is unpinned in the dev extras so
  CI installs the newest release; a version that widens its defaults would
  otherwise turn the build red without a line of this project changing. That is
  not hypothetical -- the currently installed ruff reports 550+ style findings
  across untouched files under the exact command CI runs.

### Fixed
- A title on a circular/equal-aspect figure could overlap a tip label pointing
  toward a pole. `Scene.bounds()` only sees a label's anchor point, not how far
  its rotated glyphs actually reach -- a label pointing straight up extends its
  own text length past the anchor, an amount no fixed data-unit padding can
  anticipate (it depends on font size and string length, not tree geometry).
  The renderer now measures the actual rendered extent of every label and
  expands the axes limits to contain it before placing the title.



### Added
- **Annotated NEXUS input** (`Tree.read(path, fmt="beast")`, also `"mrbayes"`,
  or `pt.read_beast`). BEAST and MrBayes write their per-node estimates into
  NEXUS comments -- node ages, 95% HPD intervals, posterior clade
  probabilities, per-branch rates -- and a plain NEXUS read keeps the topology
  and discards all of it. The new reader parses those `[&key=value,...]` blocks
  onto `node.data`, applies the TRANSLATE table, and flattens `{lower, upper}`
  intervals to `<name>_lower` / `<name>_upper`, which are exactly the keys
  `node_bars()` reads by default -- so a dated Bayesian tree now plots straight
  from the file. `pt.parse_annotation` exposes the comment parser on its own,
  and `tree_index=` selects from a file holding a posterior sample.

### Fixed
- Nested `collapse_clade()` treated an already-collapsed inner clade as a
  single tip: the outer summary undercounted the hidden tips and its triangle
  stopped well short of the real farthest leaf.
- A collapsed clade whose hidden leaves sit at zero distance (a cladogram, or
  zero-length branches) drew a zero-size, invisible triangle.
- `DensiTreeFigure(layout="circular")` scaled only the x coordinate when
  rescaling a tree onto the reference's depth. Depth is the radius on a polar
  layout, so the overlay smeared into an ellipse reaching far outside the
  reference tree; both coordinates are now scaled.
- A collapsed clade's triangle reaches out to the hidden clade's farthest leaf,
  which can be well beyond the collapsed tree's own depth. `max_x` did not
  account for it, so everything keyed to it cut through the triangle: it was
  clipped off the figure, rings were drawn on top of it, and aligned tip labels
  landed inside it. The layout now includes the collapsed extent, in the units
  it draws in (branch length, or edges on a cladogram).
- `node_bars()` and `time_axis()` each defaulted `present` to 0 independently,
  so setting it on the axis alone silently shifted every bar off the scale it
  is read against. `node_bars()` now follows the figure's time axis whatever
  order the two were added in; an explicit `present=` still wins.
- A comment following a branch length (where BEAST writes per-branch rates)
  stopped the Newick parser dead at the opening bracket.

## [0.3.0] — 2026-07-18

### Added
- **Tanglegrams** for comparing two trees of the same taxa -- e.g. a tree built
  from genomic data against one built from transcriptomic data, or two
  inference methods on one alignment:
  - `TangleFigure(left, right)` draws the two trees facing each other and links
    their shared tips. Each side is an ordinary `TreeFigure` (`fig.left` /
    `fig.right`), so every existing element, layout and colour scale works on
    either tree; ready-made `TreeFigure`s can be passed in directly. Trees with
    only partly overlapping taxa are supported -- unmatched tips are drawn but
    left unlinked.
  - `.untangle()` rotates clades to minimise crossing links (greedy hill-climb
    over single rotations; `fix="left"`/`"right"`/`None`). Rotation reorders
    children only, so topology and branch lengths are untouched -- untangling
    changes how the trees read, never what they say.
  - `.connect(...)` styles the links: a literal colour, a data column from the
    left tree's tips (with legend), or `highlight_discordant=True` to colour
    every link that crosses another.
  - `treeops.crossing_number(t1, t2)` counts crossing links (inversions between
    the two tip orders, O(n log n)) and `treeops.untangle(t1, t2)` exposes the
    rotation search on its own. Both are documented as *display* discordance:
    zero crossings does not imply identical trees, so read `robinson_foulds`
    alongside.
  - Both trees are labelled by default (`tip_labels="both"`, also `"left"`,
    `"right"` or `False`). The middle band that carries the labels and links is
    sized from the actual rendered text width, and the figure widens for long
    taxon names, so species names fit instead of colliding across the middle;
    `gap=` and `connect(inset=...)` override the estimate.
  - New `docs/tutorials/tanglegram.md` and `examples/tanglegram_demo.py` (the
    demo compares neighbour joining, UPGMA and parsimony on the bundled 16S
    alignment and shows both the discordant and the deceptive-agreement case).
  - New bundled dataset `examples/data/big16S*` -- 106 taxa across 25
    prokaryotic phyla (91 Bacteria, 15 Archaea) fetched from NCBI by
    `examples/data/fetch_large_16S.py`, for demos that need a large tree.

- **More drawing styles**, filling the gaps against iTOL / ggtree / FigTree:
  - **Collapsed clades.** `treeops.collapse_clade(tree, node)` compresses a
    clade to a single tip and `TreeFigure.collapsed_clades()` draws it as a
    triangle whose two sides reach the clade's nearest and farthest hidden
    leaf (iTOL's convention), so the wedge shows how deep and how ragged the
    hidden group is. Tip labels are offset past the triangle. Works on
    rectangular and circular layouts.
  - **Node interval bars.** `TreeFigure.node_bars(lower=, upper=)` draws the
    95% HPD age interval across each node -- the standard annotation on a
    dated Bayesian tree (FigTree's "node bars", ggtree's `geom_range`). Read
    as ages on the same scale as `time_axis()`, or as raw x with `as_age=False`.
  - **Connections.** `TreeFigure.connections(pairs)` draws curved links between
    arbitrary tips for horizontal gene transfer, gene sharing or co-occurrence
    (iTOL's `DATASET_CONNECTION`). On a circular layout the curves bend toward
    the centre; on a rectangular one they bow out past the tips. Accepts
    `(a, b)` / `(a, b, value)` tuples or a DataFrame, and `color="value"`.
  - **DensiTree.** `DensiTreeFigure(trees)` overlays a whole set of trees
    translucently so topological uncertainty is visible instead of hidden
    behind one summary tree. Trees are first rotated onto a shared tip order
    via `untangle`, which changes only how they read.
  - **Scale bar.** `TreeFigure.scale_bar()` -- a compact branch-length scale
    (ggtree's `geom_treescale`) that, unlike `time_axis()`, assumes nothing
    about branch lengths being time and works on any layout.
  - New `docs/tutorials/tree_styles.md` and `examples/tree_styles_demo.py`.
- **Grey the default state, colour the exceptions.** `baseline=` on `ring()`,
  `heatmap()`, `tip_points()` (and `build_color_scale`) renders the named
  level(s) neutral grey. Baseline levels no longer consume a palette slot, so
  the remaining levels keep the strongest hues. When one level covers most of
  the tree, colouring it as loudly as the rare ones spends the figure's ink on
  its least informative part and buries the exceptions.
- `order=` sets the categorical legend order explicitly; levels were otherwise
  sorted alphabetically, which rarely matches a meaningful progression.
- Legend keys now match the mark they stand for: filled layers (rings,
  heatmaps) get square swatches instead of dots.
- `ring(leaders=True)` draws a faint dotted guide from each tip out to the
  first ring. On a phylogram the tips sit at very different radii, so most stop
  well short of the rings and it stops being obvious which sector belongs to
  which tip. (Dropping branch lengths does *not* fix this -- a cladogram still
  places tips at different depths -- and stretching tips to a common radius
  would misrepresent the branch lengths, so a guide line is the honest fix.)

### Fixed
- A column read by two elements no longer emits the legend twice. Colouring
  tip points and a ring by the same `phylum` column stacked two identical
  legends; `RenderContext.add_scale` now ignores a key it has already
  registered.
- `color="some_column"` where the column was never joined onto the tree used
  to sail through as a literal colour and fail much later inside matplotlib as
  `Invalid RGBA argument: 'phylum'`. It now raises immediately, naming the
  columns that *are* available and pointing at `tree.join_data(df, on="name")`.
- `ring()` and `heatmap()` no longer break up into slivers on large trees. Both
  drew a fixed hairline separator around every cell; once a tree passes a few
  hundred tips that stroke is as wide as the cell itself, so a metadata ring
  rendered as a comb of thin white-gapped stripes instead of solid colour
  bands, and blocks of shared values became unreadable. Past ~150 tips the
  separator is now dropped and each cell is stroked in its own fill, so
  neighbouring cells meet with no anti-aliased seam. Force either behaviour
  with the new `separators=True/False` argument.
- `ring(pad_angle=...)` is an absolute angle, so on a large tree it could
  exceed a whole sector and produce inverted (negative-width) wedges; it is now
  clamped to leave a sliver of every sector standing.

### Changed
- Extracted the shared `draw`/`save`/`show` plumbing from `TreeFigure` into an
  internal `_Renderable` base so `TangleFigure` gets identical backend dispatch
  and export behaviour (including editable-text SVG) rather than a second copy.

## [0.2.2] — 2026-07-17

### Fixed
- `TreeFigure.branches(color=, size=)` now replaces the tree's skeleton layer
  instead of stacking a second one on top of it. `TreeFigure(tree)` already
  draws a default skeleton (`skeleton=True`), so calling `.branches(size=...)`
  again to change branch width globally used to draw a second, differently
  sized line directly over the first -- e.g. requesting a thinner line left a
  visible fringe of the original, thicker default line peeking out from
  underneath. `.branches(...)` is now guaranteed to be a single, clean, global
  override of branch color/width.

## [0.2.1] — 2026-07-15

### Changed
- SVG export (`TreeFigure.save("...svg")`) now keeps every label as a real
  `<text>` element instead of outlining glyphs to vector paths (sets
  matplotlib's `svg.fonttype="none"` for `.svg` only). The figure stays fully
  editable after importing into PowerPoint (Insert → Picture → Convert to
  Shape), Illustrator, or Inkscape -- labels can be recoloured, moved, and
  re-typed. PDF/PNG/HTML output is unchanged.

## [0.2.0] — 2026-07-15

### Added
- New `layout="circular_slanted"` (aliases `slanted_circular` / `fan_slanted`)
  -- the polar counterpart of the `slanted` layout. Each edge is a single
  straight diagonal line drawn directly from parent to child, instead of the
  ordinary circular tree's radial-spoke-plus-arc elbow, giving a cleaner
  "starburst" look that reads better on many circular trees.

### Changed
- Refreshed the default plotting colours for a more restrained, publication-
  ready look. The categorical default is now `CURATED_PALETTE` -- eight muted,
  colourblind-safe hues in a fixed order, replacing the old over-saturated
  evenly-spaced HCL hue wheel (which read as a "default plot" and, more
  seriously, collapsed green/yellow to ΔE 5.6 under protanopia -- indistinct
  for red-green colourblind readers). The new order was verified against the
  Machado-2009 CVD model (worst adjacent ΔE >= 12; >= 11 all-pairs at the full
  eight). Category counts above eight extend with a *muted* hue wheel so the
  extra colours stay in the same register. The raw wheel is still available as
  `palette="hue"`, and the named ColorBrewer palettes (`set2`/`dark2`/`tab10`)
  are unchanged. The default continuous ramp is now a single-hue blue running
  light (low, recedes) -> deep (high, salient) -- the conventional direction --
  instead of the old dark -> washed-out-light gradient. Neutral quantitative
  bar/ring fills changed from a muddy tan to a calm slate. No API changes;
  only default colours, so existing figures re-render with the new palette.

### Fixed
- Stacked continuous colorbars overlapped: with two or more continuous
  heatmap/ring columns, each colorbar's title was drawn as a rotated
  side-label, so adjacent titles ran together (e.g. "lifespanbody_mass")
  and clipped off the top edge. Titles are now horizontal labels placed
  above each bar (matching the categorical legend titles), with headroom so
  the first never clips.
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
- `reconstruct_ancestral_mutations()` (`phytreon/infer/lineage.py`): traces
  back *which* mutation/scar arose on *which* branch under the same
  Camin-Sokal model `camin_sokal_score()` minimizes, writing
  `node.data["mutations_acquired"]` for every node -- the piece that turns
  a `lineage_tree()` topology into an actual reconstructed process rather
  than just a set of relationships. `sankoff_score()`'s postorder DP loop
  is now a shared `_sankoff_dp()` helper (zero behavior change there) so
  both functions reuse the identical computation.
- Real-data validation for the general (non-CRISPR) lineage-tracing path:
  `examples/mutation_lineage_demo.py` reconstructs a clonal cell tree from
  Hou et al. 2012's real single-cell exome-sequencing mutation calls (Cell
  148:873-885, 18 genes/58 cells, via the SCITE package), reporting
  reconstruction cost and ancestral mutation acquisitions on real, noisy
  data.

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
