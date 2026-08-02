# Changelog

All notable changes to phytreon are documented here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Added
- **`highlight(by="column")`** shades every group in a joined column behind the
  branches — each in its own colour, with a key — instead of one
  `highlight(node=...)` call per clade and no way to tell the shades apart.
  Works on rectangular and circular layouts (a band and a wedge respectively).
  A group whose taxa are **not** monophyletic comes out as several bands, one
  per run of adjacent tips, rather than one band over their common ancestor:
  that ancestor reaches down over other groups' taxa, so a single band there
  would colour in tips that do not belong to the group — tidier, and false.
  The shades stay pale enough to read a name on (measured: the darkest leaves
  black text at 12.5:1).
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

### Fixed
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
