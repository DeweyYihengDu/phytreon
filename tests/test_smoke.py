"""Smoke + correctness tests. Run: pytest -q"""
import matplotlib
matplotlib.use("Agg")

import phytreon as pt


def test_io_roundtrip():
    tr = pt.datasets.primates()
    assert tr.n_leaves == 7
    nwk = tr.write()
    tr2 = pt.Tree.from_newick(nwk)
    assert sorted(tr2.leaf_names()) == sorted(tr.leaf_names())


def test_newick_quotes_special_characters_in_names():
    # Regression: names containing reserved Newick punctuation (parens,
    # commas, colons, semicolons) were written completely unquoted, so a
    # single taxon like "weird(name),here" silently split into three
    # separate leaves on round-trip, with no error raised at all.
    names = ["Homo sapiens", "weird(name),here", "quote'inside",
            "colon:test", "plain_name", "semi;colon"]
    root = pt.Node()
    for nm in names:
        root.add_child(pt.Node(name=nm, length=1.0))
    tr = pt.Tree(root=root)
    tr2 = pt.Tree.from_newick(tr.write())
    assert sorted(tr2.leaf_names()) == sorted(names)


def test_newick_leaves_simple_names_unquoted():
    tr = pt.datasets.primates()
    nwk = tr.write()
    assert "'" not in nwk


def test_mrca():
    tr = pt.datasets.primates()
    node = tr.get_mrca(["Human", "Chimp"])
    assert set(node.leaf_names()) == {"Human", "Chimp"}


def test_mrca_rejects_missing_taxon():
    # Regression: get_mrca(["Human", "NotExisting"]) used to silently return
    # the MRCA of whichever names *were* found (i.e. just "Human"'s own leaf
    # node), giving no indication that "NotExisting" isn't in the tree.
    import pytest
    tr = pt.datasets.primates()
    with pytest.raises(ValueError, match="NotExisting"):
        tr.get_mrca(["Human", "NotExisting"])
    # strict=False keeps the old lenient behaviour for callers that want it
    node = tr.get_mrca(["Human", "NotExisting"], strict=False)
    assert node.name == "Human"


def test_layouts_render_both_backends(tmp_path):
    tr = pt.datasets.primates()
    for layout in ("rectangular", "circular"):
        p = pt.TreeFigure(tr, layout=layout).tip_labels()
        fig = p.draw(backend="mpl")
        assert fig is not None
        p.save(str(tmp_path / f"{layout}.png"))
        p.save(str(tmp_path / f"{layout}.html"))
    assert (tmp_path / "rectangular.png").exists()
    assert (tmp_path / "circular.html").exists()


def test_svg_export_keeps_editable_text(tmp_path):
    # SVG must emit real <text> elements (svg.fonttype="none"), not outlined
    # glyph paths -- so labels stay editable after importing into PowerPoint
    # or a vector editor.
    import re
    tr = pt.datasets.primates()
    out = tmp_path / "tree.svg"
    pt.TreeFigure(tr).tip_labels().save(str(out))
    svg = out.read_text(encoding="utf-8")
    texts = re.findall(r"<text[^>]*>(.*?)</text>", svg, re.S)
    assert texts, "SVG has no <text> elements -- labels were outlined to paths"
    assert "Human" in " ".join(texts)     # a real tip label is editable text


def test_branches_size_is_a_single_global_override():
    # __init__(skeleton=True) already draws one branches layer; calling
    # .branches(size=...) again must replace it, not stack a second layer
    # on top (a second, differently-sized layer would leave a fringe of the
    # first layer's colour/width showing through).
    from phytreon.plot.elements import _Branches
    tr = pt.datasets.primates()
    p = pt.TreeFigure(tr).branches(size=0.3, color="red")
    layers = [e for e in p._elements if isinstance(e, _Branches)]
    assert len(layers) == 1
    assert layers[0].size == 0.3 and layers[0].color == "red"


def test_neighbor_joining():
    names = ["A", "B", "C", "D"]
    d = [[0, 5, 9, 9], [5, 0, 10, 10], [9, 10, 0, 8], [9, 10, 8, 0]]
    tr = pt.neighbor_joining(names, d)
    assert sorted(tr.leaf_names()) == names


def test_ace_parsimony_clean_clade():
    tr = pt.datasets.primates()
    clade = set(tr.get_mrca(["Human", "Orangutan"]).leaf_names())
    trait = {n: ("ape" if n in clade else "other") for n in tr.leaf_names()}
    res = pt.ace_parsimony(tr, trait)
    assert res["score"] == 1            # single clade-defining change


def test_ace_parsimony_polytomy_counts_all_steps():
    # Regression: a 3-way polytomy with 3 disjoint states used to score 1
    # (binary-Fitch all-children intersection), but correctly needs 2 state
    # changes to explain 3 distinct states from one ancestor.
    tr = pt.Tree.from_newick("(A:1,B:1,C:1);")
    res = pt.ace_parsimony(tr, {"A": "x", "B": "y", "C": "z"})
    assert res["score"] == 2


def test_ace_ml_recovers_clade():
    tr = pt.datasets.primates()
    clade = set(tr.get_mrca(["Human", "Orangutan"]).leaf_names())
    trait = {n: ("ape" if n in clade else "other") for n in tr.leaf_names()}
    pt.ace_ml(tr, trait)
    great_apes = tr.get_mrca(["Human", "Gorilla"])
    assert great_apes.data["ace_state"] == "ape"
    assert great_apes.data["ace_probs"]["ape"] > 0.9


def test_stochastic_map_paints_and_matches_ace(tmp_path):
    tr = pt.datasets.primates()
    clade = set(tr.get_mrca(["Human", "Orangutan"]).leaf_names())
    trait = {n: ("ape" if n in clade else "other") for n in tr.leaf_names()}
    pt.stochastic_map(tr, trait, n=150, seed=0)
    ga = tr.get_mrca(["Human", "Gorilla"])
    assert ga.data["ace_probs"]["ape"] > 0.8          # confident, matches ace
    assert all("paint_segments" in n.data
               for n in tr.traverse() if not n.is_root)
    p = pt.TreeFigure(tr).painted_branches().tip_labels()
    p.save(str(tmp_path / "painted.png"))
    assert (tmp_path / "painted.png").exists()


def test_stochastic_map_records_all_ard_rates():
    # Regression: only params[0] was recorded as "rate", losing information
    # for SYM/ARD models which fit more than one rate parameter.
    tr = pt.datasets.primates()
    clade = set(tr.get_mrca(["Human", "Orangutan"]).leaf_names())
    trait = {n: ("ape" if n in clade else "other") for n in tr.leaf_names()}
    pt.stochastic_map(tr, trait, n=50, model="ARD", seed=0)
    info = tr.data["stochastic_map"]
    assert info["model"] == "ARD"
    assert len(info["rates"]) == 2          # 2 states -> 2 asymmetric ARD rates


def test_continuous_ace_bounds():
    tr = pt.datasets.primates()
    trait = {n: float(i) for i, n in enumerate(tr.leaf_names())}
    vals = pt.ace_continuous(tr, trait)
    lo, hi = min(trait.values()), max(trait.values())
    assert all(lo <= v <= hi for v in vals.values())


def test_time_axis_with_geo(tmp_path):
    tr = pt.datasets.random_tree(10, seed=1)
    pt.scale_clade(tr, tr.root, 100.0 / max(n.depth() for n in tr.leaves()))
    p = pt.TreeFigure(tr).time_axis(geo=True, gridlines=True).tip_labels()
    ctx = p._build()
    assert ctx.scene.polygons               # geological bands drawn
    assert any(lb.text == "Cretaceous" for lb in ctx.scene.labels)
    p.save(str(tmp_path / "ts.png"))
    assert (tmp_path / "ts.png").exists()


def test_daylight_layout_valid():
    import math
    from phytreon.layout import get_layout
    tr = pt.datasets.random_tree(30, seed=4)
    get_layout("daylight").apply(tr)
    pts = [(t.x, t.y) for t in tr.leaves()]
    assert all(math.isfinite(x) and math.isfinite(y) for x, y in pts)
    assert len({(round(x, 6), round(y, 6)) for x, y in pts}) == len(pts)  # distinct
    # daylight must not collapse the tree to a point (it spreads tips out)
    xs = [x for x, _ in pts]
    ys = [y for _, y in pts]
    assert (max(xs) - min(xs)) > 0 and (max(ys) - min(ys)) > 0


def test_new_layouts_render(tmp_path):
    tr = pt.datasets.random_tree(15, seed=2)
    for layout in ("daylight", "dendrogram", "inward_circular", "slanted",
                   "circular_slanted"):
        p = pt.TreeFigure(tr, layout=layout).tip_labels()
        ctx = p._build()
        assert ctx.scene.paths              # branches were emitted
        out = tmp_path / f"{layout}.png"
        p.save(str(out))
        assert out.exists() and out.stat().st_size > 1000


def test_circular_slanted_draws_straight_edges_no_arcs():
    # the polar counterpart of the slanted layout: every edge is a single
    # straight 2-point line parent->child, and there is no arc connector.
    from phytreon.layout import get_layout
    tr = pt.datasets.random_tree(12, seed=1)
    lay = get_layout("circular_slanted").apply(tr)
    for node in tr.traverse():
        if node.is_root:
            continue
        bp = lay.branch_path(node)
        assert len(bp) == 2                         # straight line, not a spoke
        assert bp[0] == (node.parent.x, node.parent.y)
        assert bp[1] == (node.x, node.y)
        assert lay.child_connector(node) == []      # no arc


def test_continuous_scale_makes_colorbar(tmp_path):
    import pandas as pd
    tr = pt.datasets.primates()
    mat = pd.DataFrame({"name": tr.leaf_names(),
                        "mass": [62, 45, 160, 75, 8, 11, 25]}).set_index("name")
    p = pt.TreeFigure(tr).tip_labels().heatmap(mat)
    ctx = p._build()
    assert ctx.scene.colorbars and ctx.scene.colorbars[0][0] == "mass"
    assert not any(t == "mass" for t, _ in ctx.scene.legends)   # not a swatch legend
    p.save(str(tmp_path / "cb.png"))
    assert (tmp_path / "cb.png").exists()


def test_heatmap_accepts_name_column_not_indexed():
    import pandas as pd
    tr = pt.datasets.primates()
    # a 'name' column but NOT set as the index -- heatmap() must still match
    # rows to tips by that column, the same as ring()/bar_track() already do
    mat = pd.DataFrame({"name": tr.leaf_names(),
                        "mass": [62, 45, 160, 75, 8, 11, 25]})
    ctx = pt.TreeFigure(tr).tip_labels().heatmap(mat)._build()
    assert len(ctx.scene.polygons) == tr.n_leaves


def test_node_pies(tmp_path):
    tr = pt.datasets.primates()
    trait = {n: ("ape" if n in tr.get_mrca(["Human", "Orangutan"]).leaf_names()
                 else "other") for n in tr.leaf_names()}
    pt.ace_ml(tr, trait)
    p = pt.TreeFigure(tr).node_pies().tip_labels()
    ctx = p._build()
    assert ctx.scene.polygons                       # pie wedges drawn
    p.save(str(tmp_path / "pie.png"))
    assert (tmp_path / "pie.png").exists()


def test_numeric_scale_handles_numpy_ints():
    import numpy as np
    from phytreon.plot.figure import build_color_scale
    sc = build_color_scale("x", [np.int64(1), np.int64(5), np.int64(9)])
    assert sc.continuous          # numpy ints must read as continuous, not categorical


def test_point_size_handles_numpy_ints():
    # Regression: _resolve_size() used isinstance(v, (int, float)), which is
    # False for numpy.int64/float64 (the dtype pandas/numpy normally produce),
    # silently falling back to a constant marker size.
    import numpy as np
    tr = pt.datasets.primates()
    for tip, v in zip(tr.leaves(), [10, 20, 30, 40, 50, 60, 70]):
        tip.data["mass"] = np.int64(v)
    ctx = pt.TreeFigure(tr).tip_points(size="mass")._build()
    sizes = {m.size for m in ctx.scene.markers}
    assert len(sizes) > 1          # sizes must vary, not all fall back to 6.0


def test_metadata_rings(tmp_path):
    import pandas as pd
    tr = pt.datasets.primates()
    labs = tr.leaf_names()
    meta = pd.DataFrame({
        "name": labs,
        "habitat": ["urban", "forest", "forest", "forest", "forest", "savanna", "savanna"],
        "mass": [62, 45, 160, 75, 8, 11, 25],   # numpy int64 -> continuous ring
    }).set_index("name")
    p = (pt.TreeFigure(tr, layout="circular")
         .ring(meta, columns=["habitat", "mass"])
         .tip_labels())
    ctx = p._build()
    assert ctx.outer_radius > ctx.ring_base            # rings reserved radial space
    legend_titles = [t for t, _ in ctx.scene.legends]
    cbar_titles = [t for t, _, _, _ in ctx.scene.colorbars]
    assert "habitat" in legend_titles                  # categorical -> legend
    assert "mass" in cbar_titles                        # continuous -> colorbar
    p.save(str(tmp_path / "rings.png"))
    assert (tmp_path / "rings.png").exists()


def test_group_clade_colours_branches():
    tr = pt.datasets.primates()
    apes = tr.get_mrca(["Human", "Gibbon"])
    pt.group_clade(tr, {apes: "apes"}, key="lineage", default="other")
    assert tr.get_node("Human").data["lineage"] == "apes"
    assert tr.get_node("Macaque").data["lineage"] == "other"


def test_shape_and_bar_ring(tmp_path):
    import pandas as pd
    tr = pt.datasets.random_tree(20, seed=1)
    tips = tr.leaf_names()
    import random
    rng = random.Random(0)
    meta = pd.DataFrame({
        "name": tips,
        "city": [rng.choice(["A", "B", "C"]) for _ in tips],
        "score": [rng.random() for _ in tips],
    }).set_index("name")
    tr.join_data(meta.reset_index(), on="name")
    p = (pt.TreeFigure(tr, layout="circular")
         .tip_points(shape="city", size=5)
         .ring(meta, columns=["score"], geom="bar"))
    ctx = p._build()
    markers = {m.marker for m in ctx.scene.markers}
    assert len(markers) == 3                      # three distinct shapes for city
    titles = [t for t, _ in ctx.scene.legends]
    assert "city" in titles                        # shape legend present
    p.save(str(tmp_path / "bar.png"))
    assert (tmp_path / "bar.png").exists()

    # regression: the plotly backend must render the actual marker shapes,
    # not just the legend swatches
    plotly_fig = p.draw(backend="plotly")
    real_marker_traces = [t for t in plotly_fig.data
                          if t.mode == "markers" and t.x and len(t.x) > 1]
    symbols = {s for t in real_marker_traces for s in (t.marker.symbol or [])}
    assert len(symbols) == 3


def test_rectangular_tracks_stack(tmp_path):
    import pandas as pd
    tr = pt.datasets.primates()
    labs = tr.leaf_names()
    meta = pd.DataFrame({
        "name": labs,
        "grp": ["a", "a", "b", "b", "b", "c", "c"],
        "val": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
    }).set_index("name")
    tr.join_data(meta.reset_index(), on="name")
    p = (pt.TreeFigure(tr).tip_labels()
         .heatmap(meta[["grp"]], width=0.05)
         .bar_track(meta, "val", width=0.3))
    ctx = p._build()
    assert ctx.track_cursor > ctx.layout.max_x      # tracks stacked rightward
    p.save(str(tmp_path / "tracks.png"))
    assert (tmp_path / "tracks.png").exists()


def test_alignment_track_raster(tmp_path):
    tr = pt.datasets.primates()
    aln = pt.Alignment(tr.leaf_names(), ["ACGT-ACGT"] * tr.n_leaves)
    p = pt.TreeFigure(tr).tip_labels().alignment(aln)
    ctx = p._build()
    assert len(ctx.scene.rasters) == 1
    r = ctx.scene.rasters[0]
    assert r.codes.shape[0] == tr.n_leaves and r.codes.shape[1] == aln.ncol
    p.save(str(tmp_path / "aln.png"))
    assert (tmp_path / "aln.png").exists()


def test_plotly_shifts_aligned_paths():
    # clade_label() emits an aligned Path (the bracket bar) that must be
    # pushed past the tip labels in the plotly backend, exactly like the
    # aligned Polygon/Label/Raster primitives already are.
    tr = pt.datasets.primates()
    apes = tr.get_mrca(["Human", "Gibbon"])
    fig = pt.TreeFigure(tr).tip_labels().clade_label("Apes", node=apes)
    ctx = fig._build()
    aligned = [p for p in ctx.scene.paths if p.align]
    assert aligned, "expected clade_label() to emit an aligned Path"
    orig_x = aligned[0].points[0][0]

    plotly_fig = fig.draw(backend="plotly")
    shifted = False
    for trace in plotly_fig.data:
        xs = [x for x in (trace.x or []) if x is not None]
        if len(xs) >= 2 and len({round(x, 6) for x in xs}) == 1:
            if abs(xs[0] - orig_x) > 1e-9:
                shifted = True
    assert shifted, "aligned Path was not shifted in the plotly backend"


def _ring_cells(fig):
    # ring/heatmap cells are the polygons carrying a per-tip hover label
    return [p for p in fig._build().scene.polygons if p.label and "|" in p.label]


def test_dense_rings_drop_the_sector_separator():
    # a hairline separator is wider than the sector itself once a tree has a
    # few hundred tips, which turned the ring into a comb of white slivers
    # instead of solid colour bands
    import pandas as pd
    from phytreon.plot.elements import _RING_DENSE_TIPS

    def cells(n):
        tr = pt.datasets.random_tree(n, seed=5)
        df = pd.DataFrame({"name": tr.leaf_names(),
                           "grp": ["a", "b"] * (n // 2) + ["a"] * (n % 2)})
        return _ring_cells(pt.TreeFigure(tr, layout="circular").ring(df))

    sparse = cells(40)
    assert sparse and all(c.edgecolor == "white" for c in sparse)

    dense = cells(_RING_DENSE_TIPS + 60)
    assert dense
    # every dense cell is stroked in its own fill, so abutting sectors leave
    # no anti-aliased hairline between them
    assert all(c.edgecolor == c.facecolor for c in dense)


def test_ring_separators_can_be_forced_either_way():
    import pandas as pd
    tr = pt.datasets.random_tree(400, seed=6)
    df = pd.DataFrame({"name": tr.leaf_names(), "grp": ["a"] * 400})
    on = _ring_cells(pt.TreeFigure(tr, layout="circular").ring(df, separators=True))
    assert all(c.edgecolor == "white" for c in on)      # kept despite 400 tips

    small = pt.datasets.random_tree(30, seed=6)
    df2 = pd.DataFrame({"name": small.leaf_names(), "grp": ["a"] * 30})
    off = _ring_cells(pt.TreeFigure(small, layout="circular")
                      .ring(df2, separators=False))
    assert all(c.edgecolor == c.facecolor for c in off)


def test_large_pad_angle_cannot_invert_a_sector():
    # pad_angle is an absolute angle; on a big tree it can exceed the whole
    # sector, which used to produce negative-width (inverted) wedges
    import pandas as pd
    tr = pt.datasets.random_tree(500, seed=7)
    df = pd.DataFrame({"name": tr.leaf_names(), "grp": ["a"] * 500})
    fig = pt.TreeFigure(tr, layout="circular").ring(df, pad_angle=5.0)
    for cell in _ring_cells(fig):
        assert len(cell.points) >= 3          # still a real polygon


def test_dense_heatmap_drops_the_cell_separator():
    import pandas as pd
    from phytreon.plot.elements import _RING_DENSE_TIPS
    n = _RING_DENSE_TIPS + 60
    tr = pt.datasets.random_tree(n, seed=8)
    df = pd.DataFrame({"name": tr.leaf_names(), "v": list(range(n))})
    cells = _ring_cells(pt.TreeFigure(tr).heatmap(df))
    assert cells and all(c.edgecolor == c.facecolor for c in cells)


def test_one_column_shared_by_two_elements_yields_one_legend():
    # tip points and a ring both coloured by "phylum" used to emit two
    # identical legends stacked on top of each other
    import pandas as pd
    tr = pt.datasets.random_tree(30, seed=9)
    df = pd.DataFrame({"name": tr.leaf_names(),
                       "phylum": ["a", "b", "c"] * 10})
    tr.join_data(df, on="name")
    ctx = (pt.TreeFigure(tr, layout="circular")
           .tip_points(color="phylum")
           .ring(df, columns=["phylum"]))._build()
    assert [t for t, _ in ctx.scene.legends] == ["phylum"]


def test_unknown_colour_column_says_so_instead_of_failing_in_the_backend():
    # a column that was never joined onto the tree used to sail through as a
    # literal colour and blow up much later as "Invalid RGBA argument"
    import pytest
    tr = pt.datasets.primates()
    with pytest.raises(ValueError, match="join_data"):
        pt.TreeFigure(tr).tip_points(color="phylum")._build()
    # real colours still work untouched
    pt.TreeFigure(tr).tip_points(color="steelblue")._build()
    pt.TreeFigure(tr).tip_points(color="#ff8800")._build()


def test_baseline_category_is_greyed_and_frees_a_palette_slot():
    # when one level covers most of the tree, colouring it as loudly as the
    # rare ones spends the figure's ink on the least informative category
    from phytreon.plot.figure import BASELINE_GREY, build_color_scale
    vals = ["ref"] * 50 + ["odd", "rare"]
    plain = build_color_scale("ctx", vals)
    muted = build_color_scale("ctx", vals, baseline="ref")

    assert dict(muted.legend)["ref"] == BASELINE_GREY
    assert dict(plain.legend)["ref"] != BASELINE_GREY
    # the greyed level consumes no palette slot: the remaining levels get
    # exactly the colours they would have had if it were not in the data
    without = build_color_scale("ctx", ["odd", "rare"])
    assert dict(muted.legend)["odd"] == dict(without.legend)["odd"]
    assert dict(muted.legend)["rare"] == dict(without.legend)["rare"]
    assert len(muted.legend) == 3                      # still fully documented


def test_baseline_accepts_several_levels():
    from phytreon.plot.figure import BASELINE_GREY, build_color_scale
    scale = build_color_scale("c", ["a", "b", "c"], baseline=["a", "b"])
    lut = dict(scale.legend)
    assert lut["a"] == lut["b"] == BASELINE_GREY
    assert lut["c"] != BASELINE_GREY


def test_explicit_category_order_beats_alphabetical():
    from phytreon.plot.figure import build_color_scale
    vals = ["low", "high", "mid"]
    assert [c for c, _ in build_color_scale("c", vals).legend] == \
        ["high", "low", "mid"]                          # alphabetical default
    ordered = build_color_scale("c", vals, order=["low", "mid", "high"])
    assert [c for c, _ in ordered.legend] == ["low", "mid", "high"]
    # a level missing from order still shows up, after the listed ones
    partial = build_color_scale("c", vals, order=["high"])
    assert [c for c, _ in partial.legend] == ["high", "low", "mid"]


def test_filled_layers_get_swatch_legends_and_markers_get_dots():
    import pandas as pd
    tr = pt.datasets.random_tree(20, seed=11)
    df = pd.DataFrame({"name": tr.leaf_names(), "grp": ["a", "b"] * 10,
                       "other": ["x", "y"] * 10})
    tr.join_data(df, on="name")
    ctx = (pt.TreeFigure(tr, layout="circular")
           .tip_points(color="grp")
           .ring(df, columns=["other"]))._build()
    # a ring is a filled area -> swatch; tip points are markers -> dot
    assert ctx.scene.legend_swatch["other"] == "patch"
    assert ctx.scene.legend_swatch["grp"] == "point"


def test_ring_leaders_connect_every_tip_to_the_ring():
    # on a phylogram most tips stop well short of the rings; the dotted guide
    # is what tells the eye which sector belongs to which tip
    import pandas as pd
    tr = pt.datasets.random_tree(40, seed=10)
    df = pd.DataFrame({"name": tr.leaf_names(), "g": ["x"] * 40})
    ctx = pt.TreeFigure(tr, layout="circular").ring(df, leaders=True)._build()
    leads = [p for p in ctx.scene.paths if p.dash == "dot" and p.zorder == 0.3]
    assert len(leads) == tr.n_leaves
    radius = lambda p: (p[0] ** 2 + p[1] ** 2) ** 0.5      # noqa: E731
    for path in leads:
        inner, outer = path.points
        assert radius(outer) > radius(inner)               # points outward
    # and they stop at the ring, not short of it
    ring_inner = min(radius(pt_) for poly in ctx.scene.polygons if poly.label
                     for pt_ in poly.points)
    assert max(radius(p.points[1]) for p in leads) <= ring_inner + 1e-9
