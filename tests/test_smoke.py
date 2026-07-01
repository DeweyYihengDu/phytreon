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


def test_mrca():
    tr = pt.datasets.primates()
    node = tr.get_mrca(["Human", "Chimp"])
    assert set(node.leaf_names()) == {"Human", "Chimp"}


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
    for layout in ("daylight", "dendrogram", "inward_circular", "slanted"):
        p = pt.TreeFigure(tr, layout=layout).tip_labels()
        ctx = p._build()
        assert ctx.scene.paths              # branches were emitted
        out = tmp_path / f"{layout}.png"
        p.save(str(out))
        assert out.exists() and out.stat().st_size > 1000


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
