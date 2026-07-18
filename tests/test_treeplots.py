"""Drawing styles beyond the plain tree: collapsed clades, node interval
bars, node-to-node connections, scale bars and DensiTree clouds."""
import matplotlib
matplotlib.use("Agg")

import pytest

import phytreon as pt


def _polys(fig):
    return fig._build().scene.polygons


def _paths(fig):
    return fig._build().scene.paths


# --------------------------------------------------------------------------
# collapsed clades
# --------------------------------------------------------------------------
def test_collapse_clade_records_the_hidden_shape():
    tr = pt.Tree.from_newick("((A:1,B:3):1,(C:1,D:1):2);")
    clade = tr.get_mrca(["A", "B"])
    pt.collapse_clade(tr, clade, name="AB")
    assert tr.n_leaves == 3                     # A+B replaced by one tip
    info = clade.data["_collapsed"]
    assert info["n"] == 2
    assert info["near"] == pytest.approx(1.0)   # to A
    assert info["far"] == pytest.approx(3.0)    # to B
    assert sorted(info["leaves"]) == ["A", "B"]
    assert clade.is_leaf and clade.name == "AB"


def test_collapse_clade_autonames_with_the_hidden_count():
    tr = pt.Tree.from_newick("((A:1,B:1,E:1):1,(C:1,D:1):2);")
    node = tr.get_mrca(["A", "B", "E"])
    pt.collapse_clade(tr, node)
    assert "+2" in node.name                    # first tip plus 2 more


def test_collapse_clade_rejects_a_leaf():
    tr = pt.datasets.primates()
    with pytest.raises(ValueError, match="leaf"):
        pt.collapse_clade(tr, tr.get_node("Human"))


def test_collapsed_clades_draws_one_triangle_per_clade():
    tr = pt.Tree.from_newick("((A:1,B:3):1,(C:1,D:1):2);")
    node = tr.get_mrca(["A", "B"])
    pt.collapse_clade(tr, node)
    tris = [p for p in _polys(pt.TreeFigure(tr).collapsed_clades())
            if len(p.points) == 3]
    assert len(tris) == 1
    apex, near_pt, far_pt = tris[0].points
    assert apex[0] == pytest.approx(node.x)                 # apex at the node
    assert near_pt[0] == pytest.approx(node.x + 1.0)        # side to nearest leaf
    assert far_pt[0] == pytest.approx(node.x + 3.0)         # side to farthest


def test_collapsed_clade_label_clears_the_triangle():
    # the label used to sit on top of the triangle it belongs to
    tr = pt.Tree.from_newick("((A:1,B:3):1,(C:1,D:1):2);")
    node = tr.get_mrca(["A", "B"])
    pt.collapse_clade(tr, node, name="AB")
    ctx = pt.TreeFigure(tr).collapsed_clades().tip_labels()._build()
    lab = next(lb for lb in ctx.scene.labels if lb.text == "AB")
    assert lab.x > node.x + 3.0                 # past the farthest triangle corner


def test_collapsed_clades_renders_on_a_circular_tree(tmp_path):
    tr = pt.datasets.random_tree(40, seed=3)
    node = next(n for n in tr.traverse()
                if not n.is_leaf and not n.is_root and len(n.get_leaves()) > 4)
    pt.collapse_clade(tr, node)
    out = tmp_path / "c.png"
    pt.TreeFigure(tr, layout="circular").collapsed_clades().tip_labels().save(str(out))
    assert out.exists() and out.stat().st_size > 1000


# --------------------------------------------------------------------------
# node interval bars (95% HPD)
# --------------------------------------------------------------------------
def _dated_tree():
    tr = pt.Tree.from_newick("((A:1,B:1):2,(C:2,D:2):1);")
    for node in tr.traverse():
        if node.is_leaf:
            continue
        node.data["height_95_lower"] = 1.0
        node.data["height_95_upper"] = 2.0
    return tr


def test_node_bars_place_the_interval_on_the_time_axis():
    tr = _dated_tree()
    fig = pt.TreeFigure(tr).node_bars()
    ctx = fig._build()
    maxx = ctx.layout.max_x
    bars = [p for p in ctx.scene.paths if p.width == 3.0]
    assert len(bars) == 3                        # every internal node
    for bar in bars:
        (x0, _), (x1, _) = bar.points
        # ages run backwards from the present at max_x
        assert x0 == pytest.approx(maxx - 2.0)
        assert x1 == pytest.approx(maxx - 1.0)


def test_node_bars_can_take_raw_x_instead_of_ages():
    tr = _dated_tree()
    ctx = pt.TreeFigure(tr).node_bars(as_age=False)._build()
    bar = next(p for p in ctx.scene.paths if p.width == 3.0)
    assert bar.points[0][0] == pytest.approx(1.0)
    assert bar.points[1][0] == pytest.approx(2.0)


def test_node_bars_complain_when_no_node_carries_the_interval():
    tr = pt.datasets.primates()
    with pytest.raises(ValueError, match="height_95_lower"):
        pt.TreeFigure(tr).node_bars()._build()


def test_node_bars_reject_circular_layouts():
    tr = _dated_tree()
    with pytest.raises(NotImplementedError, match="rectangular"):
        pt.TreeFigure(tr, layout="circular").node_bars()._build()


# --------------------------------------------------------------------------
# connections
# --------------------------------------------------------------------------
def test_connections_draw_one_curve_per_pair():
    tr = pt.datasets.primates()
    pairs = [("Human", "Baboon"), ("Chimp", "Macaque")]
    ctx = pt.TreeFigure(tr).connections(pairs)._build()
    curves = [p for p in ctx.scene.paths if p.zorder == 0.8]
    assert len(curves) == 2
    for curve, (a, b) in zip(curves, pairs):
        assert len(curve.points) > 8                       # really a curve
        na, nb = tr.get_node(a), tr.get_node(b)
        assert curve.points[0] == pytest.approx((na.x, na.y))
        assert curve.points[-1] == pytest.approx((nb.x, nb.y))


def test_connections_reject_names_not_in_the_tree():
    tr = pt.datasets.primates()
    with pytest.raises(ValueError, match="Nessie"):
        pt.TreeFigure(tr).connections([("Human", "Nessie")])._build()


def test_connections_colour_by_value_adds_a_scale():
    tr = pt.datasets.primates()
    pairs = [("Human", "Baboon", 0.1), ("Chimp", "Macaque", 0.9)]
    ctx = pt.TreeFigure(tr).connections(pairs, color="value")._build()
    curves = [p for p in ctx.scene.paths if p.zorder == 0.8]
    assert len({c.color for c in curves}) == 2             # really recoloured
    assert ctx.scene.colorbars or ctx.scene.legends


def test_connections_accept_a_dataframe():
    import pandas as pd
    tr = pt.datasets.primates()
    df = pd.DataFrame({"a": ["Human"], "b": ["Baboon"], "w": [0.4]})
    ctx = pt.TreeFigure(tr).connections(df)._build()
    assert len([p for p in ctx.scene.paths if p.zorder == 0.8]) == 1


def test_connections_bend_toward_the_centre_on_a_circular_tree():
    # iTOL's CENTER_CURVES look: the curve's midpoint must fall well inside
    # the straight chord, i.e. nearer the centre of the circle
    import math
    tr = pt.datasets.primates()
    ctx = pt.TreeFigure(tr, layout="circular").connections(
        [("Human", "Baboon")])._build()
    curve = next(p for p in ctx.scene.paths if p.zorder == 0.8)
    mid = curve.points[len(curve.points) // 2]
    ends = (curve.points[0], curve.points[-1])
    chord_mid = ((ends[0][0] + ends[1][0]) / 2, (ends[0][1] + ends[1][1]) / 2)
    r = lambda p: math.hypot(*p)                          # noqa: E731
    assert r(mid) < r(chord_mid)


# --------------------------------------------------------------------------
# scale bar
# --------------------------------------------------------------------------
def test_scale_bar_picks_a_round_length():
    from phytreon.plot.elements import _ScaleBar
    assert _ScaleBar._nice(0.0834) == pytest.approx(0.05)
    assert _ScaleBar._nice(0.3) == pytest.approx(0.2)
    assert _ScaleBar._nice(7.0) == pytest.approx(5.0)
    assert _ScaleBar._nice(1.0) == pytest.approx(1.0)


def test_scale_bar_draws_a_bar_of_the_requested_length():
    tr = pt.datasets.primates()
    ctx = pt.TreeFigure(tr).scale_bar(length=0.05)._build()
    bars = [p for p in ctx.scene.paths if p.width == 1.4]
    horiz = [p for p in bars if p.points[0][1] == p.points[1][1]]
    assert horiz
    x0, x1 = horiz[0].points[0][0], horiz[0].points[1][0]
    assert (x1 - x0) == pytest.approx(0.05)
    assert any(lb.text == "0.05" for lb in ctx.scene.labels)


# --------------------------------------------------------------------------
# DensiTree
# --------------------------------------------------------------------------
def _same_taxa_trees(n=5, ntip=8):
    ref = pt.datasets.random_tree(ntip, seed=1)
    names = ref.leaf_names()
    trees = [ref]
    for s in range(2, n + 1):
        t = pt.datasets.random_tree(ntip, seed=s)
        for leaf, nm in zip(t.leaves(), names):
            leaf.name = nm
        trees.append(t)
    return trees


def test_densitree_overlays_every_tree():
    trees = _same_taxa_trees(5)
    ctx = pt.DensiTreeFigure(trees, tip_labels=False)._build()
    edges = [p for p in ctx.scene.paths if p.zorder == 1]
    # each tree contributes its own edges, all translucent
    assert len(edges) > 5 * (trees[0].n_leaves - 1)
    assert all(0 < p.opacity < 1 for p in edges)


def test_densitree_aligns_tip_order_to_the_reference():
    trees = _same_taxa_trees(4)
    ref_before = list(trees[0].leaf_names())
    fig = pt.DensiTreeFigure(trees, tip_labels=False)
    fig._build()
    # untangling rotates the others toward the reference, never the reference
    assert trees[0].leaf_names() == ref_before
    after = [pt.crossing_number(trees[0], t) for t in trees[1:]]
    fig2 = pt.DensiTreeFigure(_same_taxa_trees(4), align=False, tip_labels=False)
    fig2._build()
    raw = [pt.crossing_number(fig2.trees[0], t) for t in fig2.trees[1:]]
    assert sum(after) <= sum(raw)


def test_densitree_opacity_thins_out_for_bigger_samples():
    few = pt.DensiTreeFigure(_same_taxa_trees(3), tip_labels=False)
    many = pt.DensiTreeFigure(_same_taxa_trees(40, ntip=6), tip_labels=False)
    assert many.alpha < few.alpha


def test_densitree_needs_at_least_one_tree():
    with pytest.raises(ValueError, match="at least one"):
        pt.DensiTreeFigure([])


def test_densitree_renders(tmp_path):
    out = tmp_path / "densi.png"
    pt.DensiTreeFigure(_same_taxa_trees(6)).titled("cloud").save(str(out))
    assert out.exists() and out.stat().st_size > 1000
