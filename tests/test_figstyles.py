"""Figure styles taken from the comparative-genomics literature: ribbon
tanglegrams, multi-panel grids, domain-architecture tracks, stacked support
values, and split networks."""
import matplotlib
matplotlib.use("Agg")

import pytest

import phytreon as pt
from phytreon.plot.splitnet import conflicting, splits_from_tree


# --------------------------------------------------------------------------
# ribbons
# --------------------------------------------------------------------------
def _grouped_pair():
    a = pt.Tree.from_newick("(((A1,A2),(A3,A4)),((B1,B2),(C1,C2)));")
    b = pt.Tree.from_newick("(((A1,A2),(C1,C2)),((A3,A4),(B1,B2)));")
    return a, b, {n: n[0] for n in a.leaf_names()}


def _ribbons(ctx):
    return [p for p in ctx.scene.polygons if p.label]


def test_ribbons_draw_one_band_per_group():
    a, b, groups = _grouped_pair()
    ctx = pt.TangleFigure(a, b).ribbons(groups)._build()
    bands = _ribbons(ctx)
    assert sorted(p.label for p in bands) == ["A (4)", "B (2)", "C (2)"]
    assert len({p.facecolor for p in bands}) == 3
    for band in bands:
        assert len(band.points) > 8          # eased, not a straight quad


def test_a_ribbon_spans_all_of_its_group_on_both_sides():
    a, b, groups = _grouped_pair()
    fig = pt.TangleFigure(a, b).ribbons(groups, width=0.0)
    ctx = fig._build()
    band = next(p for p in _ribbons(ctx) if p.label.startswith("A "))
    ys = [y for _, y in band.points]
    a_rows = [t.y for t in a.leaves() if t.name.startswith("A")]
    assert min(ys) <= min(a_rows) + 1e-9
    assert max(ys) >= max(a_rows) - 1e-9


def test_ribbons_can_read_a_column_off_the_left_tree():
    a, b, groups = _grouped_pair()
    for tip in a.leaves():
        tip.data["clade"] = groups[tip.name]
    ctx = pt.TangleFigure(a, b).ribbons("clade")._build()
    assert len(_ribbons(ctx)) == 3
    assert [t for t, _ in ctx.scene.legends] == ["clade"]


def test_ribbons_complain_about_a_column_that_was_never_joined():
    a, b, _ = _grouped_pair()
    with pytest.raises(ValueError, match="join_data"):
        pt.TangleFigure(a, b).ribbons("nope")


def test_ribbons_replace_the_per_tip_links():
    a, b, groups = _grouped_pair()
    plain = pt.TangleFigure(a, b)._build()
    ribbon = pt.TangleFigure(a, b).ribbons(groups)._build()
    links = [p for p in plain.scene.paths if p.zorder == 0.4]
    still = [p for p in ribbon.scene.paths if p.zorder == 0.4]
    assert links and not still


# --------------------------------------------------------------------------
# panel grids
# --------------------------------------------------------------------------
def test_panels_lay_out_a_grid_and_drop_the_spare_axes():
    figs = [pt.TreeFigure(pt.datasets.random_tree(6, seed=i)) for i in range(5)]
    grid = pt.panels(figs, ncols=3)
    assert (grid.nrows, grid.ncols) == (2, 3)
    fig = grid.draw()
    assert len(fig.axes) >= 6
    # the sixth cell has no panel, so its frame must be off
    assert not fig.axes[5].axison


def test_panels_show_a_shared_key_once_instead_of_per_panel():
    import pandas as pd
    figs = []
    for i in range(4):
        tr = pt.datasets.random_tree(6, seed=i)
        df = pd.DataFrame({"name": tr.leaf_names(), "g": ["x", "y"] * 3})
        tr.join_data(df, on="name")
        figs.append(pt.TreeFigure(tr).tip_points(color="g"))

    shared = pt.panels(figs, ncols=2, share_legend=True)
    fig = shared.draw()
    assert getattr(fig, "_phytreon_extra_artists", None)      # one shared key
    # each panel's own key was suppressed
    assert all(not ax.get_legend() for ax in fig.axes)

    apart = pt.panels(figs, ncols=2, share_legend=False)
    per_panel = apart.draw()
    assert getattr(per_panel, "_phytreon_extra_artists", None) is None


def test_panels_reject_a_title_count_that_does_not_match():
    figs = [pt.TreeFigure(pt.datasets.random_tree(5, seed=1))]
    with pytest.raises(ValueError, match="titles"):
        pt.panels(figs, titles=["a", "b"])
    with pytest.raises(ValueError, match="at least one"):
        pt.panels([])


def test_panels_accept_every_figure_type(tmp_path):
    trees = [pt.datasets.primates() for _ in range(3)]
    mixed = [
        pt.TreeFigure(pt.datasets.primates()).titled("tree"),
        pt.TangleFigure(pt.datasets.primates(), pt.datasets.primates(),
                        tip_labels=False).titled("tanglegram"),
        pt.DensiTreeFigure(trees, tip_labels=False).titled("densitree"),
        pt.SequenceNetwork.from_pairs([("a", "b", 0.9)]).titled("network"),
    ]
    out = tmp_path / "mixed.png"
    pt.panels(mixed, ncols=2).save(str(out))
    assert out.exists() and out.stat().st_size > 1000


# --------------------------------------------------------------------------
# domain architecture
# --------------------------------------------------------------------------
def _arch_tree():
    tr = pt.Tree.from_newick("((P1:.1,P2:.1):.1,P3:.2);")
    arch = {"P1": [("A", 100), ("B", 200)],
            "P2": [("A", 100), ("B", 200), ("C", 50)],
            "P3": [("B", 200)]}
    return tr, arch


def test_domain_track_draws_one_block_per_domain():
    tr, arch = _arch_tree()
    ctx = pt.TreeFigure(tr).domains(arch)._build()
    blocks = [p for p in ctx.scene.polygons if p.label and "|" in p.label]
    assert len(blocks) == 6                       # 2 + 3 + 1
    assert [t for t, _ in ctx.scene.legends] == ["domain"]
    # one colour per distinct domain name, shared across proteins
    by_domain = {}
    for b in blocks:
        by_domain.setdefault(b.label.split("| ")[1], set()).add(b.facecolor)
    assert all(len(v) == 1 for v in by_domain.values())
    assert len(by_domain) == 3


def test_domains_to_scale_makes_a_longer_protein_draw_longer():
    tr, arch = _arch_tree()
    def width_of(protein, to_scale):
        ctx = pt.TreeFigure(tr).domains(arch, to_scale=to_scale)._build()
        xs = [x for p in ctx.scene.polygons if p.label
              and p.label.startswith(protein) for x, _ in p.points]
        return max(xs) - min(xs)

    assert width_of("P2", True) > width_of("P3", True)      # 350aa vs 200aa
    assert width_of("P2", False) == pytest.approx(width_of("P3", False), rel=0.05)


def test_domain_arrows_flip_for_a_negative_length():
    tr = pt.Tree.from_newick("(X:.1,Y:.1);")
    ctx = pt.TreeFigure(tr).domains(
        {"X": [("g", 100)], "Y": [("g", -100)]}, arrows=True)._build()
    shapes = {p.label.split("|")[0].strip(): p.points
              for p in ctx.scene.polygons if p.label and "|" in p.label}
    assert len(shapes["X"]) == 5 and len(shapes["Y"]) == 5   # block arrows

    def tip_x(pts):                        # the arrowhead's point
        xs = [x for x, _ in pts]
        return max(xs), min(xs)
    fwd_max, fwd_min = tip_x(shapes["X"])
    rev_max, rev_min = tip_x(shapes["Y"])
    # the forward arrow's apex is its rightmost x, the reverse arrow's is left
    fwd_apex = [x for x, y in shapes["X"] if x == fwd_max]
    rev_apex = [x for x, y in shapes["Y"] if x == rev_min]
    assert len(fwd_apex) == 1 and len(rev_apex) == 1


def test_domains_reject_names_that_match_no_tip():
    tr, _ = _arch_tree()
    with pytest.raises(ValueError, match="no tip name matches"):
        pt.TreeFigure(tr).domains({"nobody": ["A"]})._build()


def test_domains_reject_a_circular_layout():
    tr, arch = _arch_tree()
    with pytest.raises(NotImplementedError, match="rectangular"):
        pt.TreeFigure(tr, layout="circular").domains(arch)._build()


# --------------------------------------------------------------------------
# stacked support values
# --------------------------------------------------------------------------
def test_support_values_can_stack_with_prefixes():
    tr = pt.Tree.from_newick("((A:1,B:1):2,(C:2,D:2):1);")
    for node in tr.traverse():
        if not node.is_leaf:
            node.data.update(p=1.0, b=100, n=0.98)
    ctx = pt.TreeFigure(tr).support_labels(
        attr=["p", "b", "n"], stack=True, prefixes=["p", "b", "n"])._build()
    texts = [lb.text for lb in ctx.scene.labels]
    assert texts and all(t == "p:1\nb:100\nn:0.98" for t in texts)


def test_prefix_count_must_match_the_keys():
    tr = pt.datasets.primates()
    with pytest.raises(ValueError, match="prefixes"):
        pt.TreeFigure(tr).support_labels(attr=["a", "b"], prefixes=["only"])


# --------------------------------------------------------------------------
# split networks
# --------------------------------------------------------------------------
def test_two_conflicting_splits_make_exactly_one_box():
    names = ["A", "B", "C", "D"]
    splits = [(frozenset(["A", "B"]), 1.0), (frozenset(["B", "C"]), 1.0)]
    net = pt.SplitNetwork(names, splits)
    assert len(net.conflicts()) == 1
    verts, edges = net._median_network()
    assert len(verts) == 4 and len(edges) == 4        # a 4-cycle: the box


def test_compatible_splits_stay_a_tree():
    names = ["A", "B", "C", "D"]
    splits = [(frozenset(["A", "B"]), 1.0), (frozenset(["A", "B", "C"]), 1.0)]
    net = pt.SplitNetwork(names, splits)
    assert net.conflicts() == []
    verts, edges = net._median_network()
    assert len(edges) == len(verts) - 1               # acyclic


def test_conflicting_predicate_matches_the_definition():
    U = frozenset("ABCD")
    assert conflicting(frozenset("AB"), frozenset("BC"), U)      # all four meet
    assert not conflicting(frozenset("AB"), frozenset("ABC"), U)  # nested
    assert not conflicting(frozenset("AB"), frozenset("CD"), U)   # disjoint


def test_selection_keeps_conflicts_that_a_top_n_cut_would_discard():
    # splits present in over half the trees are the majority consensus, which
    # is compatible by construction -- ranking by weight and truncating would
    # throw away every conflict, i.e. exactly what the network exists to show
    names = list("ABCDEF")
    strong = [(frozenset(names[:i + 2]), 0.9 - 0.01 * i) for i in range(4)]
    weak = [(frozenset(["B", "C"]), 0.3)]        # conflicts with the nested set
    net = pt.SplitNetwork(names, strong + weak, max_splits=5)
    assert len(net.conflicts()) >= 1
    assert any(s == frozenset(["B", "C"]) for s, _ in net.splits)


def test_from_trees_weights_splits_by_how_often_they_appear():
    a = pt.Tree.from_newick("((A,B),(C,D));")
    b = pt.Tree.from_newick("((A,B),(C,D));")
    c = pt.Tree.from_newick("((A,C),(B,D));")
    net = pt.SplitNetwork.from_trees([a, b, c])
    weights = {s: w for s, w in net.splits}
    ab = next(w for s, w in weights.items() if s == frozenset(["A", "B"]))
    assert ab == pytest.approx(2 / 3)
    with pytest.raises(ValueError, match="at least one tree"):
        pt.SplitNetwork.from_trees([])


def test_splits_from_tree_skips_trivial_edges():
    tr = pt.Tree.from_newick("((A:1,B:1):2,(C:1,D:1):2);")
    got = splits_from_tree(tr, tr.leaf_names())
    sides = [s for s, _ in got]
    assert all(1 < len(s) < 4 for s in sides)     # no single tips, no whole set


def test_taxa_sharing_a_vertex_are_still_drawn_apart():
    # taxa no kept split separates land on one vertex; drawn there they would
    # sit exactly on top of one another
    names = ["A", "B", "C", "D"]
    splits = [(frozenset(["A", "B"]), 1.0)]        # C and D are not separated
    net = pt.SplitNetwork(names, splits)
    pos = net.positions
    assert len({(round(x, 6), round(y, 6)) for x, y in pos.values()}) == 4


def test_split_network_renders(tmp_path):
    tr = pt.datasets.random_tree(10, seed=1)
    trees = [tr, pt.datasets.random_tree(10, seed=2)]
    for t in trees[1:]:
        for leaf, nm in zip(t.leaves(), tr.leaf_names()):
            leaf.name = nm
    out = tmp_path / "net.png"
    pt.SplitNetwork.from_trees(trees).titled("splits").save(str(out))
    assert out.exists() and out.stat().st_size > 1000
