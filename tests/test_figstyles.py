"""Figure styles taken from the comparative-genomics literature: ribbon
tanglegrams, multi-panel grids, domain-architecture tracks, stacked support
values, and split networks."""
import itertools
import random

import matplotlib
matplotlib.use("Agg")

import pytest

import phytreon as pt
from phytreon.plot.splitnet import (FIT_MAX_TAXA, circular_ordering,
                                    circular_split_weights, circular_splits,
                                    conflicting, is_circular,
                                    neighbornet_ordering, splits_from_tree)


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
    verts, edges = net._network()
    assert len(verts) == 4 and len(edges) == 4        # a 4-cycle: the box


def test_compatible_splits_stay_a_tree():
    names = ["A", "B", "C", "D"]
    splits = [(frozenset(["A", "B"]), 1.0), (frozenset(["A", "B", "C"]), 1.0)]
    net = pt.SplitNetwork(names, splits)
    assert net.conflicts() == []
    verts, edges = net._network()
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


def _crossings(net):
    """Pairs of network edges that cross away from any shared node.

    Sharing a node is not a crossing, and two edges of the same split lie
    parallel, so only proper interior intersections are counted.
    """
    verts, edges = net._network()
    xy = net._vertex_coords(verts)

    def side(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    hits = 0
    for (i1, j1, _), (i2, j2, _) in itertools.combinations(edges, 2):
        if {i1, j1} & {i2, j2}:
            continue
        p, q, r, s = xy[i1], xy[j1], xy[i2], xy[j2]
        d = [side(r, s, p), side(r, s, q), side(p, q, r), side(p, q, s)]
        if any(abs(v) < 1e-12 for v in d):
            continue
        if (d[0] > 0) != (d[1] > 0) and (d[2] > 0) != (d[3] > 0):
            hits += 1
    return hits


def test_the_drawing_is_planar():
    # the whole point of the circular ordering: a split that is one arc is one
    # chord, and chords can only meet by opening a box, never by crossing
    rng = random.Random(7)
    for _ in range(25):
        n = rng.randint(6, 20)
        names = [str(i) for i in range(n)]
        sides = set()
        while len(sides) < rng.randint(5, 30):
            a, ln = rng.randrange(n), rng.randint(1, n - 1)
            sides.add(frozenset(names[(a + t) % n] for t in range(ln)))
        net = pt.SplitNetwork(names, [(s, rng.random() + 0.05) for s in sides],
                              order=names, max_splits=99)
        assert _crossings(net) == 0


def test_every_conflict_opens_exactly_one_box():
    # for a circular split system the network's independent cycles and its
    # conflicting split pairs are the same count -- no conflict is swallowed
    # and no box is invented
    rng = random.Random(3)
    for _ in range(15):
        n = rng.randint(6, 16)
        names = [str(i) for i in range(n)]
        sides = set()
        while len(sides) < rng.randint(4, 14):
            a, ln = rng.randrange(n), rng.randint(1, n - 1)
            sides.add(frozenset(names[(a + t) % n] for t in range(ln)))
        net = pt.SplitNetwork(names, [(s, 1.0) for s in sides],
                              order=names, max_splits=99)
        verts, edges = net._network()
        assert len(edges) - len(verts) + 1 == len(net.conflicts())


def test_three_mutual_conflicts_draw_as_three_rhombi():
    # the case that catches a median closure: it returns the whole 3-cube,
    # eight nodes, which in the plane can only be drawn with crossings. The
    # chord arrangement returns seven cells -- a hexagon of three rhombi
    names = [str(i) for i in range(6)]
    splits = [(frozenset(["0", "1", "2"]), 1.0),
              (frozenset(["2", "3", "4"]), 1.0),
              (frozenset(["4", "5", "0"]), 1.0)]
    net = pt.SplitNetwork(names, splits, order=names)
    assert len(net.conflicts()) == 3
    verts, edges = net._network()
    assert (len(verts), len(edges)) == (7, 9)
    assert _crossings(net) == 0


def test_circular_ordering_makes_the_hierarchy_contiguous():
    # every compatible split must come out as one arc, or it cannot be drawn
    names = list("ABCDEFGH")
    splits = [(frozenset("ABCD"), 1.0), (frozenset("AB"), 0.9),
              (frozenset("CD"), 0.8), (frozenset("EF"), 0.7),
              (frozenset("BC"), 0.4)]                    # the one conflict
    order = circular_ordering(names, splits)
    assert sorted(order) == sorted(names)
    for side, _ in splits[:4]:
        assert is_circular(side, order)


def test_splits_that_are_not_arcs_are_set_aside_not_drawn_wrong():
    names = list("ABCD")
    # AC and BD cannot both be arcs of any ordering of four taxa
    net = pt.SplitNetwork(names, [(frozenset("AB"), 1.0), (frozenset("AC"), 1.0),
                                  (frozenset("AD"), 1.0)], max_splits=99)
    assert len(net.dropped) + len(net.splits) == 3
    for side, _ in net.splits:
        assert is_circular(side, net.order)


def _four_point_distances():
    """Distances with two groupings in them at once.

    Built as the sum of two conflicting splits -- AB|CD at 1.0 and BC|AD at
    0.4 -- plus terminal lengths, so the answer is known: any method worth
    using has to report both, and the weaker one is the box.
    """
    names = list("ABCD")
    splits = [(frozenset("AB"), 1.0), (frozenset("BC"), 0.4)]
    terminal = 0.5
    mat = [[0.0] * 4 for _ in range(4)]
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i == j:
                continue
            d = 2 * terminal
            d += sum(w for s, w in splits if (a in s) != (b in s))
            mat[i][j] = d
    return names, mat


def test_a_distance_matrix_alone_can_draw_boxes():
    # splits read off one tree are compatible with each other by construction,
    # so taking them and stopping draws a tree however conflicted the data is.
    # Fitting every circular split against the distances is what finds it.
    names, mat = _four_point_distances()
    fitted = pt.SplitNetwork.from_distances(names, mat, estimate=True)
    assert fitted.estimated
    assert len(fitted.conflicts()) == 1

    from_tree = pt.SplitNetwork.from_distances(names, mat, estimate=False)
    assert not from_tree.estimated
    assert from_tree.conflicts() == []


def test_the_fit_recovers_the_weights_it_was_built_from():
    names, mat = _four_point_distances()
    net = pt.SplitNetwork.from_distances(names, mat)
    universe = frozenset(names)

    def weight(*side):
        # a split and its complement are one split, so match either side
        want = frozenset(side)
        return next(w for s, w in net.splits if s in (want, universe - want))

    assert weight("A", "B") == pytest.approx(1.0, abs=1e-6)
    assert weight("B", "C") == pytest.approx(0.4, abs=1e-6)
    for taxon in names:                      # the terminal lengths, too
        assert weight(taxon) == pytest.approx(0.5, abs=1e-6)
    assert len(net.splits) == 6              # and nothing invented


def test_circular_splits_are_every_arc_and_nothing_else():
    order = list("ABCDE")
    got = circular_splits(order)
    assert len(got) == 5 * 4 // 2            # one per pair of taxa
    for side in got:
        assert is_circular(side, order)
    assert len({frozenset(s) for s in got}) == len(got)


def test_fitting_too_many_taxa_is_refused_rather_than_hung():
    names = [str(i) for i in range(FIT_MAX_TAXA + 1)]
    mat = [[0.0 if i == j else 1.0 for j in names] for i in names]
    with pytest.raises(ValueError, match="estimate=False"):
        circular_split_weights(names, mat, names)


def _circular_distances(n, nsplits, rng):
    """Distances that are exactly the sum of a random circular split system.

    The ordering it was built on is the answer the agglomeration has to find:
    every generating split must come back as one arc, or it cannot be drawn.
    """
    names = ["t%d" % i for i in range(n)]
    truth = names[:]
    rng.shuffle(truth)
    sides = set()
    while len(sides) < nsplits:
        start, size = rng.randrange(n), rng.randint(1, n - 1)
        sides.add(frozenset(truth[(start + t) % n] for t in range(size)))
    weights = {s: rng.uniform(0.2, 1.0) for s in sides}
    for nm in names:
        weights.setdefault(frozenset([nm]), rng.uniform(0.2, 1.0))
    mat = [[0.0] * n for _ in range(n)]
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i != j:
                mat[i][j] = sum(w for s, w in weights.items()
                                if (a in s) != (b in s))
    return names, mat, weights


def test_agglomeration_finds_an_ordering_that_draws_every_split():
    # this is the property Neighbor-Net exists for: given distances that really
    # are a circular split system, the ordering it returns has to make every
    # one of those splits a single arc, or the split cannot be drawn at all
    rng = random.Random(4)
    for _ in range(20):
        n = rng.randint(6, 14)
        names, mat, weights = _circular_distances(n, rng.randint(n, 3 * n), rng)
        order = neighbornet_ordering(names, mat)
        assert sorted(order) == sorted(names)
        undrawable = [s for s in weights if not is_circular(s, order)]
        assert undrawable == []


def test_the_agglomerative_ordering_beats_the_tree_leaf_order():
    # a tree's leaf order can only respect the tree's own splits, so it loses
    # the conflicting ones -- which are exactly the splits worth drawing
    rng = random.Random(4)
    agglomerative = leaf_order = 0.0
    for _ in range(12):
        n = rng.randint(8, 14)
        names, mat, weights = _circular_distances(n, 3 * n, rng)
        total = sum(weights.values())

        def drawable(order, weights=weights, total=total):
            return sum(w for s, w in weights.items()
                       if is_circular(s, order)) / total

        agglomerative += drawable(neighbornet_ordering(names, mat))
        tree = pt.neighbor_joining(names, mat)
        leaf_order += drawable(circular_ordering(
            names, splits_from_tree(tree, names, trivial=True)))
    assert agglomerative / 12 == pytest.approx(1.0, abs=1e-9)
    assert leaf_order / 12 < 0.95


def test_neighbor_net_is_reachable_by_name():
    names, mat, _ = _circular_distances(9, 20, random.Random(1))
    net = pt.neighbor_net(names, mat)
    assert net.estimated and net.dropped == []
    assert _crossings(net) == 0
    with pytest.raises(ValueError, match="neighbornet"):
        pt.SplitNetwork.from_distances(names, mat, ordering="whatever")


def test_split_network_renders(tmp_path):
    tr = pt.datasets.random_tree(10, seed=1)
    trees = [tr, pt.datasets.random_tree(10, seed=2)]
    for t in trees[1:]:
        for leaf, nm in zip(t.leaves(), tr.leaf_names()):
            leaf.name = nm
    out = tmp_path / "net.png"
    pt.SplitNetwork.from_trees(trees).titled("splits").save(str(out))
    assert out.exists() and out.stat().st_size > 1000
