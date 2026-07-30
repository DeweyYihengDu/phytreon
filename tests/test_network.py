"""Sequence-similarity networks (the CLANS-style cluster map)."""
import math

import matplotlib
matplotlib.use("Agg")

import pytest

import phytreon as pt
from phytreon.plot.network import fruchterman_reingold


def _three_families(seed=1):
    """Three tight families plus two weak bridges -- a miniature of the
    picture CLANS is used to produce."""
    import random
    rng = random.Random(seed)
    names, groups = [], {}
    for fam, n in (("A", 8), ("B", 8), ("C", 8)):
        for i in range(n):
            nm = f"{fam}{i}"
            names.append(nm)
            groups[nm] = fam
    pairs = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if groups[a] == groups[b]:
                pairs.append((a, b, rng.uniform(0.6, 0.9)))
    pairs.append(("A0", "B0", 0.42))          # the only inter-family links
    pairs.append(("B0", "C0", 0.42))
    return names, groups, pairs


# --------------------------------------------------------------------------
# the layout
# --------------------------------------------------------------------------
def test_force_layout_returns_one_finite_point_per_node():
    names = [f"n{i}" for i in range(12)]
    edges = [(i, i + 1, 1.0) for i in range(11)]
    pos = fruchterman_reingold(names, edges, iterations=60, seed=0)
    assert len(pos) == len(names)
    assert all(math.isfinite(x) and math.isfinite(y) for x, y in pos)
    # nodes must not all collapse onto one another
    assert len({(round(x, 6), round(y, 6)) for x, y in pos}) == len(names)


def test_force_layout_handles_degenerate_inputs():
    assert fruchterman_reingold([], []) == []
    assert fruchterman_reingold(["only"], []) == [(0.0, 0.0)]
    # no edges at all: pure repulsion, still finite and distinct
    pos = fruchterman_reingold(["a", "b", "c"], [], iterations=30)
    assert all(math.isfinite(x) and math.isfinite(y) for x, y in pos)


def test_force_layout_is_deterministic_for_a_seed():
    names = [f"n{i}" for i in range(10)]
    edges = [(0, i, 1.0) for i in range(1, 10)]
    a = fruchterman_reingold(names, edges, iterations=50, seed=7)
    b = fruchterman_reingold(names, edges, iterations=50, seed=7)
    assert a == b
    assert fruchterman_reingold(names, edges, iterations=50, seed=8) != a


def test_connected_nodes_end_up_closer_than_unconnected_ones():
    # the whole point of the layout: similarity pulls, so a linked pair should
    # sit nearer than a pair with nothing between them
    names, groups, pairs = _three_families()
    net = pt.SequenceNetwork.from_pairs(pairs, names=names, seed=2)
    pos = dict(zip(net.names, net.positions))

    def dist(a, b):
        return math.hypot(pos[a][0] - pos[b][0], pos[a][1] - pos[b][1])

    within = [dist(f"A{i}", f"A{j}") for i in range(8) for j in range(i + 1, 8)]
    between = [dist(f"A{i}", f"C{j}") for i in range(8) for j in range(8)]
    assert max(within) < min(between)


# --------------------------------------------------------------------------
# constructors
# --------------------------------------------------------------------------
def test_from_distances_keeps_only_edges_above_the_cutoff():
    names = ["a", "b", "c"]
    #      a     b     c
    mat = [[0.0, 0.10, 0.90],      # a-b similar (sim 0.9), a-c not (sim 0.1)
           [0.10, 0.0, 0.95],
           [0.90, 0.95, 0.0]]
    net = pt.SequenceNetwork.from_distances(names, mat, cutoff=0.5)
    assert len(net.edges) == 1
    i, j, w = net.edges[0]
    assert {names[i], names[j]} == {"a", "b"}
    assert w == pytest.approx(0.9)


def test_from_pairs_infers_names_and_rejects_unknown_ones():
    net = pt.SequenceNetwork.from_pairs([("x", "y", 0.8), ("y", "z", 0.7)])
    assert sorted(net.names) == ["x", "y", "z"]
    assert len(net.edges) == 2
    with pytest.raises(ValueError, match="not in `names`"):
        pt.SequenceNetwork.from_pairs([("x", "nope", 0.8)], names=["x", "y"])


def test_from_pairs_defaults_missing_weights_to_one():
    net = pt.SequenceNetwork.from_pairs([("x", "y")])
    assert net.edges[0][2] == pytest.approx(1.0)


def test_from_alignment_builds_edges_from_pairwise_identity():
    aln = pt.Alignment(
        names=["s1", "s2", "s3"],
        seqs=["AAAACCCC",
              "AAAACCCG",     # nearly identical to s1
              "TTTTGGGG"],    # shares nothing
    )
    net = pt.SequenceNetwork.from_alignment(aln, cutoff=0.5)
    assert len(net.edges) == 1
    i, j, _ = net.edges[0]
    assert {net.names[i], net.names[j]} == {"s1", "s2"}


# --------------------------------------------------------------------------
# clusters and rendering
# --------------------------------------------------------------------------
def test_components_finds_the_real_clusters():
    # drop the bridges and the three families become three components
    names, groups, pairs = _three_families()
    pairs = [p for p in pairs if groups[p[0]] == groups[p[1]]]
    net = pt.SequenceNetwork.from_pairs(pairs, names=names)
    comps = net.components()
    assert len(comps) == 3
    assert sorted(len(c) for c in comps) == [8, 8, 8]
    for comp in comps:                       # each is one family, not a mix
        assert len({groups[n] for n in comp}) == 1


def test_components_lists_an_isolated_node_on_its_own():
    net = pt.SequenceNetwork.from_pairs([("a", "b", 0.9)], names=["a", "b", "lonely"])
    comps = net.components()
    assert ["lonely"] in comps


def test_scene_has_one_marker_per_sequence_and_one_path_per_edge():
    names, groups, pairs = _three_families()
    net = pt.SequenceNetwork.from_pairs(pairs, names=names)
    scene = net._build().scene
    assert len(scene.markers) == len(names)
    assert len(scene.paths) == len(net.edges)


def test_color_by_recolours_nodes_and_adds_one_legend():
    names, groups, pairs = _three_families()
    net = pt.SequenceNetwork.from_pairs(pairs, names=names).color_by(groups,
                                                                    title="family")
    ctx = net._build()
    assert [t for t, _ in ctx.scene.legends] == ["family"]
    assert len({m.color for m in ctx.scene.markers}) == 3


def test_baseline_greys_out_a_group_in_the_network_too():
    from phytreon.plot.figure import BASELINE_GREY
    names, groups, pairs = _three_families()
    net = (pt.SequenceNetwork.from_pairs(pairs, names=names)
           .color_by(groups, title="family", baseline="A"))
    ctx = net._build()
    a_nodes = [m for m, nm in zip(ctx.scene.markers, net.names)
               if groups[nm] == "A"]
    assert all(m.color == BASELINE_GREY for m in a_nodes)


def test_cluster_labels_sit_outside_their_own_cluster():
    names, groups, pairs = _three_families()
    net = (pt.SequenceNetwork.from_pairs(pairs, names=names, seed=2)
           .label_clusters({fam: [n for n in names if groups[n] == fam]
                            for fam in ("A", "B", "C")}))
    ctx = net._build()
    pos = dict(zip(net.names, net.positions))
    labels = {lb.text: (lb.x, lb.y) for lb in ctx.scene.labels}
    assert set(labels) == {"A", "B", "C"}
    for fam, (lx, ly) in labels.items():
        members = [pos[n] for n in names if groups[n] == fam]
        cx = sum(p[0] for p in members) / len(members)
        cy = sum(p[1] for p in members) / len(members)
        reach = max(math.hypot(p[0] - cx, p[1] - cy) for p in members)
        assert math.hypot(lx - cx, ly - cy) > reach     # clear of the nodes


def test_stronger_edges_are_drawn_more_firmly():
    net = pt.SequenceNetwork.from_pairs(
        [("a", "b", 0.95), ("c", "d", 0.15)], names=["a", "b", "c", "d"])
    paths = net._build().scene.paths
    strong = max(paths, key=lambda p: p.width)
    weak = min(paths, key=lambda p: p.width)
    assert strong.width > weak.width
    assert strong.opacity > weak.opacity


def test_network_renders_to_file(tmp_path):
    names, groups, pairs = _three_families()
    net = (pt.SequenceNetwork.from_pairs(pairs, names=names)
           .color_by(groups)
           .titled("cluster map"))
    for ext in ("png", "svg"):
        out = tmp_path / f"net.{ext}"
        net.save(str(out))
        assert out.exists() and out.stat().st_size > 1000
