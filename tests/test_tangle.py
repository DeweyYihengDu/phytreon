import matplotlib
matplotlib.use("Agg")

import pytest

import phytreon as pt
from phytreon.treeops import _inversions


# --------------------------------------------------------------------------
# crossing counting / untangling
# --------------------------------------------------------------------------
def test_inversions_hand_cases():
    assert _inversions([]) == 0
    assert _inversions([0, 1, 2, 3]) == 0          # sorted -> nothing crosses
    assert _inversions([0, 2, 1, 3]) == 1          # one adjacent swap
    assert _inversions([3, 2, 1, 0]) == 6          # reversed -> every pair


def test_crossing_number_counts_swapped_tips():
    a = pt.Tree.from_newick("(((A,B),C),D);")
    b = pt.Tree.from_newick("(((B,A),C),D);")      # A/B swapped -> 1 crossing
    assert pt.crossing_number(a, b) == 1
    assert pt.crossing_number(a, a) == 0


def test_untangle_resolves_a_rotation_only_difference():
    # the two trees are the same topology drawn with clades rotated; a
    # tanglegram of them should untangle to zero crossings
    a = pt.Tree.from_newick("(((A,B),(C,D)),(E,F));")
    b = pt.Tree.from_newick("(((B,A),(D,C)),(F,E));")
    assert pt.crossing_number(a, b) > 0
    assert pt.untangle(a, b) == 0
    assert pt.crossing_number(a, b) == 0


def test_untangle_preserves_topology_and_branch_lengths():
    # rotation reorders children only -- the tree itself must not change
    a = pt.datasets.primates()
    b = pt.datasets.random_tree(12, seed=7)
    before_rf = pt.robinson_foulds(b, b)
    total_before = b.total_branch_length
    names_before = set(b.leaf_names())
    pt.untangle(a, b, fix="left")
    assert set(b.leaf_names()) == names_before     # same taxa
    assert b.total_branch_length == pytest.approx(total_before)
    assert pt.robinson_foulds(b, b) == before_rf   # still a valid tree


def test_untangle_never_increases_crossings():
    a = pt.datasets.random_tree(14, seed=1)
    b = pt.datasets.random_tree(14, seed=2)
    # give both trees the same taxa so every tip links
    b = pt.Tree.from_newick(b.write())
    for leaf, name in zip(b.leaves(), a.leaf_names()):
        leaf.name = name
    before = pt.crossing_number(a, b)
    after = pt.untangle(a, b, fix="left")
    assert after <= before
    assert after == pt.crossing_number(a, b)       # reported == actual


def test_untangle_rejects_bad_fix():
    a = pt.datasets.primates()
    b = pt.datasets.primates()
    with pytest.raises(ValueError):
        pt.untangle(a, b, fix="middle")


def test_crossing_number_ignores_unshared_tips():
    a = pt.Tree.from_newick("((A,B),(C,D));")
    b = pt.Tree.from_newick("((A,B),(C,Z));")      # D vs Z not shared
    assert pt.crossing_number(a, b) == 0           # A,B,C agree


# --------------------------------------------------------------------------
# the figure
# --------------------------------------------------------------------------
def _links(ctx):
    """Link primitives are the only paths at the reserved link zorder."""
    return [p for p in ctx.scene.paths if p.zorder == 0.4]


def test_tanglegram_links_every_shared_tip(tmp_path):
    a = pt.datasets.primates()
    b = pt.Tree.from_newick(
        "(((Human,Chimp),(Gorilla,Gibbon)),(Orangutan,(Baboon,Macaque)));")
    fig = pt.TangleFigure(a, b, titles=("genome", "transcriptome"))
    ctx = fig._build()
    assert len(_links(ctx)) == 7                   # all seven taxa are shared
    out = tmp_path / "tangle.png"
    fig.save(str(out))
    assert out.exists() and out.stat().st_size > 1000


def test_tanglegram_only_links_shared_tips():
    a = pt.datasets.primates()                     # 7 taxa
    b = pt.Tree.from_newick("(((Human,Chimp),Gorilla),(Gibbon,Marmoset));")
    fig = pt.TangleFigure(a, b)
    # Human/Chimp/Gorilla/Gibbon are shared; Marmoset is not in the left tree
    assert len(_links(fig._build())) == 4


def test_right_tree_is_mirrored_and_labels_flip():
    a = pt.datasets.primates()
    b = pt.datasets.primates()
    fig = pt.TangleFigure(a, b, tip_labels="both")
    ctx = fig._build()
    tips = [lb for lb in ctx.scene.labels if lb.role == "tiplab"]
    # left tree's labels read outward to the right, the mirrored tree's to
    # the left -- otherwise the right tree's names would run back over it
    assert {lb.ha for lb in tips} == {"left", "right"}
    lefts = [lb.x for lb in tips if lb.ha == "left"]
    rights = [lb.x for lb in tips if lb.ha == "right"]
    assert max(lefts) < min(rights)                # the two label bands are disjoint


def test_tip_labels_default_names_both_sides():
    a = pt.datasets.primates()
    b = pt.datasets.primates()
    ctx = pt.TangleFigure(a, b)._build()           # default tip_labels="both"
    tips = [lb for lb in ctx.scene.labels if lb.role == "tiplab"]
    assert len(tips) == 2 * a.n_leaves             # every taxon named on each side
    assert len(set(lb.text for lb in tips)) == a.n_leaves


def test_tip_labels_can_name_one_side_only():
    a = pt.datasets.primates()
    ctx = pt.TangleFigure(a, pt.datasets.primates(), tip_labels="left")._build()
    tips = [lb for lb in ctx.scene.labels if lb.role == "tiplab"]
    assert len(tips) == a.n_leaves
    assert {lb.ha for lb in tips} == {"left"}      # only the unmirrored side


def test_middle_band_widens_for_longer_taxon_names():
    # the band is sized from the labels, so long names must not overrun the
    # links; a tree with long names gets a proportionally wider band
    short = pt.Tree.from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);")
    long_ = pt.Tree.from_newick(
        "((Methanocaldococcus_jannaschii:0.1,Bradyrhizobium_japonicum:0.1):0.1,"
        "(Verrucomicrobium_spinosum:0.1,Parachlamydia_acanthamoebae:0.1):0.1);")

    def band(tree):
        fig = pt.TangleFigure(tree, pt.Tree.from_newick(tree.write()))
        ctx = fig._build()
        depth = fig.left._build().layout.max_x
        return ctx.layout.max_x - 2 * depth        # total - both trees
    assert band(long_) > 2 * band(short)


def test_tip_labels_can_be_switched_off():
    a = pt.datasets.primates()
    ctx = pt.TangleFigure(a, pt.datasets.primates(), tip_labels=False)._build()
    assert not [lb for lb in ctx.scene.labels if lb.role == "tiplab"]
    with pytest.raises(ValueError):
        pt.TangleFigure(a, pt.datasets.primates(), tip_labels="sideways")


def test_highlight_discordant_marks_exactly_the_crossing_links():
    from phytreon.plot.tangle import DISCORDANT_COLOR
    a = pt.Tree.from_newick("(((A,B),(C,D)),E);")
    b = pt.Tree.from_newick("(((A,B),(D,C)),E);")   # only C/D disagree
    fig = pt.TangleFigure(a, b).connect(highlight_discordant=True)
    ctx = fig._build()
    flagged = [p for p in _links(ctx) if p.color == DISCORDANT_COLOR]
    assert len(flagged) == 2                       # the C and D links, nothing else


def test_link_colour_from_a_data_column_adds_one_legend():
    a = pt.datasets.primates()
    a.join_data(pt.datasets.primates_metadata().reset_index(), on="name")
    fig = pt.TangleFigure(a, pt.datasets.primates()).connect(color="habitat")
    ctx = fig._build()
    assert [t for t, _ in ctx.scene.legends] == ["habitat"]
    assert len({p.color for p in _links(ctx)}) > 1  # links really are recoloured


def test_untangle_method_chains_and_reduces_crossings():
    a = pt.Tree.from_newick("(((A,B),(C,D)),(E,F));")
    b = pt.Tree.from_newick("(((B,A),(D,C)),(F,E));")
    fig = pt.TangleFigure(a, b)
    assert fig.crossings() > 0
    assert fig.untangle() is fig                   # chains
    assert fig.crossings() == 0


def test_tanglegram_rejects_layouts_it_cannot_mirror():
    # reflecting a circular tree turns it inside out rather than facing it,
    # so this has to fail loudly instead of drawing nonsense
    a = pt.datasets.primates()
    fig = pt.TangleFigure(a, pt.datasets.primates(), layout="circular")
    with pytest.raises(ValueError, match="depth runs along x"):
        fig._build()
    # slanted keeps depth along x, so it stays allowed
    pt.TangleFigure(a, pt.datasets.primates(), layout="slanted")._build()


def test_tanglegram_renders_both_backends(tmp_path):
    a = pt.datasets.primates()
    b = pt.datasets.primates()
    fig = pt.TangleFigure(a, b).untangle()
    fig.left.tip_points(color="black", size=5)
    for ext in ("png", "svg", "html"):
        out = tmp_path / f"tangle.{ext}"
        fig.save(str(out))
        assert out.exists() and out.stat().st_size > 500
