"""Annotated NEXUS (BEAST / MrBayes) input."""
import matplotlib
matplotlib.use("Agg")

import pytest

import phytreon as pt
from phytreon.core.nexus import parse_annotation

BEAST = """#NEXUS
BEGIN TREES;
\tTRANSLATE
\t\t1 Homo_sapiens,
\t\t2 Pan_troglodytes,
\t\t3 Gorilla_gorilla
\t\t;
\tTREE tree_1 = [&R] ((1[&height=0.0]:6.5,2[&height=0.0]:6.5)\
[&height=6.5,height_95%_HPD={5.1,8.2},posterior=1.0]:2.4[&rate=0.91],\
3[&height=0.0]:8.9)[&height=8.9,height_95%_HPD={7.0,11.3},posterior=0.97];
END;
"""


@pytest.fixture()
def beast_file(tmp_path):
    path = tmp_path / "summary.tre"
    path.write_text(BEAST, encoding="utf-8")
    return str(path)


# --------------------------------------------------------------------------
# the annotation payload
# --------------------------------------------------------------------------
def test_parse_annotation_reads_numbers_strings_and_intervals():
    got = parse_annotation("&height=6.5,label=\"clade A\",range={1.5,3.5}")
    assert got["height"] == pytest.approx(6.5)
    assert got["label"] == "clade A"
    assert got["range"] == [1.5, 3.5]


def test_parse_annotation_flattens_hpd_to_lower_upper():
    got = parse_annotation("&height_95%_HPD={5.1,8.2}")
    assert got["height_95%_HPD"] == [5.1, 8.2]          # raw value kept
    assert got["height_95_lower"] == pytest.approx(5.1)  # what node_bars reads
    assert got["height_95_upper"] == pytest.approx(8.2)


def test_parse_annotation_orders_a_reversed_interval():
    got = parse_annotation("&age_95%_HPD={9.0,2.0}")
    assert got["age_95_lower"] == pytest.approx(2.0)
    assert got["age_95_upper"] == pytest.approx(9.0)


def test_parse_annotation_ignores_bare_flags_and_strips_prefixes():
    got = parse_annotation("&R")
    assert got == {}
    assert parse_annotation("&!color=#ff0000")["color"] == "#ff0000"


def test_commas_inside_braces_do_not_split_the_annotation():
    got = parse_annotation("&set={1.0,2.0,3.0},posterior=0.9")
    assert got["set"] == [1.0, 2.0, 3.0]
    assert got["posterior"] == pytest.approx(0.9)


# --------------------------------------------------------------------------
# reading a whole file
# --------------------------------------------------------------------------
def test_reading_a_beast_tree_keeps_the_node_estimates(beast_file):
    tree = pt.Tree.read(beast_file, fmt="beast")
    assert tree.leaf_names() == ["Homo_sapiens", "Pan_troglodytes",
                                 "Gorilla_gorilla"]        # TRANSLATE applied
    assert tree.root.data["posterior"] == pytest.approx(0.97)
    assert tree.root.data["height_95_lower"] == pytest.approx(7.0)
    node = tree.get_mrca(["Homo_sapiens", "Pan_troglodytes"])
    assert node.data["height"] == pytest.approx(6.5)
    assert node.data["height_95_upper"] == pytest.approx(8.2)


def test_a_comment_after_the_branch_length_is_captured(beast_file):
    # BEAST writes per-branch rates there; the parser used to stop dead at the
    # opening bracket and mangle the rest of the tree
    tree = pt.Tree.read(beast_file, fmt="beast")
    node = tree.get_mrca(["Homo_sapiens", "Pan_troglodytes"])
    assert node.data["rate"] == pytest.approx(0.91)
    assert node.length == pytest.approx(2.4)               # length still read


def test_plain_nexus_reader_still_ignores_annotations(beast_file):
    # fmt="nexus" goes through Biopython and keeps only the topology; the
    # point of fmt="beast" is that it does not
    plain = pt.Tree.read(beast_file, fmt="nexus")
    assert plain.n_leaves == 3
    assert not any("posterior" in n.data for n in plain.traverse())


def test_node_bars_work_straight_off_a_beast_file(beast_file):
    # the whole point: the defaults line up with what the reader writes
    tree = pt.Tree.read(beast_file, fmt="beast")
    ctx = pt.TreeFigure(tree).node_bars().time_axis()._build()
    bars = [p for p in ctx.scene.paths if p.width == 3.0]
    assert len(bars) == 2                                  # both internal nodes


def test_read_beast_helper_matches_tree_read(beast_file):
    assert (pt.read_beast(beast_file).leaf_names()
            == pt.Tree.read(beast_file, fmt="beast").leaf_names())


def test_helpful_errors_on_bad_input(tmp_path):
    junk = tmp_path / "x.tre"
    junk.write_text("not a nexus file", encoding="utf-8")
    with pytest.raises(ValueError, match="NEXUS"):
        pt.read_beast(str(junk))

    empty = tmp_path / "e.tre"
    empty.write_text("#NEXUS\nBEGIN TREES;\nEND;\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no TREE statement"):
        pt.read_beast(str(empty))


def test_tree_index_selects_from_a_posterior_sample(tmp_path):
    path = tmp_path / "post.tre"
    path.write_text(
        "#NEXUS\nBEGIN TREES;\n"
        "\tTREE g1 = [&R] ((A:1,B:1)[&posterior=0.5]:1,C:2);\n"
        "\tTREE g2 = [&R] ((A:1,C:1)[&posterior=0.9]:1,B:2);\n"
        "END;\n", encoding="utf-8")
    first = pt.read_beast(str(path))
    second = pt.read_beast(str(path), tree_index=1)
    assert first.get_mrca(["A", "B"]).data["posterior"] == pytest.approx(0.5)
    assert second.get_mrca(["A", "C"]).data["posterior"] == pytest.approx(0.9)
    with pytest.raises(IndexError, match="out of range"):
        pt.read_beast(str(path), tree_index=5)
