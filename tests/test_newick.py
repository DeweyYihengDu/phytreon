"""The built-in Newick parser and Node traversal on awkward input."""
import pytest

import phytreon as pt
from phytreon.core.io import parse_newick, to_newick


@pytest.mark.parametrize("text", [
    "(A, (B,C));",
    "(A,\n  (B,C)\n);",
    "( (B,C) ,A );",
])
def test_whitespace_between_tokens_adds_no_phantom_tip(text):
    tr = parse_newick(text)
    assert sorted(tr.leaf_names()) == ["A", "B", "C"]
    assert len(tr.root.children) == 2


@pytest.mark.parametrize("text", ["((A,B);", "(A,B", "((A,B),(C,D);"])
def test_unclosed_paren_raises_instead_of_hanging(text):
    with pytest.raises(ValueError, match="malformed Newick"):
        parse_newick(text)


def test_stray_close_paren_raises_instead_of_truncating():
    with pytest.raises(ValueError, match="unbalanced"):
        parse_newick("(A,B),C);")


def test_well_formed_round_trip_unchanged():
    text = "((A:0.1,B:0.2)90:0.3,'C d':0.4);"
    assert to_newick(parse_newick(text)) == text


def test_traversal_orders():
    tr = parse_newick("((A,B)x,(C,D)y)r;")
    assert [n.name for n in tr.traverse("preorder")] == ["r", "x", "A", "B", "y", "C", "D"]
    assert [n.name for n in tr.traverse("postorder")] == ["A", "B", "x", "C", "D", "y", "r"]
    assert [n.name for n in tr.traverse("levelorder")] == ["r", "x", "y", "A", "B", "C", "D"]


def test_traversal_of_deep_tree_does_not_hit_recursion_limit():
    # deeper than Python's recursion limit: out of reach of nested generators
    import sys
    depth = sys.getrecursionlimit() + 200
    root = node = pt.Node("r")
    for i in range(depth):
        node.add_child(pt.Node(f"t{i}"))
        node = node.add_child(pt.Node())
    tr = pt.Tree(root)
    assert len(tr.nodes("preorder")) == len(tr.nodes("postorder")) == 2 * depth + 1
    assert tr.n_leaves == depth + 1


def test_non_numeric_branch_length_raises():
    with pytest.raises(ValueError, match="branch length 'x'"):
        parse_newick("(A:x,B:1);")


def test_empty_branch_length_means_no_length():
    tr = parse_newick("(A:,B:1);")
    assert [n.length for n in tr.leaves()] == [None, 1.0]


def test_read_fasta_accepts_pathlib_path(tmp_path):
    path = tmp_path / "x.fasta"
    path.write_text(">a\nACGT\n>b\nAC\nGT\n")
    assert pt.read_fasta(path) == [("a", "ACGT"), ("b", "ACGT")]
    with pytest.raises(FileNotFoundError):
        pt.read_fasta(tmp_path / "missing.fasta")


def test_save_accepts_pathlib_path_and_closes_figure_on_failure(tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = pt.TreeFigure(pt.datasets.primates()).tip_labels()
    out = tmp_path / "tree.png"
    assert fig.save(out) == str(out) and out.exists()
    before = len(plt.get_fignums())
    with pytest.raises(Exception):
        fig.save(tmp_path / "no_such_dir" / "tree.png")
    assert len(plt.get_fignums()) == before
