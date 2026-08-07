"""Taxonomy-constrained tree building: constraint_tree() for a constrained
ML search, and constrained_nj() which forces the grouping outright."""
import os
import shutil
from types import SimpleNamespace

import pytest

import phytreon as pt
from phytreon.infer.constraint import constraint_tree
from phytreon.infer.distance import constrained_nj
from phytreon.infer.align import Alignment


# --------------------------------------------------------------------------
# constraint_tree: the polytomy file for IQ-TREE's -g / RAxML-NG's
# --tree-constraint
# --------------------------------------------------------------------------
def test_constraint_tree_groups_present_tips_into_polytomies():
    groups = {"a1": "A", "a2": "A", "a3": "A", "b1": "B", "b2": "B"}
    ct = constraint_tree(groups)
    a_clade = ct.get_mrca(["a1", "a2", "a3"])
    assert set(a_clade.leaf_names()) == {"a1", "a2", "a3"}
    b_clade = ct.get_mrca(["b1", "b2"])
    assert set(b_clade.leaf_names()) == {"b1", "b2"}
    assert set(ct.leaf_names()) == set(groups)


def test_constraint_tree_leaves_out_missing_and_none_labels():
    groups = {"a1": "A", "a2": "A", "free1": None}      # "free2" simply absent
    ct = constraint_tree(groups)
    assert set(ct.leaf_names()) == {"a1", "a2"}


def test_constraint_tree_singleton_group_is_a_bare_leaf_not_a_polytomy():
    groups = {"a1": "A", "a2": "A", "solo": "S"}
    ct = constraint_tree(groups)
    solo = next(n for n in ct.leaves() if n.name == "solo")
    assert solo.parent is ct.root          # sits directly under the root


def test_constraint_tree_rejects_a_grouping_with_nothing_left():
    with pytest.raises(ValueError, match="nothing"):
        constraint_tree({"a1": None, "a2": None})


def test_constraint_tree_reads_a_table_and_column():
    import pandas as pd
    df = pd.DataFrame({"name": ["a1", "a2", "b1"], "genus": ["A", "A", "B"]})
    ct = constraint_tree(df, column="genus")
    assert set(ct.get_mrca(["a1", "a2"]).leaf_names()) == {"a1", "a2"}

    with pytest.raises(TypeError, match="column"):
        constraint_tree(df)


def test_constraint_tree_writes_and_reads_back_as_plain_newick():
    ct = constraint_tree({"a1": "A", "a2": "A", "b1": "B"})
    nwk = ct.write()
    assert ":" not in nwk                   # topology only, no branch lengths
    back = pt.Tree.from_newick(nwk)
    assert set(back.get_mrca(["a1", "a2"]).leaf_names()) == {"a1", "a2"}


# --------------------------------------------------------------------------
# constrained_nj: forces the grouping outright, no search involved
# --------------------------------------------------------------------------
def _matrix(names, d):
    """``d`` maps a frozenset({a, b}) -> distance; diagonal is 0."""
    return [[0.0 if i == j else d[frozenset((names[i], names[j]))]
            for j in range(len(names))] for i in range(len(names))]


def test_constrained_nj_forces_a_group_the_data_alone_would_split():
    # b1 is genuinely closer to the A's than to b2 -- plain NJ would not put
    # the B's together -- but the constraint must win anyway
    names = ["a1", "a2", "b1", "b2"]
    d = {frozenset(("a1", "a2")): 1.0, frozenset(("a1", "b1")): 1.0,
        frozenset(("a1", "b2")): 9.0, frozenset(("a2", "b1")): 1.0,
        frozenset(("a2", "b2")): 9.0, frozenset(("b1", "b2")): 8.0}
    D = _matrix(names, d)

    plain = pt.neighbor_joining(names, D)
    assert set(plain.get_mrca(["b1", "b2"], strict=False).leaf_names()) != {"b1", "b2"}

    forced = constrained_nj(names, D, {"a1": "A", "a2": "A", "b1": "B", "b2": "B"})
    assert set(forced.leaf_names()) == set(names)
    assert set(forced.get_mrca(["b1", "b2"]).leaf_names()) == {"b1", "b2"}
    assert set(forced.get_mrca(["a1", "a2"]).leaf_names()) == {"a1", "a2"}


def test_constrained_nj_leaves_unlabelled_tips_free():
    names = ["a1", "a2", "a3", "x"]
    d = {frozenset(("a1", "a2")): 1.0, frozenset(("a1", "a3")): 1.0,
        frozenset(("a2", "a3")): 1.0, frozenset(("a1", "x")): 5.0,
        frozenset(("a2", "x")): 5.0, frozenset(("a3", "x")): 5.0}
    D = _matrix(names, d)
    tree = constrained_nj(names, D, {"a1": "A", "a2": "A", "a3": "A"})  # x absent
    assert set(tree.leaf_names()) == set(names)
    assert set(tree.get_mrca(["a1", "a2", "a3"]).leaf_names()) == {"a1", "a2", "a3"}


def test_constrained_nj_handles_exactly_two_groups():
    # NJ itself needs >= 3 taxa; the group-level backbone has only 2 "taxa"
    # (the two groups) and must not go through neighbor_joining at all
    names = ["a1", "a2", "b1", "b2", "b3"]
    d = {}
    for i, x in enumerate(names):
        for y in names[i + 1:]:
            same = x[0] == y[0]
            d[frozenset((x, y))] = 1.0 if same else 6.0
    D = _matrix(names, d)
    tree = constrained_nj(names, D, {n: n[0].upper() for n in names})
    assert set(tree.leaf_names()) == set(names)
    assert set(tree.get_mrca(["a1", "a2"]).leaf_names()) == {"a1", "a2"}
    assert set(tree.get_mrca(["b1", "b2", "b3"]).leaf_names()) == {"b1", "b2", "b3"}


def test_constrained_nj_with_one_group_is_plain_nj():
    names = ["a1", "a2", "a3", "a4"]
    d = {frozenset((x, y)): abs(i - j) + 1.0
        for i, x in enumerate(names) for j, y in enumerate(names) if i < j}
    D = _matrix(names, d)
    plain = pt.neighbor_joining(names, D)
    same = constrained_nj(names, D, {n: "ONLY" for n in names})
    assert pt.robinson_foulds(plain, same) == 0


def test_constrained_nj_two_tip_group_is_a_simple_cherry():
    names = ["a1", "a2", "b1", "b2", "b3"]
    d = {}
    for i, x in enumerate(names):
        for y in names[i + 1:]:
            d[frozenset((x, y))] = 2.0 if x[0] == y[0] else 7.0
    D = _matrix(names, d)
    tree = constrained_nj(names, D, {"a1": "A", "a2": "A",
                                     "b1": "B", "b2": "B", "b3": "B"})
    a_clade = tree.get_mrca(["a1", "a2"])
    assert set(a_clade.leaf_names()) == {"a1", "a2"}
    assert len(a_clade.get_leaves()) == 2
    kid_lengths = sorted(c.length for c in a_clade.children)
    assert kid_lengths == pytest.approx([1.0, 1.0])   # d/2 each, d = 2.0


# --------------------------------------------------------------------------
# build_tree(..., constraint=...) end to end
# --------------------------------------------------------------------------
SEQS_TWO_GENERA = [
    ("Aa1", "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"),
    ("Aa2", "ATGGCCATTGTTATGGGCCGCTGAAAGGGTGCCCGATAG"),
    ("Aa3", "ATGGCCATTGTAATGGGCCGCTGTAAGGGTGCCGATAG"),
    ("Bb1", "ATGTCGATTCTAATGAACCGCTGAAAGCGTGACCTTTAG"),
    ("Bb2", "ATGTCGATTCTAATGAACCGCTGTAAGCGTGACCTTTAG"),
    ("Bb3", "ATGTCGATTCTAATGAACCGGCTGAAAGCGTGACCTTTAG"),
]
GENUS = {"Aa1": "Aus", "Aa2": "Aus", "Aa3": "Aus",
        "Bb1": "Bus", "Bb2": "Bus", "Bb3": "Bus"}


def test_build_tree_constrained_nj_forces_the_grouping():
    tree = pt.build_tree(SEQS_TWO_GENERA, method="nj", constraint=GENUS)
    assert set(tree.get_mrca(["Aa1", "Aa2", "Aa3"]).leaf_names()) == \
        {"Aa1", "Aa2", "Aa3"}
    assert set(tree.get_mrca(["Bb1", "Bb2", "Bb3"]).leaf_names()) == \
        {"Bb1", "Bb2", "Bb3"}


def test_build_tree_rejects_constraint_on_native_ml():
    with pytest.raises(ValueError, match="native"):
        pt.build_tree(SEQS_TWO_GENERA, method="ml", ml_engine="native",
                      constraint=GENUS)


def test_build_tree_rejects_constraint_on_upgma():
    with pytest.raises(ValueError, match="UPGMA"):
        pt.build_tree(SEQS_TWO_GENERA, method="upgma", constraint=GENUS)


def test_build_tree_rejects_constraint_on_parsimony():
    with pytest.raises(ValueError, match="parsimony"):
        pt.build_tree(SEQS_TWO_GENERA, method="parsimony", constraint=GENUS)


def test_build_tree_nj_constraint_must_be_a_mapping_not_a_tree():
    ct = constraint_tree(GENUS)
    with pytest.raises(TypeError, match="mapping"):
        pt.build_tree(SEQS_TWO_GENERA, method="nj", constraint=ct)


# --------------------------------------------------------------------------
# IQ-TREE wiring: the -g flag, verified without needing the binary installed
# --------------------------------------------------------------------------
def test_infer_iqtree_passes_the_constraint_as_dash_g(monkeypatch, tmp_path):
    from phytreon.infer import ml as ml_mod

    monkeypatch.setattr(ml_mod.shutil, "which",
                        lambda tool: None if tool == "iqtree2" else "/fake/iqtree")

    calls = []
    seen_constraint = {}

    def fake_run(cmd, **kw):
        calls.append(cmd)
        # the constraint file only exists inside infer_iqtree's own
        # TemporaryDirectory, so it has to be read here, before that
        # context manager tears it down on return
        if "-g" in cmd:
            gfile = cmd[cmd.index("-g") + 1]
            seen_constraint["exists"] = os.path.exists(gfile)
            seen_constraint["tree"] = pt.Tree.read(gfile)
        with open(cmd[2] + ".treefile", "w") as f:
            f.write("(a1,a2,(b1,b2));")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(ml_mod.subprocess, "run", fake_run)

    aln = Alignment(["a1", "a2", "b1", "b2"],
                    ["ACGT", "ACGA", "TTGT", "TTGA"])
    ct = constraint_tree({"a1": "A", "a2": "A", "b1": "B", "b2": "B"})
    tree = ml_mod.infer_iqtree(aln, constraint=ct)

    assert set(tree.leaf_names()) == {"a1", "a2", "b1", "b2"}
    assert "-g" in calls[0]
    assert seen_constraint["exists"]
    assert set(seen_constraint["tree"].get_mrca(["a1", "a2"]).leaf_names()) \
        == {"a1", "a2"}


def test_infer_iqtree_accepts_a_constraint_file_path_directly(monkeypatch, tmp_path):
    from phytreon.infer import ml as ml_mod

    monkeypatch.setattr(ml_mod.shutil, "which",
                        lambda tool: None if tool == "iqtree2" else "/fake/iqtree")
    calls = []

    def fake_run(cmd, **kw):
        calls.append(cmd)
        with open(cmd[2] + ".treefile", "w") as f:
            f.write("(a1,a2,b1);")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(ml_mod.subprocess, "run", fake_run)

    gfile = tmp_path / "backbone.tre"
    gfile.write_text("(a1,a2);")
    aln = Alignment(["a1", "a2", "b1"], ["ACGT", "ACGA", "TTGT"])
    ml_mod.infer_iqtree(aln, constraint=str(gfile))
    cmd = calls[0]
    assert cmd[cmd.index("-g") + 1] == str(gfile)


@pytest.mark.skipif(not (shutil.which("iqtree2") or shutil.which("iqtree")),
                    reason="iqtree not installed")
def test_build_tree_constrained_ml_via_real_iqtree():
    tree = pt.build_tree(SEQS_TWO_GENERA, method="ml", ml_engine="iqtree",
                         constraint=GENUS)
    assert set(tree.get_mrca(["Aa1", "Aa2", "Aa3"]).leaf_names()) == \
        {"Aa1", "Aa2", "Aa3"}
    assert set(tree.get_mrca(["Bb1", "Bb2", "Bb3"]).leaf_names()) == \
        {"Bb1", "Bb2", "Bb3"}
