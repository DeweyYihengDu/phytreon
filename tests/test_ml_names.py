"""External engines get engine-safe tip names; the tree comes back with the
real ones. Verified with a fake engine, so no binary needs to be installed."""
import os
import re
from types import SimpleNamespace

import phytreon as pt
from phytreon.infer import ml as ml_mod
from phytreon.infer.align import Alignment, read_fasta
from phytreon.infer.constraint import constraint_tree

NAMES = ["sp|P1|HBA_HUMAN", "E. coli K-12", "Clade(A):x", "plain_name"]
ALN = Alignment(NAMES, ["ACGT", "ACGA", "TTGT", "TTGA"])


def _fake_iqtree(monkeypatch, seen):
    monkeypatch.setattr(ml_mod.shutil, "which",
                        lambda tool: None if tool == "iqtree2" else "/fake/iqtree")

    def fake_run(cmd, **kw):
        infile = cmd[2]
        names = [n for n, _ in read_fasta(infile)]
        seen["names"] = names
        if "-g" in cmd:
            with open(cmd[cmd.index("-g") + 1]) as f:
                seen["constraint"] = f.read()
        # like IQ-TREE: anything outside [A-Za-z0-9_.-] becomes '_'
        safe = [re.sub(r"[^A-Za-z0-9_.\-]", "_", n) for n in names]
        with open(infile + ".treefile", "w") as f:
            f.write(f"(({safe[0]},{safe[1]}),{safe[2]},{safe[3]});")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(ml_mod.subprocess, "run", fake_run)


def test_unsafe_names_are_hidden_from_the_engine_and_restored(monkeypatch):
    seen = {}
    _fake_iqtree(monkeypatch, seen)
    tree = ml_mod.infer_iqtree(ALN)
    assert all(re.fullmatch(r"[A-Za-z0-9_.\-]+", n) for n in seen["names"])
    assert set(tree.leaf_names()) == set(NAMES)
    assert set(tree.get_mrca(NAMES[:2]).leaf_names()) == set(NAMES[:2])


def test_safe_names_are_passed_through_untouched(monkeypatch):
    seen = {}
    _fake_iqtree(monkeypatch, seen)
    aln = Alignment(["a1", "a2", "b1", "b2"], ALN.seqs)
    tree = ml_mod.infer_iqtree(aln)
    assert seen["names"] == ["a1", "a2", "b1", "b2"]
    assert set(tree.leaf_names()) == {"a1", "a2", "b1", "b2"}


def test_constraint_tree_is_renamed_to_match(monkeypatch):
    seen = {}
    _fake_iqtree(monkeypatch, seen)
    ct = constraint_tree({NAMES[0]: "A", NAMES[1]: "A",
                          NAMES[2]: "B", NAMES[3]: "B"})
    before = sorted(ct.leaf_names())
    ml_mod.infer_iqtree(ALN, constraint=ct)
    sent = pt.Tree.from_newick(seen["constraint"])
    assert sorted(sent.leaf_names()) == sorted(seen["names"])
    assert sorted(ct.leaf_names()) == before        # caller's tree untouched


def test_constraint_file_path_is_renamed_to_match(monkeypatch, tmp_path):
    seen = {}
    _fake_iqtree(monkeypatch, seen)
    gfile = tmp_path / "backbone.tre"
    gfile.write_text(pt.core.io.to_newick(
        constraint_tree({NAMES[0]: "A", NAMES[1]: "A"})))
    ml_mod.infer_iqtree(ALN, constraint=str(gfile))
    sent = pt.Tree.from_newick(seen["constraint"])
    assert set(sent.leaf_names()) <= set(seen["names"])
    assert os.path.exists(gfile)
