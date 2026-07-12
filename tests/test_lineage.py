"""Single-cell CRISPR lineage-tracing tree reconstruction.

Purely additive alongside the existing (reversible) Fitch parsimony -- see
phytreon/infer/lineage.py's module docstring for the biological rationale.
Uses small, fast, offline synthetic tables; the real bundled demo dataset
(examples/data/lineage_alleletable.txt, examples/lineage_demo.py) is
exercised separately as a real-world validation script, not gated here --
there is no crisp pass/fail threshold for how closely a from-scratch NNI
parsimony search should match Cassiopeia's own specialized solver.
"""
from __future__ import annotations

import pandas as pd
import pytest

from phytreon.core.tree import Node, Tree
from phytreon.infer.align import Alignment
from phytreon.infer.parsimony import parsimony_score
from phytreon.infer.lineage import (
    read_allele_table, sankoff_score, camin_sokal_score, lineage_tree,
)


# --------------------------------------------------------------------------
# read_allele_table
# --------------------------------------------------------------------------
def test_read_allele_table_basic_recoding():
    df = pd.DataFrame([
        {"cellBC": "c1", "intBC": "i1", "r1": "AAAAA[10:2D]TTTTT", "r2": "AAAAA[None]TTTTT"},
        {"cellBC": "c2", "intBC": "i1", "r1": "AAAAA[10:2D]TTTTT", "r2": "AAAAA[20:1I]TTTTT"},
        {"cellBC": "c3", "intBC": "i1", "r1": "AAAAA[15:3D]TTTTT", "r2": "AAAAA[None]TTTTT"},
    ])
    aln = read_allele_table(df, site_cols=("r1", "r2"))
    assert sorted(aln.names) == ["c1", "c2", "c3"]
    by_name = dict(zip(aln.names, aln.seqs))
    # c1 and c2 share the same r1 edit -> same code; c3's r1 edit is distinct
    assert by_name["c1"][0] == by_name["c2"][0] != by_name["c3"][0]
    # c1 and c3's r2 is ancestral ("[None]") -> both code "0"
    assert by_name["c1"][1] == by_name["c3"][1] == "0"


def test_read_allele_table_saturated_site_keeps_ancestral_at_zero():
    # r1 is 100% edited (no cell shows "[None]") -- this occurs in the real
    # bundled dataset and, without the phantom-ancestral-row fix, would
    # make read_character_matrix recode some real edit as if ancestral.
    df = pd.DataFrame([
        {"cellBC": "c1", "intBC": "i1", "r1": "AAAAA[10:2D]TTTTT", "r2": "AAAAA[None]TTTTT"},
        {"cellBC": "c2", "intBC": "i1", "r1": "AAAAA[10:2D]TTTTT", "r2": "AAAAA[20:1I]TTTTT"},
        {"cellBC": "c3", "intBC": "i1", "r1": "AAAAA[15:3D]TTTTT", "r2": "AAAAA[None]TTTTT"},
    ])
    aln = read_allele_table(df, site_cols=("r1", "r2"))
    r1_codes = [s[0] for s in aln.seqs]
    assert "0" not in r1_codes


def test_read_allele_table_dropout_reindexed_not_omitted():
    # c2 has no row at all for intBC i2 (row-level dropout, distinct from an
    # in-row "NC" marker) -- must show up as missing, not vanish silently.
    df = pd.DataFrame([
        {"cellBC": "c1", "intBC": "i1", "r1": "AAAAA[10:2D]TTTTT"},
        {"cellBC": "c2", "intBC": "i1", "r1": "AAAAA[10:2D]TTTTT"},
        {"cellBC": "c1", "intBC": "i2", "r1": "AAAAA[None]TTTTT"},
    ])
    aln = read_allele_table(df, site_cols=("r1",))
    by_name = dict(zip(aln.names, aln.seqs))
    assert by_name["c2"][1] == "?"


def test_read_allele_table_nc_marker_is_missing():
    df = pd.DataFrame([
        {"cellBC": "c1", "intBC": "i1", "r1": "AAAAA[10:2D]TTTTT"},
        {"cellBC": "c2", "intBC": "i1", "r1": "NC"},
    ])
    aln = read_allele_table(df, site_cols=("r1",))
    by_name = dict(zip(aln.names, aln.seqs))
    assert by_name["c2"][0] == "?"


def test_read_allele_table_duplicate_rows_disagreeing_raises():
    df = pd.DataFrame([
        {"cellBC": "c1", "intBC": "i1", "r1": "AAAAA[10:2D]TTTTT"},
        {"cellBC": "c1", "intBC": "i1", "r1": "AAAAA[15:3D]TTTTT"},   # disagrees
    ])
    with pytest.raises(ValueError, match="disagree"):
        read_allele_table(df, site_cols=("r1",))


# --------------------------------------------------------------------------
# Camin-Sokal vs. Fitch: the irreversible model must actually differ
# --------------------------------------------------------------------------
def _three_taxon_tree() -> Tree:
    root = Node(name="root")
    x = Node(name="X", length=1.0)
    root.add_child(x)
    root.add_child(Node(name="leaf3", length=1.0))
    x.add_child(Node(name="leaf1", length=1.0))
    x.add_child(Node(name="leaf2", length=1.0))
    return Tree(root=root)


def test_camin_sokal_forbids_reversion_fitch_does_not():
    # ((leaf1,leaf2),leaf3), leaf1="A"(code 1), leaf2="0"(ancestral), leaf3="A"
    tree = _three_taxon_tree()
    aln = Alignment(["leaf1", "leaf2", "leaf3"], ["1", "0", "1"])
    assert parsimony_score(tree, aln) == 1.0
    assert camin_sokal_score(tree, aln) == 2.0


def test_missing_data_is_free_under_camin_sokal():
    # leaf1 and leaf2 are sisters (share stem node X); leaf3 stays ancestral
    # throughout. Forcing leaf2 to a *different* derived state than leaf1
    # ("2" vs "1") makes them conflict, so the stem X can't be shared and
    # each needs its own origin: cost 2. Marking leaf2 missing instead
    # removes that conflict -- X can freely be "1", leaf2 is compatible
    # with any state at zero cost, so the true minimum drops to 1.
    tree = _three_taxon_tree()
    conflicting = Alignment(["leaf1", "leaf2", "leaf3"], ["1", "2", "0"])
    missing = Alignment(["leaf1", "leaf2", "leaf3"], ["1", "?", "0"])
    assert camin_sokal_score(tree, conflicting) == 2.0
    assert camin_sokal_score(tree, missing) == 1.0


def test_sankoff_score_rejects_wrong_shaped_cost_matrix():
    import numpy as np
    tree = _three_taxon_tree()
    aln = Alignment(["leaf1", "leaf2", "leaf3"], ["1", "0", "1"])
    with pytest.raises(ValueError, match="cost_matrix"):
        sankoff_score(tree, aln, np.zeros((5, 5)))


# --------------------------------------------------------------------------
# lineage_tree end-to-end (small synthetic data)
# --------------------------------------------------------------------------
def test_lineage_tree_recovers_clades():
    aln = Alignment(
        names=["A1", "A2", "A3", "B1", "B2", "B3"],
        seqs=["1", "1", "1", "2", "2", "2"],
    )
    tree = lineage_tree(aln, search=True)
    leaf_sets = [frozenset(n.leaf_names()) for n in tree.traverse() if not n.is_leaf]
    assert frozenset({"A1", "A2", "A3"}) in leaf_sets
    assert frozenset({"B1", "B2", "B3"}) in leaf_sets
    assert tree.data["camin_sokal_score"] == 2.0
    assert tree.data["min_possible_score"] == 2.0
    assert tree.data["excess_origins"] == 0.0


def test_lineage_tree_on_real_data_subset():
    # small, fast, offline slice of the real bundled demo dataset -- checks
    # the full read_allele_table -> lineage_tree pipeline runs end-to-end on
    # genuine data, not a specific RF-distance-to-published-tree threshold
    # (see examples/lineage_demo.py for that honest, non-gated comparison).
    df = pd.read_csv("examples/data/lineage_alleletable.txt", sep="\t").head(60)
    aln = read_allele_table(df)
    tree = lineage_tree(aln, search=False)
    assert set(tree.leaf_names()) == set(aln.names)
    assert tree.data["camin_sokal_score"] >= tree.data["min_possible_score"] >= 0
