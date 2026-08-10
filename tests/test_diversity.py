"""Phylogenetic diversity: Faith's PD and UniFrac."""
import matplotlib
matplotlib.use("Agg")

import pytest

import phytreon as pt


def _total_branch_length(tree):
    return sum(n.length or 0.0 for n in tree.traverse() if not n.is_root)


# --------------------------------------------------------------------------
# Faith's PD
# --------------------------------------------------------------------------
def test_faiths_pd_of_every_taxon_is_the_whole_trees_branch_length():
    tr = pt.datasets.primates()
    assert pt.faiths_pd(tr, tr.leaf_names()) == pytest.approx(_total_branch_length(tr))


def test_faiths_pd_of_one_taxon_is_its_own_depth():
    tr = pt.datasets.primates()
    name = tr.leaf_names()[0]
    _, V = pt.phylo_vcv(tr, [name])
    assert pt.faiths_pd(tr, name) == pytest.approx(V[0, 0])


def test_faiths_pd_is_monotonic_in_the_taxon_set():
    tr = pt.datasets.primates()
    names = tr.leaf_names()
    small = pt.faiths_pd(tr, names[:3])
    big = pt.faiths_pd(tr, names)
    assert small <= big


def test_faiths_pd_of_nothing_is_zero():
    tr = pt.datasets.primates()
    assert pt.faiths_pd(tr, []) == 0.0


def test_faiths_pd_rejects_an_unknown_taxon():
    tr = pt.datasets.primates()
    with pytest.raises(ValueError, match="not found"):
        pt.faiths_pd(tr, ["NotReal"])


def test_faiths_pd_table_matches_the_single_sample_function():
    import pandas as pd
    tr = pt.datasets.primates()
    names = tr.leaf_names()
    table = pd.DataFrame(
        {"sampleA": [1, 1, 0, 0, 0, 0, 0], "sampleB": [0, 0, 1, 1, 1, 0, 0]},
        index=names).T
    out = pt.faiths_pd_table(tr, table)
    assert out["sampleA"] == pytest.approx(pt.faiths_pd(tr, names[:2]))
    assert out["sampleB"] == pytest.approx(pt.faiths_pd(tr, names[2:5]))


# --------------------------------------------------------------------------
# UniFrac
# --------------------------------------------------------------------------
def _split(tr):
    apes = tr.get_mrca(["Human", "Chimp", "Gorilla", "Orangutan", "Gibbon"]).leaf_names()
    monkeys = tr.get_mrca(["Macaque", "Baboon"]).leaf_names()
    return apes, monkeys


def test_unweighted_unifrac_of_a_sample_against_itself_is_zero():
    tr = pt.datasets.primates()
    apes, _ = _split(tr)
    assert pt.unweighted_unifrac(tr, apes, apes) == 0.0


def test_unweighted_unifrac_of_a_clean_bipartition_is_one():
    # apes and Old World monkeys share no branch at all in this topology --
    # every branch belongs to exactly one side, so nothing is shared
    tr = pt.datasets.primates()
    apes, monkeys = _split(tr)
    assert set(apes) | set(monkeys) == set(tr.leaf_names())
    assert pt.unweighted_unifrac(tr, apes, monkeys) == pytest.approx(1.0)


def test_unweighted_unifrac_is_symmetric():
    tr = pt.datasets.primates()
    names = tr.leaf_names()
    a, b = names[:2], names[2:5]
    assert pt.unweighted_unifrac(tr, a, b) == pytest.approx(pt.unweighted_unifrac(tr, b, a))


def test_unweighted_unifrac_rejects_an_unknown_taxon():
    tr = pt.datasets.primates()
    with pytest.raises(ValueError, match="not found"):
        pt.unweighted_unifrac(tr, ["NotReal"], ["Human"])


def test_weighted_unifrac_of_a_clean_bipartition_is_one_normalized():
    tr = pt.datasets.primates()
    apes, monkeys = _split(tr)
    a = {n: 1.0 for n in apes}
    b = {n: 3.0 for n in monkeys}       # abundances need not match scale
    assert pt.weighted_unifrac(tr, a, b, normalized=True) == pytest.approx(1.0)


def test_weighted_unifrac_matches_unweighted_for_single_tip_samples():
    # with exactly one taxon per side, "what fraction of this sample's
    # abundance sits below this branch" is 1 below the tip and 0 everywhere
    # else regardless of the actual abundance value -- precisely the
    # 0/1 membership unweighted UniFrac already computes, so the two must
    # agree exactly, in both the raw and the normalized form
    tr = pt.datasets.primates()
    a_name, b_name = "Human", "Macaque"
    unweighted = pt.unweighted_unifrac(tr, [a_name], [b_name])
    for abundance_a, abundance_b in ((1.0, 1.0), (5.0, 3.0), (0.2, 100.0)):
        a, b = {a_name: abundance_a}, {b_name: abundance_b}
        assert pt.weighted_unifrac(tr, a, b, normalized=True) == \
            pytest.approx(unweighted)


def test_weighted_unifrac_is_invariant_to_one_sides_overall_abundance_scale():
    # each side is normalised to its own total internally (a raw count
    # table and a relative-abundance table for the same sample must give
    # the same answer), so multiplying every one of a side's abundances by
    # the same constant must not change the result, raw or normalized
    tr = pt.datasets.primates()
    apes, monkeys = _split(tr)
    b = {n: 1.0 for n in monkeys}
    small = {n: 1.0 for n in apes}
    scaled = {n: 10.0 for n in apes}
    assert pt.weighted_unifrac(tr, small, b, normalized=False) == \
        pytest.approx(pt.weighted_unifrac(tr, scaled, b, normalized=False))
    assert pt.weighted_unifrac(tr, small, b) == \
        pytest.approx(pt.weighted_unifrac(tr, scaled, b))


def test_unifrac_matrix_diagonal_is_zero_and_matches_pairwise_calls():
    import pandas as pd
    tr = pt.datasets.primates()
    names = tr.leaf_names()
    table = pd.DataFrame(
        [[1, 1, 0, 0, 0, 0, 0],
        [0, 0, 1, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, 1, 1]],
        columns=names, index=["s1", "s2", "s3"])
    mat = pt.unifrac_matrix(tr, table)
    assert (mat.to_numpy().diagonal() == 0).all()
    expected = pt.unweighted_unifrac(tr, ["Human", "Chimp"], ["Gorilla", "Orangutan"])
    assert mat.loc["s1", "s2"] == pytest.approx(expected)
    assert mat.loc["s1", "s2"] == pytest.approx(mat.loc["s2", "s1"])


def test_unifrac_matrix_weighted_matches_pairwise_calls():
    import pandas as pd
    tr = pt.datasets.primates()
    names = tr.leaf_names()
    table = pd.DataFrame(
        [[1, 2, 0, 0, 0, 0, 0],
        [0, 0, 3, 1, 0, 0, 0]],
        columns=names, index=["s1", "s2"])
    mat = pt.unifrac_matrix(tr, table, weighted=True)
    expected = pt.weighted_unifrac(
        tr, {"Human": 1, "Chimp": 2}, {"Gorilla": 3, "Orangutan": 1})
    assert mat.loc["s1", "s2"] == pytest.approx(expected)


def test_unifrac_matrix_matches_the_pairwise_functions_on_random_tables():
    # unifrac_matrix computes every pair as one numpy expression over a
    # samples x edges array, while the pairwise functions walk dicts per call --
    # two genuinely different code paths for the same quantity, so this pins
    # them together across all three modes on trees and tables it has not been
    # tuned for. Guards the vectorised path against a rewrite quietly changing
    # the numbers rather than only the runtime.
    import numpy as np
    import pandas as pd
    rng = np.random.default_rng(0)
    for trial in range(6):
        n_tips = int(rng.integers(6, 40))
        tr = pt.datasets.random_tree(n_tips, seed=trial)
        tips = tr.leaf_names()
        n_s = int(rng.integers(3, 7))
        rows = []
        for _ in range(n_s):
            v = np.zeros(n_tips)
            k = int(rng.integers(1, n_tips + 1))
            v[rng.choice(n_tips, k, replace=False)] = rng.integers(1, 200, k)
            rows.append(v)
        table = pd.DataFrame(rows, columns=tips,
                            index=[f"s{i}" for i in range(n_s)])
        mats = {
            "unweighted": pt.unifrac_matrix(tr, table),
            "weighted": pt.unifrac_matrix(tr, table, weighted=True),
            "raw": pt.unifrac_matrix(tr, table, weighted=True, normalized=False),
        }
        for m in mats.values():
            assert np.allclose(np.diag(m), 0.0)
            assert np.allclose(m, m.T)
        for i in range(n_s):
            for j in range(i + 1, n_s):
                a, b = table.iloc[i], table.iloc[j]
                pa = list(table.columns[a.to_numpy() > 0])
                pb = list(table.columns[b.to_numpy() > 0])
                assert mats["unweighted"].iat[i, j] == pytest.approx(
                    pt.unweighted_unifrac(tr, pa, pb))
                assert mats["weighted"].iat[i, j] == pytest.approx(
                    pt.weighted_unifrac(tr, a.to_dict(), b.to_dict()))
                assert mats["raw"].iat[i, j] == pytest.approx(
                    pt.weighted_unifrac(tr, a.to_dict(), b.to_dict(),
                                        normalized=False))


# --------------------------------------------------------------------------
# Table columns that are not tips of the tree. The realistic case, not an
# exotic one: tree building drops sequences (too short, failed alignment,
# chimeras filtered after the table was counted), so an ASV table normally has
# more ASVs than the tree has tips.
# --------------------------------------------------------------------------
def _table_with_an_extra_taxon():
    import pandas as pd
    tr = pt.Tree.from_newick("((A:1,B:1):1,(C:1,D:1):1);")
    table = pd.DataFrame(
        {"A": [10, 0], "B": [0, 10], "C": [10, 0], "D": [0, 10],
         "NotInTree": [90, 0]},
        index=["s1", "s2"])
    return tr, table


def test_table_functions_reject_columns_that_are_not_tips_of_the_tree():
    # weighted unifrac_matrix used to accept this silently and return a wrong
    # number rather than no number: the extra taxon's abundance still landed in
    # each sample's total, so every real taxon's fraction came out too small
    # (0.846 where the answer is 0.5, on this exact table). Unweighted raised a
    # bare KeyError, and the single-pair functions a clear ValueError -- three
    # behaviours for one mistake.
    tr, table = _table_with_an_extra_taxon()
    for call in (lambda: pt.unifrac_matrix(tr, table),
                 lambda: pt.unifrac_matrix(tr, table, weighted=True),
                 lambda: pt.faiths_pd_table(tr, table)):
        with pytest.raises(ValueError, match="not tips of the tree"):
            call()


def test_the_rejection_message_says_how_many_and_names_them():
    tr, table = _table_with_an_extra_taxon()
    try:
        pt.unifrac_matrix(tr, table, weighted=True)
    except ValueError as exc:
        assert "1 of 5" in str(exc)
        assert "NotInTree" in str(exc)


def test_subsetting_the_table_to_the_tree_gives_the_pairwise_answer():
    # the fix the error message tells you to apply, and the value it lands on
    tr, table = _table_with_an_extra_taxon()
    kept = [c for c in table.columns if c in set(tr.leaf_names())]
    mat = pt.unifrac_matrix(tr, table[kept], weighted=True)
    assert mat.loc["s1", "s2"] == pytest.approx(pt.weighted_unifrac(
        tr, table.loc["s1", kept].to_dict(), table.loc["s2", kept].to_dict()))
    assert mat.loc["s1", "s2"] == pytest.approx(0.5)


def test_an_all_zero_extra_column_is_rejected_too():
    # checked over the columns up front, not per row as a side effect of which
    # taxa happened to be nonzero there -- otherwise the same table raises or
    # not depending on the data in it
    import pandas as pd
    tr = pt.Tree.from_newick("((A:1,B:1):1,(C:1,D:1):1);")
    table = pd.DataFrame({"A": [1, 0], "B": [0, 1], "C": [1, 0], "D": [0, 1],
                          "NotInTree": [0, 0]}, index=["s1", "s2"])
    with pytest.raises(ValueError, match="not tips of the tree"):
        pt.faiths_pd_table(tr, table)
