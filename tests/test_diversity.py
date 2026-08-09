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
