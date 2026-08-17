"""Phylogenetic factorization: which edge of a tree best explains a covariate.

Validated by construction: abundances are simulated so a KNOWN clade's
response to a covariate is the only real signal in the data, and the test
checks that phylofactor's first factor names that exact clade (or its
complement -- an edge and its bipartition's other side are the same split).
A second construction plants two independent clades at different effect
strengths and checks both are recovered, strongest first, which is what the
greedy "residual bin" structure is supposed to guarantee.

Also measured directly rather than assumed: the reported p-value is the
winner of many candidate-edge tests, and reporting it as if it were one
calibrated test inflates the false-positive rate sharply -- confirmed on data
with no real association at all, and kept as a permanent test so the module
docstring's specific number does not quietly go stale.
"""
import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

import phytreon as pt


def _clade_of_size(tree, n):
    return set(next(nd.leaf_names() for nd in tree.traverse("postorder")
                    if not nd.is_leaf and len(nd.leaf_names()) == n))


def _lognormal_table(tips, n_samples, seed):
    rng = np.random.default_rng(seed)
    return rng.lognormal(2.0, 0.3, size=(n_samples, len(tips))), rng


def _side_matches(found: set, target: set, all_taxa: set) -> bool:
    """An edge and its bipartition's other side are the same split."""
    return found == target or found == (all_taxa - target)


TREE = pt.datasets.random_tree(30, seed=7)
TIPS = TREE.leaf_names()
ALL_TAXA = set(TIPS)
TARGET_CLADE = _clade_of_size(TREE, 6)


def test_recovers_a_single_planted_clade_exactly():
    baseline, rng = _lognormal_table(TIPS, 60, seed=0)
    covariate = rng.uniform(-1, 1, 60)
    effect = np.array([3.0 if t in TARGET_CLADE else 0.0 for t in TIPS])
    abundance = baseline * np.exp(np.outer(covariate, effect))
    table = pd.DataFrame(abundance, columns=TIPS, index=[f"s{i}" for i in range(60)])

    res = pt.phylofactor(TREE, table, covariate, n_factors=1)
    found = set(res["factors"].iloc[0]["side1"])
    assert _side_matches(found, TARGET_CLADE, ALL_TAXA)
    assert res["factors"].iloc[0]["p"] < 1e-10
    assert res["factors"].iloc[0]["r2"] > 0.9
    assert res["n_factors_found"] == 1
    assert list(res["balances"].columns) == ["factor1"]
    assert len(res["balances"]) == 60


def test_recovers_two_planted_clades_strongest_first():
    weak_clade = _clade_of_size(TREE, 4) - TARGET_CLADE
    # guard against an unlucky choice of tree/seed where the two clades overlap
    assert not (weak_clade & TARGET_CLADE)

    baseline, rng = _lognormal_table(TIPS, 60, seed=0)
    covariate = rng.uniform(-1, 1, 60)
    eff_strong = np.array([4.0 if t in TARGET_CLADE else 0.0 for t in TIPS])
    eff_weak = np.array([1.5 if t in weak_clade else 0.0 for t in TIPS])
    abundance = baseline * np.exp(np.outer(covariate, eff_strong)
                                  + np.outer(covariate, eff_weak))
    table = pd.DataFrame(abundance, columns=TIPS, index=[f"s{i}" for i in range(60)])

    res = pt.phylofactor(TREE, table, covariate, n_factors=2)
    factor1 = set(res["factors"].iloc[0]["side1"])
    factor2 = set(res["factors"].iloc[1]["side1"])
    assert _side_matches(factor1, TARGET_CLADE, ALL_TAXA)
    assert _side_matches(factor2, weak_clade, ALL_TAXA - TARGET_CLADE)
    assert res["factors"].iloc[0]["F"] > res["factors"].iloc[1]["F"]
    assert res["n_factors_found"] == 2


def test_categorical_covariate_recovers_the_same_planted_clade():
    baseline, rng = _lognormal_table(TIPS, 60, seed=0)
    groups = rng.choice(["treatment", "control"], size=60)
    effect = np.array([2.5 if t in TARGET_CLADE else 0.0 for t in TIPS])
    shift = (groups == "treatment").astype(float)
    abundance = baseline * np.exp(np.outer(shift, effect))
    table = pd.DataFrame(abundance, columns=TIPS, index=[f"s{i}" for i in range(60)])

    res = pt.phylofactor(TREE, table, groups, n_factors=1, categorical=True)
    found = set(res["factors"].iloc[0]["side1"])
    assert _side_matches(found, TARGET_CLADE, ALL_TAXA)
    assert res["factors"].iloc[0]["r2"] is None    # not defined for the ANOVA path


def test_covariate_can_be_pulled_from_a_table_column_by_name():
    baseline, rng = _lognormal_table(TIPS, 60, seed=0)
    covariate = rng.uniform(-1, 1, 60)
    effect = np.array([3.0 if t in TARGET_CLADE else 0.0 for t in TIPS])
    abundance = baseline * np.exp(np.outer(covariate, effect))
    table = pd.DataFrame(abundance, columns=TIPS, index=[f"s{i}" for i in range(60)])
    table["env"] = covariate

    res = pt.phylofactor(TREE, table, "env", n_factors=1)
    found = set(res["factors"].iloc[0]["side1"])
    assert _side_matches(found, TARGET_CLADE, ALL_TAXA)
    assert "env" not in res["factors"].iloc[0]["side1"] + res["factors"].iloc[0]["side2"]


# --------------------------------------------------------------------------
# The honesty check: the reported p is the winner of many tests
# --------------------------------------------------------------------------
def test_top_factor_p_value_is_badly_miscalibrated_on_pure_noise():
    # documented prominently in the module's own docstring for exactly this
    # reason -- kept as a test so that number cannot go stale silently
    tr = pt.datasets.random_tree(20, seed=2)
    tips = tr.leaf_names()
    n_reps, hits = 40, 0
    for rep in range(n_reps):
        rng = np.random.default_rng(rep)
        covariate = rng.uniform(-1, 1, 40)
        abundance = rng.lognormal(2.0, 0.3, size=(40, 20))   # no real relationship
        table = pd.DataFrame(abundance, columns=tips, index=[f"s{i}" for i in range(40)])
        res = pt.phylofactor(tr, table, covariate, n_factors=1)
        hits += res["factors"].iloc[0]["p"] < 0.05
    rate = hits / n_reps
    assert rate > 0.3, (
        f"expected pronounced multiple-testing inflation on pure noise, got "
        f"only {rate:.2f} significant at nominal 0.05 -- if this dropped, "
        f"something about the search or scoring changed"
    )


# --------------------------------------------------------------------------
# Mechanics and errors
# --------------------------------------------------------------------------
def test_rejects_table_columns_that_are_not_tips_of_the_tree():
    baseline, rng = _lognormal_table(TIPS, 20, seed=1)
    table = pd.DataFrame(baseline, columns=TIPS, index=[f"s{i}" for i in range(20)])
    table["Ghost"] = 1.0
    with pytest.raises(ValueError, match="not tips of the tree"):
        pt.phylofactor(TREE, table, rng.uniform(-1, 1, 20), n_factors=1)


def test_rejects_a_covariate_of_the_wrong_length():
    baseline, rng = _lognormal_table(TIPS, 20, seed=1)
    table = pd.DataFrame(baseline, columns=TIPS, index=[f"s{i}" for i in range(20)])
    with pytest.raises(ValueError, match="15 entries for"):
        pt.phylofactor(TREE, table, rng.uniform(-1, 1, 15), n_factors=1)


def test_rejects_an_unknown_covariate_column_name():
    baseline, rng = _lognormal_table(TIPS, 20, seed=1)
    table = pd.DataFrame(baseline, columns=TIPS, index=[f"s{i}" for i in range(20)])
    with pytest.raises(ValueError, match="not a column"):
        pt.phylofactor(TREE, table, "not_a_real_column", n_factors=1)


def test_rejects_more_factors_than_the_taxa_can_support():
    small_tree = pt.Tree.from_newick("((A:1,B:1):1,(C:1,D:1):1);")
    table = pd.DataFrame(np.ones((5, 4)), columns=["A", "B", "C", "D"])
    with pytest.raises(ValueError, match="n_factors must be >= 1"):
        pt.phylofactor(small_tree, table, [0, 1, 0, 1, 0], n_factors=0)
    with pytest.raises(ValueError, match="cannot support"):
        pt.phylofactor(small_tree, table, [0, 1, 0, 1, 0], n_factors=3)
