"""Phylogenetic signal (Blomberg's K, Pagel's lambda) and PGLS.

Mostly Monte Carlo checks against known statistical properties rather than
hand-computed numbers: these formulas have no simple closed-form answer to
assert against for an arbitrary tree, but each has a well-defined *expected*
behaviour over many simulated datasets (K's own defining property is that it
averages to 1 under the Brownian motion its own tree implies; PGLS's is the
Felsenstein 1985 demonstration that it, unlike OLS, does not inflate Type I
error for two independently-evolved traits on the same tree). Confirmed by
hand first against a real simulation run before writing these -- see the
commit message for the numbers -- these thresholds are set well inside the
observed margins, not right at them.
"""
import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest

import phytreon as pt


# --------------------------------------------------------------------------
# phylo_vcv: the shared building block
# --------------------------------------------------------------------------
def test_phylo_vcv_diagonal_is_root_to_tip_depth():
    tr = pt.datasets.primates()
    names, V = pt.phylo_vcv(tr)
    for i, name in enumerate(names):
        leaf = next(lf for lf in tr.leaves() if lf.name == name)
        depth = 0.0
        node = leaf
        while node.parent is not None:
            depth += node.length or 0.0
            node = node.parent
        assert V[i, i] == pytest.approx(depth)


def test_phylo_vcv_offdiagonal_is_shared_depth_to_the_mrca():
    tr = pt.datasets.primates()
    names, V = pt.phylo_vcv(tr)
    i, j = names.index("Human"), names.index("Chimp")
    mrca = tr.get_mrca(["Human", "Chimp"])
    shared = 0.0
    node = mrca
    while node.parent is not None:
        shared += node.length or 0.0
        node = node.parent
    assert V[i, j] == pytest.approx(shared)
    assert V[i, j] == V[j, i]


def test_phylo_vcv_respects_a_requested_taxon_order():
    tr = pt.datasets.primates()
    names = list(reversed(tr.leaf_names()))
    got, V = pt.phylo_vcv(tr, names)
    assert got == names
    assert V.shape == (len(names), len(names))


def test_phylo_vcv_rejects_an_unknown_taxon():
    tr = pt.datasets.primates()
    with pytest.raises(ValueError, match="not found"):
        pt.phylo_vcv(tr, ["NotReal"])


# --------------------------------------------------------------------------
# Blomberg's K
# --------------------------------------------------------------------------
def test_blomberg_k_averages_to_one_under_the_brownian_motion_it_defines():
    # K's own definition scales the observed statistic by its expectation
    # under BM on this exact tree, so simulating real BM data on that same
    # tree and averaging K over many replicates has to land near 1 -- this
    # is checking the formula, not just that it runs (measured: 300 sims on
    # a 30-tip random tree gave mean 0.98-1.05 across several seeds)
    tr = pt.datasets.random_tree(30, seed=7)
    names, V = pt.phylo_vcv(tr)
    rng = np.random.default_rng(0)
    ks = [pt.blomberg_k(tr, dict(zip(names, rng.multivariate_normal(
        np.zeros(len(names)), V))))["K"] for _ in range(300)]
    assert 0.8 < np.mean(ks) < 1.25


def test_blomberg_k_is_lower_on_average_for_signal_free_data():
    tr = pt.datasets.random_tree(30, seed=7)
    names, _ = pt.phylo_vcv(tr)
    rng = np.random.default_rng(1)
    ks_iid = [pt.blomberg_k(tr, dict(zip(names, rng.normal(size=len(names)))))["K"]
             for _ in range(300)]
    assert np.mean(ks_iid) < 0.85     # well under the ~1.0 real-BM baseline above


def test_blomberg_k_permutation_test_runs_and_bounds_p_in_zero_one():
    tr = pt.datasets.random_tree(15, seed=2)
    names, V = pt.phylo_vcv(tr)
    rng = np.random.default_rng(3)
    trait = dict(zip(names, rng.multivariate_normal(np.zeros(len(names)), V)))
    result = pt.blomberg_k(tr, trait, n_perm=200, seed=0)
    assert 0.0 <= result["p"] <= 1.0
    assert result["n_perm"] == 200


def test_blomberg_k_needs_at_least_three_taxa():
    tr = pt.datasets.random_tree(10, seed=1)
    names = tr.leaf_names()[:2]
    with pytest.raises(ValueError, match="at least 3"):
        pt.blomberg_k(tr, {n: 1.0 for n in names})


def test_blomberg_k_rejects_a_trait_naming_an_unknown_tip():
    # via phylo_vcv's own check -- blomberg_k does not duplicate it
    tr = pt.datasets.primates()
    trait = {n: 1.0 for n in tr.leaf_names()}
    trait["NotReal"] = 2.0
    with pytest.raises(ValueError, match="not found"):
        pt.blomberg_k(tr, trait)


# --------------------------------------------------------------------------
# Pagel's lambda
# --------------------------------------------------------------------------
def test_pagels_lambda_recovers_one_for_real_brownian_motion_on_average():
    # small trees give a genuinely bimodal lambda estimate (a documented
    # property of the estimator itself, not a bug -- Freckleton et al. 2002),
    # so this needs enough tips for that to have settled down; measured at
    # n=60 tips, mean 0.98 / median 1.0 across 200 sims
    tr = pt.datasets.random_tree(60, seed=5)
    names, V = pt.phylo_vcv(tr)
    rng = np.random.default_rng(0)
    lams = [pt.pagels_lambda(tr, dict(zip(names, rng.multivariate_normal(
        np.zeros(len(names)), V))))["lambda"] for _ in range(150)]
    assert np.median(lams) > 0.8


def test_pagels_lambda_recovers_zero_for_signal_free_data_on_average():
    tr = pt.datasets.random_tree(60, seed=5)
    names, _ = pt.phylo_vcv(tr)
    rng = np.random.default_rng(1)
    lams = [pt.pagels_lambda(tr, dict(zip(names, rng.normal(size=len(names)))))["lambda"]
           for _ in range(150)]
    assert np.median(lams) < 0.2


def test_pagels_lambda_is_bounded_and_reports_a_likelihood_ratio_test():
    tr = pt.datasets.random_tree(20, seed=9)
    names, V = pt.phylo_vcv(tr)
    rng = np.random.default_rng(2)
    trait = dict(zip(names, rng.multivariate_normal(np.zeros(len(names)), V)))
    result = pt.pagels_lambda(tr, trait)
    assert 0.0 <= result["lambda"] <= 1.0
    assert result["LR"] >= 0.0          # a likelihood-ratio statistic
    assert 0.0 <= result["p"] <= 1.0
    assert result["logLik"] >= result["logLik0"] - 1e-6  # lambda=ML fits >= lambda=0


def test_pagels_lambda_needs_at_least_three_taxa():
    tr = pt.datasets.random_tree(10, seed=1)
    names = tr.leaf_names()[:2]
    with pytest.raises(ValueError, match="at least 3"):
        pt.pagels_lambda(tr, {n: 1.0 for n in names})


# --------------------------------------------------------------------------
# PGLS: the Felsenstein (1985) demonstration -- OLS on two independently-
# evolved traits sharing one tree inflates Type I error; GLS against the
# same tree's covariance structure must not.
# --------------------------------------------------------------------------
def test_pgls_corrects_the_type_i_error_ols_gets_wrong_on_shared_ancestry():
    tr = pt.datasets.random_tree(40, seed=11)
    names, V = pt.phylo_vcv(tr)
    rng = np.random.default_rng(4)
    n_sims = 150
    ols_sig = pgls_sig = 0
    for _ in range(n_sims):
        xv = rng.multivariate_normal(np.zeros(len(names)), V)
        yv = rng.multivariate_normal(np.zeros(len(names)), V)   # independent of x
        xm, ym = xv - xv.mean(), yv - yv.mean()
        b = (xm @ ym) / (xm @ xm)
        resid = ym - b * xm
        n = len(xv)
        se_b = np.sqrt((resid @ resid) / (n - 2) / (xm @ xm))
        from scipy.stats import t as t_dist
        p_ols = 2 * t_dist.sf(abs(b / se_b), df=n - 2)
        ols_sig += p_ols < 0.05
        res = pt.pgls(tr, dict(zip(names, yv)), dict(zip(names, xv)))
        pgls_sig += res["p"]["x"] < 0.05
    # measured on this exact setup: OLS ~0.35-0.45, PGLS ~0.03-0.07
    assert ols_sig / n_sims > 0.20, "OLS should show its well-known inflation here"
    assert pgls_sig / n_sims < 0.15, "PGLS should not inflate the way OLS just did"


def test_pgls_recovers_a_known_slope_from_a_direct_causal_relationship():
    tr = pt.datasets.random_tree(30, seed=6)
    names, V = pt.phylo_vcv(tr)
    rng = np.random.default_rng(5)
    xv = rng.multivariate_normal(np.zeros(len(names)), V)
    yv = 3.0 * xv + rng.normal(scale=0.05, size=len(names))    # y = 3x + tiny noise
    res = pt.pgls(tr, dict(zip(names, yv)), dict(zip(names, xv)))
    assert res["coefficients"]["x"] == pytest.approx(3.0, abs=0.1)
    assert res["p"]["x"] < 0.001
    assert res["r2"] > 0.98


def test_pgls_accepts_a_dataframe_of_several_predictors():
    import pandas as pd
    tr = pt.datasets.random_tree(30, seed=6)
    names, V = pt.phylo_vcv(tr)
    rng = np.random.default_rng(5)
    x1 = rng.multivariate_normal(np.zeros(len(names)), V)
    x2 = rng.normal(size=len(names))
    yv = 2.0 * x1 - 1.0 * x2 + rng.normal(scale=0.05, size=len(names))
    X = pd.DataFrame({"x1": x1, "x2": x2}, index=names)
    res = pt.pgls(tr, dict(zip(names, yv)), X)
    assert res["coefficients"]["x1"] == pytest.approx(2.0, abs=0.15)
    assert res["coefficients"]["x2"] == pytest.approx(-1.0, abs=0.15)
    assert set(res["coefficients"]) == {"Intercept", "x1", "x2"}


def test_pgls_fixed_lambda_zero_matches_ordinary_least_squares_on_an_ultrametric_tree():
    # lambda=0 leaves each tip's own variance (the diagonal of phylo_vcv, its
    # root-to-tip depth) untouched and only zeroes the *shared* history off
    # the diagonal -- a star phylogeny only if every tip is equally deep to
    # begin with (an ultrametric tree), which is when this GLS fit actually
    # reduces to an equal-weight, ordinary least-squares fit; a tree with
    # unequal tip depths makes lambda=0 a *weighted* fit instead (still
    # correct, just not the comparison this test is making), which is why
    # random_tree() (not ultrametric) is not used here -- UPGMA always is.
    rng = np.random.default_rng(6)
    names = [f"t{i}" for i in range(20)]
    n = len(names)
    coords = rng.normal(size=(n, 2))
    D = [[float(np.linalg.norm(coords[i] - coords[j])) for j in range(n)] for i in range(n)]
    tr = pt.upgma(names, D)
    _, V = pt.phylo_vcv(tr, names)
    assert np.allclose(np.diag(V), np.diag(V)[0]), "UPGMA must be ultrametric"

    xv = rng.normal(size=n)
    yv = 2.0 * xv + rng.normal(scale=0.3, size=n)
    res = pt.pgls(tr, dict(zip(names, yv)), dict(zip(names, xv)), lambda_=0.0)
    ols_slope = np.polyfit(xv, yv, 1)[0]
    assert res["coefficients"]["x"] == pytest.approx(ols_slope, rel=1e-6)


def test_pgls_rejects_too_few_taxa_for_the_predictors_given():
    # one predictor + an intercept needs at least 3 taxa for a meaningful
    # fit (2 parameters plus at least one residual degree of freedom); this
    # tree has only 2
    tr = pt.Tree.from_newick("(A:1,B:1);")
    y = {"A": 1.0, "B": 2.0}
    x = {"A": 1.0, "B": 2.0}
    with pytest.raises(ValueError, match="needs more taxa"):
        pt.pgls(tr, y, x)
