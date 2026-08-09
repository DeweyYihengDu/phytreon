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


def test_pagels_lambda_reports_an_actual_log_likelihood():
    # valid absolutely, not just up to the additive constant that cancels
    # inside the likelihood ratio -- checked against the multivariate normal
    # density at the fitted parameters, computed independently by scipy
    from scipy.stats import multivariate_normal
    tr = pt.datasets.random_tree(15, seed=12)
    names, V = pt.phylo_vcv(tr)
    rng = np.random.default_rng(4)
    x = np.linalg.cholesky(V) @ rng.normal(size=15)
    trait = dict(zip(names, x))
    res = pt.pagels_lambda(tr, trait)

    Vl = V * res["lambda"]
    np.fill_diagonal(Vl, np.diag(V))
    Vi = np.linalg.inv(Vl)
    ones = np.ones(len(x))
    a_hat = float(ones @ Vi @ x / (ones @ Vi @ ones))
    resid = x - a_hat
    sigma2 = float(resid @ Vi @ resid) / len(x)
    expected = multivariate_normal.logpdf(x, mean=a_hat * ones, cov=sigma2 * Vl)
    assert res["logLik"] == pytest.approx(expected)


def test_pagels_lambda_significance_errs_conservative_not_permissive():
    # trait drawn iid across tips, so there is no signal and every p < 0.05 is
    # a false positive. lambda = 0 sits on the boundary of [0, 1], which makes
    # the LR's true null a 50:50 chi2_0/chi2_1 mixture; testing against plain
    # chi2_1 demands a larger LR than that and so under-rejects, which is the
    # safe direction. Measured at 800 replicates: 0.5-2.0% against a nominal 5%.
    tr = pt.datasets.random_tree(20, seed=120)
    names = tr.leaf_names()
    rng = np.random.default_rng(7)
    rate = np.mean([pt.pagels_lambda(tr, dict(zip(names, rng.normal(size=20))))["p"]
                    < 0.05 for _ in range(300)])
    assert rate < 0.05


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
    # measured on this exact 40-tip setup over 300 sims at three seeds:
    # OLS 0.307-0.373, PGLS 0.053-0.057. The OLS figure is specific to this
    # tree size -- the inflation comes from shared ancestry, so it grows with
    # the tree rather than being one constant (16% at 10 tips, 43% at 80), and
    # the sweep behind those numbers lives in the README.
    assert ols_sig / n_sims > 0.20, "OLS should show its well-known inflation here"
    # tight enough to catch a real regression: PGLS sits near the nominal 0.05
    # here, so anything approaching OLS territory should fail rather than pass
    assert pgls_sig / n_sims < 0.10, "PGLS should not inflate the way OLS just did"


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


# --------------------------------------------------------------------------
# How lambda is estimated, and what that costs -- the small-sample
# false-positive story. Measured over a calibration sweep (tree shape x tree
# size x true lambda) before these were written: pooled over the cells with
# 10-20 taxa, at a nominal 5%, lambda_="ML" rejected 8.1% of true nulls,
# lambda_="REML" 7.1%, lambda_=1.0 8.5% (600 replicates a cell). A separate
# and larger run put ML at 7.9% and REML at 6.8% (6400 replicates) -- the same
# quantities, and the spread between the two runs is why the README quotes the
# REML figure as "~7%" rather than implying more precision than there is.
# --------------------------------------------------------------------------
def test_pgls_estimates_lambda_by_reml_by_default():
    tr = pt.datasets.random_tree(20, seed=3)
    names, V = pt.phylo_vcv(tr)
    rng = np.random.default_rng(0)
    xv = np.linalg.cholesky(V) @ rng.normal(size=20)
    yv = 2.0 * xv + rng.normal(scale=0.5, size=20)
    res = pt.pgls(tr, dict(zip(names, yv)), dict(zip(names, xv)))
    assert res["lambda_method"] == "REML"
    assert pt.pgls(tr, dict(zip(names, yv)), dict(zip(names, xv)),
                   lambda_="ML")["lambda_method"] == "ML"
    assert pt.pgls(tr, dict(zip(names, yv)), dict(zip(names, xv)),
                   lambda_=1.0)["lambda_method"] == "fixed"


def test_reml_lambda_is_less_downward_biased_than_ml_in_small_samples():
    # data simulated on the untransformed tree, so the true lambda is 1. Both
    # estimators fall short of it at 10 tips; ML falls much further, and a
    # too-small lambda is what understates the tips' dependence and so
    # understates the standard errors. Measured: ML ~0.60, REML ~0.78.
    tr = pt.datasets.random_tree(10, seed=110)
    names, V = pt.phylo_vcv(tr)
    L = np.linalg.cholesky(V)
    rng = np.random.default_rng(0)
    lam_ml, lam_reml = [], []
    for _ in range(150):
        yv = L @ rng.normal(size=10)
        xv = L @ rng.normal(size=10)
        yd, xd = dict(zip(names, yv)), dict(zip(names, xv))
        lam_ml.append(pt.pgls(tr, yd, xd, lambda_="ML")["lambda"])
        lam_reml.append(pt.pgls(tr, yd, xd, lambda_="REML")["lambda"])
    assert np.mean(lam_ml) < np.mean(lam_reml) - 0.1
    assert np.mean(lam_reml) < 1.0          # still biased low, just less so


def test_ml_lambda_claims_significance_reml_does_not_and_never_the_reverse():
    # paired -- both estimators see the same simulated data, and x and y are
    # independent so every rejection is a false positive. The extra rejections
    # ML produces are strictly one-directional: measured over 400 replicates
    # at three seeds, "REML rejected but ML did not" happened 0 times, while
    # "ML rejected but REML did not" happened 5, 12 and 16 times.
    tr = pt.datasets.random_tree(10, seed=110)
    names, V = pt.phylo_vcv(tr)
    L = np.linalg.cholesky(V)
    rng = np.random.default_rng(2)
    ml_only = reml_only = 0
    for _ in range(400):
        xv, yv = L @ rng.normal(size=10), L @ rng.normal(size=10)
        yd, xd = dict(zip(names, yv)), dict(zip(names, xv))
        sig_ml = pt.pgls(tr, yd, xd, lambda_="ML")["p"]["x"] < 0.05
        sig_reml = pt.pgls(tr, yd, xd, lambda_="REML")["p"]["x"] < 0.05
        ml_only += sig_ml and not sig_reml
        reml_only += sig_reml and not sig_ml
    assert reml_only == 0
    assert ml_only >= 4


def test_estimating_lambda_costs_a_degree_of_freedom():
    tr = pt.datasets.random_tree(20, seed=3)
    names, V = pt.phylo_vcv(tr)
    rng = np.random.default_rng(0)
    xv = np.linalg.cholesky(V) @ rng.normal(size=20)
    yv = 2.0 * xv + rng.normal(scale=0.5, size=20)
    yd, xd = dict(zip(names, yv)), dict(zip(names, xv))
    given = pt.pgls(tr, yd, xd, lambda_=1.0)
    estimated = pt.pgls(tr, yd, xd)
    assert given["dof"] == 20 - 2                 # intercept + one predictor
    assert estimated["dof"] == given["dof"] - 1   # plus lambda itself


def test_pgls_rejects_an_unrecognised_lambda_keyword():
    tr = pt.datasets.random_tree(10, seed=1)
    names = tr.leaf_names()
    d = {n: float(i) for i, n in enumerate(names)}
    with pytest.raises(ValueError, match="'REML', 'ML', or a number"):
        pt.pgls(tr, d, d, lambda_="bogus")


def test_pgls_refuses_a_fit_with_no_residual_degrees_of_freedom_left():
    # 4 taxa, intercept + 2 predictors + an estimated lambda = 4 parameters
    import pandas as pd
    tr = pt.datasets.random_tree(4, seed=1)
    names = tr.leaf_names()
    rng = np.random.default_rng(0)
    y = dict(zip(names, rng.normal(size=4)))
    X = pd.DataFrame({"x1": rng.normal(size=4), "x2": rng.normal(size=4)},
                     index=names)
    with pytest.raises(ValueError, match="no residual degrees of freedom"):
        pt.pgls(tr, y, X)
    # handing lambda over instead frees that degree of freedom back up
    assert pt.pgls(tr, y, X, lambda_=1.0)["dof"] == 1


# --------------------------------------------------------------------------
# The parametric bootstrap p-value
# --------------------------------------------------------------------------
def test_pgls_bootstrap_covers_the_predictors_only_and_is_reproducible():
    tr = pt.datasets.random_tree(15, seed=4)
    names, V = pt.phylo_vcv(tr)
    rng = np.random.default_rng(0)
    xv = np.linalg.cholesky(V) @ rng.normal(size=15)
    yv = xv + rng.normal(scale=1.0, size=15)
    yd, xd = dict(zip(names, yv)), dict(zip(names, xv))
    a = pt.pgls(tr, yd, xd, n_boot=100, seed=7)
    b = pt.pgls(tr, yd, xd, n_boot=100, seed=7)
    assert a["p_boot"] == b["p_boot"]              # same seed, same answer
    assert set(a["p_boot"]) == {"x"}               # no null for the intercept
    assert a["n_boot"] == 100
    assert 1 / 101 <= a["p_boot"]["x"] <= 1.0
    assert "p_boot" not in pt.pgls(tr, yd, xd)     # opt-in


def test_pgls_bootstrap_p_cannot_go_below_its_own_resolution():
    # an overwhelming real effect: the t-based p is astronomically small, but a
    # p-value counted out of 200 null draws can only ever report 1/201
    tr = pt.datasets.random_tree(15, seed=4)
    names, V = pt.phylo_vcv(tr)
    rng = np.random.default_rng(0)
    xv = np.linalg.cholesky(V) @ rng.normal(size=15)
    yv = 5.0 * xv + rng.normal(scale=0.01, size=15)
    res = pt.pgls(tr, dict(zip(names, yv)), dict(zip(names, xv)),
                  n_boot=200, seed=1)
    assert res["p"]["x"] < 1e-10
    assert res["p_boot"]["x"] == pytest.approx(1 / 201)


def test_pgls_bootstrap_tests_each_predictor_against_a_reduced_null_model():
    # x1 drives y, x2 does not; each predictor's null keeps the other one, so
    # x2's non-effect has to survive x1 being in the model
    import pandas as pd
    tr = pt.datasets.random_tree(25, seed=8)
    names, V = pt.phylo_vcv(tr)
    rng = np.random.default_rng(3)
    x1 = np.linalg.cholesky(V) @ rng.normal(size=25)
    x2 = rng.normal(size=25)
    yv = 3.0 * x1 + rng.normal(scale=0.3, size=25)
    X = pd.DataFrame({"x1": x1, "x2": x2}, index=names)
    res = pt.pgls(tr, dict(zip(names, yv)), X, n_boot=150, seed=2)
    assert set(res["p_boot"]) == {"x1", "x2"}
    assert res["p_boot"]["x1"] < 0.02      # real effect
    assert res["p_boot"]["x2"] > 0.10      # none


def test_pgls_rejects_too_few_taxa_for_the_predictors_given():
    # one predictor + an intercept needs at least 3 taxa for a meaningful
    # fit (2 parameters plus at least one residual degree of freedom); this
    # tree has only 2
    tr = pt.Tree.from_newick("(A:1,B:1);")
    y = {"A": 1.0, "B": 2.0}
    x = {"A": 1.0, "B": 2.0}
    with pytest.raises(ValueError, match="needs more taxa"):
        pt.pgls(tr, y, x)
