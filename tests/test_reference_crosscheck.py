"""Comparative methods against something other than themselves.

The rest of the comparative tests check each statistic against the property it
is *defined* to have -- Blomberg's K averaging to 1 under Brownian motion, PGLS
not inflating Type I error, and so on. Strong, but all of it circular in one
respect: the same implementation produces both the number and the behaviour
being checked. These tests bring in outside arithmetic instead.

Two sources, both pure Python, since a Python package that needed R to verify
itself would be a contradiction:

* ``statsmodels.GLS`` -- PGLS at a fixed lambda *is* generalised least squares
  with a known error covariance, so a mature independent implementation of that
  exact computation can be held against ours coefficient by coefficient.
* ``sympy`` -- exact rational arithmetic on a tree small enough that the
  expected values can also be read off it by hand, which removes floating point
  and indexing from the question entirely.

Neither is a runtime dependency (both live in the ``dev`` extra), so each test
skips rather than fails where they are missing.
"""
import numpy as np
import pytest

import phytreon as pt


# A tree small enough to read the answers off by hand:
#
#         root
#        /    \
#     n1:1    n2:2
#     /  \     /  \
#   A:1  B:3 C:1  D:2
#
# root-to-tip depths  A=2, B=4, C=3, D=4
# shared depth        A&B = 1 (n1), C&D = 2 (n2), anything across the root = 0
# total branch length 1+3+1+1+2+2 = 10
HAND_TREE = "((A:1,B:3):1,(C:1,D:2):2);"
HAND_VCV = [[2, 1, 0, 0],
            [1, 4, 0, 0],
            [0, 0, 3, 2],
            [0, 0, 2, 4]]


def test_phylo_vcv_matches_a_matrix_read_off_the_tree_by_hand():
    tr = pt.Tree.from_newick(HAND_TREE)
    names, V = pt.phylo_vcv(tr, ["A", "B", "C", "D"])
    assert names == ["A", "B", "C", "D"]
    assert np.allclose(V, np.array(HAND_VCV, dtype=float))


def test_faiths_pd_and_unifrac_match_hand_arithmetic():
    tr = pt.Tree.from_newick(HAND_TREE)
    # every taxon -> the tree's whole branch length
    assert pt.faiths_pd(tr, ["A", "B", "C", "D"]) == pytest.approx(10.0)
    # A alone -> its own branch plus n1 above it
    assert pt.faiths_pd(tr, ["A"]) == pytest.approx(2.0)
    assert pt.faiths_pd(tr, ["B"]) == pytest.approx(4.0)
    # {A,B} and {C,D} share no branch at all
    assert pt.unweighted_unifrac(tr, ["A", "B"], ["C", "D"]) == pytest.approx(1.0)
    # {A} vs {A,B}: edges {A,n1} against {A,B,n1}; B (3) unshared of 1+3+1 = 5
    assert pt.unweighted_unifrac(tr, ["A"], ["A", "B"]) == pytest.approx(0.6)
    # weighted, all abundance on one tip each side, is the unweighted answer
    assert pt.weighted_unifrac(tr, {"A": 5.0}, {"C": 9.0}) == pytest.approx(
        pt.unweighted_unifrac(tr, ["A"], ["C"]))


def test_blomberg_k_matches_exact_rational_arithmetic():
    sp = pytest.importorskip("sympy")
    tr = pt.Tree.from_newick(HAND_TREE)
    values = [1, 4, 2, 7]

    # the definition, in exact rationals, with nothing of phytreon's involved
    V = sp.Matrix(HAND_VCV)
    x = sp.Matrix([sp.Rational(v) for v in values])
    n = 4
    V_inv = V.inv()
    one = sp.ones(n, 1)
    a_hat = (one.T * V_inv * x)[0] / (one.T * V_inv * one)[0]
    resid = x - a_hat * one
    mse0 = (resid.T * resid)[0] / (n - 1)
    mse = (resid.T * V_inv * resid)[0] / (n - 1)
    expected_ratio = (V.trace() - sp.Rational(n) / (one.T * V_inv * one)[0]) / (n - 1)
    k_exact = (mse0 / mse) / expected_ratio

    got = pt.blomberg_k(tr, dict(zip(["A", "B", "C", "D"], map(float, values))))["K"]
    assert got == pytest.approx(float(k_exact), rel=1e-12)
    # and it really is an exact rational, not a float that happened to agree
    assert k_exact == sp.Rational(70014, 93775)


@pytest.mark.parametrize("n_tips,lam,n_pred", [
    (20, 1.0, 1),     # untransformed Brownian motion
    (20, 0.0, 1),     # star phylogeny -- the ordinary least squares end
    (30, 0.5, 1),     # partly transformed
    (40, 1.0, 2),     # several predictors
    (25, 0.3, 3),
])
def test_pgls_matches_statsmodels_gls_at_a_fixed_lambda(n_tips, lam, n_pred):
    sm = pytest.importorskip("statsmodels.api")
    import pandas as pd

    tr = pt.datasets.random_tree(n_tips, seed=n_tips)
    names, V = pt.phylo_vcv(tr)
    rng = np.random.default_rng(n_tips * 10 + n_pred)
    chol = np.linalg.cholesky(V)
    x_cols = {f"x{k}": chol @ rng.normal(size=n_tips) for k in range(n_pred)}
    y = sum(x_cols.values()) + rng.normal(scale=0.5, size=n_tips)

    x_arg = (pd.DataFrame(x_cols, index=names) if n_pred > 1
             else dict(zip(names, x_cols["x0"])))
    mine = pt.pgls(tr, dict(zip(names, y)), x_arg, lambda_=lam)

    # pgls sorts the taxa it uses, so line statsmodels up the same way
    taxa = sorted(names)
    order = [names.index(t) for t in taxa]
    _, V_sorted = pt.phylo_vcv(tr, taxa)
    sigma = V_sorted * lam
    np.fill_diagonal(sigma, np.diag(V_sorted))
    design = sm.add_constant(
        np.column_stack([x_cols[f"x{k}"][order] for k in range(n_pred)]))
    ref = sm.GLS(y[order], design, sigma=sigma).fit()

    keys = ["Intercept"] + (["x"] if n_pred == 1
                            else [f"x{k}" for k in range(n_pred)])
    assert mine["dof"] == ref.df_resid
    for label, got, want in (("coefficients", mine["coefficients"], ref.params),
                             ("se", mine["se"], ref.bse),
                             ("t", mine["t"], ref.tvalues),
                             ("p", mine["p"], ref.pvalues)):
        assert np.allclose([got[k] for k in keys], np.asarray(want),
                           rtol=1e-10, atol=1e-12), label
