"""Models of continuous trait evolution (BM/OU/EB/lambda/white), phylogenetic
PCA, and Fritz & Purvis' D for binary traits.

Each is checked against something with a known answer rather than against
itself:

* OU and EB both contain Brownian motion as a limit (``alpha -> 0``, ``b -> 0``),
  so their likelihoods have to converge on the BM likelihood there. That tests
  the covariance formulas, which were taken from the literature and verified
  symbolically rather than recalled.
* Model choice is tested by simulating under each model in turn and checking
  AICc recovers the one that generated the data -- the property the whole
  exercise exists for.
* Phylogenetic PCA has an exact reference case: on a star tree there is no
  phylogeny to correct for, so it must reduce to ordinary PCA.
* D is *defined* by two calibration points -- 0 for a Brownian threshold trait,
  1 for a randomly scattered one -- so simulating both and checking the means
  land there tests the scaling itself, not merely that it runs.
"""
import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest

import phytreon as pt
from phytreon.comparative.models import _log_lik, _shape


def _bm_setup(n=40, seed=2):
    tr = pt.datasets.random_tree(n, seed=seed)
    names, V = pt.phylo_vcv(tr)
    depths = np.diag(V).copy()
    patristic = depths[:, None] + depths[None, :] - 2.0 * V
    return tr, names, V, depths, patristic


def _simulate(names, C, seed):
    L = np.linalg.cholesky(C + 1e-10 * np.eye(len(names)))
    rng = np.random.default_rng(seed)
    return dict(zip(names, L @ rng.normal(size=len(names))))


# --------------------------------------------------------------------------
# The covariance formulas, via the Brownian limits they must contain
# --------------------------------------------------------------------------
def test_ou_likelihood_converges_on_brownian_motion_as_alpha_goes_to_zero():
    tr, names, V, depths, patristic = _bm_setup()
    trait = _simulate(names, V, 0)
    y = np.array([trait[n] for n in names])
    bm = pt.fit_continuous(tr, trait, "BM")["logLik"]
    close = _log_lik(_shape("OU", 1e-8, V, depths, patristic), y)
    far = _log_lik(_shape("OU", 1e-4, V, depths, patristic), y)
    assert close == pytest.approx(bm, abs=1e-7)
    # and it is actually approaching, not coincidentally equal at one point
    assert abs(close - bm) < abs(far - bm)


def test_eb_likelihood_converges_on_brownian_motion_as_b_goes_to_zero():
    tr, names, V, depths, patristic = _bm_setup()
    trait = _simulate(names, V, 0)
    y = np.array([trait[n] for n in names])
    bm = pt.fit_continuous(tr, trait, "BM")["logLik"]
    assert _log_lik(_shape("EB", -1e-9, V, depths, patristic),
                    y) == pytest.approx(bm, abs=1e-6)


def test_the_ou_covariance_is_symmetric_on_a_non_ultrametric_tree():
    # the form usually quoted is written with T_i alone, which is fine on an
    # ultrametric tree and silently asymmetric when tips differ in depth --
    # which random_tree's do. The symmetric form is used instead.
    tr, names, V, depths, patristic = _bm_setup()
    assert not np.allclose(depths, depths[0]), "this tree should not be ultrametric"
    C = _shape("OU", 2.0, V, depths, patristic)
    assert np.allclose(C, C.T)
    # its diagonal is the known OU variance sigma^2/(2a) (1 - exp(-2 a T))
    alpha = 2.0
    assert np.allclose(np.diag(C),
                       (1.0 - np.exp(-2.0 * alpha * depths)) / (2.0 * alpha))


@pytest.mark.parametrize("n,seed", [(30, 1), (50, 4), (20, 7)])
def test_two_implementations_of_pagels_lambda_agree(n, seed):
    # pagels_lambda (signal.py) and fit_continuous(model="lambda") (models.py)
    # estimate the same quantity by ML through separate code -- separate
    # covariance construction, separate likelihood, separate optimiser call.
    # Nothing forced them to agree, so this is a real cross-check rather than a
    # tautology, and it is the guard against the two drifting apart later.
    tr = pt.datasets.random_tree(n, seed=seed)
    names, V = pt.phylo_vcv(tr)
    rng = np.random.default_rng(seed)
    trait = dict(zip(names, np.linalg.cholesky(V) @ rng.normal(size=n)))
    a = pt.pagels_lambda(tr, trait)
    b = pt.fit_continuous(tr, trait, "lambda")
    assert b["lambda"] == pytest.approx(a["lambda"], abs=1e-6)
    assert b["logLik"] == pytest.approx(a["logLik"], abs=1e-8)


def test_brownian_log_likelihood_agrees_across_the_two_modules():
    # fit_continuous's BM fit and signal's own profile log-likelihood at
    # lambda = 1 are the same model reached two ways, constants included -- which
    # is what makes the AIC in one comparable with a likelihood ratio in the other
    from phytreon.comparative.signal import _profile_log_lik
    tr, names, V, _, _ = _bm_setup(n=30, seed=1)
    trait = _simulate(names, V, 1)
    y = np.array([trait[n] for n in names])
    direct = _profile_log_lik(V, np.ones((len(names), 1)), y, 1.0)
    assert pt.fit_continuous(tr, trait, "BM")["logLik"] == pytest.approx(direct)


def test_shape_rejects_an_unknown_model():
    tr, names, V, depths, patristic = _bm_setup(n=10)
    with pytest.raises(ValueError, match="unknown model"):
        _shape("nope", 1.0, V, depths, patristic)
    with pytest.raises(ValueError, match="unknown model"):
        pt.fit_continuous(tr, _simulate(names, V, 0), "nope")


# --------------------------------------------------------------------------
# Model choice
# --------------------------------------------------------------------------
@pytest.mark.parametrize("model,param", [
    ("BM", 0.0),
    ("OU", 3.0),
    ("EB", -3.0),
    ("white", 0.0),
])
def test_aicc_recovers_the_model_that_generated_the_data(model, param):
    # the property the whole module exists for. Not every replicate: a parameter
    # value near the Brownian boundary is genuinely indistinguishable from BM, so
    # OU and EB legitimately lose some replicates *to BM specifically*. Measured
    # 18-23 of 25 for each model; the threshold is set well inside that.
    tr, names, V, depths, patristic = _bm_setup()
    scale = param / depths.max() if model == "EB" else param
    wins = 0
    for s in range(12):
        trait = _simulate(names, _shape(model, scale, V, depths, patristic), s)
        if pt.compare_continuous_models(tr, trait).index[0] == model:
            wins += 1
    assert wins >= 7, f"{model} recovered only {wins}/12 times"


def test_akaike_weights_sum_to_one_and_the_table_is_ranked():
    tr, names, V, depths, patristic = _bm_setup()
    df = pt.compare_continuous_models(tr, _simulate(names, V, 0))
    assert set(df.index) == set(pt.CONTINUOUS_MODELS)
    assert df["weight"].sum() == pytest.approx(1.0)
    assert df["AICc"].is_monotonic_increasing
    assert df["delta"].iloc[0] == 0.0
    assert (df["delta"] >= 0).all()


def _ultrametric_tree(n=30, seed=0):
    """A UPGMA tree, which is ultrametric by construction."""
    rng = np.random.default_rng(seed)
    names = [f"t{i}" for i in range(n)]
    coords = rng.normal(size=(n, 3))
    D = [[float(np.linalg.norm(coords[i] - coords[j])) for j in range(n)]
         for i in range(n)]
    tr = pt.upgma(names, D)
    _, V = pt.phylo_vcv(tr, names)
    assert np.allclose(np.diag(V), np.diag(V)[0]), "UPGMA must be ultrametric"
    return tr, names


def test_white_noise_wins_when_the_trait_ignores_the_tree():
    # A trait drawn iid across tips has no phylogenetic structure, so 'white'
    # should win -- but only on an ultrametric tree, which is why one is built
    # here. See the next test for why that qualifier is not pedantry.
    tr, names = _ultrametric_tree()
    rng = np.random.default_rng(1)
    trait = dict(zip(names, rng.normal(size=len(names))))
    df = pt.compare_continuous_models(tr, trait)
    assert df.index[0] == "white"
    assert df.loc["white", "delta"] == 0.0
    # lambda=0 and alpha->infinity both reach white noise on an ultrametric
    # tree, so those two tie with it on likelihood and lose only on parameters
    assert df.loc["lambda", "logLik"] == pytest.approx(df.loc["white", "logLik"])


def test_lambda_zero_is_not_white_noise_on_a_non_ultrametric_tree():
    # lambda=0 zeroes the *shared* history but keeps each tip's own variance --
    # the original diagonal of phylo_vcv, i.e. its root-to-tip depth. That is
    # white noise only if every tip is equally deep. On a tree whose depths span
    # 12x it is a *weighted* model instead, and on iid data it genuinely fits
    # better than equal-variance white noise and earns its extra parameter.
    #
    # The same fact bit pgls(lambda_=0.0) earlier and is worth pinning in both
    # places: it is a property of the parameterisation, not a bug, and the only
    # way a reader learns it is if something says so.
    tr, names, V, depths, patristic = _bm_setup()
    assert depths.max() / depths.min() > 5.0, "tree should be strongly non-ultrametric"
    rng = np.random.default_rng(1)
    trait = dict(zip(names, rng.normal(size=len(names))))
    y = np.array([trait[n] for n in names])
    ll_white = _log_lik(_shape("white", 0.0, V, depths, patristic), y)
    ll_lambda0 = _log_lik(_shape("lambda", 0.0, V, depths, patristic), y)
    assert ll_lambda0 != pytest.approx(ll_white)
    assert ll_lambda0 > ll_white
    # ... whereas on an ultrametric tree the two coincide exactly
    utr, unames = _ultrametric_tree()
    _, uV = pt.phylo_vcv(utr, unames)
    ud = np.diag(uV).copy()
    upat = ud[:, None] + ud[None, :] - 2.0 * uV
    uy = rng.normal(size=len(unames))
    assert _log_lik(_shape("lambda", 0.0, uV, ud, upat), uy) == pytest.approx(
        _log_lik(_shape("white", 0.0, uV, ud, upat), uy))


def test_ou_half_life_is_log_two_over_alpha():
    # the interpretable form of alpha, in the tree's own branch-length units
    tr, names, V, depths, patristic = _bm_setup()
    trait = _simulate(names, _shape("OU", 5.0, V, depths, patristic), 0)
    fit = pt.fit_continuous(tr, trait, "OU")
    assert fit["half_life"] == pytest.approx(np.log(2.0) / fit["alpha"])


def test_aicc_penalises_more_than_aic_and_converges_on_it_with_more_taxa():
    trait_small = None
    gaps = []
    for n in (12, 120):
        tr, names, V, depths, patristic = _bm_setup(n=n, seed=3)
        trait_small = _simulate(names, V, 0)
        fit = pt.fit_continuous(tr, trait_small, "OU")
        assert fit["AICc"] > fit["AIC"]
        gaps.append(fit["AICc"] - fit["AIC"])
    assert gaps[1] < gaps[0]


ZERO_LENGTH_TREE = "((A:0.0,B:0.0):1.0,(C:1.0,D:1.0):1.0,E:2.0);"
ZERO_LENGTH_TRAIT = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 1.5, "E": 2.5}


def test_the_singular_tree_guard_covers_every_path_that_inverts_the_covariance():
    # blomberg_k, pagels_lambda and pgls got this guard earlier; the models added
    # afterwards bypassed it and fell back to numpy's bare "Singular matrix",
    # which says nothing about which tips or what to do.
    tr = pt.Tree.from_newick(ZERO_LENGTH_TREE)
    import pandas as pd
    frame = pd.DataFrame({"a": list(ZERO_LENGTH_TRAIT.values()),
                          "b": [1.0, 2.0, 3.0, 4.0, 5.0]},
                         index=list(ZERO_LENGTH_TRAIT))
    for call in (lambda: pt.fit_continuous(tr, ZERO_LENGTH_TRAIT, "BM"),
                 lambda: pt.fit_continuous(tr, ZERO_LENGTH_TRAIT, "OU"),
                 lambda: pt.fit_continuous(tr, ZERO_LENGTH_TRAIT, "EB"),
                 lambda: pt.fit_continuous(tr, ZERO_LENGTH_TRAIT, "lambda"),
                 lambda: pt.compare_continuous_models(tr, ZERO_LENGTH_TRAIT),
                 lambda: pt.phylo_pca(tr, frame)):
        with pytest.raises(ValueError, match="zero distance from each other"):
            call()
        try:
            call()
        except ValueError as exc:
            assert "'A'" in str(exc) and "'B'" in str(exc)


def test_the_guard_is_not_applied_where_the_covariance_is_never_inverted():
    # The guard belongs where the matrix is INVERTED, not merely where it is
    # used. White noise's covariance is the identity -- it ignores the tree, so a
    # singular tree cannot hurt it. fritz_purvis_d only *simulates from* the
    # covariance, which is well defined when singular: two tips at zero distance
    # simply come out perfectly correlated, which is what the tree says of them.
    # Guarding either would reject a case that is genuinely computable.
    tr = pt.Tree.from_newick(ZERO_LENGTH_TREE)
    white = pt.fit_continuous(tr, ZERO_LENGTH_TRAIT, "white")
    assert np.isfinite(white["logLik"])
    assert np.isfinite(white["AIC"])

    binary = {"A": 1, "B": 1, "C": 0, "D": 0, "E": 1}
    res = pt.fritz_purvis_d(tr, binary, n_sim=99, seed=0)
    assert np.isfinite(res["D"])

    # and the perfect correlation that makes that legitimate
    names, V = pt.phylo_vcv(tr, list(ZERO_LENGTH_TRAIT))
    chol = np.linalg.cholesky(V + 1e-12 * np.eye(len(names)))
    rng = np.random.default_rng(0)
    draws = np.array([chol @ rng.normal(size=len(names)) for _ in range(200)])
    ia, ib = names.index("A"), names.index("B")
    assert np.corrcoef(draws[:, ia], draws[:, ib])[0, 1] == pytest.approx(1.0)


def test_fit_continuous_needs_at_least_three_taxa_and_a_known_criterion():
    tr = pt.Tree.from_newick("(A:1,B:1);")
    with pytest.raises(ValueError, match="at least 3 taxa"):
        pt.fit_continuous(tr, {"A": 1.0, "B": 2.0}, "BM")
    tr2, names, V, _, _ = _bm_setup(n=10)
    with pytest.raises(ValueError, match="AIC"):
        pt.compare_continuous_models(tr2, _simulate(names, V, 0), criterion="BIC")


# --------------------------------------------------------------------------
# Phylogenetic PCA
# --------------------------------------------------------------------------
def test_phylo_pca_reduces_to_ordinary_pca_on_a_star_tree():
    # the exact reference case: a star tree has no shared history, so there is
    # nothing to correct for and the phylogenetic estimate must coincide with
    # the ordinary one -- mean, loadings and explained variance alike
    import pandas as pd
    n, p = 30, 4
    names = [f"t{i}" for i in range(n)]
    star = pt.Tree.from_newick("(" + ",".join(f"{t}:1.0" for t in names) + ");")
    rng = np.random.default_rng(0)
    X = pd.DataFrame(rng.normal(size=(n, p)), index=names,
                     columns=[f"tr{j}" for j in range(p)])
    res = pt.phylo_pca(star, X)

    values = X.to_numpy()
    centred = values - values.mean(axis=0)
    S = centred.T @ centred / (n - 1)
    ev, evec = np.linalg.eigh(S)
    order = np.argsort(ev)[::-1]
    ev = ev[order]
    assert np.allclose(res["phylo_mean"].to_numpy(), values.mean(axis=0))
    assert np.allclose(res["explained"].to_numpy(), ev / ev.sum())
    assert np.allclose(np.abs(res["loadings"].to_numpy()),
                       np.abs(evec[:, order]))


def test_phylo_pca_differs_from_ordinary_pca_when_the_tree_has_structure():
    # the reason the function exists: with shared history, the GLS mean is not
    # the column mean, so the centring -- and everything after it -- changes
    import pandas as pd
    n, p = 30, 4
    tr = pt.datasets.random_tree(n, seed=5)
    tips = tr.leaf_names()
    _, V = pt.phylo_vcv(tr, tips)
    chol = np.linalg.cholesky(V + 1e-12 * np.eye(n))
    rng = np.random.default_rng(0)
    X = pd.DataFrame(np.column_stack([chol @ rng.normal(size=n) for _ in range(p)]),
                     index=tips, columns=[f"tr{j}" for j in range(p)])
    res = pt.phylo_pca(tr, X)
    assert not np.allclose(res["phylo_mean"].to_numpy(),
                           X.to_numpy().mean(axis=0), atol=1e-6)


def test_phylo_pca_structural_properties():
    import pandas as pd
    n, p = 30, 4
    tr = pt.datasets.random_tree(n, seed=5)
    tips = tr.leaf_names()
    rng = np.random.default_rng(2)
    X = pd.DataFrame(rng.normal(size=(n, p)), index=tips,
                     columns=[f"tr{j}" for j in range(p)])
    for mode in ("cov", "corr"):
        res = pt.phylo_pca(tr, X, mode=mode)
        assert res["explained"].sum() == pytest.approx(1.0)
        assert np.all(np.diff(res["eigenvalues"].to_numpy()) <= 1e-12)
        loadings = res["loadings"].to_numpy()
        assert np.allclose(loadings.T @ loadings, np.eye(p))
        assert list(res["scores"].index) == tips
        assert list(res["loadings"].index) == list(X.columns)
        assert res["mode"] == mode


def test_phylo_pca_axis_signs_are_stable_across_runs():
    # eigenvectors are only defined up to sign, so without a convention the same
    # data can produce mirrored axes on different runs and confuse a biplot
    import pandas as pd
    tr = pt.datasets.random_tree(25, seed=6)
    tips = tr.leaf_names()
    rng = np.random.default_rng(3)
    X = pd.DataFrame(rng.normal(size=(25, 3)), index=tips, columns=list("abc"))
    a = pt.phylo_pca(tr, X)["loadings"].to_numpy()
    b = pt.phylo_pca(tr, X)["loadings"].to_numpy()
    assert np.allclose(a, b)
    for j in range(a.shape[1]):
        assert a[np.argmax(np.abs(a[:, j])), j] > 0


def test_phylo_pca_rejects_a_bad_mode_and_too_few_taxa():
    import pandas as pd
    tr = pt.datasets.random_tree(10, seed=1)
    tips = tr.leaf_names()
    X = pd.DataFrame(np.zeros((10, 2)), index=tips, columns=["a", "b"])
    with pytest.raises(ValueError, match="mode must be"):
        pt.phylo_pca(tr, X, mode="nope")
    small = pt.Tree.from_newick("(A:1,B:1);")
    with pytest.raises(ValueError, match="at least 3 taxa"):
        pt.phylo_pca(small, pd.DataFrame(np.zeros((2, 2)), index=["A", "B"],
                                         columns=["a", "b"]))


# --------------------------------------------------------------------------
# Fritz & Purvis' D
# --------------------------------------------------------------------------
def test_sister_clade_differences_of_a_basal_split_is_exactly_one():
    # the reference case that fixes the statistic's scale: nodal values are 1
    # throughout one clade, 0 throughout the other and 0.5 at the root, so only
    # the root's two edges contribute, 0.5 each
    from phytreon.comparative.signal import _sister_clade_differences
    tr = pt.Tree.from_newick(
        "(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);")
    trait = {t: (1.0 if t in "ABCD" else 0.0) for t in "ABCDEFGH"}
    assert _sister_clade_differences(tr, trait) == pytest.approx(1.0)


def test_d_is_near_zero_for_a_brownian_threshold_trait():
    # D's own definition: 0 is "as clumped as a Brownian trait pushed through a
    # threshold", so traits built exactly that way must average there. This
    # tests the scaling, not just that the function runs.
    tr = pt.datasets.random_tree(60, seed=4)
    names, V = pt.phylo_vcv(tr)
    chol = np.linalg.cholesky(V + 1e-12 * np.eye(60))
    rng = np.random.default_rng(0)
    ds = []
    for s in range(6):
        cont = chol @ rng.normal(size=60)
        cut = np.sort(cont)[::-1][19]              # top 20 tips get state 1
        trait = {n: float(v >= cut) for n, v in zip(names, cont)}
        ds.append(pt.fritz_purvis_d(tr, trait, n_sim=149, seed=s)["D"])
    assert abs(float(np.mean(ds))) < 0.4, f"mean D {np.mean(ds):.3f}"


def test_d_is_near_one_for_a_randomly_scattered_trait():
    tr = pt.datasets.random_tree(60, seed=4)
    names = tr.leaf_names()
    rng = np.random.default_rng(1)
    ds = []
    for s in range(6):
        vals = np.array([1.0] * 20 + [0.0] * 40)
        rng.shuffle(vals)
        ds.append(pt.fritz_purvis_d(tr, dict(zip(names, vals)),
                                    n_sim=149, seed=s)["D"])
    assert abs(float(np.mean(ds)) - 1.0) < 0.4, f"mean D {np.mean(ds):.3f}"


def test_d_goes_below_zero_for_a_trait_confined_to_one_clade():
    tr = pt.datasets.random_tree(60, seed=4)
    names = tr.leaf_names()
    clade = next(nd.leaf_names() for nd in tr.traverse("postorder")
                 if not nd.is_leaf and 15 <= len(nd.leaf_names()) <= 25)
    trait = {n: float(n in set(clade)) for n in names}
    res = pt.fritz_purvis_d(tr, trait, n_sim=299, seed=0)
    assert res["D"] < 0.5
    assert res["p_random"] < 0.05        # clearly not randomly scattered
    assert res["n_ones"] == len(clade)
    assert res["mean_random"] > res["mean_brownian"]


def test_d_reports_both_references_separately():
    # "different from random" and "different from Brownian" are separate claims
    # and a trait is often the first without being the second
    tr = pt.datasets.random_tree(40, seed=7)
    names, V = pt.phylo_vcv(tr)
    chol = np.linalg.cholesky(V + 1e-12 * np.eye(40))
    rng = np.random.default_rng(0)
    cont = chol @ rng.normal(size=40)
    cut = np.sort(cont)[::-1][14]
    trait = {n: float(v >= cut) for n, v in zip(names, cont)}
    res = pt.fritz_purvis_d(tr, trait, n_sim=199, seed=0)
    assert 0.0 < res["p_random"] <= 1.0
    assert 0.0 < res["p_brownian"] <= 1.0
    assert set(res) >= {"D", "observed", "mean_random", "mean_brownian",
                        "p_random", "p_brownian", "n", "n_ones", "n_sim"}


def test_d_rejects_data_it_cannot_use():
    tr = pt.datasets.random_tree(12, seed=1)
    names = tr.leaf_names()
    ok = {n: float(i % 2) for i, n in enumerate(names)}
    # every tip needs a value -- an unlabelled tip has no state to be clumped
    with pytest.raises(ValueError, match="value for every tip"):
        pt.fritz_purvis_d(tr, {k: v for k, v in list(ok.items())[:-1]}, n_sim=9)
    with pytest.raises(ValueError, match="not tips of the tree"):
        pt.fritz_purvis_d(tr, {**ok, "Ghost": 1.0}, n_sim=9)
    with pytest.raises(ValueError, match="binary"):
        pt.fritz_purvis_d(tr, {n: float(i) for i, n in enumerate(names)}, n_sim=9)
    with pytest.raises(ValueError, match="both states"):
        pt.fritz_purvis_d(tr, {n: 1.0 for n in names}, n_sim=9)


def test_d_accepts_booleans_and_ints():
    tr = pt.datasets.random_tree(20, seed=2)
    names = tr.leaf_names()
    as_bool = {n: (i % 3 == 0) for i, n in enumerate(names)}
    as_int = {n: int(i % 3 == 0) for i, n in enumerate(names)}
    a = pt.fritz_purvis_d(tr, as_bool, n_sim=99, seed=0)
    b = pt.fritz_purvis_d(tr, as_int, n_sim=99, seed=0)
    assert a["D"] == pytest.approx(b["D"])
    assert a["observed"] == pytest.approx(b["observed"])
