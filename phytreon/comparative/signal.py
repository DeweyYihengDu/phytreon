"""Phylogenetic signal (Blomberg's K, Pagel's lambda) and PGLS.

All three start from the same object: the Brownian-motion expected
variance-covariance matrix of the tips (:func:`phylo_vcv`, ape's
``vcv.phylo``) -- tip i's variance is its root-to-tip depth, and two tips'
covariance is the depth of their shared history (root to MRCA). Phylogenetic
signal asks how well a trait's *actual* covariance across tips matches that
expectation; PGLS uses the same matrix as the error structure for a
regression, so that two traits' correlation is not inflated by the fact
that close relatives are not independent data points.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple, Union

from ..core.tree import Node, Tree


def phylo_vcv(tree: Tree, taxa: Optional[Sequence[str]] = None
             ) -> Tuple[List[str], "np.ndarray"]:  # noqa: F821
    """Brownian-motion expected variance-covariance matrix of the tips.

    ``V[i, i]`` is tip i's root-to-tip depth (its variance under BM, run for
    that long from a fixed root value); ``V[i, j]`` is the depth of the tips'
    shared history -- the path length from the root to their MRCA, since
    that is the portion of each tip's random walk the two share and
    everything after their split is independent. Returns ``(names, V)`` with
    ``names`` in the same order as ``V``'s rows/columns; pass ``taxa`` to fix
    that order (and to work on a subset) instead of the tree's own leaf order.
    """
    import numpy as np
    leaves = tree.leaves()
    if taxa is not None:
        want = list(taxa)
        by_name = {lf.name: lf for lf in leaves}
        missing = set(want) - set(by_name)
        if missing:
            raise ValueError(f"taxa not found in tree: {sorted(missing)}")
        leaves = [by_name[n] for n in want]
    names = [lf.name for lf in leaves]
    n = len(names)
    # root -> leaf node lists, so two tips' shared depth is just their
    # paths' common prefix length
    paths: List[List[Node]] = []
    for lf in leaves:
        path, node = [], lf
        while node is not None:
            path.append(node)
            node = node.parent
        paths.append(list(reversed(path)))
    depths = [sum(node.length or 0.0 for node in path[1:]) for path in paths]
    V = np.zeros((n, n))
    for i in range(n):
        V[i, i] = depths[i]
        for j in range(i + 1, n):
            shared = 0.0
            for a, b in zip(paths[i][1:], paths[j][1:]):
                if a is not b:
                    break
                shared += a.length or 0.0
            V[i, j] = V[j, i] = shared
    return names, V


def _gls_mean(V_inv, x):
    import numpy as np
    ones = np.ones(len(x))
    denom = ones @ V_inv @ ones
    return float(ones @ V_inv @ x / denom), float(denom)


def _check_invertible(names, V, caller: str) -> None:
    """Refuse a singular covariance matrix, saying which tips caused it.

    Not an exotic input: any two tips at zero patristic distance -- which is
    what a zero-length terminal branch produces, and IQ-TREE and RAxML emit
    those routinely -- have identical rows in ``V``, so it cannot be inverted.
    Nothing downstream can separate such tips, and left to itself each of these
    functions failed differently and unhelpfully: a bare ``LinAlgError`` from
    numpy, or (worse) an optimiser quietly escaping to ``lambda = 0``, where
    the off-diagonals vanish and the matrix inverts again -- reporting "no
    phylogenetic signal, p = 1.0" for what is really "not computable on this
    tree".
    """
    import numpy as np
    n = len(names)
    try:
        np.linalg.cholesky(V)
        return          # positive definite: nothing to report, and this is the
                        # common case, so it is the cheap check that runs first
    except np.linalg.LinAlgError:
        pass            # fall through and work out what to tell the caller
    flat = [i for i in range(n) if V[i, i] <= 0.0]
    if flat:
        raise ValueError(
            f"{caller}: {sorted(names[i] for i in flat)} sit at zero distance "
            f"from the root, so the tree implies no variance for them"
        )
    # patristic distance from the covariance: d_ij = V_ii + V_jj - 2 V_ij
    groups, seen = [], set()
    for i in range(n):
        if i in seen:
            continue
        tied = [j for j in range(i + 1, n)
                if abs(V[i, i] + V[j, j] - 2.0 * V[i, j]) <= 0.0]
        if tied:
            seen.update(tied)
            groups.append([names[i]] + [names[j] for j in tied])
    if groups:
        raise ValueError(
            f"{caller}: the tree has tips at zero distance from each other, so "
            f"its covariance matrix is singular and nothing can tell them "
            f"apart: {groups}. Give the branches between them a length, or "
            f"drop all but one tip from each group."
        )
    raise ValueError(
        f"{caller}: the tree's covariance matrix is singular (rank "
        f"{np.linalg.matrix_rank(V)} of {n}), so it cannot be inverted"
    )


def blomberg_k(tree: Tree, trait: Dict[str, float], n_perm: int = 0,
              seed: Optional[int] = None) -> Dict[str, object]:
    """Blomberg et al. (2003)'s K: how strongly a continuous trait's
    similarity across tips tracks the tree, relative to what Brownian motion
    predicts for this exact topology and these exact branch lengths.

    ``K = 1`` matches Brownian motion; ``K > 1`` means close relatives
    resemble each other *more* than that (e.g. strong stabilising selection
    within a clade); ``K < 1`` means less (convergence, measurement error,
    or simply a trait that is not tracking the tree). Computed as the ratio
    of two mean squared errors around the GLS-estimated phylogenetic mean --
    one plain, one weighted by the inverse of :func:`phylo_vcv` -- scaled by
    that ratio's own expected value under Brownian motion on this tree,
    which is what makes K comparable across different trees and trait sets
    rather than only meaningful relative to a null distribution (Blomberg,
    Garland & Ives 2003, *Evolution* 57:717-745, eq. 1-3).

    ``n_perm > 0`` adds a permutation test (shuffling the trait across tips,
    the original paper's own significance procedure): ``p`` is the fraction
    of shuffles whose K is at least as large as the observed one.
    """
    import numpy as np
    # phylo_vcv(tree, list(trait)) itself rejects any name in trait that is
    # not a real tip, so nothing further is needed here to catch that case.
    names, V = phylo_vcv(tree, list(trait))
    x = np.array([trait[n] for n in names], dtype=float)
    n = len(names)
    if n < 3:
        raise ValueError("blomberg_k needs at least 3 taxa with trait values")
    _check_invertible(names, V, "blomberg_k")
    V_inv = np.linalg.inv(V)
    ones = np.ones(n)

    def k_stat(values):
        a_hat, denom = _gls_mean(V_inv, values)
        resid = values - a_hat
        mse0 = float(resid @ resid) / (n - 1)
        mse = float(resid @ V_inv @ resid) / (n - 1)
        return mse0 / mse if mse > 0 else float("inf")

    observed = k_stat(x)
    expected_ratio = (float(np.trace(V)) - n / float(ones @ V_inv @ ones)) / (n - 1)
    k = observed / expected_ratio

    result: Dict[str, object] = {"K": k, "n": n}
    if n_perm > 0:
        rng = np.random.default_rng(seed)
        raw_ratios = np.array([k_stat(rng.permutation(x)) for _ in range(n_perm)])
        # compared on the raw ratio rather than K: the two differ only by
        # expected_ratio, which is a property of the tree and so identical for
        # every permutation, and dividing both sides by it changes nothing
        result["p"] = float(np.mean(raw_ratios >= observed))
        result["n_perm"] = n_perm
    return result


def _lambda_transform(V, lam: float):
    import numpy as np
    Vl = V * lam
    np.fill_diagonal(Vl, np.diag(V))
    return Vl


def _gls_profile_crit(V, X, Y, lam: float, reml: bool) -> float:
    """Profile (concentrated) criterion for a GLS fit at a fixed ``lam``, with
    both ``beta`` and ``sigma^2`` substituted by their conditional optima so
    only ``lam`` is left to optimise over. Constants are dropped -- this is for
    minimising, not for reporting; see :func:`_profile_log_lik`.

    ``reml=False`` is the ordinary ML criterion. ``reml=True`` is the
    restricted one: it divides the residual sum of squares by ``n - p`` instead
    of ``n`` and adds ``log|X' V^-1 X|``, correcting for the fact that the mean
    structure was estimated from the same data -- without which ML pulls the
    variance component (here ``lam``) systematically towards zero in small
    samples, and a too-small ``lam`` understates how non-independent the tips
    are and so understates the standard errors.
    """
    import numpy as np
    n, p = X.shape
    Vl = _lambda_transform(V, lam)
    try:
        V_inv = np.linalg.inv(Vl)
    except np.linalg.LinAlgError:
        return 1e18
    sign, logdet = np.linalg.slogdet(Vl)
    if sign <= 0:
        return 1e18
    XtVinv = X.T @ V_inv
    A = XtVinv @ X
    try:
        beta = np.linalg.solve(A, XtVinv @ Y)
    except np.linalg.LinAlgError:
        return 1e18
    resid = Y - X @ beta
    ssq = float(resid @ V_inv @ resid)
    if ssq <= 0:
        return 1e18
    if reml:
        sign_a, logdet_a = np.linalg.slogdet(A)
        if sign_a <= 0:
            return 1e18
        return 0.5 * ((n - p) * np.log(ssq / (n - p)) + logdet + logdet_a)
    return 0.5 * (n * np.log(ssq / n) + logdet)


def _profile_log_lik(V, X, Y, lam: float) -> float:
    """The ML criterion above turned back into an actual log-likelihood, by
    restoring the constants dropped for optimisation, so the reported value is
    comparable with other software's (and usable for AIC) rather than only
    valid up to an additive constant that cancels inside a likelihood ratio.
    """
    import numpy as np
    n = X.shape[0]
    crit = _gls_profile_crit(V, X, Y, lam, reml=False)
    if crit >= 1e17:
        return -float("inf")
    return -crit - 0.5 * n * (np.log(2.0 * np.pi) + 1.0)


def pagels_lambda(tree: Tree, trait: Dict[str, float]) -> Dict[str, object]:
    """Pagel (1999)'s lambda: the multiplier on the tree's internal (shared)
    branch lengths that best fits a continuous trait's covariance pattern,
    found by maximum likelihood.

    ``lambda = 1`` leaves the tree as given (full Brownian motion);
    ``lambda = 0`` shrinks every shared branch to nothing, i.e. a star
    phylogeny with no phylogenetic signal at all; in between, relatives
    resemble each other less than the tree alone would predict. Reported
    alongside a likelihood-ratio test against ``lambda = 0`` (`chi2`, 1 df) --
    more robust than :func:`blomberg_k` to polytomies and uncertain branch
    lengths (Freckleton, Harvey & Pagel 2002), so prefer it when either is a
    concern for this particular tree.
    """
    import numpy as np
    from scipy.optimize import minimize_scalar
    from scipy.stats import chi2

    # phylo_vcv(tree, list(trait)) itself rejects any name in trait that is
    # not a real tip, so nothing further is needed here to catch that case.
    names, V = phylo_vcv(tree, list(trait))
    x = np.array([trait[n] for n in names], dtype=float)
    n = len(names)
    if n < 3:
        raise ValueError("pagels_lambda needs at least 3 taxa with trait values")
    _check_invertible(names, V, "pagels_lambda")

    X = np.ones((n, 1))          # intercept only: the phylogenetic mean
    fit = minimize_scalar(lambda lam: _gls_profile_crit(V, X, x, lam, reml=False),
                          bounds=(0.0, 1.0), method="bounded")
    lam = float(fit.x)
    ll_lam = _profile_log_lik(V, X, x, lam)
    ll_0 = _profile_log_lik(V, X, x, 0.0)
    lr = 2.0 * (ll_lam - ll_0)
    # lambda = 0 sits on the boundary of [0, 1], so the null distribution of LR
    # is really a 50:50 mixture of chi2_0 and chi2_1 (Self & Liang 1987) rather
    # than chi2_1. Using plain chi2_1 needs a *larger* LR to reject, i.e. it
    # errs conservative, which is the safe direction for a signal test -- kept
    # deliberately, and matching what phytools' phylosig reports.
    p = float(chi2.sf(max(lr, 0.0), df=1))
    return {"lambda": lam, "logLik": ll_lam, "logLik0": ll_0,
           "LR": max(lr, 0.0), "p": p, "n": n}


def _fit_gls(V, X, Y, lambda_):
    """One GLS fit at a given design matrix, returning everything downstream
    needs. Split out so the bootstrap's null model (this same design minus the
    column being tested) is fitted by the identical estimator rather than a
    parallel reimplementation of it.
    """
    import numpy as np
    from scipy.optimize import minimize_scalar
    n, p = X.shape
    if isinstance(lambda_, str):
        method = lambda_.upper()
        if method not in ("REML", "ML"):
            raise ValueError(
                f"lambda_ must be 'REML', 'ML', or a number, not {lambda_!r}"
            )
        lam = float(minimize_scalar(
            lambda L: _gls_profile_crit(V, X, Y, L, reml=method == "REML"),
            bounds=(0.0, 1.0), method="bounded").x)
        # lambda came out of the same data as the coefficients, so it costs a
        # degree of freedom like any other estimated parameter; not charging it
        # is a large part of why the small-sample false-positive rate was high
        estimated = 1
    else:
        method = "fixed"
        lam = float(lambda_)
        estimated = 0

    V_inv = np.linalg.inv(_lambda_transform(V, lam))
    XtVinv = X.T @ V_inv
    cov = np.linalg.inv(XtVinv @ X)
    beta = cov @ XtVinv @ Y
    resid = Y - X @ beta
    dof = n - p - estimated
    if dof < 1:
        raise ValueError(
            f"pgls has no residual degrees of freedom left ({n} taxa, "
            f"{p} coefficients"
            f"{' and an estimated lambda' if estimated else ''})"
        )
    ss_res = float(resid @ V_inv @ resid)
    sigma2 = ss_res / dof
    se = np.sqrt(np.diag(cov) * sigma2)
    return {"lam": lam, "method": method, "beta": beta, "se": se,
            "t": beta / se, "dof": dof, "sigma2": sigma2, "ss_res": ss_res,
            "resid": resid, "V_inv": V_inv, "fitted": X @ beta}


def pgls(tree: Tree, y: Dict[str, float], x: Union[Dict[str, float], "pd.DataFrame"],  # noqa: F821
        lambda_: Union[float, str] = "REML", n_boot: int = 0,
        seed: Optional[int] = None) -> Dict[str, object]:
    """Phylogenetic generalised least squares: regress ``y`` on one or more
    predictors, using the tree's own covariance structure as the error term
    instead of assuming every tip is an independent observation.

    ``x`` is a single ``{tip: value}`` mapping (one predictor) or a
    :class:`pandas.DataFrame` indexed by tip name, one column per predictor
    (several at once); an intercept is always added.

    ``lambda_`` controls the error covariance, i.e. how much of the tree's
    shared history the residuals are assumed to carry (see
    :func:`pagels_lambda`):

    * ``"REML"`` (default) estimates it by restricted maximum likelihood,
      jointly with the regression.
    * ``"ML"`` estimates it by plain maximum likelihood -- for reproducing
      software that does it that way (``caper::pgls``). Prefer ``"REML"``:
      ML pulls ``lambda`` towards zero in small samples, which understates how
      dependent the tips are and so understates the standard errors (measured
      at 10 tips on data whose real lambda was 1.0: ML averaged 0.60, REML
      0.78, and the false-positive rate over trees of 10-20 taxa was 7.9%
      under ML against 6.8% under REML on the same 6400 datasets, at a
      nominal 5%).
    * a fixed number, skipping estimation entirely -- ``1.0`` for
      untransformed Brownian-motion GLS, ``0.0`` for an ordinary
      non-phylogenetic least-squares fit. Only do this if the value is
      genuinely known in advance: fixing ``lambda`` at the wrong value
      misspecifies the error structure and invalidates the p-values no matter
      how many taxa there are (measured: 12% false positives at 80 tips for
      ``lambda_=1.0`` on data whose real lambda was 0.5, where the estimated
      default gave 5.5%).

    ``n_boot`` adds a parametric bootstrap p-value (``"p_boot"``) for each
    predictor, alongside the t-based one. Worth the cost below roughly 20 taxa:
    the t-test treats the estimated ``lambda`` as if it were known exactly, and
    when it happens to come out too low the standard errors come out too small,
    which leaves a residual ~7% false-positive rate over trees of 10-20 taxa
    that REML alone does not remove. The bootstrap simulates from the fitted
    reduced model and re-estimates ``lambda`` on every replicate, pricing that
    uncertainty in rather than conditioning on one point estimate of it
    (measured over 3200 replicates: 7.3% -> 5.5%, i.e. from six standard errors
    above the nominal 5% to within one and a half of it). Costs ``n_boot`` extra
    fits per predictor; ``seed`` makes it reproducible. Being a count of null draws, ``p_boot`` cannot go below
    ``1 / (n_boot + 1)``, so read the t-based ``p`` instead for effects far
    below that floor -- the bootstrap earns its keep near the decision
    threshold, not out in the tail. Reported for the predictors only, not the
    intercept, whose null (a regression through the origin) is not the
    hypothesis anyone is asking about here.

    Returns coefficients, standard errors, t-values and p-values (one row
    per predictor plus the intercept), R^2, the fitted ``lambda`` and how it
    was obtained, the number of taxa used, and the residual degrees of freedom
    (one lower when ``lambda`` was estimated rather than given, since it is
    then one more parameter read out of the same data).
    """
    import numpy as np
    from scipy.stats import t as t_dist

    if hasattr(x, "columns"):
        x_names = list(x.columns)
        x_by_tip = {tip: [float(x.loc[tip, c]) for c in x_names] for tip in x.index}
    else:
        x_names = ["x"]
        x_by_tip = {tip: [float(v)] for tip, v in x.items()}

    taxa = sorted(set(y) & set(x_by_tip))
    if len(taxa) < len(x_names) + 2:
        raise ValueError(
            f"pgls needs more taxa ({len(taxa)}) than predictors "
            f"({len(x_names)}) plus two for a meaningful fit"
        )
    names, V = phylo_vcv(tree, taxa)
    if isinstance(lambda_, str):
        # only when lambda is being estimated. That is the path where a singular
        # V is dangerous rather than merely fatal: the optimiser can escape to
        # lambda = 0, where the off-diagonals vanish and the matrix inverts
        # again, and report a confident fit that used none of the tree. A
        # lambda_ handed over as a number is the caller asserting the error
        # structure -- including lambda_=0.0, which ignores the tree entirely
        # and is perfectly well defined on a tree this check would reject -- so
        # leave that to np.linalg.inv, which raises loudly if it really is
        # singular at the value asked for.
        _check_invertible(names, V, "pgls")
    n = len(names)
    Y = np.array([y[t] for t in names], dtype=float)
    X = np.column_stack([np.ones(n)] + [
        [x_by_tip[t][i] for t in names] for i in range(len(x_names))
    ])

    full = _fit_gls(V, X, Y, lambda_)
    beta, se, t_vals, dof, V_inv = (full["beta"], full["se"], full["t"],
                                    full["dof"], full["V_inv"])
    p_vals = 2.0 * t_dist.sf(np.abs(t_vals), df=dof)

    y_mean_gls, _ = _gls_mean(V_inv, Y)
    ss_tot = float((Y - y_mean_gls) @ V_inv @ (Y - y_mean_gls))
    r2 = 1.0 - full["ss_res"] / ss_tot if ss_tot > 0 else float("nan")

    coef_names = ["Intercept"] + x_names
    result: Dict[str, object] = {
        "coefficients": dict(zip(coef_names, beta.tolist())),
        "se": dict(zip(coef_names, se.tolist())),
        "t": dict(zip(coef_names, t_vals.tolist())),
        "p": dict(zip(coef_names, p_vals.tolist())),
        "lambda": full["lam"],
        "lambda_method": full["method"],
        "r2": r2,
        "n": n,
        "dof": dof,
    }

    if n_boot > 0:
        rng = np.random.default_rng(seed)
        p_boot: Dict[str, float] = {}
        for j, cname in enumerate(coef_names):
            if j == 0:
                continue        # a no-intercept null is not a hypothesis anyone
                                # is asking about here; t-based p stays for it
            # simulate under "this predictor has no effect" -- the other
            # predictors are kept and refitted, so the null is the reduced
            # model rather than the full one with one coefficient blanked out
            X0 = np.delete(X, j, axis=1)
            null = _fit_gls(V, X0, Y, lambda_)
            chol = np.linalg.cholesky(
                null["sigma2"] * _lambda_transform(V, null["lam"])
                + 1e-12 * np.eye(n)
            )
            t_obs = abs(float(t_vals[j]))
            at_least = 0
            for _ in range(n_boot):
                y_star = null["fitted"] + chol @ rng.normal(size=n)
                at_least += abs(float(_fit_gls(V, X, y_star, lambda_)["t"][j])) >= t_obs
            # the +1s count the observed data as one of its own null draws,
            # which is what keeps the p-value from ever being exactly 0
            p_boot[cname] = (1 + at_least) / (n_boot + 1)
        result["p_boot"] = p_boot
        result["n_boot"] = n_boot

    return result
