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
        obs_raw = k_stat(x)
        result["p"] = float(np.mean(raw_ratios >= obs_raw))
        result["n_perm"] = n_perm
    return result


def _lambda_transform(V, lam: float):
    import numpy as np
    Vl = V * lam
    np.fill_diagonal(Vl, np.diag(V))
    return Vl


def _profile_neg_ll(V, x, lam: float) -> float:
    import numpy as np
    n = len(x)
    Vl = _lambda_transform(V, lam)
    try:
        V_inv = np.linalg.inv(Vl)
    except np.linalg.LinAlgError:
        return 1e18
    sign, logdet = np.linalg.slogdet(Vl)
    if sign <= 0:
        return 1e18
    a_hat, _ = _gls_mean(V_inv, x)
    resid = x - a_hat
    ssq = float(resid @ V_inv @ resid)
    if ssq <= 0:
        return 1e18
    sigma2 = ssq / n
    return 0.5 * (n * np.log(sigma2) + logdet)   # + const, dropped for optimisation


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

    fit = minimize_scalar(lambda lam: _profile_neg_ll(V, x, lam),
                          bounds=(0.0, 1.0), method="bounded")
    lam = float(fit.x)
    ll_lam = -float(fit.fun)
    ll_0 = -_profile_neg_ll(V, x, 0.0)
    lr = 2.0 * (ll_lam - ll_0)
    p = float(chi2.sf(max(lr, 0.0), df=1))
    return {"lambda": lam, "logLik": ll_lam, "logLik0": ll_0,
           "LR": max(lr, 0.0), "p": p, "n": n}


def pgls(tree: Tree, y: Dict[str, float], x: Union[Dict[str, float], "pd.DataFrame"],  # noqa: F821
        lambda_: Union[float, str] = "ML") -> Dict[str, object]:
    """Phylogenetic generalised least squares: regress ``y`` on one or more
    predictors, using the tree's own covariance structure as the error term
    instead of assuming every tip is an independent observation.

    ``x`` is a single ``{tip: value}`` mapping (one predictor) or a
    :class:`pandas.DataFrame` indexed by tip name, one column per predictor
    (several at once); an intercept is always added. ``lambda_="ML"``
    (default) estimates :func:`pagels_lambda` on the *residuals'* structure
    jointly with the regression, so the error covariance is not forced to be
    a pure, untransformed Brownian-motion tree when the data do not support
    that; pass a fixed number (``1.0`` = untransformed BM-GLS, ``0.0`` =
    an ordinary, non-phylogenetic least-squares fit) to skip that estimation.

    Returns coefficients, standard errors, t-values and p-values (one row
    per predictor plus the intercept), R^2, the fitted ``lambda``, and the
    number of taxa used.
    """
    import numpy as np
    from scipy.optimize import minimize_scalar
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
    n = len(names)
    Y = np.array([y[t] for t in names], dtype=float)
    X = np.column_stack([np.ones(n)] + [
        [x_by_tip[t][i] for t in names] for i in range(len(x_names))
    ])

    def fit_at(lam: float):
        Vl = _lambda_transform(V, lam) if lam != 1.0 else V
        V_inv = np.linalg.inv(Vl)
        XtVinv = X.T @ V_inv
        cov = np.linalg.inv(XtVinv @ X)
        beta = cov @ XtVinv @ Y
        resid = Y - X @ beta
        dof = n - X.shape[1]
        sigma2 = float(resid @ V_inv @ resid) / dof
        se = np.sqrt(np.diag(cov) * sigma2)
        return beta, se, sigma2, dof, resid, V_inv

    if lambda_ == "ML":
        def neg_ll(lam):
            Vl = _lambda_transform(V, lam)
            try:
                V_inv = np.linalg.inv(Vl)
            except np.linalg.LinAlgError:
                return 1e18
            sign, logdet = np.linalg.slogdet(Vl)
            if sign <= 0:
                return 1e18
            XtVinv = X.T @ V_inv
            beta = np.linalg.solve(XtVinv @ X, XtVinv @ Y)
            resid = Y - X @ beta
            ssq = float(resid @ V_inv @ resid)
            if ssq <= 0:
                return 1e18
            sigma2 = ssq / n
            return 0.5 * (n * np.log(sigma2) + logdet)
        lam = float(minimize_scalar(neg_ll, bounds=(0.0, 1.0), method="bounded").x)
    else:
        lam = float(lambda_)

    beta, se, sigma2, dof, resid, V_inv = fit_at(lam)
    t_vals = beta / se
    p_vals = 2.0 * t_dist.sf(np.abs(t_vals), df=dof)

    y_mean_gls, _ = _gls_mean(V_inv, Y)
    ss_res = float(resid @ V_inv @ resid)
    ss_tot = float((Y - y_mean_gls) @ V_inv @ (Y - y_mean_gls))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    coef_names = ["Intercept"] + x_names
    return {
        "coefficients": dict(zip(coef_names, beta.tolist())),
        "se": dict(zip(coef_names, se.tolist())),
        "t": dict(zip(coef_names, t_vals.tolist())),
        "p": dict(zip(coef_names, p_vals.tolist())),
        "lambda": lam,
        "r2": r2,
        "n": n,
        "dof": dof,
    }
