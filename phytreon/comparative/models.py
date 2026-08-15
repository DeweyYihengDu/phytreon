"""Models of continuous trait evolution beyond Brownian motion, and
phylogenetic PCA.

:func:`~phytreon.comparative.signal.blomberg_k` and
:func:`~phytreon.comparative.signal.pagels_lambda` ask how well a trait fits
Brownian motion. That is a useful question but a narrow one: a trait can depart
from BM in several distinguishable ways, and "it does not look like BM" does not
say which. This module fits the alternatives and lets the data choose between
them by likelihood.

=========  ==================================================================
model      what it says happened
=========  ==================================================================
``BM``     variance accumulates steadily; a pure random walk, no pull, no trend
``OU``     the trait is pulled back towards an optimum with strength ``alpha``
           -- stabilising selection or a constraint. BM is the ``alpha -> 0``
           limit
``EB``     the rate itself changes exponentially through time at rate ``b``;
           ``b < 0`` is an early burst (fast radiation, then slowdown), ``b > 0``
           a late acceleration. BM is the ``b -> 0`` limit
``lambda`` the tree's shared history is discounted by a factor
           (see :func:`~phytreon.comparative.signal.pagels_lambda`)
``white``  no phylogenetic structure at all -- the null worth keeping in the
           comparison, because a trait that fits it best is one the tree says
           nothing about
=========  ==================================================================

``lambda`` at 0 and ``white`` are **not** the same model unless the tree is
ultrametric, which is worth knowing before reading a comparison table:
``lambda = 0`` removes the shared history but keeps each tip's own variance,
which is its root-to-tip depth, so on a tree whose tips differ in depth it is a
*weighted* model rather than an equal-variance one. On such a tree it can
legitimately beat ``white`` on data with no phylogenetic signal at all, because
the unequal variances are themselves a better description of the data than
equal ones. They coincide exactly when every tip is equally deep.

All five are fitted as the same generalised least squares problem with a
different expected covariance matrix, so they are directly comparable by AIC.
The covariance forms are those of Hansen (1997) / Butler & King (2004) for OU
and Blomberg et al. (2003) / Harmon et al. (2010) for EB.
"""
from __future__ import annotations

from typing import Dict, Sequence

from ..core.tree import Tree
from .signal import _check_invertible, phylo_vcv


MODELS = ("BM", "OU", "EB", "lambda", "white")

# free parameters per model: the phylogenetic mean and sigma^2 are always
# estimated, so BM and white cost 2 and the rest cost one more for their own
_N_PARAMS = {"BM": 2, "OU": 3, "EB": 3, "lambda": 3, "white": 2}


def _shape(model: str, param: float, shared, depths, patristic):
    """The expected tip covariance up to a scale factor.

    ``sigma^2`` is profiled out analytically, so only the *shape* of the matrix
    depends on the model's own parameter and this returns that shape. ``shared``
    is the root-to-MRCA time matrix (:func:`phylo_vcv`, whose diagonal is each
    tip's ``depths``), and ``patristic`` the tip-to-tip distance
    ``T_i + T_j - 2 s_ij``.
    """
    import numpy as np
    if model == "BM":
        return shared
    if model == "white":
        return np.eye(len(depths))
    if model == "lambda":
        out = shared * param
        np.fill_diagonal(out, depths)
        return out
    if model == "OU":
        alpha = param
        if alpha <= 1e-12:
            return shared                     # the BM limit, exactly
        # Symmetric in i and j by construction. Harmon writes this with T_i
        # alone, which is the same thing on an ultrametric tree but is *not*
        # symmetric when tips differ in depth -- and phytreon's trees usually do.
        return (np.exp(-alpha * patristic)
                * (1.0 - np.exp(-2.0 * alpha * shared)) / (2.0 * alpha))
    if model == "EB":
        b = param
        if abs(b) < 1e-12:
            return shared                     # the BM limit, exactly
        return (np.exp(b * shared) - 1.0) / b
    raise ValueError(f"unknown model {model!r}; choose from {MODELS}")


def _log_lik(C, y):
    """Log-likelihood of ``y`` under covariance ``sigma^2 * C``, with the
    phylogenetic mean and ``sigma^2`` at their conditional optima.

    A real log-likelihood, constants included, so it can go into an AIC rather
    than only into a ratio against another fit of the same shape.
    """
    import numpy as np
    n = len(y)
    try:
        C_inv = np.linalg.inv(C)
    except np.linalg.LinAlgError:
        return -np.inf
    sign, logdet = np.linalg.slogdet(C)
    if sign <= 0 or not np.isfinite(logdet):
        return -np.inf
    ones = np.ones(n)
    denom = ones @ C_inv @ ones
    if denom <= 0:
        return -np.inf
    mean = float(ones @ C_inv @ y / denom)
    resid = y - mean
    ssq = float(resid @ C_inv @ resid)
    if ssq <= 0:
        return -np.inf
    sigma2 = ssq / n
    return float(-0.5 * (n * np.log(2.0 * np.pi) + n * np.log(sigma2)
                         + logdet + n))


def fit_continuous(tree: Tree, trait: Dict[str, float], model: str = "BM"
                  ) -> Dict[str, object]:
    """Fit one model of continuous trait evolution by maximum likelihood.

    ``model`` is one of ``"BM"``, ``"OU"``, ``"EB"``, ``"lambda"`` or
    ``"white"`` -- see the module docstring for what each claims. Returns the
    log-likelihood, AIC, small-sample-corrected AICc, the number of free
    parameters, the fitted phylogenetic mean and ``sigma2``, and the model's own
    parameter where it has one (``alpha`` for OU, ``b`` for EB, ``lambda`` for
    lambda).

    Fitted by profiling: the mean and ``sigma2`` have closed-form conditional
    optima, so only the model's own parameter is searched over numerically,
    which is both faster and far less likely to stall than optimising all three
    at once.
    """
    import numpy as np
    from scipy.optimize import minimize_scalar
    if model not in MODELS:
        raise ValueError(f"unknown model {model!r}; choose from {MODELS}")
    names, shared = phylo_vcv(tree, list(trait))
    y = np.array([trait[n] for n in names], dtype=float)
    n = len(names)
    if n < 3:
        raise ValueError("fit_continuous needs at least 3 taxa with trait values")
    if model != "white":
        # Every model except white builds its covariance out of the tree's, so a
        # singular tree covariance (tips at zero distance -- a zero-length
        # terminal branch) makes them unfittable. Checked here rather than left
        # to numpy, which reports only "Singular matrix" from somewhere inside
        # the optimiser. White noise is exempt because its covariance is the
        # identity: it ignores the tree entirely and is well defined regardless.
        _check_invertible(names, shared, f"fit_continuous({model!r})")
    depths = np.diag(shared).copy()
    patristic = depths[:, None] + depths[None, :] - 2.0 * shared

    def at(param):
        return _log_lik(_shape(model, param, shared, depths, patristic), y)

    if model in ("BM", "white"):
        param, ll = None, at(0.0)
    elif model == "lambda":
        fit = minimize_scalar(lambda p: -at(p), bounds=(0.0, 1.0),
                              method="bounded")
        param, ll = float(fit.x), -float(fit.fun)
    elif model == "OU":
        # alpha searched on a log scale: it is a rate of pull whose plausible
        # range spans orders of magnitude, and a linear search wastes almost all
        # its evaluations at the strong-pull end
        hi = max(depths.max(), 1e-9)
        fit = minimize_scalar(lambda lg: -at(np.exp(lg)),
                              bounds=(np.log(1e-6 / hi), np.log(1e4 / hi)),
                              method="bounded")
        param, ll = float(np.exp(fit.x)), -float(fit.fun)
    else:                                     # EB
        hi = max(depths.max(), 1e-9)
        fit = minimize_scalar(lambda p: -at(p), bounds=(-10.0 / hi, 1e-9),
                              method="bounded")
        param, ll = float(fit.x), -float(fit.fun)

    C = _shape(model, 0.0 if param is None else param, shared, depths, patristic)
    C_inv = np.linalg.inv(C)
    ones = np.ones(n)
    mean = float(ones @ C_inv @ y / (ones @ C_inv @ ones))
    resid = y - mean
    sigma2 = float(resid @ C_inv @ resid) / n

    k = _N_PARAMS[model]
    aic = -2.0 * ll + 2.0 * k
    # AICc: the small-sample correction, which matters here because comparative
    # data sets are usually small enough for AIC to under-penalise
    aicc = (aic + (2.0 * k * (k + 1)) / (n - k - 1) if n - k - 1 > 0
            else float("inf"))
    out: Dict[str, object] = {"model": model, "logLik": ll, "AIC": aic,
                              "AICc": aicc, "k": k, "n": n,
                              "sigma2": sigma2, "mean": mean}
    if model == "OU":
        out["alpha"] = param
        # the time for half the pull towards the optimum -- the interpretable
        # form of alpha, in the same units as the tree's branch lengths
        out["half_life"] = float(np.log(2.0) / param) if param > 0 else float("inf")
    elif model == "EB":
        out["b"] = param
    elif model == "lambda":
        out["lambda"] = param
    return out


def compare_continuous_models(tree: Tree, trait: Dict[str, float],
                              models: Sequence[str] = MODELS,
                              criterion: str = "AICc") -> "pd.DataFrame":  # noqa: F821
    """Fit several models to one trait and rank them.

    Returns a DataFrame sorted best-first, with each model's log-likelihood,
    parameter count, AIC and AICc, the difference from the best
    (``delta``), and Akaike weights -- the relative support for each model,
    ``exp(-delta/2)`` normalised to sum to 1.

    ``criterion`` selects which of ``"AIC"``/``"AICc"`` the ranking, ``delta``
    and weights are based on; AICc by default, since comparative data sets are
    usually small enough for the correction to matter and it converges on AIC
    when they are not.

    Read the weights, not just the winner: three models within a couple of AIC
    units of each other means the data does not separate them, which is a real
    result about the data rather than a licence to report the top row.
    """
    import numpy as np
    import pandas as pd
    if criterion not in ("AIC", "AICc"):
        raise ValueError(f"criterion must be 'AIC' or 'AICc', not {criterion!r}")
    rows = [fit_continuous(tree, trait, m) for m in models]
    df = pd.DataFrame(rows).set_index("model")
    score = df[criterion].to_numpy()
    df["delta"] = score - np.nanmin(score)
    w = np.exp(-0.5 * df["delta"].to_numpy())
    df["weight"] = w / np.nansum(w)
    keep = [c for c in ("logLik", "k", "AIC", "AICc", "delta", "weight",
                        "alpha", "half_life", "b", "lambda", "sigma2", "mean",
                        "n") if c in df.columns]
    return df[keep].sort_values(criterion)


# --------------------------------------------------------------------------
# Phylogenetic PCA
# --------------------------------------------------------------------------
def phylo_pca(tree: Tree, data, mode: str = "cov") -> Dict[str, object]:
    """Phylogenetic principal components analysis (Revell 2009).

    Ordinary PCA estimates the trait covariance matrix as if every species were
    an independent sample, which for species on a tree they are not -- close
    relatives resemble each other whatever the traits are doing, so the leading
    axes of a plain PCA can describe the phylogeny rather than the traits. This
    estimates the same covariance matrix with the tree's expected covariance
    divided out first, and projects the species onto its eigenvectors.

    ``data`` is a :class:`pandas.DataFrame` of traits indexed by tip name, one
    column per trait. ``mode="corr"`` standardises the traits first, which is
    what you want when they are in different units and the largest-variance
    trait would otherwise dominate every axis.

    Returns the eigenvalues, the fraction of variance each axis explains, the
    loadings (traits x axes), the species scores (tips x axes), and the
    phylogenetic mean each trait was centred on -- which is a GLS estimate, not
    a plain column mean, and is the step that makes this phylogenetic rather
    than ordinary.
    """
    import numpy as np
    import pandas as pd
    if mode not in ("cov", "corr"):
        raise ValueError(f"mode must be 'cov' or 'corr', not {mode!r}")
    taxa = [str(i) for i in data.index]
    names, C = phylo_vcv(tree, taxa)
    X = data.to_numpy(dtype=float)
    n, p = X.shape
    if n < 3:
        raise ValueError("phylo_pca needs at least 3 taxa")
    # the tree's covariance is inverted below, so a singular one is fatal --
    # reported by which tips caused it rather than as a bare "Singular matrix"
    _check_invertible(names, C, "phylo_pca")
    C_inv = np.linalg.inv(C)
    ones = np.ones((n, 1))

    # the phylogenetic mean of each trait: a GLS estimate, which weights species
    # by how much independent information they carry rather than equally
    a = (np.linalg.inv(ones.T @ C_inv @ ones) @ (ones.T @ C_inv @ X)).ravel()
    centred = X - a
    # the trait covariance with the tree's covariance divided out
    R = (centred.T @ C_inv @ centred) / (n - 1)
    if mode == "corr":
        sd = np.sqrt(np.diag(R))
        R = R / np.outer(sd, sd)
        centred = centred / sd

    values, vectors = np.linalg.eigh(R)
    order = np.argsort(values)[::-1]          # eigh returns ascending
    values, vectors = values[order], vectors[:, order]
    # sign convention: make each axis's largest-magnitude loading positive, so
    # repeated runs do not flip axes arbitrarily
    for j in range(vectors.shape[1]):
        if vectors[np.argmax(np.abs(vectors[:, j])), j] < 0:
            vectors[:, j] *= -1.0

    axes = [f"PC{j + 1}" for j in range(p)]
    total = values.sum()
    return {
        "eigenvalues": pd.Series(values, index=axes),
        "explained": pd.Series(values / total if total > 0 else values * np.nan,
                               index=axes),
        "loadings": pd.DataFrame(vectors, index=list(data.columns), columns=axes),
        "scores": pd.DataFrame(centred @ vectors, index=names, columns=axes),
        "phylo_mean": pd.Series(a, index=list(data.columns)),
        "mode": mode,
    }
