"""Pure-Python maximum-likelihood phylogenetics (no external engine).

Implements the standard ML pipeline for nucleotide data:

* Felsenstein's pruning likelihood, vectorised over compressed site
  patterns (numpy), with rescaling to avoid underflow.
* Time-reversible substitution models: JC69, K80, HKY85, GTR (with
  empirical base frequencies; ti/tv and GTR exchangeabilities estimated
  by ML).  Eigendecomposition gives a fast P(t)=exp(Qt).
* ML branch-length optimisation (Brent per edge) and model-parameter
  optimisation, alternated to convergence.
* NNI hill-climbing topology search from a start tree (default NJ).

This is exact ML for small/medium problems.  It is pure Python, so it is
far slower than RAxML/IQ-TREE -- use those (via :mod:`phytreon.infer.ml`)
for large alignments; this engine is for convenience and self-containment.
"""
from __future__ import annotations

from typing import Dict, List, Optional

from ..core.tree import Tree
from .align import Alignment

_STATES = "ACGT"
_S2I = {c: i for i, c in enumerate(_STATES)}
_S2I["U"] = 3


# --------------------------------------------------------------------------
# data: compress alignment to unique site patterns
# --------------------------------------------------------------------------
def _encode(aln: Alignment):
    import numpy as np
    names = list(aln.names)
    ncol = aln.ncol
    # column tuple -> weight; and per-tip state index per unique column
    patterns: Dict[tuple, int] = {}
    order: List[tuple] = []
    for j in range(ncol):
        col = tuple(s[j].upper() for s in aln.seqs)
        if col not in patterns:
            patterns[col] = 0
            order.append(col)
        patterns[col] += 1
    npat = len(order)
    states = np.full((len(names), npat), -1, dtype=np.int8)
    for p, col in enumerate(order):
        for i, ch in enumerate(col):
            states[i, p] = _S2I.get(ch, -1)
    weights = np.array([patterns[c] for c in order], dtype=float)
    freqs = _empirical_freqs(states)
    return names, states, weights, freqs


def _empirical_freqs(states):
    import numpy as np
    counts = np.zeros(4)
    for s in range(4):
        counts[s] = (states == s).sum()
    if counts.sum() == 0:
        return np.full(4, 0.25)
    f = counts / counts.sum()
    f = np.clip(f, 1e-6, None)
    return f / f.sum()


# --------------------------------------------------------------------------
# substitution models -> normalised Q and its eigendecomposition
# --------------------------------------------------------------------------
def _exchange_matrix(model: str, params) -> "list":
    import numpy as np
    R = np.ones((4, 4))
    if model in ("JC69", "JC"):
        pass
    elif model in ("K80", "K2P", "HKY85", "HKY"):
        kappa = params[0]
        R = np.ones((4, 4))
        R[0, 2] = R[2, 0] = kappa      # A<->G transition
        R[1, 3] = R[3, 1] = kappa      # C<->T transition
    elif model == "GTR":
        ac, ag, at, cg, ct = params      # gt fixed to 1
        R = np.array([[0, ac, ag, at],
                      [ac, 0, cg, ct],
                      [ag, cg, 0, 1.0],
                      [at, ct, 1.0, 0.0]], dtype=float)
    else:
        raise ValueError(f"unknown model {model!r}")
    np.fill_diagonal(R, 0.0)
    return R


def _build_Q(model, params, pi):
    import numpy as np
    R = _exchange_matrix(model, params)
    Q = R * pi[None, :]
    for i in range(4):
        Q[i, i] = -Q[i].sum()
    scale = -(pi * np.diag(Q)).sum()
    return Q / scale


class _Model:
    """Holds the model + cached eigendecomposition for fast P(t)."""

    def __init__(self, name, pi, params, ncat=1, shape=0.5):
        self.name = name
        self.pi = pi
        self.params = list(params)
        self.ncat = ncat              # discrete-gamma rate categories (1 = off)
        self.shape = shape            # gamma shape alpha
        self._decompose()

    def rate_categories(self):
        """Return (rates, weights) of the discrete-gamma model (mean=1)."""
        import numpy as np
        if self.ncat <= 1:
            return np.array([1.0]), np.array([1.0])
        from scipy.stats import gamma
        from scipy.special import gammainc
        a, k = self.shape, self.ncat
        bounds = [0.0] + [gamma.ppf(i / k, a=a, scale=1.0 / a)
                          for i in range(1, k)] + [np.inf]
        rates = []
        for i in range(k):
            gl = gammainc(a + 1, a * bounds[i])
            gh = 1.0 if not np.isfinite(bounds[i + 1]) else gammainc(a + 1, a * bounds[i + 1])
            rates.append((gh - gl) * k)
        return np.array(rates), np.full(k, 1.0 / k)

    def set_shape(self, shape):
        self.shape = max(shape, 1e-3)

    def _decompose(self):
        import numpy as np
        Q = _build_Q(self.name, self.params, self.pi)
        vals, vecs = np.linalg.eig(Q)
        self.vals = vals.real
        self.vecs = vecs.real
        self.vinv = np.linalg.inv(self.vecs)

    def set_params(self, params):
        self.params = list(params)
        self._decompose()

    def P(self, t):
        import numpy as np
        return (self.vecs * np.exp(self.vals * t)) @ self.vinv

    @property
    def nparams(self):
        return {"JC69": 0, "K80": 1, "HKY85": 1, "GTR": 5}[_canon(self.name)]


def _canon(name):
    return {"JC": "JC69", "K2P": "K80", "HKY": "HKY85"}.get(name, name)


# --------------------------------------------------------------------------
# likelihood (Felsenstein pruning, vectorised over patterns, rescaled)
# --------------------------------------------------------------------------
def _site_logliks(tree: Tree, model: _Model, names, states, rate=1.0):
    """Per-pattern log site-likelihoods at one relative rate (Felsenstein)."""
    import numpy as np
    idx = {n: i for i, n in enumerate(names)}
    npat = states.shape[1]
    cache_L: Dict[int, "np.ndarray"] = {}
    cache_s: Dict[int, "np.ndarray"] = {}

    for node in tree.traverse("postorder"):
        if node.is_leaf:
            L = np.ones((npat, 4))
            srow = states[idx[node.name]]
            known = srow >= 0
            L[known] = 0.0
            L[known, srow[known]] = 1.0
            cache_L[id(node)] = L
            cache_s[id(node)] = np.zeros(npat)
        else:
            L = np.ones((npat, 4))
            scal = np.zeros(npat)
            for c in node.children:
                P = model.P(max(c.length or 0.0, 1e-9) * rate)
                L = L * (cache_L[id(c)] @ P.T)
                scal = scal + cache_s[id(c)]
            m = L.max(axis=1)
            m[m <= 0] = 1.0
            L = L / m[:, None]
            scal = scal + np.log(m)
            cache_L[id(node)] = L
            cache_s[id(node)] = scal

    root = tree.root
    site = (cache_L[id(root)] * model.pi[None, :]).sum(axis=1)
    return np.log(site) + cache_s[id(root)]


def _log_likelihood(tree: Tree, model: _Model, names, states, weights) -> float:
    import numpy as np
    rates, wts = model.rate_categories()
    if len(rates) == 1:
        return float((weights * _site_logliks(tree, model, names, states)).sum())
    # average site likelihood over discrete-gamma rate categories
    from scipy.special import logsumexp
    stacked = np.stack([_site_logliks(tree, model, names, states, r) + np.log(w)
                        for r, w in zip(rates, wts)])
    return float((weights * logsumexp(stacked, axis=0)).sum())


# --------------------------------------------------------------------------
# optimisation
# --------------------------------------------------------------------------
def _optimize_branches(tree, model, data, rounds=3) -> float:
    from scipy.optimize import minimize_scalar
    names, states, weights, _ = data
    edges = [n for n in tree.traverse() if not n.is_root]
    best = _log_likelihood(tree, model, names, states, weights)
    for _ in range(rounds):
        improved = False
        for node in edges:
            old = node.length or 0.1

            def neg(t, _node=node):
                _node.length = t
                return -_log_likelihood(tree, model, names, states, weights)

            res = minimize_scalar(neg, bounds=(1e-6, 5.0), method="bounded")
            if -res.fun > best + 1e-6:
                node.length = res.x
                best = -res.fun
                improved = True
            else:
                node.length = old
        if not improved:
            break
    return best


def _optimize_model(tree, model, data) -> float:
    from scipy.optimize import minimize
    names, states, weights, _ = data
    np_ = model.nparams
    use_gamma = model.ncat > 1
    if np_ == 0 and not use_gamma:
        return _log_likelihood(tree, model, names, states, weights)

    # parameter vector = model exchangeabilities (+ gamma shape if enabled)
    x0 = list(model.params) + ([model.shape] if use_gamma else [])

    def neg(x):
        if any(v <= 0 for v in x):
            return 1e18
        if np_:
            model.set_params(x[:np_])
        if use_gamma:
            model.set_shape(x[np_])
            model._decompose() if not np_ else None
        return -_log_likelihood(tree, model, names, states, weights)

    res = minimize(neg, x0, method="Nelder-Mead",
                   options={"xatol": 1e-3, "fatol": 1e-3, "maxiter": 300})
    if np_:
        model.set_params(res.x[:np_])
    if use_gamma:
        model.set_shape(res.x[np_])
    return -res.fun


# --------------------------------------------------------------------------
# NNI topology search (shared helpers live in infer/_search.py)
# --------------------------------------------------------------------------
def _nni_search(tree, model, data, max_sweeps=20) -> float:
    from ._search import internal_edges, nni_neighbors
    best = _optimize_branches(tree, model, data, rounds=2)
    for _ in range(max_sweeps):
        improved = False
        for node in internal_edges(tree):
            for swap in list(nni_neighbors(node)):
                swap()                                  # apply NNI
                ll = _optimize_branches(tree, model, data, rounds=1)
                if ll > best + 1e-4:
                    best = ll
                    improved = True
                    break                               # keep move, next edge
                swap()                                  # involution -> undo
                _optimize_branches(tree, model, data, rounds=1)
        if not improved:
            break
    return best


# --------------------------------------------------------------------------
# public entry point
# --------------------------------------------------------------------------
def _new_model(name, freqs, gamma):
    import numpy as np
    pi = freqs if name in ("HKY85", "GTR") else np.full(4, 0.25)
    init = {"JC69": [], "K80": [2.0], "HKY85": [2.0], "GTR": [1, 1, 1, 1, 1]}[name]
    return _Model(name, pi, init, ncat=max(1, gamma), shape=0.5)


def ml_tree(alignment: Alignment, model: str = "HKY85", search: bool = True,
            gamma: int = 0, start: Optional[Tree] = None, max_sweeps: int = 20) -> Tree:
    """Maximum-likelihood tree (pure Python).

    ``model`` is one of ``JC69`` / ``K80`` / ``HKY85`` / ``GTR``.  ``gamma`` is
    the number of discrete +G rate categories (0 = off; 4 is typical) -- this
    models among-site rate variation and usually fits real data much better.
    With ``search=True`` the topology is refined by NNI hill-climbing from the
    start tree (NJ by default).  Result carries ``tree.data['logL']``,
    ``['model']``, ``['model_params']``, ``['gamma_shape']``.
    """
    from .bootstrap import nj_builder
    from ..treeops import midpoint_root

    data = _encode(alignment)
    names, states, weights, freqs = data
    tree = start or midpoint_root(nj_builder(alignment))
    name = _canon(model)
    mdl = _new_model(name, freqs, gamma)

    prev = -1e18
    trajectory = []
    converged = False
    for _ in range(8):
        if search:
            ll = _nni_search(tree, mdl, data, max_sweeps=max_sweeps)
        else:
            ll = _optimize_branches(tree, mdl, data, rounds=4)
        ll = _optimize_model(tree, mdl, data)
        ll = _optimize_branches(tree, mdl, data, rounds=2)
        trajectory.append(round(ll, 4))
        if ll - prev < 1e-3:
            converged = True
            break
        prev = ll

    import math
    ntips = len(names)
    nsites = float(weights.sum())
    freq_params = 3 if name in ("HKY85", "GTR") else 0      # empirical base freqs
    k = (2 * ntips - 3) + mdl.nparams + freq_params + (1 if mdl.ncat > 1 else 0)
    tree.data["logL"] = ll
    tree.data["model"] = name + ("+G" if mdl.ncat > 1 else "")
    tree.data["model_params"] = list(mdl.params)
    tree.data["base_freqs"] = list(mdl.pi)
    tree.data["gamma_shape"] = mdl.shape if mdl.ncat > 1 else None
    tree.data["free_params"] = k
    tree.data["AIC"] = 2 * k - 2 * ll
    tree.data["BIC"] = k * math.log(nsites) - 2 * ll
    tree.data["logL_trajectory"] = trajectory       # per-iteration logL
    tree.data["converged"] = converged
    return tree


def model_finder(alignment: Alignment, models=("JC69", "K80", "HKY85", "GTR"),
                 gammas=(0, 4), start: Optional[Tree] = None):
    """Rank substitution models by AIC on a fixed (NJ) tree -- a lightweight
    ModelFinder.  Returns a list of dicts sorted by AIC (best first)."""
    from .bootstrap import nj_builder
    from ..treeops import midpoint_root
    base = start or midpoint_root(nj_builder(alignment))
    rows = []
    for m in models:
        for g in gammas:
            t = base  # branch lengths re-optimised per model inside ml_tree
            from copy import deepcopy
            t = deepcopy(base)
            t = ml_tree(alignment, model=m, gamma=g, search=False, start=t)
            rows.append({"model": t.data["model"], "logL": round(t.data["logL"], 2),
                         "k": t.data["free_params"], "AIC": round(t.data["AIC"], 2),
                         "BIC": round(t.data["BIC"], 2)})
    rows.sort(key=lambda r: r["AIC"])
    return rows


def log_likelihood(tree: Tree, alignment: Alignment, model: str = "HKY85",
                   params=None, gamma: int = 0, shape: float = 0.5) -> float:
    """Felsenstein log-likelihood of ``tree`` for ``alignment`` under a model
    (optionally with ``gamma`` discrete-rate categories)."""
    names, states, weights, freqs = _encode(alignment)
    name = _canon(model)
    mdl = _new_model(name, freqs, gamma)
    if params is not None:
        mdl.set_params(params)
    if gamma > 1:
        mdl.set_shape(shape)
        mdl._decompose()
    return _log_likelihood(tree, mdl, names, states, weights)
