"""Ancestral character estimation.

Two discrete-trait reconstructions are provided:

* :func:`ace_parsimony` -- Fitch parsimony (fast, model-free).
* :func:`ace_ml` -- marginal maximum likelihood under an equal-rates
  Mk model (Felsenstein pruning for the likelihood, an up/down pass for
  the marginal posterior at every internal node, single rate estimated
  by ML).

Both write results onto ``node.data`` so the plotting layer can map them
to aesthetics (e.g. pie charts / coloured node points at ancestors).
Continuous-trait reconstruction (PIC / REML) is a documented extension
point -- see :func:`ace_continuous`.
"""
from __future__ import annotations

from typing import Dict, List, Optional

from ..core.tree import Node, Tree


# --------------------------------------------------------------------------
# Fitch parsimony
# --------------------------------------------------------------------------
def ace_parsimony(tree: Tree, trait: Dict[str, str]) -> Dict[str, object]:
    """Fitch parsimony reconstruction of a discrete trait.

    ``trait`` maps tip name -> state.  Tips missing from the dict are
    treated as fully ambiguous.  Each node gets ``data['ace_state']``;
    the function returns ``{"score": int, "states": [...]}``.
    """
    states = sorted({v for v in trait.values()})
    full = set(states)

    # down-pass: state sets + score
    score = 0
    for node in tree.traverse("postorder"):
        if node.is_leaf:
            s = trait.get(node.name)
            node.data["_pset"] = {s} if s is not None else set(full)
        else:
            sets = [c.data["_pset"] for c in node.children]
            inter = set.intersection(*sets) if sets else set(full)
            if inter:
                node.data["_pset"] = inter
            else:
                node.data["_pset"] = set.union(*sets)
                score += 1

    # up-pass: resolve, preferring the parent's chosen state
    for node in tree.traverse("preorder"):
        pset = node.data.pop("_pset")
        if node.is_root:
            chosen = sorted(pset)[0]
        else:
            parent_state = node.parent.data["ace_state"]
            chosen = parent_state if parent_state in pset else sorted(pset)[0]
        node.data["ace_state"] = chosen

    return {"score": score, "states": states}


# --------------------------------------------------------------------------
# Marginal ML under an equal-rates Mk model
# --------------------------------------------------------------------------
def mk_q_builder(k: int, model: str = "ER"):
    """Return ``(build_Q, n_params, init)`` for an Mk rate matrix.

    ``model``: ``"ER"`` (equal rates, 1 param), ``"SYM"`` (symmetric,
    k(k-1)/2 params), ``"ARD"`` (all rates different, k(k-1) params).
    """
    import numpy as np
    if model == "ER":
        def build(p):
            Q = np.full((k, k), p[0], dtype=float)
            np.fill_diagonal(Q, 0.0)
            for i in range(k):
                Q[i, i] = -Q[i].sum()
            return Q
        return build, 1, [1.0]
    if model == "SYM":
        ut = [(i, j) for i in range(k) for j in range(i + 1, k)]

        def build(p):
            Q = np.zeros((k, k))
            for v, (i, j) in zip(p, ut):
                Q[i, j] = Q[j, i] = v
            for i in range(k):
                Q[i, i] = -Q[i].sum()
            return Q
        return build, len(ut), [1.0] * len(ut)
    if model == "ARD":
        pairs = [(i, j) for i in range(k) for j in range(k) if i != j]

        def build(p):
            Q = np.zeros((k, k))
            for v, (i, j) in zip(p, pairs):
                Q[i, j] = v
            for i in range(k):
                Q[i, i] = -Q[i].sum()
            return Q
        return build, len(pairs), [1.0] * len(pairs)
    raise ValueError(f"unknown Mk model {model!r}; use ER/SYM/ARD")


def ace_ml(tree: Tree, trait: Dict[str, str], model: str = "ER",
           default_length: float = 1.0) -> Dict[Node, Dict[str, float]]:
    """Marginal ML ancestral states under an Mk model (``ER``/``SYM``/``ARD``).

    Returns ``{node: {state: posterior_prob}}`` for internal nodes and writes
    ``data['ace_probs']`` / ``data['ace_state']`` onto every node.  Fitted
    rates are stored on ``tree.root.data['_ace_model']``.
    """
    import numpy as np
    from scipy.linalg import expm
    from scipy.optimize import minimize, minimize_scalar

    states = sorted({v for v in trait.values()})
    k = len(states)
    idx = {s: i for i, s in enumerate(states)}
    pi = np.full(k, 1.0 / k)
    nodes = tree.nodes("postorder")
    build_Q, nparams, init = mk_q_builder(k, model)

    def branch_p(params, t):
        return expm(build_Q(params) * t)

    def down_partials(params) -> Dict[Node, "np.ndarray"]:
        D: Dict[Node, np.ndarray] = {}
        for node in nodes:
            if node.is_leaf:
                vec = np.zeros(k)
                s = trait.get(node.name)
                if s is None:
                    vec[:] = 1.0
                else:
                    vec[idx[s]] = 1.0
                D[node] = vec
            else:
                vec = np.ones(k)
                for c in node.children:
                    vec = vec * (branch_p(params, c.length or default_length) @ D[c])
                D[node] = vec
        return D

    def neg_log_lik(params) -> float:
        if np.any(np.asarray(params) <= 0):
            return 1e18
        D = down_partials(params)
        return -np.log(float(pi @ D[tree.root]) + 1e-300)

    if nparams == 1:
        r = float(minimize_scalar(lambda x: neg_log_lik([x]),
                                  bounds=(1e-4, 100.0), method="bounded").x)
        params = [r]
    else:
        params = list(minimize(neg_log_lik, init, method="Nelder-Mead",
                               options={"xatol": 1e-3, "fatol": 1e-3,
                                        "maxiter": 400}).x)

    D = down_partials(params)
    U: Dict[Node, np.ndarray] = {tree.root: pi.copy()}
    for node in tree.nodes("preorder"):
        if node.is_root:
            continue
        p = node.parent
        msg = U[p].copy()
        for sib in p.children:
            if sib is node:
                continue
            msg = msg * (branch_p(params, sib.length or default_length) @ D[sib])
        U[node] = branch_p(params, node.length or default_length).T @ msg

    result: Dict[Node, Dict[str, float]] = {}
    for node in tree.traverse():
        post = D[node] * U[node]
        total = post.sum()
        post = post / total if total > 0 else np.full(k, 1.0 / k)
        probs = {states[i]: float(post[i]) for i in range(k)}
        node.data["ace_probs"] = probs
        node.data["ace_state"] = max(probs, key=probs.get)
        if not node.is_leaf:
            result[node] = probs

    tree.root.data["_ace_model"] = {
        "model": model, "rates": [float(p) for p in params], "states": states,
        "rate": float(params[0]) if nparams == 1 else None,
    }
    return result


# --------------------------------------------------------------------------
# Continuous traits -- extension point
# --------------------------------------------------------------------------
def ace_continuous(tree: Tree, trait: Dict[str, float]) -> Dict[Node, float]:
    """Ancestral states for a continuous trait via Felsenstein's
    independent-contrasts weighted averaging (Brownian motion).

    Returns ``{node: estimated_value}`` and writes ``data['ace_value']``.
    """
    # postorder: weighted average of children, weights = 1 / (branch len + accumulated var)
    var: Dict[Node, float] = {}
    val: Dict[Node, float] = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf:
            if node.name not in trait:
                raise KeyError(f"missing trait value for tip {node.name!r}")
            val[node] = float(trait[node.name])
            var[node] = 0.0
        else:
            ws, xs, vs = [], [], []
            for c in node.children:
                v = (c.length or 1.0) + var[c]
                ws.append(1.0 / v)
                xs.append(val[c])
                vs.append(v)
            wsum = sum(ws)
            val[node] = sum(w * x for w, x in zip(ws, xs)) / wsum
            var[node] = 1.0 / wsum
    for node in tree.traverse():
        node.data["ace_value"] = val[node]
    return val
