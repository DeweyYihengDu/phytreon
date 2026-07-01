"""Stochastic character mapping.

Simulates discrete-trait histories along the branches under an Mk model,
conditioned on the observed tip states (Nielsen 2002; Huelsenbeck 2003):

1. fit the rate by ML and compute the down-pass conditional likelihoods,
2. for each replicate, sample internal node states (preorder, conditioned)
   then simulate the CTMC path along every branch (rejection sampling),
3. summarise over replicates: per-node posterior state probabilities and
   per-branch average dwell-time fractions.

Results are written to ``node.data`` (``ace_probs``, ``ace_state``,
``paint_segments``) so :meth:`TreeFigure.painted_branches` can paint the
branches and :meth:`TreeFigure.node_pies` can chart node posteriors.
"""
from __future__ import annotations

from typing import Dict, Optional

from ..core.tree import Tree


def stochastic_map(tree: Tree, trait: Dict[str, str], n: int = 200,
                   model: str = "ER", rate: Optional[float] = None,
                   seed: Optional[int] = None) -> Tree:
    """Stochastic character map a discrete ``trait`` ({tip: state}) over ``tree``.

    ``model`` is the Mk model (``ER``/``SYM``/``ARD``); ``n`` is the number of
    simulated histories (default 200 -- raise for smoother posteriors).
    """
    import numpy as np
    from scipy.linalg import expm
    from scipy.optimize import minimize, minimize_scalar
    from .ace import mk_q_builder

    states = sorted({v for v in trait.values()})
    k = len(states)
    idx = {s: i for i, s in enumerate(states)}
    pi = np.full(k, 1.0 / k)
    post = tree.nodes("postorder")
    build_Q, nparams, init = mk_q_builder(k, model)

    def down(params):
        Q = build_Q(params)
        D = {}
        for node in post:
            if node.is_leaf:
                v = np.zeros(k)
                s = trait.get(node.name)
                if s is None:
                    v[:] = 1.0
                else:
                    v[idx[s]] = 1.0
                D[id(node)] = v
            else:
                v = np.ones(k)
                for c in node.children:
                    v = v * (expm(Q * max(c.length or 0.0, 1e-6)) @ D[id(c)])
                D[id(node)] = v
        return D

    if rate is not None:
        params = [rate] if nparams == 1 else [rate] * nparams
    else:
        def negll(params):
            if np.any(np.asarray(params) <= 0):
                return 1e18
            return -np.log(pi @ down(params)[id(tree.root)] + 1e-300)
        if nparams == 1:
            params = [float(minimize_scalar(lambda r: negll([r]),
                                            bounds=(1e-4, 100.0),
                                            method="bounded").x)]
        else:
            params = list(minimize(negll, init, method="Nelder-Mead",
                                   options={"maxiter": 400}).x)

    D = down(params)
    Q = build_Q(params)
    rate = float(params[0])
    P = {id(n): expm(Q * max(n.length or 0.0, 1e-6))
         for n in tree.traverse() if not n.is_root}
    rng = np.random.default_rng(seed)

    counts = {id(n): np.zeros(k) for n in tree.traverse()}
    dwell = {id(n): np.zeros(k) for n in tree.traverse() if not n.is_root}

    def sim_branch(a, b, t):
        for _ in range(50):
            segs, cur, rem = [], a, t
            while True:
                q = -Q[cur, cur]
                dt = rng.exponential(1.0 / q) if q > 0 else rem
                if dt >= rem:
                    segs.append((cur, rem))
                    break
                segs.append((cur, dt))
                rem -= dt
                pr = Q[cur].copy()
                pr[cur] = 0.0
                cur = int(rng.choice(k, p=pr / pr.sum()))
            if cur == b:
                return segs
        return [(a, t / 2), (b, t / 2)] if a != b else [(a, t)]

    pre = tree.nodes("preorder")
    for _ in range(n):
        s = {}
        pr = pi * D[id(tree.root)]
        s[id(tree.root)] = int(rng.choice(k, p=pr / pr.sum()))
        counts[id(tree.root)][s[id(tree.root)]] += 1
        for node in pre:
            if node.is_root:
                continue
            sp = s[id(node.parent)]
            pr = P[id(node)][sp] * D[id(node)]
            si = int(rng.choice(k, p=pr / pr.sum()))
            s[id(node)] = si
            counts[id(node)][si] += 1
            for st, dt in sim_branch(sp, si, max(node.length or 0.0, 1e-6)):
                dwell[id(node)][st] += dt

    for node in tree.traverse():
        c = counts[id(node)]
        node.data["ace_probs"] = {states[i]: float(c[i] / n) for i in range(k)}
        node.data["ace_state"] = states[int(c.argmax())]
    for node in tree.traverse():
        if node.is_root:
            continue
        d = dwell[id(node)]
        tot = d.sum() or 1.0
        node.data["paint_segments"] = [(states[i], float(d[i] / tot))
                                       for i in range(k) if d[i] > 0]
    tree.data["stochastic_map"] = {"rate": rate, "states": states, "n": n}
    return tree
