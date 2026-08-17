"""PACo: does one phylogeny's structure depend on another's?

:mod:`community` asks whether a *single* community is phylogenetically
clustered. This asks a different question about *two* trees at once -- given
a table of which host lineage associates with which symbiont lineage (or
which predator with which prey, which phage with which bacterium), is the
symbiont phylogeny's shape predictable from the host phylogeny's, more than
chance association would produce? That is the cophylogeny question, and it is
what "does taxon group A's phylogenetic structure track taxon group B's"
actually means as a testable claim, as opposed to the community-level
clustering :func:`~phytreon.comparative.community.ses_mpd` and relatives
measure.

Implements PACo (Procrustean Approach to Cophylogeny; Balbuena,
Miguez-Lozano & Blasco-Costa 2013, *PLoS ONE* 8(4):e61048): each phylogeny's
patristic distances become a configuration of points via Principal
Coordinates Analysis (PCoA), the host-symbiont association table expands
both configurations to one row per observed link, and Procrustes
superimposition asks how well the symbiont configuration can be rotated,
reflected and scaled onto the host one. Significance comes from permuting
which host each link is attributed to, not from a parametric null -- the
same logic as every permutation test in this package, applied here to
"is the observed host assignment special, or would a random one fit about as
well".

Not implemented: the Cailliez/Lingoes corrections PACo's original R
implementation offers for negative PCoA eigenvalues (this drops
non-positive-eigenvalue axes instead, which is that implementation's own
default, `correction="none"`); and the full jackknife-with-confidence-interval
treatment of individual links (the original paper's "pseudovalues" -- here
each link's contribution is its own squared residual from the fitted
superimposition, simpler and exact, but without the jackknife's bias
correction or interval).
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple, Union

from ..core.tree import Tree


def _classical_pcoa(D: "np.ndarray") -> "np.ndarray":  # noqa: F821
    """Principal coordinates of a distance matrix: double-centre ``-0.5 D^2``
    and eigendecompose, keeping only positive-eigenvalue axes (``correction=
    "none"``, matching PACo's own reference implementation's default -- see
    the module docstring for what is not implemented here instead).

    Exact for a genuinely Euclidean input (verified: recovers a random point
    cloud's own pairwise distances from its PCoA embedding to floating-point
    precision), and the standard embedding used for non-Euclidean distances
    such as a tree's patristic distances otherwise.
    """
    import numpy as np
    n = D.shape[0]
    D2 = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    keep = vals > 1e-8
    return vecs[:, keep] * np.sqrt(vals[keep])


def _as_distances(x, taxa_kind: str) -> Tuple["list", "np.ndarray"]:  # noqa: F821
    if isinstance(x, Tree):
        from .community import patristic_distances
        return patristic_distances(x)
    if isinstance(x, tuple) and len(x) == 2:
        names, D = x
        return list(names), D
    raise TypeError(
        f"paco: {taxa_kind} must be a Tree or an (names, distance_matrix) "
        f"tuple (e.g. from patristic_distances), not {type(x).__name__}"
    )


def _procrustes_m2(Xh: "np.ndarray", Xs: "np.ndarray",  # noqa: F821
                   host_rows: "np.ndarray", sym_rows: "np.ndarray") -> float:  # noqa: F821
    import numpy as np
    from scipy.spatial import procrustes
    X, Y = Xh[host_rows], Xs[sym_rows]
    dh, ds = X.shape[1], Y.shape[1]
    if dh < ds:
        X = np.column_stack([X, np.zeros((len(host_rows), ds - dh))])
    elif ds < dh:
        Y = np.column_stack([Y, np.zeros((len(host_rows), dh - ds))])
    return float(procrustes(X, Y)[2])


def paco(host: Union[Tree, Tuple], symbiont: Union[Tree, Tuple], links,
        n_perm: int = 999, seed: Optional[int] = None) -> Dict[str, object]:
    """PACo cophylogeny test: is ``symbiont``'s phylogenetic structure
    predictable from ``host``'s, given which lineages are observed together?

    ``host`` and ``symbiont`` are each a :class:`~phytreon.Tree` or an
    ``(names, distance_matrix)`` pair (e.g. from
    :func:`~phytreon.patristic_distances`, worth passing directly if calling
    this repeatedly on the same tree). ``links`` is a host x symbiont table
    (a :class:`pandas.DataFrame` indexed by host name with symbiont names as
    columns; any nonzero entry counts as a link) -- one row of it per
    observed association, however that association was decided (a literal
    host-parasite record, or lineages thresholded from a co-occurrence
    matrix across samples).

    Returns:

    ``m2``
        the observed sum of squared Procrustes residuals -- **smaller means
        more congruent**, the opposite direction from a correlation
        coefficient, so do not read it like one.
    ``p``
        permutation p-value: the fraction of random reassignments of host
        identity to the observed links whose ``m2`` is at least as small as
        the real one. Small ``p`` says the observed host-symbiont pairing
        fits distinctly better than a random one would -- i.e. genuine
        cophylogenetic structure, not merely that the trees each have
        structure on their own.
    ``link_residuals``
        a DataFrame, one row per link, with its own squared residual from the
        fitted superimposition -- sums to ``m2`` and ranks which specific
        associations drive (or resist) the overall congruence.
    ``n_links``, ``host_axes``, ``symbiont_axes``
        how many associations were used and how many positive-eigenvalue PCoA
        axes each side kept, for judging how much the fit rests on.
    """
    import numpy as np
    import pandas as pd

    host_names, Dh = _as_distances(host, "host")
    sym_names, Ds = _as_distances(symbiont, "symbiont")

    if not hasattr(links, "index") or not hasattr(links, "columns"):
        raise TypeError("paco: links must be a pandas DataFrame (hosts x symbionts)")
    unknown_hosts = sorted(set(links.index) - set(host_names))
    unknown_syms = sorted(set(links.columns) - set(sym_names))
    if unknown_hosts or unknown_syms:
        raise ValueError(
            f"paco: links has names not present in the corresponding tree/distances -- "
            f"hosts: {unknown_hosts[:10]}, symbionts: {unknown_syms[:10]}"
        )

    host_idx = {n: i for i, n in enumerate(host_names)}
    sym_idx = {n: i for i, n in enumerate(sym_names)}
    arr = links.to_numpy()
    hi_link, si_link = np.nonzero(arr)
    if len(hi_link) < 3:
        raise ValueError(
            f"paco: only {len(hi_link)} nonzero links in `links`; need at least 3 "
            f"for a meaningful superimposition"
        )
    host_rows = np.array([host_idx[links.index[i]] for i in hi_link])
    sym_rows = np.array([sym_idx[links.columns[j]] for j in si_link])

    Xh, Xs = _classical_pcoa(Dh), _classical_pcoa(Ds)
    observed = _procrustes_m2(Xh, Xs, host_rows, sym_rows)

    rng = np.random.default_rng(seed)
    n_hosts = len(host_names)
    ge = 0
    for _ in range(n_perm):
        perm = rng.permutation(n_hosts)
        if _procrustes_m2(Xh, Xs, perm[host_rows], sym_rows) <= observed:
            ge += 1
    p = (ge + 1) / (n_perm + 1)

    # per-link squared residuals from the FITTED (observed) superimposition,
    # which sum exactly to m2 -- ranks links by contribution without the
    # jackknife machinery (see module docstring)
    from scipy.spatial import procrustes as _procrustes_fn
    X, Y = Xh[host_rows], Xs[sym_rows]
    dh, ds = X.shape[1], Y.shape[1]
    if dh < ds:
        X = np.column_stack([X, np.zeros((len(host_rows), ds - dh))])
    elif ds < dh:
        Y = np.column_stack([Y, np.zeros((len(host_rows), dh - ds))])
    mtx1, mtx2, _ = _procrustes_fn(X, Y)
    per_link = ((mtx1 - mtx2) ** 2).sum(axis=1)

    link_table = pd.DataFrame({
        "host": [links.index[i] for i in hi_link],
        "symbiont": [links.columns[j] for j in si_link],
        "squared_residual": per_link,
    }).sort_values("squared_residual", ascending=False).reset_index(drop=True)

    return {
        "m2": observed,
        "p": float(p),
        "n_perm": n_perm,
        "n_links": len(hi_link),
        "host_axes": int(Xh.shape[1]),
        "symbiont_axes": int(Xs.shape[1]),
        "link_residuals": link_table,
    }
