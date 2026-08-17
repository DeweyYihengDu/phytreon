"""Phylogenetic factorization: which edges of a tree best explain variation
in a covariate, found automatically rather than checked one clade at a time.

:func:`~phytreon.comparative.community.ses_mpd` and relatives ask "is *this*
pre-chosen sample or clade phylogenetically clustered". Phylofactorization
(Washburne et al. 2017, *PeerJ* 5:e2969; graph-partitioning generalisation
2019, *Ecological Monographs* 89:e01353) asks a different, more automatic
question: given a samples x taxa abundance table and a covariate over
samples, **which edge of the tree** -- i.e. which split of taxa into two
clades -- best explains variation in that covariate, with no clade chosen in
advance?

The mechanism is an isometric log-ratio (ILR) balance: for a candidate edge
splitting a set of taxa into two sides, the balance is
``sqrt(n1*n2/(n1+n2)) * log(geometric_mean(side1) / geometric_mean(side2))``
per sample -- the standard compositional-data coordinate for "how much more
of side 1 is here relative to side 2, scale-invariant to how you counted".
Scoring every edge's balance against the covariate and keeping the best one
is one *factor*; :func:`phylofactor` repeats this on the resulting two-piece
partition, greedily and recursively (Washburne's "graph partitioning"
framing: cutting the tree at the winning edge produces two components, and
future edges are searched for within one of those components rather than
across the whole tree again), so each successive factor explains the
strongest remaining association after the previous ones are accounted for.

Not implemented: the original paper's stopping rule (a Kolmogorov-Smirnov
test on the sequence of explained-variance values, deciding automatically how
many factors are real); ``n_factors`` is given explicitly instead, the
simpler and more transparent choice for a first cut. Also not implemented:
Washburne's iterative regression-based multiple-covariate GLM extension --
only a single covariate (continuous, via linear regression, or categorical,
via one-way ANOVA) is scored per edge.
"""
from __future__ import annotations

from typing import Dict, Sequence, Union

from ..core.tree import Tree


def _ilr_balance(values: "np.ndarray", idx1: "np.ndarray",  # noqa: F821
                 idx2: "np.ndarray", pseudocount: float) -> "np.ndarray":  # noqa: F821
    import numpy as np
    n1, n2 = len(idx1), len(idx2)
    g1 = np.exp(np.log(values[:, idx1] + pseudocount).mean(axis=1))
    g2 = np.exp(np.log(values[:, idx2] + pseudocount).mean(axis=1))
    return float(np.sqrt(n1 * n2 / (n1 + n2))) * np.log(g1 / g2)


def _score(balance: "np.ndarray", y, categorical: bool):  # noqa: F821
    """F-statistic, p-value and (continuous case only) r-squared of
    ``balance`` against the covariate ``y``."""
    import numpy as np
    if categorical:
        from scipy import stats as _stats
        groups = [balance[np.asarray(y) == g] for g in sorted(set(y))]
        groups = [g for g in groups if len(g) > 0]
        if len(groups) < 2:
            return 0.0, 1.0, None
        f, p = _stats.f_oneway(*groups)
        return float(f), float(p), None
    y = np.asarray(y, dtype=float)
    n = len(balance)
    x = y - y.mean()
    denom = float(x @ x)
    if denom <= 0:
        return 0.0, 1.0, None
    b = float(x @ balance) / denom
    resid = balance - b * x - balance.mean()
    ss_res = float(resid @ resid)
    ss_tot = float(((balance - balance.mean()) ** 2).sum())
    if ss_tot <= 0 or ss_res <= 0 or n <= 2:
        return 0.0, 1.0, None
    r2 = 1.0 - ss_res / ss_tot
    if r2 >= 1.0:
        return float("inf"), 0.0, 1.0
    f = r2 / (1.0 - r2) * (n - 2)
    from scipy import stats as _stats
    p = float(_stats.f.sf(f, 1, n - 2))
    return float(f), p, r2


def _candidate_edges(tree: Tree, bin_taxa: set):
    """Every edge of ``tree`` whose bipartition splits ``bin_taxa`` (and only
    ``bin_taxa``) into two nonempty groups, as ``(side1_taxa, side2_taxa)``.

    Restricting candidates to edges internal to one current partition bin
    is the "graph partitioning" step: once an edge has been factored out, the
    two pieces it created are explored independently, so a later factor
    cannot re-cross a split already made.
    """
    out = []
    seen = set()
    for node in tree.traverse("postorder"):
        if node.parent is None:
            continue
        side1 = frozenset(node.leaf_names()) & bin_taxa
        if not side1 or side1 == bin_taxa:
            continue
        if side1 in seen:
            continue
        seen.add(side1)
        out.append((side1, bin_taxa - side1))
    return out


def phylofactor(tree: Tree, table, covariate: Union[Sequence, str],
                n_factors: int = 1, categorical: bool = False,
                pseudocount: float = 1e-6) -> Dict[str, object]:
    """Greedy phylogenetic factorization: the ``n_factors`` tree edges whose
    ILR balance best explains ``covariate``, found in order.

    ``table`` is a samples x taxa abundance :class:`~pandas.DataFrame` (every
    column must be a tip of ``tree`` -- see
    :func:`~phytreon.comparative.diversity.unifrac_matrix` for why this
    matters: an ASV table usually has more columns than the tree has tips,
    and the extras have to be dropped explicitly rather than silently
    counted). ``covariate`` is either a length-``nsamples`` sequence aligned
    to ``table``'s rows, or a column name already in ``table`` to pull out
    and exclude from the taxa being factored. ``categorical=True`` scores
    each edge by a one-way ANOVA F-statistic on the balance across
    covariate groups instead of linear regression -- for a treatment/control
    label rather than a continuous measurement.

    Returns a DataFrame, one row per factor in the order found, with the
    winning edge's two sides (``side1``/``side2``, as taxon lists),
    ``side1_size``, the F-statistic and p-value, and the fraction of the
    *current bin's* balance variance the edge explains (``r2`` for continuous
    covariates; omitted for categorical, where F itself is the natural scale).
    Also returns ``"balances"``, the per-sample ILR balance of each factor
    (a samples x n_factors DataFrame) -- the actual explanatory coordinate,
    ready to correlate against anything else or plot.

    **Read the p-value as "how strong", not as a calibrated significance
    test.** Each factor is the single best-scoring edge out of every
    candidate in its bin -- up to roughly ``2 x n_taxa`` of them for the
    first factor -- and reporting the winner of many tests as if it were one
    test inflates the false-positive rate sharply: measured on data with
    *no* real association at all, the top factor's p was below 0.05 in 62%
    of replicates, not 5%. Use the p-values to rank factors and to sanity
    check the overall search (near 1 says nothing in the tree explains the
    covariate at all), not as an unadjusted claim of significance for the
    winning edge; correct for the search (e.g. a label-permutation null on
    the whole procedure) before reporting one.
    """
    import numpy as np
    import pandas as pd

    if isinstance(covariate, str):
        if covariate not in table.columns:
            raise ValueError(
                f"phylofactor: covariate {covariate!r} is not a column of table"
            )
        y = table[covariate].to_numpy()
        table = table.drop(columns=[covariate])
    else:
        y = np.asarray(list(covariate))
        if len(y) != len(table):
            raise ValueError(
                f"phylofactor: covariate has {len(y)} entries for "
                f"{len(table)} samples"
            )

    tree_taxa = set(tree.leaf_names())
    unknown = [str(c) for c in table.columns if c not in tree_taxa]
    if unknown:
        raise ValueError(
            f"phylofactor: {len(unknown)} of {len(table.columns)} table "
            f"columns are not tips of the tree: {sorted(unknown)[:10]}"
            f"{' ...' if len(unknown) > 10 else ''}"
        )
    if n_factors < 1:
        raise ValueError(f"phylofactor: n_factors must be >= 1, got {n_factors}")
    if table.shape[1] < 2 ** n_factors:
        raise ValueError(
            f"phylofactor: {table.shape[1]} taxa cannot support {n_factors} "
            f"factors (each factor needs its bin to have >= 2 taxa)"
        )

    taxa = list(table.columns)
    col_idx = {t: i for i, t in enumerate(taxa)}
    values = table.to_numpy(dtype=float)
    bins = [frozenset(taxa)]

    rows = []
    balances = {}
    for k in range(n_factors):
        best = None    # (f, p, r2, side1, side2, bin_index, balance)
        for bi, bin_taxa in enumerate(bins):
            if len(bin_taxa) < 2:
                continue
            for side1, side2 in _candidate_edges(tree, set(bin_taxa)):
                i1 = np.array([col_idx[t] for t in side1])
                i2 = np.array([col_idx[t] for t in side2])
                bal = _ilr_balance(values, i1, i2, pseudocount)
                f, p, r2 = _score(bal, y, categorical)
                if best is None or f > best[0]:
                    best = (f, p, r2, side1, side2, bi, bal)
        if best is None:
            break
        f, p, r2, side1, side2, bi, bal = best
        rows.append({
            "side1": sorted(side1), "side2": sorted(side2),
            "side1_size": len(side1), "side2_size": len(side2),
            "F": f, "p": p, "r2": r2,
        })
        balances[f"factor{k + 1}"] = bal
        bins[bi:bi + 1] = [frozenset(side1), frozenset(side2)]

    factor_table = pd.DataFrame(rows)
    balance_table = pd.DataFrame(balances, index=table.index)
    return {
        "factors": factor_table,
        "balances": balance_table,
        "n_factors_found": len(rows),
        "remaining_bins": [sorted(b) for b in bins],
    }
