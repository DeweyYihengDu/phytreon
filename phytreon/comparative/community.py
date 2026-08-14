"""Community phylogenetics: is a community's membership phylogenetically
structured, and does that structure imply selection or chance?

:mod:`~phytreon.comparative.diversity` answers "how much evolutionary history is
in this sample" (Faith's PD) and "how different are two samples' branches"
(UniFrac). This module answers the question after that one: *given* that two
samples differ, are the taxa in a sample more closely related to each other than
a random draw from the same pool would be -- and is the turnover between samples
more phylogenetic than chance can explain?

The metrics come in pairs, a raw distance and a standardised index built on it
by comparison with a null model:

============  ==================  ==============================================
raw           standardised        reads as
============  ==================  ==============================================
:func:`mpd`   :func:`ses_mpd`     relatedness across the whole community
:func:`mntd`  :func:`ses_mntd`    relatedness among each taxon's closest kin
:func:`beta_mntd`  :func:`beta_nti`  whether turnover between samples is
                                  more phylogenetic than chance
============  ==================  ==============================================

Sign conventions differ between the two families and are a classic source of
inverted results, so they are stated explicitly on each function. NRI and NTI
carry a factor of -1 (Webb 2000; they are the negations of the standardised
effect sizes, Kembel 2009), so **positive means clustered**. betaNTI does not
(Stegen et al. 2012), so **positive means more turnover than expected**.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

from ..core.tree import Node, Tree


# --------------------------------------------------------------------------
# The primitive everything here needs: tip-to-tip distance through the tree
# --------------------------------------------------------------------------
def patristic_distances(tree: Tree, taxa: Optional[Sequence[str]] = None
                       ) -> Tuple[List[str], "np.ndarray"]:  # noqa: F821
    """Patristic (tip-to-tip) distances: the summed branch length of the path
    through the tree between every pair of tips.

    Returns ``(names, D)`` with ``names`` in the same order as ``D``'s rows and
    columns; pass ``taxa`` to fix that order or to restrict to a subset. This is
    the tree's own notion of how far apart two tips are, as opposed to
    :func:`~phytreon.infer.distance.distance_matrix`, which measures sequences
    against each other without reference to any tree.

    Independent of where the tree is rooted: the path between two tips runs
    through their MRCA either way, so moving the root moves depths but not
    differences of them.
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

    # root -> tip node lists, so a pair's shared path is their common prefix
    paths: List[List[Node]] = []
    for lf in leaves:
        path, node = [], lf
        while node is not None:
            path.append(node)
            node = node.parent
        paths.append(list(reversed(path)))
    depths = np.array([sum(nd.length or 0.0 for nd in p[1:]) for p in paths])

    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            shared = 0.0
            for a, b in zip(paths[i][1:], paths[j][1:]):
                if a is not b:
                    break
                shared += a.length or 0.0
            # depth_i + depth_j - 2 * depth_of_their_MRCA
            D[i, j] = D[j, i] = depths[i] + depths[j] - 2.0 * shared
    return names, D


def _resolve_sample(sample, names: List[str]):
    """A sample as (indices into names, relative abundances over those tips)."""
    import numpy as np
    index = {nm: i for i, nm in enumerate(names)}
    if isinstance(sample, dict):
        items = [(k, float(v)) for k, v in sample.items() if float(v) > 0.0]
    else:
        items = [(k, 1.0) for k in sample]
    missing = [k for k, _ in items if k not in index]
    if missing:
        raise ValueError(f"taxa not found in tree: {sorted(set(missing))}")
    if not items:
        return np.array([], dtype=int), np.array([])
    idx = np.array([index[k] for k, _ in items], dtype=int)
    w = np.array([v for _, v in items], dtype=float)
    return idx, w / w.sum()


# --------------------------------------------------------------------------
# Raw metrics
# --------------------------------------------------------------------------
def _mpd_at(D, idx, w, abundance_weighted: bool) -> float:
    """MPD for a sample already resolved to positions in ``D`` and weights.

    Split out because the null model works by permuting *positions* -- a tip
    relabelling is exactly ``idx -> perm[idx]`` -- so every null draw can reuse
    this without going back through tip names.
    """
    import numpy as np
    if len(idx) < 2:
        return float("nan")
    sub = D[np.ix_(idx, idx)]
    if not abundance_weighted:
        return float(sub[np.triu_indices(len(idx), k=1)].mean())
    weight = np.outer(w, w)
    np.fill_diagonal(weight, 0.0)      # a taxon is not its own partner
    total = weight.sum()
    return float((sub * weight).sum() / total) if total > 0 else float("nan")


def _mntd_at(D, idx, w, abundance_weighted: bool) -> float:
    import numpy as np
    if len(idx) < 2:
        return float("nan")
    sub = D[np.ix_(idx, idx)].astype(float, copy=True)
    np.fill_diagonal(sub, np.inf)      # exclude each taxon as its own kin
    nearest = sub.min(axis=1)
    return float(nearest @ w) if abundance_weighted else float(nearest.mean())


def _beta_mntd_at(D, ia, wa, ib, wb, abundance_weighted: bool) -> float:
    import numpy as np
    if len(ia) == 0 or len(ib) == 0:
        return float("nan")
    cross = D[np.ix_(ia, ib)]
    if not abundance_weighted:
        wa = np.full(len(ia), 1.0 / len(ia))
        wb = np.full(len(ib), 1.0 / len(ib))
    return float(0.5 * (cross.min(axis=1) @ wa + cross.min(axis=0) @ wb))


def mpd(tree_or_D, sample, abundance_weighted: bool = False,
        names: Optional[List[str]] = None) -> float:
    """Mean pairwise distance: the average patristic distance between all pairs
    of taxa in one sample.

    Reads the community's *overall* relatedness, weighting deep splits heavily,
    which makes it sensitive to structure near the root -- use :func:`mntd` for
    the same question asked only of each taxon's closest relative.

    ``sample`` is an iterable of tip names, or a ``{tip: abundance}`` mapping.
    ``abundance_weighted=True`` weights each pair by the product of the two
    taxa's relative abundances instead of counting every pair once, so a pair of
    abundant taxa matters more than a pair of rare ones. ``tree_or_D`` may be a
    :class:`~phytreon.core.tree.Tree` or an already-computed distance matrix
    with its ``names`` (worth passing when calling this repeatedly).
    """
    names, D = _as_distances(tree_or_D, names)
    idx, w = _resolve_sample(sample, names)
    return _mpd_at(D, idx, w, abundance_weighted)


def mntd(tree_or_D, sample, abundance_weighted: bool = False,
         names: Optional[List[str]] = None) -> float:
    """Mean nearest-taxon distance: for each taxon in a sample, the distance to
    its closest relative *also in that sample*, averaged.

    The terminal counterpart to :func:`mpd` -- it only ever looks at each
    taxon's nearest neighbour, so it detects clustering at the tips of the tree
    even when the community spans the whole of it, which MPD would average away.

    ``abundance_weighted=True`` weights each taxon's nearest-neighbour distance
    by its own relative abundance.
    """
    names, D = _as_distances(tree_or_D, names)
    idx, w = _resolve_sample(sample, names)
    return _mntd_at(D, idx, w, abundance_weighted)


def beta_mntd(tree_or_D, sample_a, sample_b,
              abundance_weighted: bool = True,
              names: Optional[List[str]] = None) -> float:
    """Beta mean nearest-taxon distance: how far each taxon in one sample is
    from its closest relative *in the other* sample, averaged both ways.

    The between-sample form of :func:`mntd`, and the raw quantity
    :func:`beta_nti` standardises. 0 when both samples hold the same taxa;
    large when each sample's members have no close relatives in the other.
    Symmetric by construction -- the two directions are averaged, since "how
    far is A's nearest kin in B" and the reverse are not the same number
    (Stegen et al. 2012, *ISME J* 6:1653-1664).
    """
    names, D = _as_distances(tree_or_D, names)
    ia, wa = _resolve_sample(sample_a, names)
    ib, wb = _resolve_sample(sample_b, names)
    return _beta_mntd_at(D, ia, wa, ib, wb, abundance_weighted)


def _as_distances(tree_or_D, names):
    import numpy as np
    if isinstance(tree_or_D, Tree):
        return patristic_distances(tree_or_D)
    if names is None:
        raise ValueError(
            "pass names= alongside a precomputed distance matrix, so its rows "
            "can be matched to tip names"
        )
    return list(names), np.asarray(tree_or_D, dtype=float)


# --------------------------------------------------------------------------
# Standardised indices: the raw metric against a null model
# --------------------------------------------------------------------------
def _null_indices(n_tips: int, n_null: int, rng):
    """Tip-label shuffles of the distance matrix.

    The null asks "what would this metric look like for the same number of taxa,
    with the same abundances, drawn without regard to the phylogeny" -- so the
    tips are relabelled and the community left alone, rather than the reverse.
    This is Stegen et al.'s randomisation and picante's ``taxa.labels``, and it
    keeps both community richness and the tree's shape exactly as observed,
    which a null that resampled taxa instead would not.
    """
    for _ in range(n_null):
        yield rng.permutation(n_tips)


def _standardise(observed, nulls, negate: bool):
    import numpy as np
    nulls = np.asarray([v for v in nulls if np.isfinite(v)])
    if len(nulls) < 2 or not np.isfinite(observed):
        return {"null_mean": float("nan"), "null_sd": float("nan"),
                "ses": float("nan"), "p": float("nan")}
    mean, sd = float(nulls.mean()), float(nulls.std(ddof=1))
    ses = (observed - mean) / sd if sd > 0 else float("nan")
    # two-sided quantile p, counting the observation among the draws so it can
    # never come out as exactly zero
    more_extreme = int(np.sum(np.abs(nulls - mean) >= abs(observed - mean)))
    p = (more_extreme + 1) / (len(nulls) + 1)
    return {"null_mean": mean, "null_sd": sd,
            "ses": -ses if negate else ses, "p": float(p)}


def _check_columns(table, names, caller: str) -> None:
    unknown = [str(c) for c in table.columns if c not in set(names)]
    if unknown:
        raise ValueError(
            f"{caller}: {len(unknown)} of {len(table.columns)} table columns "
            f"are not tips of the tree: {sorted(unknown)[:10]}"
            f"{' ...' if len(unknown) > 10 else ''}"
        )


def _ses_table(tree, table, metric_at, negate, n_null, seed,
               abundance_weighted, obs_name, index_name):
    import numpy as np
    import pandas as pd
    names, D = patristic_distances(tree)
    _check_columns(table, names, obs_name)
    rng = np.random.default_rng(seed)
    perms = list(_null_indices(len(names), n_null, rng))

    rows = {}
    for s in table.index:
        sample = {c: float(v) for c, v in table.loc[s].items() if float(v) > 0}
        idx, w = _resolve_sample(sample, names)
        observed = metric_at(D, idx, w, abundance_weighted)
        # a tip relabelling is just a permutation of positions, so the sample's
        # abundances stay attached to each other and only which tips they sit on
        # changes -- no name bookkeeping, and it is the same operation for every
        # metric here
        nulls = [metric_at(D, perm[idx], w, abundance_weighted)
                 for perm in perms]
        stats = _standardise(observed, nulls, negate)
        rows[s] = {obs_name: observed, "null_mean": stats["null_mean"],
                   "null_sd": stats["null_sd"], index_name: stats["ses"],
                   "p": stats["p"], "n_taxa": len(sample)}
    return pd.DataFrame.from_dict(rows, orient="index")


def ses_mpd(tree: Tree, table, abundance_weighted: bool = False,
            n_null: int = 999, seed: Optional[int] = None) -> "pd.DataFrame":  # noqa: F821
    """:func:`mpd` for every sample, standardised against a null model, reported
    with the net relatedness index (NRI).

    ``NRI = -(MPD_obs - mean(MPD_null)) / sd(MPD_null)`` -- note the leading
    minus (Webb 2000), which makes NRI the *negation* of the standardised effect
    size (Kembel 2009) and is where implementations most often go wrong:

    * ``NRI > 0`` -- co-occurring taxa are **more** closely related than chance,
      i.e. phylogenetically clustered, the pattern habitat filtering leaves.
    * ``NRI < 0`` -- **less** closely related, i.e. overdispersed, the pattern
      competitive exclusion among close relatives leaves.

    Returns one row per sample with the observed MPD, the null mean and standard
    deviation it is being judged against, ``NRI``, a two-sided p-value, and how
    many taxa the sample had. ``p`` is a quantile p-value with the observation
    counted among the draws, so it is never exactly 0 and never below
    ``1 / (n_null + 1)``.
    """
    return _ses_table(tree, table, _mpd_at, True, n_null, seed,
                      abundance_weighted, "mpd", "NRI")


def ses_mntd(tree: Tree, table, abundance_weighted: bool = False,
             n_null: int = 999, seed: Optional[int] = None) -> "pd.DataFrame":  # noqa: F821
    """:func:`mntd` for every sample, standardised, reported with the nearest
    taxon index (NTI).

    Same convention and same leading minus as :func:`ses_mpd` -- ``NTI > 0`` is
    clustering -- but asked of each taxon's nearest relative rather than of all
    pairs, so it registers clustering at the tips that MPD averages away. A
    community drawn from a few tight clusters spread across the tree can have
    ``NTI > 0`` and ``NRI`` near zero at the same time; the pair is more
    informative than either alone.
    """
    return _ses_table(tree, table, _mntd_at, True, n_null, seed,
                      abundance_weighted, "mntd", "NTI")


def beta_nti(tree: Tree, table, abundance_weighted: bool = True,
             n_null: int = 999, seed: Optional[int] = None) -> "pd.DataFrame":  # noqa: F821
    """Pairwise beta nearest taxon index: for every pair of samples, how many
    standard deviations the observed :func:`beta_mntd` sits from its null.

    ``betaNTI = (betaMNTD_obs - mean(null)) / sd(null)``, with **no** sign flip,
    unlike :func:`ses_mpd` and :func:`ses_mntd` (Stegen et al. 2012):

    * ``betaNTI > +2`` -- more phylogenetic turnover than chance explains,
      read as variable selection: the two samples' conditions differ enough to
      favour different lineages.
    * ``betaNTI < -2`` -- less turnover than chance, read as homogeneous
      selection: the same conditions in both, favouring the same lineages.
    * ``|betaNTI| < 2`` -- not distinguishable from chance, so phylogeny alone
      does not separate the two samples; the usual follow-up is a
      taxonomic-turnover null (Raup-Crick on Bray-Curtis) to split that
      remainder into dispersal limitation and drift, which this function does
      not do and does not pretend to.

    Returned as a full square DataFrame of betaNTI values with a zero diagonal,
    so it can be indexed by sample pair directly.

    A pair of samples holding exactly the same taxa comes out as ``NaN`` rather
    than as a number, and that is the right answer rather than a gap: one
    relabelling is applied to both samples of a pair, so relabelling two
    identical samples leaves them identical, every null draw gives a betaMNTD of
    exactly 0, and there is no null variation to measure the observed 0 against.
    """
    import numpy as np
    import pandas as pd
    names, D = patristic_distances(tree)
    _check_columns(table, names, "beta_nti")
    samples = list(table.index)
    resolved = [_resolve_sample(
        {c: float(v) for c, v in table.loc[s].items() if float(v) > 0}, names)
        for s in samples]
    rng = np.random.default_rng(seed)
    perms = list(_null_indices(len(names), n_null, rng))

    mat = np.zeros((len(samples), len(samples)))
    for i in range(len(samples)):
        ia, wa = resolved[i]
        for j in range(i + 1, len(samples)):
            ib, wb = resolved[j]
            observed = _beta_mntd_at(D, ia, wa, ib, wb, abundance_weighted)
            # one permutation per draw, applied to *both* samples: a taxon shared
            # by the pair has one position, so it lands on one relabelled tip for
            # both and stays shared. Permuting the two independently would erase
            # the overlap the metric is about and inflate the null.
            nulls = [_beta_mntd_at(D, perm[ia], wa, perm[ib], wb,
                                   abundance_weighted) for perm in perms]
            stats = _standardise(observed, nulls, negate=False)
            mat[i, j] = mat[j, i] = stats["ses"]
    return pd.DataFrame(mat, index=samples, columns=samples)


# --------------------------------------------------------------------------
# Hypothesis tests on a distance matrix -- the step after unifrac_matrix
# --------------------------------------------------------------------------
def permanova(distances, groups, n_perm: int = 999,
              seed: Optional[int] = None) -> Dict[str, object]:
    """PERMANOVA: do these groups of samples differ in composition?

    A permutational analysis of variance on a distance matrix (Anderson 2001) --
    the standard test to run on a :func:`~phytreon.comparative.diversity.unifrac_matrix`
    once samples carry labels. Partitions total dissimilarity into within- and
    between-group parts, giving a pseudo-F, and gets its p-value by reshuffling
    the group labels rather than from an F distribution, since distances between
    the same samples are not independent observations and no closed-form null
    applies.

    ``distances`` is a square distance matrix (DataFrame or array), ``groups`` a
    same-ordered sequence of labels (or a Series indexed by sample name).
    Returns pseudo-F, ``R2`` (the fraction of dissimilarity between groups), the
    p-value, and the group sizes.
    """
    import numpy as np
    D, labels = _square(distances)
    g = _align_groups(groups, labels)
    uniq, counts = np.unique(g, return_counts=True)
    if len(uniq) < 2:
        raise ValueError(f"permanova needs at least 2 groups, got {list(uniq)}")
    if (counts < 2).any():
        small = [str(u) for u, c in zip(uniq, counts) if c < 2]
        raise ValueError(
            f"permanova needs at least 2 samples per group; {small} have one"
        )
    n = len(g)
    D2 = D ** 2

    def pseudo_f(assign):
        total = D2[np.triu_indices(n, 1)].sum() / n
        within = 0.0
        for u in uniq:
            m = assign == u
            k = int(m.sum())
            sub = D2[np.ix_(m, m)]
            within += sub[np.triu_indices(k, 1)].sum() / k
        between = total - within
        a, r = len(uniq), n - len(uniq)
        return (between / (a - 1)) / (within / r), between / total

    f_obs, r2 = pseudo_f(g)
    rng = np.random.default_rng(seed)
    ge = sum(pseudo_f(rng.permutation(g))[0] >= f_obs for _ in range(n_perm))
    return {"pseudo_F": float(f_obs), "R2": float(r2),
            "p": (ge + 1) / (n_perm + 1), "n_perm": n_perm,
            "groups": {str(u): int(c) for u, c in zip(uniq, counts)}}


def mantel(distances_a, distances_b, method: str = "pearson",
           n_perm: int = 999, seed: Optional[int] = None) -> Dict[str, object]:
    """Mantel test: are two distance matrices over the same samples correlated?

    The test for questions of the form "does community dissimilarity track
    environmental or geographic distance" -- correlate a
    :func:`~phytreon.comparative.diversity.unifrac_matrix` against a matrix of
    differences in temperature, depth, or kilometres apart. Significance comes
    from permuting one matrix's samples, because its entries share samples and
    so are not independent, which an ordinary correlation test would assume.

    ``method`` is ``"pearson"`` or ``"spearman"`` (rank-based, for a monotone
    but non-linear relationship). Returns the correlation, its p-value, and the
    number of samples compared.
    """
    import numpy as np
    from scipy.stats import pearsonr, spearmanr
    A, labels_a = _square(distances_a)
    B, labels_b = _square(distances_b)
    if labels_a is not None and labels_b is not None and labels_a != labels_b:
        shared = [x for x in labels_a if x in set(labels_b)]
        if len(shared) < 3:
            raise ValueError(
                f"mantel needs at least 3 samples in common, found {len(shared)}"
            )
        ia = [labels_a.index(x) for x in shared]
        ib = [labels_b.index(x) for x in shared]
        A, B = A[np.ix_(ia, ia)], B[np.ix_(ib, ib)]
    if A.shape != B.shape:
        raise ValueError(f"matrices are different sizes: {A.shape} vs {B.shape}")
    n = A.shape[0]
    if n < 3:
        raise ValueError(f"mantel needs at least 3 samples, got {n}")
    iu = np.triu_indices(n, 1)
    corr = {"pearson": pearsonr, "spearman": spearmanr}.get(method)
    if corr is None:
        raise ValueError(f"method must be 'pearson' or 'spearman', not {method!r}")
    r_obs = float(corr(A[iu], B[iu])[0])
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n_perm):
        p = rng.permutation(n)
        ge += abs(float(corr(A[np.ix_(p, p)][iu], B[iu])[0])) >= abs(r_obs)
    return {"r": r_obs, "p": (ge + 1) / (n_perm + 1), "n": n,
            "method": method, "n_perm": n_perm}


def _square(d):
    """A square distance matrix and its labels, from a DataFrame or an array."""
    import numpy as np
    if hasattr(d, "to_numpy") and hasattr(d, "index"):
        labels = [str(x) for x in d.index]
        arr = d.to_numpy(dtype=float)
    else:
        labels, arr = None, np.asarray(d, dtype=float)
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
        raise ValueError(f"expected a square distance matrix, got {arr.shape}")
    return arr, labels


def _align_groups(groups, labels):
    import numpy as np
    if hasattr(groups, "reindex") and labels is not None:
        return np.asarray(groups.reindex(labels).to_numpy())
    g = np.asarray(list(groups))
    if labels is not None and len(g) != len(labels):
        raise ValueError(
            f"got {len(g)} group labels for {len(labels)} samples"
        )
    return g
