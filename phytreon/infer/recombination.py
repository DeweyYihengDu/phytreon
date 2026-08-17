"""Scanning an alignment for evidence of recombination via site incompatibility.

:func:`~phytreon.gene_tree_conflict` and
:func:`~phytreon.comparative.domains.compare_domain_trees` detect exchange
between whole domains or whole genes -- one tree per unit, compared against
another. This looks *inside* a single gene's own alignment for the same kind
of signal at finer resolution: do nearby sites disagree with each other more
than chance, in a way that a single, non-recombining tree cannot produce?

**What this is not.** The literature calls this family of methods the PHI
test (Bruen, Poss & Bryant 2006, *Genetics* 172:2665-2681), built on a
statistic called "refined incompatibility" between pairs of sites -- a real
number, not just compatible-or-not, computed from the minimum number of
mutations any tree would need to explain both sites jointly. That quantity is
not implemented here: after reading the primary paper and attempting to
verify a simplifying special case independently (whether it reduces to the
plain binary compatibility test for biallelic sites), the check did not come
back clean, and shipping a mislabelled statistic under a well-known method's
name would be worse than not having it. What follows instead is the same
overall *framework* -- a windowed scan of pairwise site incompatibility, with
significance from permuting site order -- built on the classical, unambiguous
binary four-gamete compatibility test (Hudson & Kaplan 1985) rather than
Bruen et al.'s own refined statistic. Call it a four-gamete compatibility
scan, not PHI.

Two sites (each reduced to two states -- see ``biallelic_recode``) are
*compatible* if at most three of the four possible state combinations occur
across the sequences; a fourth appearing is exactly the classical signature
that no single tree, with each site mutating only once, can explain both --
recombination or recurrent mutation, and this test does not try to tell
those apart (nothing based on incompatibility alone can).
"""
from __future__ import annotations

from typing import Dict, Optional

from .matrix import Alignment


def biallelic_recode(aln: Alignment, min_count: int = 2) -> Dict[str, object]:
    """Reduce every column to at most two states for the four-gamete test.

    The classical compatibility test is defined for biallelic sites. Real
    alignments have gaps, ambiguity codes and occasionally three or more
    observed states at one column; this keeps, at each column, the two most
    common non-gap states and recodes every sequence's entry there as 0
    (matches the more common of the two), 1 (matches the other), or missing
    (anything else -- a gap, an ambiguity code, or a third allele). A column
    is **informative** only if both states reach ``min_count`` -- a singleton
    variant can never by itself create the four-gamete pattern the test looks
    for, so keeping it would only cost power for no possible gain.

    Returns ``{"states": array (nseq, n_informative), int8, -1=missing,
    "columns": the original alignment column index of each kept site}``.
    """
    import numpy as np
    from collections import Counter
    ncol = aln.ncol
    nseq = aln.nseq
    kept_cols = []
    rows = []
    for j in range(ncol):
        col = [s[j].upper() for s in aln.seqs]
        counts = Counter(c for c in col if c not in "-.NX?")
        top = counts.most_common(2)
        if len(top) < 2 or top[0][1] < min_count or top[1][1] < min_count:
            continue
        major, minor = top[0][0], top[1][0]
        row = np.full(nseq, -1, dtype=np.int8)
        for i, c in enumerate(col):
            if c == major:
                row[i] = 0
            elif c == minor:
                row[i] = 1
        kept_cols.append(j)
        rows.append(row)
    states = np.array(rows, dtype=np.int8).T if rows else np.empty((nseq, 0), dtype=np.int8)
    return {"states": states, "columns": kept_cols}


def _incompatible(a: "np.ndarray", b: "np.ndarray") -> Optional[bool]:  # noqa: F821
    """Do sites ``a`` and ``b`` (each 0/1/-1 over the same sequences) show all
    four gamete combinations? ``None`` if too little shared, non-missing data
    to tell (fewer than 2 sequences carry a definite state at both).

    The reference implementation: a direct, easy-to-trust pairwise check, used
    to verify :func:`_incompatibility_matrix`'s vectorised version rather than
    in the scan itself (which needs this for every pair of a few hundred
    sites, hundreds of permutations over -- a plain Python loop calling this
    once per pair took 20 seconds on a 270-site alignment; the matrix version
    below computes the same thing for every pair at once and cut that to
    under a hundredth of a second, verified to agree with this function
    exactly before replacing it).
    """
    both = (a >= 0) & (b >= 0)
    if both.sum() < 2:
        return None
    pairs = {(int(x), int(y)) for x, y in zip(a[both], b[both])}
    return len(pairs) == 4


def _incompatibility_matrix(states: "np.ndarray"):  # noqa: F821
    """``(incompatible, scored)``, both ``(n_sites, n_sites)`` boolean: for
    every pair of sites at once, whether all four gamete combinations occur
    (``incompatible``) and whether at least two sequences had non-missing
    data at both sites to tell (``scored``).

    Each of the four gamete types' presence, for *every* pair of sites
    simultaneously, is a matrix product: with ``M0``/``M1`` marking which
    sequences carry state 0 / state 1 at each site, ``M0.T @ M1 > 0`` at
    entry ``(i, j)`` is true exactly when some sequence has state 0 at site
    i and state 1 at site j -- gamete (0, 1) for that pair -- for every pair
    at once, rather than a Python-level loop recomputing the same thing site
    pair by site pair.
    """
    import numpy as np
    miss = states < 0
    m0 = (states == 0).astype(np.int32)
    m1 = (states == 1).astype(np.int32)
    g00 = (m0.T @ m0) > 0
    g01 = (m0.T @ m1) > 0
    g10 = (m1.T @ m0) > 0
    g11 = (m1.T @ m1) > 0
    incompatible = g00 & g01 & g10 & g11
    both_present = (~miss).T.astype(np.int32) @ (~miss).astype(np.int32)
    scored = both_present >= 2
    np.fill_diagonal(incompatible, False)
    np.fill_diagonal(scored, False)
    return incompatible, scored


def four_gamete_scan(aln: Alignment, window: int,
                     min_count: int = 2, n_perm: int = 999,
                     seed: Optional[int] = None) -> Dict[str, object]:
    """Windowed four-gamete compatibility scan for recombination.

    ``window`` restricts scored site-pairs to those within ``window``
    *informative* sites of each other, and is required rather than optional:
    without a window, every pair of sites is scored regardless of order, so
    permuting the sites changes nothing about *which* pairs get counted, and
    the "permutation test" compares the observed statistic against copies of
    itself -- always significant at exactly the same trivial value, which is
    not a test of anything (checked directly: literal, exact equality across
    every permutation, not merely "usually similar"). A window is what makes
    the statistic depend on site order at all, which is what a permutation
    test needs to have something to say.

    Significance comes from permuting the **order of informative sites**
    ``n_perm`` times and rescoring: if there is a real recombination
    breakpoint, physically nearby sites carry more shared history than
    distant ones, and this signal is destroyed by putting sites in a random
    order (Bruen et al.'s stated rationale for exactly this permutation,
    which this scan reuses even though its own statistic is different from
    theirs) -- with no true structure, the windowed statistic should not care
    what order the sites are in.

    **Power is narrow and genuinely sensitive to `window`, checked directly
    rather than assumed.** A 30-base foreign tract inserted into an otherwise
    clonal 300-base sequence was detected in 60% of replicates at
    `window=20` with low divergence, and in under 10% at `window=5`,
    `window=10`, or higher divergence -- the same simulated event, only the
    window or the mutation rate changed. A window much larger than the
    true recombinant tract dilutes the local signal into the (often already
    substantial, from ordinary homoplasy) background rate; too small a
    window may not gather enough pairs to detect it at all. There is no
    default that works well across cases: pick `window` relative to how long
    a recombination tract would plausibly be for the data at hand, and do
    not treat a non-significant result as evidence recombination did not
    happen with a differently-scaled event.

    Returns ``mean_incompatibility`` (fraction of scored pairs showing all
    four gametes -- the statistic itself), ``p`` (permutation p-value, one-
    sided: higher incompatibility than expected), ``n_informative_sites``,
    ``n_pairs_scored``, and ``n_perm``.
    """
    import numpy as np
    if window is None:
        raise ValueError(
            "four_gamete_scan: window=None scores every pair regardless of "
            "site order, which makes the permutation test compare the "
            "observed statistic against exact copies of itself -- pick a "
            "window (see the window= docstring for what that changes)"
        )
    if window < 1:
        raise ValueError(f"four_gamete_scan: window must be >= 1, got {window}")

    recoded = biallelic_recode(aln, min_count=min_count)
    states = recoded["states"]
    n_sites = states.shape[1]
    if n_sites < 2:
        raise ValueError(
            f"four_gamete_scan: only {n_sites} informative (biallelic, "
            f"min_count={min_count}) site(s) -- need at least 2 to test any pair"
        )

    # computed once for every pair of sites; a permutation only changes which
    # PAIRS fall inside the window, not whether any given pair is itself
    # incompatible, so nothing here needs recomputing per permutation
    incompatible, scored = _incompatibility_matrix(states)
    # upper triangle only: incompatible/scored/window_mask are all symmetric,
    # so summing the full matrix would count every pair twice. The mean
    # (numerator and denominator both double) is unaffected either way, but
    # n_pairs_scored should mean what it says.
    upper = np.triu(np.ones((n_sites, n_sites), dtype=bool), k=1)
    offset = np.arange(n_sites)[None, :] - np.arange(n_sites)[:, None]
    window_mask = upper & (offset <= window)

    def score(order: "np.ndarray") -> float:  # noqa: F821
        # order[k] = which original site sits at position k; reindexing the
        # precomputed matrices by `order` asks "what would be windowed
        # together if the sites came in this order" without recomputing any
        # pairwise incompatibility
        sc = scored[np.ix_(order, order)] & window_mask
        n_scored = int(sc.sum())
        if n_scored == 0:
            return float("nan"), 0
        n_hits = int((incompatible[np.ix_(order, order)] & sc).sum())
        return n_hits / n_scored, n_scored

    identity = np.arange(n_sites)
    observed, n_scored = score(identity)
    if n_scored == 0:
        raise ValueError(
            "four_gamete_scan: no site pair had enough shared non-missing "
            "data to score -- data too sparse, or window too small"
        )
    if not np.isfinite(observed):
        raise ValueError("four_gamete_scan: could not compute an observed statistic")

    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n_perm):
        val, _ = score(rng.permutation(n_sites))
        if np.isfinite(val) and val >= observed:
            ge += 1
    p = (ge + 1) / (n_perm + 1)

    return {
        "mean_incompatibility": float(observed),
        "p": float(p),
        "n_informative_sites": n_sites,
        "n_pairs_scored": n_scored,
        "n_perm": n_perm,
        "window": window,
    }
