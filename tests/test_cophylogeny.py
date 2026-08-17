"""PACo: cophylogeny congruence between a host and a symbiont tree.

The permutation test's calibration was checked directly on the public
`pt.paco` function before anything here was written, and the first attempt
at that check was itself wrong in an instructive way, kept in the comments
below because the mistake is easy to repeat: using each side's *default*
distance-matrix row order as a stand-in for "no particular correspondence"
between two independent trees is not neutral. `random_tree`'s row order
(whatever `patristic_distances` returns it in) turned out to carry enough
shared statistical structure across independently-seeded trees that pairing
"row i with row i" was measurably, systematically better than a random
permutation -- a 90% false-positive rate instead of the nominal 5%, for a
question that had nothing to do with `paco()` itself: it never engages
whatever a caller's default row order happens to be, since `links` is always
reindexed by explicit name. Once the test used an explicit, independently
randomised host<->symbiont name correspondence instead, calibration was
correct (4%). All tests here build their link tables by explicit name for
that reason, never by relying on either tree's incidental default order.
"""
import matplotlib
matplotlib.use("Agg")

import re

import numpy as np
import pandas as pd
import pytest

import phytreon as pt


def _relabelled(tree: pt.Tree, mapping: dict) -> pt.Tree:
    """A copy of ``tree`` with every leaf name replaced via ``mapping``.

    Done as one atomic regex substitution, not sequential ``.replace()``
    calls -- sequential replacement is unsafe here and genuinely bit an
    early draft of this file's own tests. When the new names are built by
    prefixing the old ones (as "S" + host name is, below), an early
    replacement's *output* can contain a later key's *search pattern* as a
    substring: mapping ``t03 -> St08`` first turns "t03:" into "St08:" in
    the string, and the later key ``t08 -> St11`` then matches that
    newly-written "St08:" too, cascading into "SSt11:". A single regex pass
    keyed off the ORIGINAL text has no such ordering dependency.
    """
    nw = tree.write()
    pattern = re.compile("|".join(re.escape(f"{old}:") for old in mapping))
    return pt.Tree.from_newick(
        pattern.sub(lambda m: mapping[m.group(0)[:-1]] + ":", nw))


def _cospeciating_pair(n, seed):
    """A host tree and a symbiont tree of IDENTICAL topology (true
    cospeciation), with the correct one-to-one link table between them."""
    host = pt.datasets.random_tree(n, seed=seed)
    hnames = host.leaf_names()
    snames = [f"S{h}" for h in hnames]
    sym = _relabelled(host, dict(zip(hnames, snames)))
    links = pd.DataFrame(np.eye(n), index=hnames, columns=snames)
    return host, sym, links


def _independent_pair_with_random_links(n, seed):
    """Two topologically unrelated trees, linked by an explicitly randomised
    (not default-row-order) one-to-one correspondence."""
    host = pt.datasets.random_tree(n, seed=seed)
    sym_raw = pt.datasets.random_tree(n, seed=seed + 50000)
    hnames = host.leaf_names()
    snames = [f"S{i}" for i in range(n)]
    sym = _relabelled(sym_raw, dict(zip(sym_raw.leaf_names(), snames)))
    rng = np.random.default_rng(seed + 90000)
    shuffled = rng.permutation(snames)
    links = pd.DataFrame(np.eye(n), index=hnames, columns=shuffled)
    return host, sym, links


# --------------------------------------------------------------------------
# Calibration: the property that matters most, checked on the real function
# --------------------------------------------------------------------------
def test_type_i_error_is_near_nominal_for_independent_trees():
    n_reps, hits = 50, 0
    for rep in range(n_reps):
        host, sym, links = _independent_pair_with_random_links(12, seed=1000 + rep)
        res = pt.paco(host, sym, links, n_perm=149, seed=rep)
        hits += res["p"] < 0.05
    rate = hits / n_reps
    # +-2 SE at n=50, p=0.05 is about +-0.06; a generous margin around that
    assert rate < 0.16, f"false-positive rate {rate:.3f} over {n_reps} reps"


def test_power_is_high_for_true_cospeciation():
    n_reps, hits = 25, 0
    for rep in range(n_reps):
        host, sym, links = _cospeciating_pair(12, seed=2000 + rep)
        res = pt.paco(host, sym, links, n_perm=149, seed=rep)
        hits += res["p"] < 0.05
    assert hits / n_reps > 0.9


def test_m2_is_lower_for_congruent_than_incongruent_pairings():
    # directional correctness, not just "is the p-value calibrated": true
    # cospeciation must fit distinctly better than an unrelated pairing
    host, sym, links = _cospeciating_pair(12, seed=3)
    congruent = pt.paco(host, sym, links, n_perm=9, seed=0)

    unrelated_sym = pt.datasets.random_tree(12, seed=8000)
    snames = [f"S{h}" for h in host.leaf_names()]
    unrelated_sym = _relabelled(unrelated_sym, dict(zip(unrelated_sym.leaf_names(), snames)))
    incongruent = pt.paco(host, unrelated_sym, links, n_perm=9, seed=0)

    assert congruent["m2"] < incongruent["m2"]


# --------------------------------------------------------------------------
# Mechanics: link residuals, inputs, errors
# --------------------------------------------------------------------------
def test_link_residuals_sum_exactly_to_m2():
    host, sym, links = _cospeciating_pair(10, seed=4)
    res = pt.paco(host, sym, links, n_perm=19, seed=0)
    assert res["link_residuals"]["squared_residual"].sum() == pytest.approx(res["m2"])
    assert len(res["link_residuals"]) == res["n_links"] == 10


def test_link_residuals_flag_a_deliberately_swapped_pair():
    # true cospeciation, EXCEPT hosts h0 and h1 are linked to each other's
    # symbionts instead of their own -- those two links should sit near the
    # top of link_residuals, an untouched one near the bottom
    host, sym, links = _cospeciating_pair(12, seed=5)
    h0, h1 = host.leaf_names()[0], host.leaf_names()[1]
    s0, s1 = f"S{h0}", f"S{h1}"
    swapped = links.copy()
    swapped.loc[h0, [s0, s1]] = swapped.loc[h0, [s1, s0]].to_numpy()
    swapped.loc[h1, [s0, s1]] = swapped.loc[h1, [s1, s0]].to_numpy()

    res = pt.paco(host, sym, swapped, n_perm=19, seed=0)
    top_pairs = set(zip(res["link_residuals"]["host"].head(2),
                        res["link_residuals"]["symbiont"].head(2)))
    assert (h0, s1) in top_pairs or (h1, s0) in top_pairs


def test_accepts_precomputed_distance_tuples_as_well_as_trees():
    host, sym, links = _cospeciating_pair(10, seed=6)
    host_d = pt.patristic_distances(host)
    sym_d = pt.patristic_distances(sym)
    res_tree = pt.paco(host, sym, links, n_perm=19, seed=0)
    res_tuple = pt.paco(host_d, sym_d, links, n_perm=19, seed=0)
    assert res_tree["m2"] == pytest.approx(res_tuple["m2"])


def test_rejects_a_links_argument_that_is_not_a_dataframe():
    host, sym, links = _cospeciating_pair(10, seed=7)
    with pytest.raises(TypeError, match="must be a pandas DataFrame"):
        pt.paco(host, sym, links.to_numpy())


def test_rejects_link_names_absent_from_either_tree():
    host, sym, links = _cospeciating_pair(10, seed=8)
    bad = links.copy()
    bad = bad.rename(index={links.index[0]: "NotAHost"})
    with pytest.raises(ValueError, match="hosts:"):
        pt.paco(host, sym, bad)


def test_rejects_too_few_links():
    host, sym, links = _cospeciating_pair(10, seed=9)
    sparse = links.copy()
    sparse.iloc[2:, :] = 0    # keep only 2 links
    with pytest.raises(ValueError, match="at least 3"):
        pt.paco(host, sym, sparse)


def test_pcoa_recovers_exact_euclidean_distances():
    # the embedding step, checked in isolation against an exact property:
    # PCoA of a genuinely Euclidean distance matrix must recover the
    # original pairwise distances from the embedded coordinates
    from phytreon.comparative.cophylogeny import _classical_pcoa
    rng = np.random.default_rng(0)
    pts = rng.normal(size=(15, 4))
    D = np.linalg.norm(pts[:, None] - pts[None, :], axis=-1)
    coords = _classical_pcoa(D)
    D_back = np.linalg.norm(coords[:, None] - coords[None, :], axis=-1)
    assert np.allclose(D, D_back, atol=1e-8)
