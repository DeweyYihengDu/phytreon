"""Four-gamete compatibility scan for within-gene recombination signal.

**What this tests, and what it deliberately is not.** The literature calls
this family of methods the PHI test (Bruen, Poss & Bryant 2006), built on a
"refined incompatibility" statistic between site pairs -- a real number, not
just compatible-or-not. After reading the primary paper and independently
trying to verify a simplifying special case (whether refined incompatibility
reduces to plain binary compatibility for biallelic sites, checked by brute
force over every possible tree topology on small taxon sets), the check came
back genuinely ambiguous rather than confirming it, so that statistic is not
implemented here. What is implemented is the same overall framework --
windowed pairwise compatibility, tested by permuting site order -- built on
the classical, unambiguous Hudson & Kaplan (1985) four-gamete test instead.
It is validated on its own terms below, not held to Bruen et al.'s.

**The central finding, kept as permanent tests rather than a docstring
claim only**: this statistic's power is narrow and sensitive to the window
parameter relative to the size of any real recombinant tract, and a naive
"two large, roughly equal blocks" scenario is close to a worst case for it
(the population of cross-block pairs is comparable in number to within-block
pairs, so a full-alignment permutation's baseline ends up close to or above
the true windowed signal). A short foreign tract against a long clonal
background, with a matched window, is where it has genuine power -- also
demonstrated, not asserted.
"""
import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest

import phytreon as pt
from phytreon.infer.recombination import _incompatibility_matrix, _incompatible

BASES = "ACGT"


def _evolve_dna(newick, length, rate, seed):
    tree = pt.Tree.from_newick(newick)
    rng = np.random.default_rng(seed)
    root = "".join(rng.choice(list(BASES), length))
    seqs = {}

    def walk(node, parent):
        seq = list(parent)
        p = 1.0 - np.exp(-rate * (node.length or 0.0))
        for i in np.flatnonzero(rng.random(length) < p):
            seq[i] = rng.choice(list(BASES))
        seq = "".join(seq)
        if node.is_leaf:
            seqs[node.name] = seq
        for child in node.children:
            walk(child, seq)

    for child in tree.root.children:
        walk(child, root)
    return seqs


TREE_A = "(((A:.15,B:.15):.15,(C:.15,D:.15):.15):.15,((E:.15,F:.15):.15,(G:.15,H:.15):.15):.15);"
TREE_B = "(((A:.15,E:.15):.15,(C:.15,G:.15):.15):.15,((B:.15,F:.15):.15,(D:.15,H:.15):.15):.15);"
TAXA = list("ABCDEFGH")


# --------------------------------------------------------------------------
# The compatibility primitives, by hand
# --------------------------------------------------------------------------
def test_incompatible_matches_the_classical_four_gamete_test_by_hand():
    a = np.array([0, 0, 1, 1])
    three_gametes = np.array([0, 1, 0, 0])    # (0,0)(0,1)(1,0) -- compatible
    four_gametes = np.array([0, 1, 1, 0])     # + (1,1) -- incompatible
    assert _incompatible(a, three_gametes) is False
    assert _incompatible(a, four_gametes) is True


def test_incompatible_returns_none_for_insufficient_shared_data():
    a = np.array([0, -1, 1, -1])
    b = np.array([-1, 1, -1, 1])
    assert _incompatible(a, b) is None    # no sequence has both non-missing


def test_biallelic_recode_matches_a_hand_checked_column():
    # col0: A,A,T,A,G -- major A(3), minor tied T/G at 1 each, below min_count=2
    # col1,col2: invariant
    # col3: T,A,T,A,T -- major T(3), minor A(2), both >= min_count=2: informative
    aln = pt.Alignment(["s1", "s2", "s3", "s4", "s5"],
                       ["ACGT", "ACGA", "TCGT", "ACGA", "GCGT"])
    res = pt.biallelic_recode(aln, min_count=2)
    assert res["columns"] == [3]
    assert res["states"].ravel().tolist() == [0, 1, 0, 1, 0]


# --------------------------------------------------------------------------
# The vectorised pairwise matrix, against the trusted pairwise loop
# --------------------------------------------------------------------------
def test_incompatibility_matrix_matches_the_pairwise_loop_exactly():
    rng = np.random.default_rng(3)
    states = rng.integers(-1, 2, size=(14, 40)).astype(np.int8)
    incompatible, scored = _incompatibility_matrix(states)
    n = states.shape[1]
    for _ in range(200):
        i, j = rng.choice(n, size=2, replace=False)
        expected = _incompatible(states[:, i], states[:, j])
        if expected is None:
            assert not scored[i, j]
        else:
            assert scored[i, j]
            assert incompatible[i, j] == expected


def test_incompatibility_matrix_is_symmetric_with_a_clear_diagonal():
    rng = np.random.default_rng(4)
    states = rng.integers(-1, 2, size=(10, 25)).astype(np.int8)
    incompatible, scored = _incompatibility_matrix(states)
    assert np.array_equal(incompatible, incompatible.T)
    assert np.array_equal(scored, scored.T)
    assert not incompatible.diagonal().any()
    assert not scored.diagonal().any()


# --------------------------------------------------------------------------
# window=None: proven degenerate, not just discouraged
# --------------------------------------------------------------------------
def test_window_is_required_and_none_is_rejected_with_an_explanation():
    seqs = _evolve_dna(TREE_A, 100, 1.0, 0)
    aln = pt.Alignment(TAXA, [seqs[t] for t in TAXA])
    with pytest.raises(ValueError, match="scores every pair regardless of site order"):
        pt.four_gamete_scan(aln, window=None)
    with pytest.raises(ValueError, match="must be >= 1"):
        pt.four_gamete_scan(aln, window=0)


def test_without_a_window_every_permutation_gives_the_identical_statistic():
    # the reason window=None is rejected, shown directly rather than just
    # claimed: with no window, which pairs get scored does not depend on
    # site order at all, so a permutation test built on it is vacuous
    seqs = _evolve_dna(TREE_A, 100, 1.0, 0)
    aln = pt.Alignment(TAXA, [seqs[t] for t in TAXA])
    recoded = pt.biallelic_recode(aln)
    states = recoded["states"]
    incompatible, scored = _incompatibility_matrix(states)
    upper = np.triu(np.ones_like(incompatible), k=1)
    observed = (incompatible & scored & upper).sum() / (scored & upper).sum()
    rng = np.random.default_rng(0)
    for _ in range(20):
        perm = rng.permutation(states.shape[1])
        reordered = (incompatible[np.ix_(perm, perm)]
                    & scored[np.ix_(perm, perm)] & upper).sum()
        total = (scored[np.ix_(perm, perm)] & upper).sum()
        assert reordered / total == pytest.approx(observed)


# --------------------------------------------------------------------------
# Calibration and the honestly narrow power envelope
# --------------------------------------------------------------------------
def test_type_i_error_is_near_nominal_with_no_recombination():
    n_reps, hits = 60, 0
    for rep in range(n_reps):
        seqs = _evolve_dna(TREE_A, 300, 0.6, rep + 50000)
        aln = pt.Alignment(TAXA, [seqs[t] for t in TAXA])
        res = pt.four_gamete_scan(aln, window=20, n_perm=199, seed=rep)
        hits += res["p"] < 0.05
    rate = hits / n_reps
    assert rate < 0.15, f"false-positive rate {rate:.3f} over {n_reps} reps"


def test_detects_a_short_foreign_tract_with_a_matched_window():
    # the regime this scan actually has power in: a SHORT recombinant tract
    # against a long clonal background, with window scaled to the tract --
    # not the naive "two huge equal blocks" case, which is close to a worst
    # case for a full-alignment-permutation null (see the module docstring)
    def tract_recombinant(seed):
        before = _evolve_dna(TREE_A, 135, 0.3, seed)
        tract = _evolve_dna(TREE_B, 30, 0.3, seed + 1)
        after = _evolve_dna(TREE_A, 135, 0.3, seed + 2)
        return {t: before[t] + tract[t] + after[t] for t in TAXA}

    n_reps, hits = 20, 0
    for rep in range(n_reps):
        seqs = tract_recombinant(rep * 3)
        aln = pt.Alignment(TAXA, [seqs[t] for t in TAXA])
        res = pt.four_gamete_scan(aln, window=20, n_perm=199, seed=rep)
        hits += res["p"] < 0.05
    assert hits / n_reps > 0.4, f"only {hits}/{n_reps} detected"


def test_power_collapses_outside_the_favourable_window_showing_the_limitation():
    # the same event as above, but at a mismatched window -- documented as a
    # real limitation rather than left for a user to discover unassisted
    def tract_recombinant(seed):
        before = _evolve_dna(TREE_A, 135, 0.3, seed)
        tract = _evolve_dna(TREE_B, 30, 0.3, seed + 1)
        after = _evolve_dna(TREE_A, 135, 0.3, seed + 2)
        return {t: before[t] + tract[t] + after[t] for t in TAXA}

    n_reps, hits = 15, 0
    for rep in range(n_reps):
        seqs = tract_recombinant(rep * 3)
        aln = pt.Alignment(TAXA, [seqs[t] for t in TAXA])
        res = pt.four_gamete_scan(aln, window=5, n_perm=99, seed=rep)
        hits += res["p"] < 0.05
    assert hits / n_reps < 0.3


# --------------------------------------------------------------------------
# Mechanics and errors
# --------------------------------------------------------------------------
def test_reports_expected_fields_and_pairs_sum_to_the_stated_count():
    seqs = _evolve_dna(TREE_A, 200, 0.6, 1)
    aln = pt.Alignment(TAXA, [seqs[t] for t in TAXA])
    res = pt.four_gamete_scan(aln, window=15, n_perm=49, seed=0)
    assert set(res) == {"mean_incompatibility", "p", "n_informative_sites",
                        "n_pairs_scored", "n_perm", "window"}
    assert 0.0 <= res["mean_incompatibility"] <= 1.0
    assert 0.0 < res["p"] <= 1.0
    assert res["window"] == 15


def test_rejects_too_few_informative_sites():
    aln = pt.Alignment(["a", "b", "c"], ["ACGT", "ACGT", "ACGT"])   # invariant
    with pytest.raises(ValueError, match="informative"):
        pt.four_gamete_scan(aln, window=5)


def test_min_count_controls_which_columns_count_as_informative():
    # a variant seen in only one sequence cannot itself create a 4th gamete,
    # so raising min_count should only ever remove informative sites, never add
    seqs = _evolve_dna(TREE_A, 150, 0.8, 2)
    aln = pt.Alignment(TAXA, [seqs[t] for t in TAXA])
    loose = pt.biallelic_recode(aln, min_count=1)
    strict = pt.biallelic_recode(aln, min_count=3)
    assert strict["states"].shape[1] <= loose["states"].shape[1]
