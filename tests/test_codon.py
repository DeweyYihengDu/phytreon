"""Codon-based selection tests: GY94 model mechanics, and simulation
recovery for M0, free-ratio, and the corrected branch-site test.

**A real numerical bug found during validation, kept as a permanent
regression guard below.** ``_CodonModel`` originally diagonalised its
61x61 rate matrix with the general (non-symmetric) ``numpy.linalg.eig``.
At ``omega=1`` -- which every branch-site fit visits unconditionally
(site class 1 always; class 2 as well under the null it is tested
against) -- the model stops distinguishing synonymous from
nonsynonymous change, and with near-uniform codon frequencies this makes
the rate matrix's eigenvalues genuinely, not just numerically, repeated.
``eig``'s eigenvector matrix is only reliable when it can be inverted
cleanly, which repeated eigenvalues do not guarantee; in practice this
made it return transition "probabilities" as negative as -1.96 (kappa=
3.101, uniform pi), propagating to NaN log-likelihoods that silently
corrupted every optimizer step built on top of them. The fix exploits
reversibility (``pi_i * Q_ij == pi_j * Q_ji``, true here because the
GY94 rate depends only on the unordered codon pair) to symmetrise into
``S = diag(sqrt(pi)) @ Q @ diag(1/sqrt(pi))`` and diagonalise that with
``eigh`` instead, which returns an orthogonal eigenbasis regardless of
eigenvalue multiplicity.

A second, independent bug was found in :func:`~phytreon.infer.codon.
branch_site_test`'s optimiser: ``omega2`` is fit via a smooth
reparametrisation starting at ``raw=0``, and scipy's default
Nelder-Mead initial simplex perturbs a coordinate starting at exactly
0.0 by a tiny absolute step -- too tiny to move it at all against the
other four parameters' stronger initial gradients. Two unrelated
simulated datasets (including one with an enormous, unmistakable
``p2=0.5, omega2=15`` signal) both converged to ``omega2 ~= 50.5``,
the untouched starting value, regardless of what the data said. Fixed
with an explicit, deliberately generous initial simplex.

**Caveat found during validation, not a bug.** Once both fixes were in
place, the likelihood surface itself was checked directly (bypassing
the optimiser) and found to be genuinely flat for ``omega2`` beyond
roughly 30-50 on the validation tree used here: the likelihood-ratio
test correctly and strongly detects that *some* elevation above 1 is
present, but the point estimate of *how much* is only weakly identified
once the foreground branch is long enough to be near-saturated, and can
land anywhere from the true value up to the method's own upper bound
(100) with little likelihood cost either way. This matches documented
behaviour of the branch-site test elsewhere (PAML/HyPhy users routinely
see foreground omega estimates run to the software's own cap) and is
tested below as an expected property (the LRT stays informative; the
point estimate does not), not treated as something to eliminate.
"""
import numpy as np
import pytest

import phytreon as pt
from phytreon.core.tree import Tree
from phytreon.infer.align import Alignment
from phytreon.infer.codon import (
    CODON_AA, SENSE_CODONS, STOP_CODONS, _build_Q_codon, _CodonModel,
    _encode_codons, _foreground_edges,
)

# --------------------------------------------------------------------------
# A small tree with an unambiguous two-taxon foreground clade, and a
# from-scratch forward simulator built on the model's own primitives --
# the same approach :mod:`~phytreon.infer.ml_native`'s own validation uses
# for nucleotide/protein models, generalised to 61 codon states and a
# 4-class site mixture.
# --------------------------------------------------------------------------
NEWICK = "(((A:0.1,B:0.1):1.0,C:0.3):0.2,(D:0.25,(E:0.2,F:0.2):0.15):0.2);"
FOREGROUND = ["A", "B"]


def _simulate_branch_site(tree, kappa, omega0, omega2, p0, p1, n_codons, seed,
                          pi=None):
    from phytreon.infer.codon import _class_weights
    rng = np.random.default_rng(seed)
    if pi is None:
        pi = np.full(61, 1.0 / 61)
    fg_stem = _foreground_edges(tree, FOREGROUND)
    model0 = _CodonModel(kappa, omega0, pi)
    model1 = _CodonModel(kappa, 1.0, pi)
    model2 = _CodonModel(kappa, omega2, pi)
    class_bm = [
        {n: model0 for n in tree.traverse() if not n.is_root},
        {n: model1 for n in tree.traverse() if not n.is_root},
        {n: (model2 if n in fg_stem else model0) for n in tree.traverse() if not n.is_root},
        {n: (model2 if n in fg_stem else model1) for n in tree.traverse() if not n.is_root},
    ]
    weights = _class_weights(p0, p1)
    names = tree.leaf_names()
    seqs = {n: [] for n in names}
    for _ in range(n_codons):
        bm = class_bm[rng.choice(4, p=weights)]
        state = {tree.root: rng.choice(61, p=pi)}
        for node in tree.traverse("preorder"):
            if node.is_root:
                continue
            P = bm[node].P(node.length or 0.0)
            row = P[state[node.parent]]
            state[node] = rng.choice(61, p=row / row.sum())
            if node.is_leaf:
                seqs[node.name].append(SENSE_CODONS[state[node]])
    return Alignment(names=list(names), seqs=["".join(seqs[n]) for n in names])


def _simulate_free_ratio(tree, kappa, omega_bg, omega_fg, n_codons, seed):
    rng = np.random.default_rng(seed)
    pi = np.full(61, 1.0 / 61)
    fg_stem = _foreground_edges(tree, FOREGROUND)
    model_bg = _CodonModel(kappa, omega_bg, pi)
    model_fg = _CodonModel(kappa, omega_fg, pi)
    names = tree.leaf_names()
    seqs = {n: [] for n in names}
    for _ in range(n_codons):
        state = {tree.root: rng.choice(61, p=pi)}
        for node in tree.traverse("preorder"):
            if node.is_root:
                continue
            m = model_fg if node in fg_stem else model_bg
            row = m.P(node.length or 0.0)[state[node.parent]]
            state[node] = rng.choice(61, p=row / row.sum())
            if node.is_leaf:
                seqs[node.name].append(SENSE_CODONS[state[node]])
    return Alignment(names=list(names), seqs=["".join(seqs[n]) for n in names])


# --------------------------------------------------------------------------
# The genetic code and the GY94 rate matrix, hand- and cross-checked
# --------------------------------------------------------------------------
def test_genetic_code_matches_biopython_exactly():
    from Bio.Data import CodonTable
    bio = CodonTable.unambiguous_dna_by_id[1]
    for codon in CODON_AA:
        if codon in bio.stop_codons:
            assert CODON_AA[codon] == "*"
        else:
            assert CODON_AA[codon] == bio.forward_table[codon]
    assert len(SENSE_CODONS) == 61
    assert len(STOP_CODONS) == 3


def test_codon_frequencies_f3x4_matches_a_hand_computed_example():
    # two codons, positions independently countable by hand:
    # pos0: A,A -> A=1.0   pos1: T,C -> T=0.5,C=0.5   pos2: T,T -> T=1.0
    aln = Alignment(["s1", "s2"], ["ATT", "ACT"])
    pi = pt.codon_frequencies(aln)
    idx = {c: i for i, c in enumerate(SENSE_CODONS)}
    assert pi[idx["ATT"]] == pytest.approx(0.5)
    assert pi[idx["ACT"]] == pytest.approx(0.5)
    assert pi.sum() == pytest.approx(1.0)
    assert (pi >= 0).all()


def test_build_Q_codon_ratios_match_gy94_by_hand():
    pi = np.full(61, 1.0 / 61)
    kappa, omega = 2.5, 0.3
    Q = _build_Q_codon(kappa, omega, pi)
    idx = {c: i for i, c in enumerate(SENSE_CODONS)}
    # TTT(Phe)->TTC(Phe): synonymous transition (T<->C)
    # TTT(Phe)->TTA(Leu): nonsynonymous transition (T<->A is a transversion,
    # so pick a genuine nonsynonymous TRANSITION instead: TTT->CTT (Leu),
    # T<->C at position 0 is a transition
    assert CODON_AA["TTT"] == "F" and CODON_AA["TTC"] == "F"      # synonymous
    assert CODON_AA["TTT"] == "F" and CODON_AA["CTT"] == "L"      # nonsynonymous
    q_syn_ts = Q[idx["TTT"], idx["TTC"]]     # synonymous transition: kappa*pi_j*scale
    q_nonsyn_ts = Q[idx["TTT"], idx["CTT"]]  # nonsynonymous transition: kappa*omega*pi_j*scale
    assert q_nonsyn_ts / q_syn_ts == pytest.approx(omega)
    # TTT->TCT (Ser): synonymous transversion (T<->C? no: position1 T->C is
    # a transition too -- use TTT->TGT (Cys), T<->G transversion, nonsynonymous
    q_nonsyn_tv = Q[idx["TTT"], idx["TGT"]]
    assert q_nonsyn_tv / q_syn_ts == pytest.approx(omega / kappa)
    assert Q.sum(axis=1) == pytest.approx(0.0, abs=1e-9)   # rows sum to zero
    assert (pi * np.diag(Q)).sum() == pytest.approx(-1.0)  # scaled to 1 subst/unit time


# --------------------------------------------------------------------------
# The eigh fix: a permanent regression guard against the omega=1 failure
# --------------------------------------------------------------------------
def test_transition_matrix_is_valid_at_omega_one_with_near_uniform_frequencies():
    # the exact failure regime found during validation: omega=1 (used
    # unconditionally by site class 1, and by class 2 under the branch-site
    # null) combined with close-to-uniform codon frequencies. Before the
    # eigh fix this produced entries as negative as -1.96 for some kappa.
    for kappa in [1.0, 2.0, 3.101, 5.0]:
        for pi in [np.full(61, 1.0 / 61),
                  np.full(61, 1.0 / 61) + np.concatenate([[0.005], np.full(60, -0.005 / 60)])]:
            P = _CodonModel(kappa, 1.0, pi).P(0.3)
            assert np.all(np.isfinite(P)), f"kappa={kappa}: non-finite P(t)"
            assert P.min() >= -1e-9, f"kappa={kappa}: negative transition probability {P.min()}"
            assert P.sum(axis=1) == pytest.approx(1.0, abs=1e-6)


def test_transition_matrix_valid_across_a_kappa_omega_grid():
    rng = np.random.default_rng(0)
    pi = rng.dirichlet(np.full(61, 5.0))
    for kappa in [0.5, 1.0, 2.0, 4.0]:
        for omega in [0.05, 0.3, 1.0, 3.0, 20.0]:
            P = _CodonModel(kappa, omega, pi).P(0.5)
            assert np.all(np.isfinite(P))
            assert P.min() >= -1e-9
            assert P.sum(axis=1) == pytest.approx(1.0, abs=1e-6)


# --------------------------------------------------------------------------
# M0: recovers the simulating kappa/omega
# --------------------------------------------------------------------------
def test_fit_m0_recovers_simulating_parameters():
    kappa_true, omega_true = 3.0, 0.25
    errs_k, errs_w = [], []
    for seed in range(3):
        tree = Tree.from_newick(NEWICK)
        aln = _simulate_branch_site(tree, kappa_true, omega_true, 1.0,
                                    p0=1.0, p1=0.0, n_codons=600, seed=seed)
        res = pt.fit_m0(tree, aln, rounds=4)
        errs_k.append(abs(res["kappa"] - kappa_true))
        errs_w.append(abs(res["omega"] - omega_true))
        assert res["tree"] is not tree   # never mutates the caller's tree
    assert np.mean(errs_k) < 0.5
    assert np.mean(errs_w) < 0.08


# --------------------------------------------------------------------------
# Free-ratio: recovers two different omegas and flags the difference
# --------------------------------------------------------------------------
def test_fit_free_ratio_recovers_two_omegas_and_detects_the_difference():
    for seed in range(3):
        tree = Tree.from_newick(NEWICK)
        aln = _simulate_free_ratio(tree, kappa=3.0, omega_bg=0.2, omega_fg=1.8,
                                   n_codons=800, seed=seed)
        res = pt.fit_free_ratio(tree, aln, FOREGROUND, rounds=4)
        assert res["omega_background"] < 0.5
        assert res["omega_foreground"] > 1.0
        assert res["p"] < 0.01


def test_fit_free_ratio_null_case_no_false_positive():
    # same omega on both -- the LRT should not manufacture a difference
    tree = Tree.from_newick(NEWICK)
    aln = _simulate_free_ratio(tree, kappa=2.5, omega_bg=0.3, omega_fg=0.3,
                               n_codons=600, seed=42)
    res = pt.fit_free_ratio(tree, aln, FOREGROUND, rounds=4)
    assert res["p"] > 0.05


# --------------------------------------------------------------------------
# Branch-site test: detects a strong signal, has real power on a moderate
# one, stays calibrated under the null, and -- the initial-simplex bug's
# regression guard -- does not collapse to the same omega2 regardless of
# what the data say.
# --------------------------------------------------------------------------
# branch-site fits are expensive (a 4-class site mixture, each class its own
# Felsenstein pruning pass, inside a 5-parameter Nelder-Mead with a
# deliberately wide initial simplex); rounds=3 is used throughout below
# rather than the module default of 6 because a direct timing/accuracy probe
# (varying rounds from 3 to 6 on identical data) found IDENTICAL LR/p/omega2
# at every round count -- the fit already converges (its own
# ll-improvement-below-tolerance early exit) well before round 3, so extra
# rounds cost real time for zero change in the answer.
def test_branch_site_test_detects_a_strong_planted_signal():
    tree = Tree.from_newick(NEWICK)
    aln = _simulate_branch_site(tree, kappa=3.0, omega0=0.08, omega2=15.0,
                                p0=0.4, p1=0.1, n_codons=800, seed=1)
    res = pt.branch_site_test(tree, aln, FOREGROUND, rounds=3)
    assert res["p"] < 0.001
    assert res["full"]["omega2"] > 5.0     # not stuck near the null


def test_branch_site_test_has_power_on_a_moderate_signal():
    # n_codons=600 was tried here first and cut power for real (1/3
    # detected, not noise) -- n=1000 is the setting with actual track record
    # (6/6 across two independent batches during validation), so it stays,
    # even though it costs more time than the other tests in this file.
    n_reps, hits = 3, 0
    for seed in range(n_reps):
        tree = Tree.from_newick(NEWICK)
        aln = _simulate_branch_site(tree, kappa=3.0, omega0=0.1, omega2=6.0,
                                    p0=0.5, p1=0.2, n_codons=1000, seed=10 + seed)
        res = pt.branch_site_test(tree, aln, FOREGROUND, rounds=3)
        hits += res["p"] < 0.05
    assert hits >= 2, f"only {hits}/{n_reps} detected"


def test_branch_site_test_type_i_error_near_nominal():
    from phytreon.infer.codon import _CodonModel as _CM
    n_reps, hits = 6, 0
    for rep in range(n_reps):
        tree = Tree.from_newick(NEWICK)
        rng = np.random.default_rng(100 + rep)
        pi = np.full(61, 1.0 / 61)
        model0 = _CM(2.5, 0.2, pi)
        model1 = _CM(2.5, 1.0, pi)
        names = tree.leaf_names()
        seqs = {n: [] for n in names}
        for _ in range(300):
            model = model0 if rng.random() < 0.7 else model1
            state = {tree.root: rng.choice(61, p=pi)}
            for node in tree.traverse("preorder"):
                if node.is_root:
                    continue
                row = model.P(node.length or 0.0)[state[node.parent]]
                state[node] = rng.choice(61, p=row / row.sum())
                if node.is_leaf:
                    seqs[node.name].append(SENSE_CODONS[state[node]])
        aln = Alignment(names=list(names), seqs=["".join(seqs[n]) for n in names])
        res = pt.branch_site_test(tree, aln, FOREGROUND, rounds=3)
        hits += res["p"] < 0.05
    assert hits / n_reps < 0.34, f"false-positive rate {hits}/{n_reps} too high"


def test_branch_site_omega2_tracks_data_not_its_own_starting_point():
    # regression guard for the initial-simplex bug: two datasets with very
    # different true omega2 both used to converge to omega2 ~= 50.5 --
    # exactly the untouched value at the optimiser's own raw=0 starting
    # point -- regardless of which was fitted. They must no longer match.
    tree1 = Tree.from_newick(NEWICK)
    aln1 = _simulate_branch_site(tree1, kappa=3.0, omega0=0.08, omega2=15.0,
                                 p0=0.4, p1=0.1, n_codons=700, seed=1)
    res1 = pt.branch_site_test(tree1, aln1, FOREGROUND, rounds=3)

    tree2 = Tree.from_newick(NEWICK)
    aln2 = _simulate_branch_site(tree2, kappa=3.0, omega0=0.3, omega2=4.0,
                                 p0=0.6, p1=0.2, n_codons=400, seed=99)
    res2 = pt.branch_site_test(tree2, aln2, FOREGROUND, rounds=3)

    assert res1["full"]["omega2"] != pytest.approx(50.5, abs=1.0)
    assert res2["full"]["omega2"] != pytest.approx(50.5, abs=1.0)


# --------------------------------------------------------------------------
# Mechanics and errors
# --------------------------------------------------------------------------
def test_foreground_must_form_a_clade():
    # A and C are both under the ((A,B),C) node, but so is B -- {A, C} is
    # not that node's (or any node's) exact leaf set, so this must be
    # rejected even though get_mrca(["A", "C"]) itself finds a real node
    tree = Tree.from_newick(NEWICK)
    with pytest.raises(ValueError, match="do not form a clade"):
        _foreground_edges(tree, ["A", "C"])


def test_foreground_mrca_cannot_be_the_root():
    tree = Tree.from_newick(NEWICK)
    with pytest.raises(ValueError, match="root"):
        _foreground_edges(tree, tree.leaf_names())


def test_encode_codons_rejects_length_not_a_multiple_of_three():
    aln = Alignment(["a", "b"], ["ATGC", "ATGA"])
    with pytest.raises(ValueError, match="not a multiple of 3"):
        _encode_codons(aln, ["a", "b"])


def test_rejects_mismatched_taxa_between_tree_and_alignment():
    tree = Tree.from_newick(NEWICK)
    aln = Alignment(["A", "B", "C", "D", "E", "X"], ["ATG" * 5] * 6)
    with pytest.raises(ValueError, match="same taxa"):
        pt.fit_m0(tree, aln)
