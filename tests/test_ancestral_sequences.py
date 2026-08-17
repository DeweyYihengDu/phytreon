"""Marginal ML ancestral protein sequence reconstruction.

Three kinds of check, matching how everything statistical in this codebase
gets validated -- never trust a novel pruning/up-pass derivation on its own
say-so:

* the down-pass, in isolation, against `ml_native._site_logliks_aa` -- a
  second, already-published-model-validated (against IQ-TREE2) implementation
  of the identical quantity, computed a structurally different way
  (per-column here, pattern-compressed there).
* the full pipeline's own total log-likelihood against the branch-length
  optimiser's internally reported one -- same tree, same model, same data,
  two different code paths, on the actual public function rather than an
  internal helper.
* recovery on simulated data with a KNOWN true ancestor: reconstruction
  accuracy and reported confidence should closely track each other at low
  divergence, and should both degrade together as data saturates -- the
  behaviour that makes "confidence" mean something rather than just being a
  number in the output.
"""
import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest

import phytreon as pt
from phytreon.comparative.ancestral_sequences import _down_pass

AA = "ARNDCQEGHILKMFPSTWYV"

BALANCED_TREE = (
    "(((A:.05,B:.05):.05,(C:.05,D:.05):.05):.05,"
    "((E:.05,F:.05):.05,(G:.05,H:.05):.05):.05);"
)


def _evolve(newick, length, rate, seed):
    """Evolve a protein down a tree from a random root; independent sites,
    uniform replacement -- same simulator style used throughout this repo's
    other ML-validation tests."""
    tree = pt.Tree.from_newick(newick)
    rng = np.random.default_rng(seed)
    root = "".join(rng.choice(list(AA), length))
    seqs = {}

    def walk(node, parent):
        seq = list(parent)
        p = 1.0 - np.exp(-rate * (node.length or 0.0))
        for i in np.flatnonzero(rng.random(length) < p):
            seq[i] = rng.choice(list(AA))
        seq = "".join(seq)
        if node.is_leaf:
            seqs[node.name] = seq
        for child in node.children:
            walk(child, seq)

    for child in tree.root.children:
        walk(child, root)
    return root, seqs


def _root_of(result):
    return next(n for n in result["tree"].traverse() if n.is_root)


# --------------------------------------------------------------------------
# The down-pass alone, against the IQ-TREE-validated implementation
# --------------------------------------------------------------------------
def test_down_pass_matches_the_iqtree_validated_pruning_exactly():
    from phytreon.infer.ml_native import _encode_aa, _new_model_aa, _site_logliks_aa

    rng = np.random.default_rng(0)
    tr = pt.datasets.random_tree(8, seed=1)
    names = tr.leaf_names()
    ncol = 40
    seqs = {n: "".join(rng.choice(list(AA), ncol)) for n in names}
    aln = pt.Alignment(names, [seqs[n] for n in names])

    model = _new_model_aa("LG", gamma=4)
    enc_names, enc_states, weights, freqs = _encode_aa(aln)
    # this test's validity depends on no two columns being identical, so the
    # pattern-compressed encoding preserves original column order 1:1
    assert enc_states.shape[1] == ncol, "columns collapsed -- pick a different seed"
    assert enc_names == names

    s2i = {c: i for i, c in enumerate(AA)}
    states = np.array([[s2i[ch] for ch in seqs[n]] for n in names])
    idx = {n: i for i, n in enumerate(names)}

    for rate in (0.3, 1.0, 2.5):
        _D, _logscale, mine = _down_pass(tr, model, idx, states, rate)
        theirs = _site_logliks_aa(tr, model, enc_names, enc_states, rate=rate)
        assert np.allclose(mine, theirs, atol=1e-9), f"mismatch at rate={rate}"


# --------------------------------------------------------------------------
# The full pipeline's log-likelihood, two ways
# --------------------------------------------------------------------------
def test_reported_loglik_matches_the_branch_length_optimiser():
    # fitted_logLik comes from the pattern-compressed optimiser that fit the
    # branch lengths; logLik comes from this module's own full per-column
    # pruning over the fitted tree. Same model, same data -- must agree.
    root, seqs = _evolve(BALANCED_TREE, 200, 3.0, 0)
    names = sorted(seqs)
    aln = pt.Alignment(names, [seqs[n] for n in names])
    res = pt.reconstruct_ancestral_sequences(
        pt.Tree.from_newick(BALANCED_TREE), aln, model="LG", gamma=4, fit_model=True)
    assert res["fitted_logLik"] == pytest.approx(res["logLik"], abs=0.05)


# --------------------------------------------------------------------------
# Recovery and calibration on simulated data with a known ancestor
# --------------------------------------------------------------------------
def test_root_recovery_is_high_at_low_divergence_and_confidence_tracks_it():
    true_root, seqs = _evolve(BALANCED_TREE, 200, 3.0, 0)
    names = sorted(seqs)
    aln = pt.Alignment(names, [seqs[n] for n in names])
    res = pt.reconstruct_ancestral_sequences(
        pt.Tree.from_newick(BALANCED_TREE), aln, gamma=4)
    root = _root_of(res)
    recon = res["sequences"][root.name]
    conf = res["confidence"][root.name]

    identity = np.mean([a == b for a, b in zip(recon, true_root)])
    assert identity > 0.85
    # confidence and identity should be close, not just both "high" --
    # that is what makes the reported number a genuine calibration
    assert abs(float(conf.mean()) - identity) < 0.1


def test_confidence_is_calibrated_high_confidence_sites_are_more_often_right():
    true_root, seqs = _evolve(BALANCED_TREE, 250, 4.0, 3)
    names = sorted(seqs)
    aln = pt.Alignment(names, [seqs[n] for n in names])
    res = pt.reconstruct_ancestral_sequences(
        pt.Tree.from_newick(BALANCED_TREE), aln, gamma=4)
    root = _root_of(res)
    recon = res["sequences"][root.name]
    conf = res["confidence"][root.name]
    correct = np.array([a == b for a, b in zip(recon, true_root)])

    high = conf >= 0.9
    low = conf < 0.7
    assert high.sum() >= 10 and low.sum() >= 5, "test data too easy/hard to be informative"
    assert correct[high].mean() > correct[low].mean()


def test_recovery_degrades_towards_equilibrium_as_divergence_saturates():
    low_root, low_seqs = _evolve(BALANCED_TREE, 200, 3.0, 0)
    high_root, high_seqs = _evolve(BALANCED_TREE, 200, 60.0, 1)
    names = sorted(low_seqs)

    def identity(root_seq, seqs):
        aln = pt.Alignment(names, [seqs[n] for n in names])
        res = pt.reconstruct_ancestral_sequences(
            pt.Tree.from_newick(BALANCED_TREE), aln, gamma=4)
        root = _root_of(res)
        recon = res["sequences"][root.name]
        return np.mean([a == b for a, b in zip(recon, root_seq)])

    low_identity = identity(low_root, low_seqs)
    high_identity = identity(high_root, high_seqs)
    assert low_identity > 0.85
    assert high_identity < low_identity - 0.3


# --------------------------------------------------------------------------
# Mechanics: masking, gap detection, non-mutation, node naming
# --------------------------------------------------------------------------
def _short_random_alignment(seed=2, ncol=30, gap_col=10):
    names = list("ABCDEFGH")
    rng = np.random.default_rng(seed)
    seqs = {n: "".join(rng.choice(list(AA), ncol)) for n in names}
    if gap_col is not None:
        seqs = {n: s[:gap_col] + "-" + s[gap_col + 1:] for n, s in seqs.items()}
    return pt.Alignment(names, [seqs[n] for n in names])


def test_all_gap_columns_are_detected_and_reported():
    aln = _short_random_alignment(gap_col=10)
    res = pt.reconstruct_ancestral_sequences(
        pt.Tree.from_newick(BALANCED_TREE), aln, gamma=0)
    assert res["all_gap_columns"] == [10]


def test_no_gap_columns_reports_an_empty_list():
    aln = _short_random_alignment(gap_col=None)
    res = pt.reconstruct_ancestral_sequences(
        pt.Tree.from_newick(BALANCED_TREE), aln, gamma=0)
    assert res["all_gap_columns"] == []


def test_mask_below_replaces_uncertain_calls_with_x():
    root, seqs = _evolve(BALANCED_TREE, 200, 3.0, 0)
    names = sorted(seqs)
    aln = pt.Alignment(names, [seqs[n] for n in names])
    unmasked = pt.reconstruct_ancestral_sequences(
        pt.Tree.from_newick(BALANCED_TREE), aln, gamma=4)
    masked = pt.reconstruct_ancestral_sequences(
        pt.Tree.from_newick(BALANCED_TREE), aln, gamma=4, mask_below=0.9)
    root_name = _root_of(masked).name
    seq_u = unmasked["sequences"][root_name]
    seq_m = masked["sequences"][root_name]
    conf = unmasked["confidence"][root_name]
    for i, (u, m, c) in enumerate(zip(seq_u, seq_m, conf)):
        if c < 0.9:
            assert m == "X"
        else:
            assert m == u


def test_input_tree_is_never_mutated_and_internal_nodes_get_stable_names():
    aln = _short_random_alignment(gap_col=None)
    caller_tree = pt.Tree.from_newick(BALANCED_TREE)
    before = {frozenset(n.leaf_names()): n.length for n in caller_tree.traverse()}

    result = pt.reconstruct_ancestral_sequences(caller_tree, aln, gamma=0, fit_model=True)

    after = {frozenset(n.leaf_names()): n.length for n in caller_tree.traverse()}
    assert before == after

    internal_names = [n.name for n in result["tree"].traverse() if not n.is_leaf]
    assert all(internal_names)                      # none empty or None
    assert len(internal_names) == len(set(internal_names))  # all distinct


def test_fit_model_false_preserves_branch_lengths_exactly():
    aln = _short_random_alignment(gap_col=None)
    caller_tree = pt.Tree.from_newick(BALANCED_TREE)
    before = {frozenset(n.leaf_names()): n.length for n in caller_tree.traverse()
             if n.length is not None}

    result = pt.reconstruct_ancestral_sequences(caller_tree, aln, gamma=4, fit_model=False)

    after = {frozenset(n.leaf_names()): n.length for n in result["tree"].traverse()
             if n.length is not None}
    assert set(before) == set(after)
    for k in before:
        assert before[k] == pytest.approx(after[k], abs=1e-12)
    assert result["fitted_logLik"] is None


def test_gamma_shape_override_is_never_re_estimated():
    aln = _short_random_alignment(gap_col=None)
    result = pt.reconstruct_ancestral_sequences(
        pt.Tree.from_newick(BALANCED_TREE), aln, gamma=4,
        fit_model=True, gamma_shape=0.77)
    assert result["gamma_shape"] == 0.77


# --------------------------------------------------------------------------
# Errors, and the FASTA-export convenience
# --------------------------------------------------------------------------
def test_rejects_an_unknown_model():
    aln = _short_random_alignment(gap_col=None)
    with pytest.raises(ValueError, match="model must be one of"):
        pt.reconstruct_ancestral_sequences(
            pt.Tree.from_newick(BALANCED_TREE), aln, model="NOPE")


def test_rejects_tree_alignment_taxon_mismatch():
    names = list("ABCDEFG") + ["Ghost"]
    rng = np.random.default_rng(9)
    seqs = ["".join(rng.choice(list(AA), 20)) for _ in names]
    aln = pt.Alignment(names, seqs)
    with pytest.raises(ValueError, match="same taxa"):
        pt.reconstruct_ancestral_sequences(pt.Tree.from_newick(BALANCED_TREE), aln)


def test_ancestral_alignment_wraps_reconstructed_sequences_for_export():
    aln = _short_random_alignment(gap_col=None)
    result = pt.reconstruct_ancestral_sequences(
        pt.Tree.from_newick(BALANCED_TREE), aln, gamma=0, fit_model=True)
    internal_names = [n.name for n in result["tree"].traverse() if not n.is_leaf]

    anc_only = pt.ancestral_alignment(result)
    assert anc_only.nseq == len(internal_names)
    assert anc_only.ncol == aln.ncol
    assert set(anc_only.names) == set(internal_names)

    with_tips = pt.ancestral_alignment(result, include_tips=True)
    assert with_tips.nseq == len(internal_names) + 8

    one = pt.ancestral_alignment(result, nodes=[internal_names[0]])
    assert one.names == [internal_names[0]]
    assert one.seqs[0] == result["sequences"][internal_names[0]]

    # round-trips through the alignment machinery's own FASTA writer
    fasta = anc_only.to_fasta()
    assert fasta.startswith(">")
    assert all(name in fasta for name in internal_names)


def test_ancestral_alignment_rejects_an_unknown_node_name():
    aln = _short_random_alignment(gap_col=None)
    result = pt.reconstruct_ancestral_sequences(
        pt.Tree.from_newick(BALANCED_TREE), aln, gamma=0)
    with pytest.raises(ValueError, match="not internal node names"):
        pt.ancestral_alignment(result, nodes=["not_a_real_node"])
