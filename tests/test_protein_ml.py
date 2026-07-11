"""Protein (amino acid) maximum-likelihood support.

This is purely additive on top of the pre-existing nucleotide-only native
ML engine (:mod:`phytreon.infer.ml_native`): every nucleotide default and
code path is unchanged (checked explicitly below), and protein support
requires opting in with an explicit model name.

There is no external ML tool (IQ-TREE/RAxML/PhyML) available in this
environment to cross-validate against, so correctness is instead checked
by two independent, install-free methods:

* the model's ``P(t) = exp(Qt)`` (computed here via eigendecomposition,
  see :class:`phytreon.infer.ml_native._ModelAA`) is compared against
  ``scipy.linalg.expm(Qt)`` -- a different numerical algorithm for the
  same matrix exponential;
* the tree log-likelihood from the pruning/dynamic-programming traversal
  (:func:`phytreon.infer.ml_native._site_logliks_aa`) is compared against
  a brute-force nested sum over ancestral states, computed independently
  of the traversal/caching code.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from phytreon.core.tree import Node, Tree
from phytreon.infer.align import Alignment
from phytreon.infer.aa_models import AA_MODELS, AA_STATES
from phytreon.infer.ml_native import (
    ml_tree, log_likelihood, model_finder,
    _ModelAA, _build_Q_aa, _new_model_aa, _site_logliks_aa, _AA_S2I,
)
from phytreon.infer.bootstrap import distance_matrix_model, bootstrap_support
from phytreon.infer.pipeline import build_tree

# Two well-separated clusters (A1-3, B1-3) plus an outgroup (C1-2) --
# realistic enough that ML should recover the two clusters as clades.
PROT_SEQS = Alignment(
    names=["A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2"],
    seqs=[
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKI",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILARVGDGTQDNLSGAEKAVQVKV",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILARVGDGTQDNLSGAEKAVQVKV",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILARVGDGTQDNLSGAEKAVQVKW",
        "MKTAYIAKQRQMSFVKSHFSRQLEERLGLIEVQAPILWRVGDGTQDNLSGAEKAVQVYV",
        "MKTAYIAKQRQMSFVKSHFSRQLEERLGLIEVQAPILWRVGDGTQDNLSGAEKAVQVYA",
    ],
)

NUC_SEQS = Alignment(
    names=["T1", "T2", "T3", "T4"],
    seqs=[
        "ACGTACGTACGTACGTACGTACGTACGTACGT",
        "ACGTACGTACGTACGTACGTACGTACGTACGT",
        "ACGAACGTACGTACGTACGTACGTACGTACGA",
        "ACGAACGTACGTACGTACGTACGTACGTACGA",
    ],
)


# --------------------------------------------------------------------------
# model data sanity
# --------------------------------------------------------------------------
def test_aa_model_matrices_are_well_formed():
    for name, (R, pi) in AA_MODELS.items():
        assert len(R) == 20 and all(len(row) == 20 for row in R)
        assert len(pi) == 20
        assert abs(sum(pi) - 1.0) < 1e-4, name
        for i in range(20):
            assert R[i][i] == 0.0
            for j in range(20):
                assert R[i][j] == R[j][i], f"{name} not symmetric at ({i},{j})"
    assert len(AA_STATES) == 20


# --------------------------------------------------------------------------
# independent check #1: eigendecomposition P(t) vs. scipy.linalg.expm
# --------------------------------------------------------------------------
@pytest.mark.parametrize("model", ["LG", "WAG", "JTT"])
@pytest.mark.parametrize("t", [0.01, 0.1, 0.5, 1.5])
def test_transition_matrix_matches_scipy_expm(model, t):
    from scipy.linalg import expm
    _, pi = AA_MODELS[model]
    Q = _build_Q_aa(model, np.array(pi))
    mdl = _new_model_aa(model, gamma=0)
    P_eig = mdl.P(t)
    P_expm = expm(Q * t)
    assert P_eig == pytest.approx(P_expm, abs=1e-6)
    # a valid transition matrix: rows are probability distributions
    assert P_eig.sum(axis=1) == pytest.approx(np.ones(20), abs=1e-6)
    assert (P_eig >= -1e-8).all()


# --------------------------------------------------------------------------
# independent check #2: pruning traversal vs. brute-force ancestral sum
# --------------------------------------------------------------------------
def _brute_force_loglik(model: _ModelAA, seqs: dict, tA, tB, tC, tX) -> float:
    """Directly compute log P(site) for ((A:tA,B:tB):tX,C:tC) by summing
    over the two unobserved ancestral states, with no traversal/caching --
    an independent path through the same math the pruning code implements."""
    PA, PB, PC, PX = model.P(tA), model.P(tB), model.P(tC), model.P(tX)
    a, b, c = (_AA_S2I[seqs[k]] for k in ("A", "B", "C"))
    total = 0.0
    for r in range(20):
        for x in range(20):
            total += model.pi[r] * PX[r, x] * PA[x, a] * PB[x, b] * PC[r, c]
    return math.log(total)


def _small_tree(tA, tB, tC, tX) -> Tree:
    root = Node(name="root")
    x = Node(name="X", length=tX)
    root.add_child(x)
    root.add_child(Node(name="C", length=tC))
    x.add_child(Node(name="A", length=tA))
    x.add_child(Node(name="B", length=tB))
    return Tree(root=root)


@pytest.mark.parametrize("model", ["LG", "WAG", "JTT"])
def test_pruning_matches_brute_force_ancestral_sum(model):
    mdl = _new_model_aa(model, gamma=0)
    tA, tB, tC, tX = 0.05, 0.12, 0.30, 0.08
    tree = _small_tree(tA, tB, tC, tX)
    obs = {"A": "L", "B": "I", "C": "V"}
    names = ["A", "B", "C"]
    states = np.array([[_AA_S2I[obs[n]]] for n in names])  # one site column
    ll_pruning = _site_logliks_aa(tree, mdl, names, states)[0]
    ll_brute = _brute_force_loglik(mdl, obs, tA, tB, tC, tX)
    assert ll_pruning == pytest.approx(ll_brute, abs=1e-8)


# --------------------------------------------------------------------------
# alphabet-mismatch guard: prevents silent nonsense results
# --------------------------------------------------------------------------
def test_ml_tree_rejects_mismatched_alphabet():
    with pytest.raises(ValueError, match="protein"):
        ml_tree(NUC_SEQS, model="LG")
    with pytest.raises(ValueError, match="nucleotide"):
        ml_tree(PROT_SEQS, model="HKY85")


def test_log_likelihood_rejects_nucleotide_model_on_protein_data():
    tree = ml_tree(PROT_SEQS, model="LG", search=False)
    with pytest.raises(ValueError, match="protein"):
        log_likelihood(tree, PROT_SEQS, model="GTR")


# --------------------------------------------------------------------------
# protein ML end-to-end
# --------------------------------------------------------------------------
def test_native_ml_recovers_groups_on_protein():
    tree = ml_tree(PROT_SEQS, model="LG", search=True)
    a_clade = {"A1", "A2", "A3"}
    b_clade = {"B1", "B2", "B3"}
    leaf_sets = [frozenset(n.leaf_names()) for n in tree.traverse() if not n.is_leaf]
    assert any(set(s) == a_clade for s in leaf_sets)
    assert any(set(s) == b_clade for s in leaf_sets)
    assert tree.data["model"].startswith("LG")
    assert math.isfinite(tree.data["logL"])


def test_ml_bootstrap_support_on_protein():
    tree = ml_tree(PROT_SEQS, model="JTT", search=False)

    def builder(a):
        return ml_tree(a, model="JTT", search=False)

    ref, support = bootstrap_support(PROT_SEQS, builder=builder, n=5, reference=tree, seed=0)
    assert all(0 <= v <= 100 for v in support.values())


# --------------------------------------------------------------------------
# additive, not modified: nucleotide defaults are unchanged
# --------------------------------------------------------------------------
def test_build_tree_ml_for_protein_needs_explicit_model():
    # No "auto" default -- build_tree's ml_model default stays "HKY85"
    # (unchanged), so calling it on protein data without an explicit
    # ml_model must fail exactly like ml_tree does directly.
    with pytest.raises(ValueError, match="protein"):
        build_tree(PROT_SEQS, aligner="none", method="ml", ml_search=False)


def test_build_tree_ml_default_still_hky85_for_nucleotide():
    tree = build_tree(NUC_SEQS, aligner="none", method="ml", ml_search=False)
    assert tree.data["model"].startswith("HKY85")


def test_build_tree_ml_for_protein_with_explicit_model():
    tree = build_tree(PROT_SEQS, aligner="none", method="ml",
                      ml_model="LG", ml_search=False, bootstrap=0)
    assert tree.data["model"].startswith("LG")


def test_distance_matrix_model_poisson_is_an_explicit_opt_in():
    # unchanged: "jc69" (the default) falls back to raw p-distance on
    # protein data -- it does NOT auto-upgrade to poisson.
    _, d_jc = distance_matrix_model(PROT_SEQS, model="jc69")
    _, d_raw = distance_matrix_model(PROT_SEQS, model="raw")
    assert d_jc == d_raw
    _, d_poisson = distance_matrix_model(PROT_SEQS, model="poisson")
    assert d_poisson != d_raw


def test_build_tree_nj_on_protein_with_explicit_poisson_distance():
    tree = build_tree(PROT_SEQS, aligner="none", method="nj", dist_model="poisson")
    assert set(tree.leaf_names()) == set(PROT_SEQS.names)


def test_model_finder_defaults_to_protein_models_on_protein_data():
    rows = model_finder(PROT_SEQS, gammas=(0,))
    assert {r["model"] for r in rows} == {"JTT", "WAG", "LG"}
    rows_nuc = model_finder(NUC_SEQS, gammas=(0,))
    assert {r["model"] for r in rows_nuc} == {"JC69", "K80", "HKY85", "GTR"}
