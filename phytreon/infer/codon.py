"""Codon-based selection tests: is a gene, or specific branches of its tree,
under positive selection?

Every other model in this package treats a DNA or protein alignment as
letters evolving independently of what they encode. This is the one place
that reads DNA as *codons* and asks the question those other models cannot:
whether amino-acid-changing (nonsynonymous) substitutions are more or less
common than silent (synonymous) ones, and specifically whether that ratio --
omega, dN/dS -- exceeds 1 on particular branches, the signature of positive
selection rather than the pervasive omega<1 of purifying selection nearly
everywhere else in a genome.

Three tests, all built on the same Goldman & Yang (1994) codon substitution
model and the same Felsenstein pruning machinery
:func:`~phytreon.infer.ml_native.ml_tree` already uses for nucleotide and
protein trees, generalised from 4/20 states to 61 (the sense codons of the
standard genetic code, stop codons excluded):

``fit_m0``
    one omega for the whole tree -- is the gene under selection at all, on
    average.
``fit_free_ratio``
    two omegas, one for a chosen set of "foreground" branches and one for
    everything else, tested against ``fit_m0`` by a likelihood-ratio test --
    does *this specific lineage* show a different selective regime.
``branch_site_test``
    the corrected branch-site test (Zhang, Nielsen & Yang 2005 -- the
    published fix for Yang & Nielsen (2002)'s original version, which that
    later paper showed gives excessive false positives): does *some* codons
    on the foreground branches, not necessarily the whole gene, show
    positive selection episodically.

**Scope, stated plainly.** This is M0 / free-ratio / the corrected
branch-site test only -- not the full PAML suite (no site models M1a/M2a/M7/
M8, no Bayes empirical Bayes site posterior identification, no clade models).
Equilibrium codon frequencies are F3x4 (the product of empirical
position-specific nucleotide frequencies, Goldman & Yang's own default) --
not F61 or F1x4.
"""
from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

from ..core.tree import Node, Tree
from .matrix import Alignment

# --------------------------------------------------------------------------
# The standard genetic code
# --------------------------------------------------------------------------
_BASES = "TCAG"
_AA_TABLE = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)
CODONS = [a + b + c for a in _BASES for b in _BASES for c in _BASES]
CODON_AA = dict(zip(CODONS, _AA_TABLE))
STOP_CODONS = frozenset(c for c, aa in CODON_AA.items() if aa == "*")
SENSE_CODONS = [c for c in CODONS if c not in STOP_CODONS]
assert len(SENSE_CODONS) == 61

_TRANSITIONS = {frozenset("AG"), frozenset("CT")}


def _nt_diff(a: str, b: str) -> List[Tuple[int, str, str]]:
    return [(i, x, y) for i, (x, y) in enumerate(zip(a, b)) if x != y]


def _is_transition(x: str, y: str) -> bool:
    return frozenset((x, y)) in _TRANSITIONS


# --------------------------------------------------------------------------
# Codon frequencies (F3x4) and the GY94 rate matrix
# --------------------------------------------------------------------------
def codon_frequencies(aln: Alignment, method: str = "F3x4") -> "np.ndarray":  # noqa: F821
    """Equilibrium frequencies of the 61 sense codons (:data:`SENSE_CODONS`
    order), estimated from the alignment.

    ``method="F3x4"`` (the only one implemented, and Goldman & Yang's own
    default): the product of the empirical nucleotide frequency at each of
    the three codon positions, taken across every codon in the alignment
    (gaps and incomplete codons excluded), renormalised over the 61 sense
    codons once stop-codon combinations are dropped -- the product of three
    independent position frequencies puts some mass on codons that would be
    stops, which cannot appear in a coding sequence and so cannot be counted
    as part of the 61-codon probability simplex.
    """
    import numpy as np
    if method != "F3x4":
        raise ValueError(f"codon_frequencies: unknown method {method!r}, only 'F3x4'")
    pos_counts = [dict.fromkeys(_BASES, 0) for _ in range(3)]
    for seq in aln.seqs:
        seq = seq.upper().replace("U", "T")
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            if codon in CODON_AA and codon not in STOP_CODONS:
                for p in range(3):
                    if codon[p] in pos_counts[p]:
                        pos_counts[p][codon[p]] += 1
    pos_freq = []
    for p in range(3):
        total = sum(pos_counts[p].values())
        if total == 0:
            raise ValueError("codon_frequencies: no complete sense codons found")
        pos_freq.append({b: pos_counts[p][b] / total for b in _BASES})
    raw = np.array([pos_freq[0][c[0]] * pos_freq[1][c[1]] * pos_freq[2][c[2]]
                    for c in SENSE_CODONS])
    return raw / raw.sum()


def _build_Q_codon(kappa: float, omega: float, pi: "np.ndarray"  # noqa: F821
                   ) -> "np.ndarray":  # noqa: F821
    """The Goldman & Yang (1994) 61x61 rate matrix.

    ``q_ij`` (codons i != j, differing at exactly one nucleotide position):
    ``pi_j`` for a synonymous transversion, ``kappa * pi_j`` for a
    synonymous transition, ``omega * pi_j`` for a nonsynonymous
    transversion, ``omega * kappa * pi_j`` for a nonsynonymous transition;
    0 for any pair differing at more than one position. Scaled so the
    expected number of substitutions per codon per unit branch length is 1,
    the same convention :func:`~phytreon.infer.ml_native._build_Q_aa` uses.
    """
    import numpy as np
    n = len(SENSE_CODONS)
    Q = np.zeros((n, n))
    for i, ci in enumerate(SENSE_CODONS):
        for j, cj in enumerate(SENSE_CODONS):
            if i == j:
                continue
            diffs = _nt_diff(ci, cj)
            if len(diffs) != 1:
                continue
            _, x, y = diffs[0]
            rate = kappa if _is_transition(x, y) else 1.0
            if CODON_AA[ci] != CODON_AA[cj]:
                rate *= omega
            Q[i, j] = rate * pi[j]
    for i in range(n):
        Q[i, i] = -Q[i].sum()
    scale = -(pi * np.diag(Q)).sum()
    if scale <= 0:
        raise ValueError("codon model: degenerate rate matrix (scale <= 0)")
    return Q / scale


class _CodonModel:
    """Eigendecomposed GY94 model, cached the same way
    :class:`~phytreon.infer.ml_native._Model` caches a nucleotide or protein
    one -- built fresh whenever ``kappa``/``omega`` change, reused across
    every branch-length evaluation in between."""

    def __init__(self, kappa: float, omega: float, pi: "np.ndarray"):  # noqa: F821
        self.kappa = kappa
        self.omega = omega
        self.pi = pi
        self._decompose()

    def _decompose(self):
        import numpy as np
        Q = _build_Q_codon(self.kappa, self.omega, self.pi)
        vals, vecs = np.linalg.eig(Q)
        self.vals = vals.real
        self.vecs = vecs.real
        self.vinv = np.linalg.inv(self.vecs)

    def set_params(self, kappa: float, omega: float):
        self.kappa = kappa
        self.omega = omega
        self._decompose()

    def P(self, t: float) -> "np.ndarray":  # noqa: F821
        import numpy as np
        return (self.vecs * np.exp(self.vals * max(t, 1e-9))) @ self.vinv


# --------------------------------------------------------------------------
# Encoding an alignment as codons, and the generalised pruning likelihood
# --------------------------------------------------------------------------
def _encode_codons(aln: Alignment, names: Sequence[str]) -> "np.ndarray":  # noqa: F821
    import numpy as np
    if aln.ncol % 3 != 0:
        raise ValueError(
            f"codon model: alignment length {aln.ncol} is not a multiple of 3"
        )
    idx_taxon = {n: i for i, n in enumerate(aln.names)}
    ncodon = aln.ncol // 3
    codon_idx = {c: i for i, c in enumerate(SENSE_CODONS)}
    states = np.full((len(names), ncodon), -1, dtype=np.int16)
    for i, name in enumerate(names):
        seq = aln.seqs[idx_taxon[name]].upper().replace("U", "T")
        for j in range(ncodon):
            states[i, j] = codon_idx.get(seq[3 * j:3 * j + 3], -1)
    return states


def _site_logliks_codon(tree: Tree, names: Sequence[str], states: "np.ndarray",  # noqa: F821
                        branch_model: Dict[Node, "_CodonModel"],  # noqa: F821
                        pi: "np.ndarray") -> "np.ndarray":  # noqa: F821
    """Felsenstein pruning with a possibly different model *per branch* --
    the generalisation :func:`~phytreon.infer.ml_native._site_logliks_aa`
    does not need (one shared model everywhere) but the free-ratio and
    branch-site models do (foreground branches use a different omega).
    ``branch_model[node]`` is the model to use for the branch *above*
    ``node``; the root's own entry is never read.
    """
    import numpy as np
    idx = {n: i for i, n in enumerate(names)}
    npat = states.shape[1]
    cache_L: Dict[int, "np.ndarray"] = {}
    cache_s: Dict[int, "np.ndarray"] = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf:
            L = np.ones((npat, 61))
            row = states[idx[node.name]]
            known = row >= 0
            L[known] = 0.0
            L[known, row[known]] = 1.0
            cache_L[id(node)] = L
            cache_s[id(node)] = np.zeros(npat)
        else:
            L = np.ones((npat, 61))
            scal = np.zeros(npat)
            for c in node.children:
                P = branch_model[c].P(c.length or 0.0)
                L = L * (cache_L[id(c)] @ P.T)
                scal = scal + cache_s[id(c)]
            m = L.max(axis=1)
            m = np.where(m > 0, m, 1.0)
            L = L / m[:, None]
            scal = scal + np.log(m)
            cache_L[id(node)] = L
            cache_s[id(node)] = scal
    root = tree.root
    site_ll = np.log((cache_L[id(root)] * pi[None, :]).sum(axis=1) + 1e-300)
    return site_ll + cache_s[id(root)]


def _foreground_edges(tree: Tree, foreground: Sequence[str]) -> set:
    """The single stem edge leading to the MRCA of ``foreground``'s taxa, as
    a one-element set of :class:`Node` (the node whose *own* branch is that
    stem edge) -- the standard, simplest branch-site labelling: does
    selection act on the branch where this lineage originates.
    """
    mrca = tree.get_mrca(list(foreground))
    if mrca is None:
        raise ValueError(
            f"codon model: {list(foreground)} do not form a clade in this tree "
            f"(no single node has exactly this leaf set)"
        )
    if mrca.is_root:
        raise ValueError(
            "codon model: foreground taxa's MRCA is the tree's root, which has "
            "no branch of its own to label as foreground"
        )
    return {mrca}


def _validate_codon_inputs(tree: Tree, aln: Alignment) -> List[str]:
    tree_taxa = set(tree.leaf_names())
    aln_taxa = set(aln.names)
    if tree_taxa != aln_taxa:
        only_tree = sorted(tree_taxa - aln_taxa)
        only_aln = sorted(aln_taxa - tree_taxa)
        raise ValueError(
            "codon model: tree and alignment must have the same taxa. In the "
            f"tree but not the alignment: {only_tree[:10]}. In the alignment "
            f"but not the tree: {only_aln[:10]}"
        )
    return tree.leaf_names()
