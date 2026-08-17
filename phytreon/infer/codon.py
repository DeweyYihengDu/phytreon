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

from typing import Dict, List, Optional, Sequence, Tuple

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
        # Q is reversible (pi_i * Q_ij == pi_j * Q_ji: the GY94 rate between
        # two codons depends only on their unordered relationship --
        # transition/transversion, synonymous/nonsynonymous -- never on
        # direction), so it is similar to a SYMMETRIC matrix via
        # S = diag(sqrt(pi)) @ Q @ diag(1/sqrt(pi)). Decomposing S with eigh
        # rather than Q with the general eig is not just faster: eigh always
        # returns an orthogonal eigenbasis, so it stays exact even when
        # eigenvalues coincide, whereas eig's eigenvector matrix can become
        # ill-conditioned or fail outright at a true degeneracy.
        #
        # That degeneracy is not a corner case here -- it is guaranteed at
        # omega=1, which every branch-site fit (:func:`branch_site_test`)
        # visits unconditionally (site class 1 always; class 2 as well under
        # the null it tests against). At omega=1 the model stops
        # distinguishing synonymous from nonsynonymous changes, so with
        # near-uniform codon frequencies it reduces to three i.i.d. copies of
        # the same per-position nucleotide process -- found via a validation
        # run whose omega2 estimate collapsed to its own starting value no
        # matter the data: the old eig-based P(t) was silently returning
        # transition "probabilities" as negative as -1.96 there (confirmed by
        # direct inspection, kappa=3.101, uniform pi), which propagated to
        # NaN log-likelihoods and corrupted the optimizer's every step.
        import numpy as np
        Q = _build_Q_codon(self.kappa, self.omega, self.pi)
        sqrt_pi = np.sqrt(self.pi)
        S = (sqrt_pi[:, None] * Q) / sqrt_pi[None, :]
        S = (S + S.T) / 2.0   # exact symmetry, clearing float round-off
        vals, U = np.linalg.eigh(S)
        self.vals = vals
        self.vecs = U / sqrt_pi[:, None]
        self.vinv = U.T * sqrt_pi[None, :]

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
    # get_mrca returns the graph-theoretic common ancestor regardless of
    # whether OTHER leaves also descend from it, so "mrca is None" alone
    # does not catch a non-clade foreground (e.g. two taxa from a 3-leaf
    # clade) -- found via a test that named such a pair and got a silent,
    # wrong stem edge back instead of the documented error. The exact-leaf-
    # set check below is what actually enforces "forms a clade".
    if mrca is None or set(mrca.leaf_names()) != set(foreground):
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


# --------------------------------------------------------------------------
# M0: one omega for the whole tree
# --------------------------------------------------------------------------
def _optimize_branches_codon(tree: Tree, branch_model: Dict[Node, "_CodonModel"],  # noqa: F821
                             names, states, pi, rounds: int = 3) -> float:
    from scipy.optimize import minimize_scalar
    edges = [n for n in tree.traverse() if not n.is_root]
    best = float(_site_logliks_codon(tree, names, states, branch_model, pi).sum())
    for _ in range(rounds):
        improved = False
        for node in edges:
            old = node.length or 0.1

            def neg(t, _node=node):
                _node.length = float(t)
                return -float(_site_logliks_codon(tree, names, states, branch_model, pi).sum())

            res = minimize_scalar(neg, bounds=(1e-6, 10.0), method="bounded")
            if -res.fun > best + 1e-6:
                node.length = float(res.x)
                best = -res.fun
                improved = True
            else:
                node.length = old
        if not improved:
            break
    return best


def fit_m0(tree: Tree, aln: Alignment, kappa: Optional[float] = None,
          omega: Optional[float] = None, fit_model: bool = True,
          rounds: int = 6) -> Dict[str, object]:
    """M0: a single kappa and omega shared by the whole tree.

    The baseline codon model -- is this gene under selection *on average*,
    with no attempt to localise it to particular branches or sites. Fits
    branch lengths, kappa and omega jointly by ML on a **copy** of ``tree``
    (the input is never mutated), unless ``kappa``/``omega`` are given, in
    which case that one is held fixed rather than estimated (branch lengths
    are still fit either way when ``fit_model=True``).

    Equilibrium codon frequencies are F3x4, estimated once from ``aln`` and
    not re-optimised. Returns the fitted tree, ``kappa``, ``omega``,
    ``logLik``, and ``codon_frequencies``.
    """
    from scipy.optimize import minimize

    _validate_codon_inputs(tree, aln)
    work = Tree.from_newick(tree.write())
    work_names = work.leaf_names()
    states = _encode_codons(aln, work_names)
    pi = codon_frequencies(aln)

    k0 = kappa if kappa is not None else 2.0
    w0 = omega if omega is not None else 0.4
    model = _CodonModel(k0, w0, pi)
    branch_model = {n: model for n in work.traverse() if not n.is_root}

    fit_kappa = kappa is None
    fit_omega = omega is None
    ll = float(_site_logliks_codon(work, work_names, states, branch_model, pi).sum())
    if fit_model:
        prev = -1e18
        for _ in range(rounds):
            ll = _optimize_branches_codon(work, branch_model, work_names, states, pi, rounds=2)
            if fit_kappa or fit_omega:
                x0 = ([model.kappa] if fit_kappa else []) + ([model.omega] if fit_omega else [])

                def neg(x):
                    xi = iter(x)
                    k = next(xi) if fit_kappa else model.kappa
                    w = next(xi) if fit_omega else model.omega
                    if k <= 0 or w <= 0:
                        return 1e18
                    model.set_params(k, w)
                    return -float(_site_logliks_codon(work, work_names, states,
                                                       branch_model, pi).sum())

                res = minimize(neg, x0, method="Nelder-Mead",
                              options={"xatol": 1e-4, "fatol": 1e-4, "maxiter": 200})
                xi = iter(res.x)
                model.set_params(next(xi) if fit_kappa else model.kappa,
                                 next(xi) if fit_omega else model.omega)
                ll = -float(res.fun)
            if ll - prev < 1e-3:
                break
            prev = ll
        ll = float(_site_logliks_codon(work, work_names, states, branch_model, pi).sum())

    return {"tree": work, "kappa": model.kappa, "omega": model.omega,
            "logLik": ll, "codon_frequencies": pi, "n_codons": states.shape[1]}


# --------------------------------------------------------------------------
# Free-ratio: a second omega for a labelled set of foreground branches
# --------------------------------------------------------------------------
def fit_free_ratio(tree: Tree, aln: Alignment, foreground: Sequence[str],
                   rounds: int = 6) -> Dict[str, object]:
    """Two omegas -- one for the branches leading to ``foreground``'s MRCA,
    one for the rest of the tree -- tested against :func:`fit_m0` by a
    likelihood-ratio test (1 df, an ordinary chi-square: unlike the
    branch-site test below, omega is not on a boundary here, so no mixture
    correction is needed).

    ``foreground`` is a set of taxon names; the single branch immediately
    ancestral to their MRCA is what gets the second omega (the standard,
    simplest labelling -- does selection differ on the branch where this
    lineage originates). kappa and branch lengths are shared across both
    omega classes.

    Returns the two-omega fit (``tree``, ``kappa``, ``omega_foreground``,
    ``omega_background``, ``logLik``) plus the comparison against M0:
    ``m0`` (that fit's own result dict), ``LR``, ``p``.
    """
    from scipy.optimize import minimize

    _validate_codon_inputs(tree, aln)
    m0 = fit_m0(tree, aln, rounds=rounds)

    work = Tree.from_newick(tree.write())
    work_names = work.leaf_names()
    states = _encode_codons(aln, work_names)
    pi = m0["codon_frequencies"]
    fg_stem = _foreground_edges(work, foreground)

    kappa = m0["kappa"]
    model_bg = _CodonModel(kappa, m0["omega"], pi)
    model_fg = _CodonModel(kappa, m0["omega"], pi)
    branch_model = {n: (model_fg if n in fg_stem else model_bg)
                    for n in work.traverse() if not n.is_root}
    # start from M0's own FITTED branch lengths, not the caller's original
    # (possibly NJ, unoptimised) ones -- already a good starting point
    for n, m0n in zip(work.traverse(), m0["tree"].traverse()):
        if not n.is_root:
            n.length = m0n.length

    def total_ll():
        return float(_site_logliks_codon(work, work_names, states, branch_model, pi).sum())

    prev = -1e18
    ll = total_ll()
    for _ in range(rounds):
        ll = _optimize_branches_codon(work, branch_model, work_names, states, pi, rounds=2)

        def neg(x):
            k, w_bg, w_fg = x
            if k <= 0 or w_bg <= 0 or w_fg <= 0:
                return 1e18
            model_bg.set_params(k, w_bg)
            model_fg.set_params(k, w_fg)
            return -total_ll()

        res = minimize(neg, [model_bg.kappa, model_bg.omega, model_fg.omega],
                       method="Nelder-Mead",
                       options={"xatol": 1e-4, "fatol": 1e-4, "maxiter": 300})
        model_bg.set_params(res.x[0], res.x[1])
        model_fg.set_params(res.x[0], res.x[2])
        ll = -float(res.fun)
        if ll - prev < 1e-3:
            break
        prev = ll

    from scipy.stats import chi2
    lr = max(2.0 * (ll - m0["logLik"]), 0.0)
    p = float(chi2.sf(lr, df=1))
    return {"tree": work, "kappa": model_bg.kappa,
            "omega_background": model_bg.omega, "omega_foreground": model_fg.omega,
            "logLik": ll, "m0": m0, "LR": lr, "p": p}


# --------------------------------------------------------------------------
# Branch-site test (Zhang, Nielsen & Yang 2005's corrected Model A)
# --------------------------------------------------------------------------
def _mixture_logliks(tree: Tree, names, states, pi,
                     class_weights: Sequence[float],
                     class_branch_models: Sequence[Dict[Node, "_CodonModel"]]  # noqa: F821
                     ) -> "np.ndarray":  # noqa: F821
    """Per-codon log-likelihood under a mixture of site classes, each with
    its own (possibly per-branch-varying) model -- one down-pass per class,
    combined by log-sum-exp so a class with negligible weight cannot
    numerically swamp one that fits far better."""
    import numpy as np
    from scipy.special import logsumexp
    per_class = np.stack([
        _site_logliks_codon(tree, names, states, bm, pi) + np.log(max(w, 1e-300))
        for w, bm in zip(class_weights, class_branch_models)
    ])
    return logsumexp(per_class, axis=0)


def _class_weights(p0: float, p1: float) -> Tuple[float, float, float, float]:
    """The four site-class proportions from Zhang, Nielsen & Yang's own
    parametrisation: class 2's total mass (1 - p0 - p1) split between the
    2a/2b sub-classes in the same ratio as p0:p1 themselves."""
    p2 = max(1.0 - p0 - p1, 0.0)
    denom = p0 + p1
    if denom <= 0:
        return p0, p1, p2 / 2, p2 / 2
    return p0, p1, p2 * p0 / denom, p2 * p1 / denom


def _fit_branch_site_model(tree: Tree, names, states, pi, fg_stem,
                           fix_omega2: bool, rounds: int) -> Dict[str, object]:
    from scipy.optimize import minimize

    import numpy as np

    kappa = 2.0
    omega0 = 0.3
    omega2_raw = 0.0             # omega2 = 1 + exp(omega2_raw), see below
    p0_raw, p1_raw = 0.5, 0.5    # stick-breaking: see props() below

    def props(p0_raw, p1_raw):
        p0 = 1.0 / (1.0 + np.exp(-p0_raw))
        p1 = (1.0 - p0) / (1.0 + np.exp(-p1_raw))
        return _class_weights(p0, p1)

    OMEGA2_MAX = 100.0   # generous: no biologically real dN/dS ratio approaches this

    def omega2_of(raw):
        # Mapped smoothly into (1, OMEGA2_MAX) via a logistic squeeze rather
        # than a hard omega2 >= 1 rejection wall OR an unbounded 1 + exp(raw):
        #
        # - a hard wall at 1 traps Nelder-Mead exactly there with nowhere
        #   feasible to step towards (found first: simulated genuine positive
        #   selection, got back exactly omega2 = 1.0, LR = 0.0 -- not noise,
        #   identical every replicate);
        # - removing the wall entirely (1 + exp(raw), unbounded above) traps
        #   it differently: when the data barely identify omega2 (few sites in
        #   the omega2-using classes, short foreground branch, or -- found
        #   second -- a joint ridge where several parameters move together and
        #   omega2 stops mattering to the likelihood at all), Nelder-Mead can
        #   drift along that flat direction to a nonsensical omega2 in the
        #   billions with no convergence check ever flagging it, since the
        #   likelihood is not improving, just failing to get worse.
        #
        # A finite upper bound closes off that drift while keeping the same
        # smooth, wall-free landscape near 1 -- and 100 costs nothing on
        # genuine signal: a likelihood surface that actually wants omega2 = 20
        # is not going to prefer 100 to get there, it is only the *directionless*
        # drift on already-flat data this stops.
        return 1.0 + (OMEGA2_MAX - 1.0) / (1.0 + np.exp(-raw))

    model0 = _CodonModel(kappa, omega0, pi)               # class 0's uniform model
    model1 = _CodonModel(kappa, 1.0, pi)                  # class 1's uniform model (omega=1)
    model2 = _CodonModel(kappa, omega2_of(omega2_raw), pi)  # foreground-only, 2a/2b

    def branch_models():
        all_bg = {n: model0 for n in tree.traverse() if not n.is_root}
        all_bg1 = {n: model1 for n in tree.traverse() if not n.is_root}
        mix2a = {n: (model2 if n in fg_stem else model0)
                for n in tree.traverse() if not n.is_root}
        mix2b = {n: (model2 if n in fg_stem else model1)
                for n in tree.traverse() if not n.is_root}
        return [all_bg, all_bg1, mix2a, mix2b]

    def total_ll(p0_raw, p1_raw):
        w0, w1, w2a, w2b = props(p0_raw, p1_raw)
        return float(_mixture_logliks(tree, names, states, pi,
                                      [w0, w1, w2a, w2b], branch_models()).sum())

    prev = -1e18
    ll = total_ll(p0_raw, p1_raw)
    for _ in range(rounds):
        # branch lengths: optimise against the CURRENT mixture, one edge at a
        # time, exactly as _optimize_branches_codon does for a single model
        from scipy.optimize import minimize_scalar
        edges = [n for n in tree.traverse() if not n.is_root]
        best = total_ll(p0_raw, p1_raw)
        for node in edges:
            old = node.length or 0.1

            def neg(t, _node=node):
                _node.length = float(t)
                return -total_ll(p0_raw, p1_raw)

            res = minimize_scalar(neg, bounds=(1e-6, 10.0), method="bounded")
            if -res.fun > best + 1e-6:
                node.length = float(res.x)
                best = -res.fun
            else:
                node.length = old
        ll = best

        # model parameters: kappa, omega0, p0_raw, p1_raw always; omega2_raw too
        # unless omega2 is fixed at 1 (the null model) -- omega0 and kappa keep
        # a hard rejection wall rather than their own reparametrisation since,
        # unlike omega2, neither ever needed to move away from its own starting
        # boundary in validation (the true values sat in the interior)
        x0 = [kappa, omega0, p0_raw, p1_raw] + ([] if fix_omega2 else [omega2_raw])

        def neg(x):
            k, w0_, pr0, pr1 = x[:4]
            w2 = omega2_of(x[4]) if not fix_omega2 else 1.0
            if k <= 0 or w0_ <= 0 or w0_ > 1.0:
                return 1e18
            model0.set_params(k, w0_)
            model1.set_params(k, 1.0)
            model2.set_params(k, w2)
            return -total_ll(pr0, pr1)

        # Explicit initial simplex: scipy's default step for a coordinate
        # starting at exactly 0.0 (omega2_raw) is tiny, and validation showed
        # it left that axis unexplored -- two very different simulated
        # datasets (one with an enormous p2=0.5, omega2=15 signal) both
        # converged to omega2 ~= 50.5, exactly omega2_of(raw=0), i.e. the
        # untouched starting point, not a fitted value. A deliberately
        # generous per-parameter step replaces scipy's own collapsing default.
        step = np.array([1.0, 0.3, 1.5, 1.5, 3.0][:len(x0)])
        simplex = np.vstack([x0] + [np.array(x0) + step * np.eye(len(x0))[i]
                                    for i in range(len(x0))])
        res = minimize(neg, x0, method="Nelder-Mead",
                       options={"xatol": 1e-4, "fatol": 1e-4, "maxiter": 800,
                               "initial_simplex": simplex})
        kappa, omega0, p0_raw, p1_raw = res.x[:4]
        omega2_raw = res.x[4] if not fix_omega2 else 0.0
        omega2 = omega2_of(omega2_raw) if not fix_omega2 else 1.0
        model0.set_params(kappa, omega0)
        model1.set_params(kappa, 1.0)
        model2.set_params(kappa, omega2)
        ll = -float(res.fun)
        if ll - prev < 1e-3:
            break
        prev = ll

    w0, w1, w2a, w2b = props(p0_raw, p1_raw)
    return {"kappa": kappa, "omega0": omega0, "omega2": omega2,
            "p0": w0, "p1": w1, "p2a": w2a, "p2b": w2b, "logLik": ll}


def branch_site_test(tree: Tree, aln: Alignment, foreground: Sequence[str],
                     rounds: int = 6) -> Dict[str, object]:
    """The corrected branch-site test for positive selection on specific
    codons along specific branches (Zhang, Nielsen & Yang 2005's Model A,
    fixing the excess false positives of Yang & Nielsen (2002)'s original).

    Four site classes, mixed: class 0 (proportion ``p0``) has ``0 < omega0 <
    1`` on every branch; class 1 (``p1``) has ``omega=1`` fixed on every
    branch; classes 2a/2b (splitting the remaining ``1 - p0 - p1`` in the
    same ratio as ``p0:p1``) share ``omega0``/``1`` respectively on
    background branches but switch to a single shared ``omega2`` on the
    foreground branches -- so evidence for positive selection is evidence
    that *some sites*, not necessarily most of them, shifted to ``omega2``
    specifically on the labelled lineage.

    ``foreground`` labels the single stem branch to its taxa's MRCA, as in
    :func:`fit_free_ratio`. The null model fixes ``omega2 = 1``; the
    likelihood-ratio test against it uses a 50:50 mixture of a point mass at
    0 and a chi-square(1 df) null distribution rather than a plain
    chi-square, because omega2 = 1 sits on the *boundary* of the alternative
    model's parameter space (omega2 >= 1) rather than in its interior --
    the same boundary-mixture logic used for
    :func:`~phytreon.comparative.pagels_lambda`'s test elsewhere in this
    package, here in the direction that matters: the plain chi-square this
    replaces would be anti-conservative, not merely imprecise, because it
    demands a smaller LR than the true null requires.

    Returns the full model's fit, the null model's fit, ``LR``, and ``p``.
    """
    _validate_codon_inputs(tree, aln)
    work = Tree.from_newick(tree.write())
    work_names = work.leaf_names()
    states = _encode_codons(aln, work_names)
    pi = codon_frequencies(aln)
    fg_stem = _foreground_edges(work, foreground)

    full = _fit_branch_site_model(work, work_names, states, pi, fg_stem,
                                  fix_omega2=False, rounds=rounds)

    work_null = Tree.from_newick(tree.write())
    for n, wn in zip(work_null.traverse(), work.traverse()):
        if not n.is_root:
            n.length = wn.length
    fg_stem_null = _foreground_edges(work_null, foreground)
    null = _fit_branch_site_model(work_null, work_names, states, pi, fg_stem_null,
                                  fix_omega2=True, rounds=rounds)

    from scipy.stats import chi2
    lr = max(2.0 * (full["logLik"] - null["logLik"]), 0.0)
    p = float(0.5 * chi2.sf(lr, df=1))
    return {"full": full, "null": null, "LR": lr, "p": p, "tree": work}
