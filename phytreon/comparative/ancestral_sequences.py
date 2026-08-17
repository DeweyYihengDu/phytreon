"""Marginal maximum-likelihood reconstruction of ancestral protein sequences.

:mod:`ace` reconstructs a discrete or continuous *trait* at internal nodes.
This reconstructs a *sequence* -- the residue at every alignment column, for
every internal node -- under an empirical amino-acid substitution model
(JTT/WAG/LG, the same models :func:`~phytreon.build_tree` and
:func:`~phytreon.infer.ml_tree` use to score a protein tree) with +Gamma
rate variation across sites.

This is the step that turns a tree into something a wet lab can act on: pick
a node (the root, or the MRCA of a clade of interest -- from
:func:`~phytreon.comparative.domains.compare_domain_trees` or
:func:`~phytreon.infer.gene_tree_conflict`, an ancestor is exactly what a
resurrection experiment needs), pull its reconstructed sequence and its
per-site confidence, and decide whether it is trustworthy enough to
synthesize.

Reuses the pruning/optimisation engine :func:`~phytreon.infer.ml_tree`
already uses for tree search (:mod:`phytreon.infer.ml_native`), which is
validated against IQ-TREE2's log-likelihood on the same data -- rather than a
second, independent implementation of amino-acid pruning that could disagree
with the first one silently.

**Known limitation, stated rather than hidden**: this reconstructs the
residue at each *existing* alignment column. It does not reconstruct
insertion/deletion history -- whether a column was present at all in a given
ancestor -- which needs a separate gap model (as FastML and similar tools
provide) and is not implemented here. A reconstructed sequence is exactly
``aln.ncol`` positions long; deciding which columns to keep for a
synthesis-ready sequence (e.g. dropping columns gapped in most descendants of
that node) is left to the caller.
"""
from __future__ import annotations

from typing import Dict, Iterable, Optional

from ..core.tree import Node, Tree
from ..infer.matrix import Alignment


def _down_pass(tree: Tree, model, idx: Dict[str, int], states, rate: float):
    """Felsenstein pruning, vectorised over alignment columns.

    Returns ``(D, logscale, site_loglik)``: ``D[node]`` is a rescaled
    ``(ncol, 20)`` array of down-partial likelihoods (rescaled per node per
    column for numerical stability, exactly as
    :func:`~phytreon.infer.ml_native._site_logliks_aa` does -- the rescaling
    is a per-column scalar shared by all 20 states, so it cancels exactly out
    of any state *ratio* computed downstream and does not bias the
    reconstruction). ``logscale`` is the accumulated log of those scalars, and
    ``site_loglik`` is the true per-column log-likelihood at this rate,
    ``logscale`` added back in.
    """
    import numpy as np
    ncol = states.shape[1]
    D: Dict[Node, "np.ndarray"] = {}
    logscale: Dict[Node, "np.ndarray"] = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf:
            row = states[idx[node.name]]
            L = np.ones((ncol, 20))
            known = row >= 0
            L[known] = 0.0
            L[known, row[known]] = 1.0
            D[node] = L
            logscale[node] = np.zeros(ncol)
        else:
            L = np.ones((ncol, 20))
            scal = np.zeros(ncol)
            for c in node.children:
                P = model.P(max(c.length or 0.0, 1e-9) * rate)
                L = L * (D[c] @ P.T)
                scal = scal + logscale[c]
            m = L.max(axis=1)
            m = np.where(m > 0, m, 1.0)
            L = L / m[:, None]
            scal = scal + np.log(m)
            D[node] = L
            logscale[node] = scal
    root = tree.root
    site_loglik = np.log((D[root] * model.pi[None, :]).sum(axis=1) + 1e-300) + logscale[root]
    return D, logscale, site_loglik


def _up_pass(tree: Tree, model, D: Dict[Node, "np.ndarray"], rate: float,  # noqa: F821
            ncol: int):
    """The complementary pass: for every node, the likelihood of everything
    *outside* its own subtree, as a function of its own state.

    Mirrors :func:`~phytreon.comparative.ace.ace_ml`'s single-site up-pass
    exactly (``U[node] = P(t).T @ msg`` in column-vector form), only batched
    over columns -- see the module's derivation in the accompanying tests for
    why the batched form of a transposed column-vector product is
    ``msg @ P`` rather than ``msg @ P.T``.
    """
    import numpy as np
    U: Dict[Node, "np.ndarray"] = {tree.root: np.tile(model.pi, (ncol, 1))}
    for node in tree.traverse("preorder"):
        if node.is_root:
            continue
        parent = node.parent
        msg = U[parent].copy()
        for sib in parent.children:
            if sib is node:
                continue
            P_sib = model.P(max(sib.length or 0.0, 1e-9) * rate)
            msg = msg * (D[sib] @ P_sib.T)
        P_self = model.P(max(node.length or 0.0, 1e-9) * rate)
        U[node] = msg @ P_self
    return U


def reconstruct_ancestral_sequences(
    tree: Tree, aln: Alignment, model: str = "LG", gamma: int = 4,
    fit_model: bool = True, gamma_shape: Optional[float] = None,
    mask_below: Optional[float] = None,
) -> Dict[str, object]:
    """Marginal ML ancestral protein sequences at every internal node of ``tree``.

    ``model`` is ``"JTT"``/``"WAG"``/``"LG"``. ``gamma`` is the number of
    discrete +G rate categories (as in :func:`~phytreon.build_tree`'s
    ``ml_gamma``); real protein sites vary enormously in how fast they evolve,
    and ignoring that does not just lose resolution -- it makes reconstructed
    conserved sites *overconfident*, which is the wrong direction to err in
    before deciding what to synthesize. Default 4, not 0.

    ``fit_model=True`` (default) refits ``tree``'s branch lengths and the
    gamma shape by ML under this exact model, on a **copy** -- the tree
    passed in is never mutated. Sensible when ``tree`` came from
    ``method="nj"`` or was built under a different model. Pass ``False`` to
    keep externally-fit branch lengths exactly as given instead -- the right
    choice if ``tree`` already came from ``build_tree(..., method="ml",
    ml_engine="iqtree", ...)`` with a site-heterogeneous model
    (LG+C60/PMSF), since refitting under phytreon's own site-homogeneous
    LG/WAG/JTT would silently pull the tree's metric structure toward the
    weaker model.

    ``gamma_shape``, if given, is used as-is and never re-estimated
    regardless of ``fit_model`` -- for reusing a shape already fit
    externally (branch lengths can still be refit around it). Without it,
    shape is fit together with branch lengths when ``fit_model=True``, or
    left at an uninformative 0.5 when ``fit_model=False``.

    ``mask_below``, if given, replaces the reconstructed residue with ``X``
    at any column whose posterior probability for its best state falls below
    this threshold, rather than reporting a confident-looking letter for an
    unconfident call.

    Returns a dict:

    ``tree``
        a copy of the input tree (branch lengths possibly refit); unnamed
        internal nodes are given stable names (``"anc0"``, ``"anc1"``, ...)
        so results can be indexed by name. Also carries
        ``node.data["asr_sequence"]`` and ``["asr_mean_confidence"]``.
    ``sequences``
        ``{name: str}`` for every node, tips (observed, as given) and
        internal (reconstructed) alike -- pass to
        ``Alignment(list(sequences), list(sequences.values()))`` and then
        ``.to_fasta()`` to hand a specific ancestor to synthesis.
    ``posterior``
        ``{internal_node_name: (ncol, 20) array}``, columns in
        :data:`~phytreon.infer.aa_models.AA_STATES` order.
    ``confidence`` / ``mean_confidence``
        per-site and per-node summaries of the posterior at the
        reconstructed state -- read before trusting a sequence.
    ``logLik``
        the model's total log-likelihood on this tree and alignment, computed
        directly from the same per-column pruning used for reconstruction
        (independent of, and checked in the test suite against, the
        pattern-compressed value ``fit_model=True`` reports internally).
    ``all_gap_columns``
        alignment columns with no data at all across every taxon in ``tree``
        -- their "reconstruction" is really just the model's own equilibrium
        frequencies reshaped by the tree, not evidence, and is reported so it
        is not mistaken for one.
    """
    import numpy as np
    from ..infer.aa_models import AA_MODELS, AA_STATES
    from ..infer.ml_native import _alphabet_mismatch, _encode_aa, _new_model_aa

    name = model.upper()
    if name not in AA_MODELS:
        raise ValueError(
            f"reconstruct_ancestral_sequences: model must be one of "
            f"{sorted(AA_MODELS)}, not {model!r}"
        )
    err = _alphabet_mismatch(aln, want_protein=True)
    if err:
        raise ValueError(f"reconstruct_ancestral_sequences: {err}")

    tree_taxa = set(tree.leaf_names())
    aln_taxa = set(aln.names)
    if tree_taxa != aln_taxa:
        only_tree = sorted(tree_taxa - aln_taxa)
        only_aln = sorted(aln_taxa - tree_taxa)
        raise ValueError(
            "reconstruct_ancestral_sequences: tree and alignment must have "
            f"the same taxa. In the tree but not the alignment: {only_tree[:10]}"
            f"{' ...' if len(only_tree) > 10 else ''}. In the alignment but "
            f"not the tree: {only_aln[:10]}{' ...' if len(only_aln) > 10 else ''}"
        )
    if aln.ncol < 1:
        raise ValueError("reconstruct_ancestral_sequences: alignment has no columns")

    work = Tree.from_newick(tree.write())
    counter = 0
    for node in work.traverse("postorder"):
        if not node.is_leaf and not node.name:
            node.name = f"anc{counter}"
            counter += 1

    model_obj = _new_model_aa(name, gamma)
    if gamma_shape is not None:
        model_obj.set_shape(gamma_shape)

    fitted_loglik = None
    if fit_model:
        data = _encode_aa(aln)
        from ..infer.ml_native import _optimize_branches, _optimize_model
        prev = -1e18
        ll = prev
        for _ in range(8):
            ll = _optimize_branches(work, model_obj, data, rounds=4)
            if gamma_shape is None:
                ll = _optimize_model(work, model_obj, data)
            ll = _optimize_branches(work, model_obj, data, rounds=2)
            if ll - prev < 1e-3:
                break
            prev = ll
        fitted_loglik = ll

    leaves = work.leaf_names()
    idx = {n: i for i, n in enumerate(leaves)}
    aln_idx = {n: i for i, n in enumerate(aln.names)}
    s2i = {c: i for i, c in enumerate(AA_STATES)}
    ncol = aln.ncol
    states = np.full((len(leaves), ncol), -1, dtype=np.int8)
    for i, leaf_name in enumerate(leaves):
        seq = aln.seqs[aln_idx[leaf_name]]
        for j, ch in enumerate(seq.upper()):
            states[i, j] = s2i.get(ch, -1)
    all_gap_columns = sorted(np.flatnonzero((states < 0).all(axis=0)).tolist())

    rates, wts = model_obj.rate_categories()
    ncat = len(rates)
    internal_nodes = [n for n in work.traverse() if not n.is_leaf]

    cat_ll = np.zeros((ncat, ncol))
    per_cat_D = []
    for ci, r in enumerate(rates):
        D, _logscale, site_ll = _down_pass(work, model_obj, idx, states, r)
        cat_ll[ci] = site_ll
        per_cat_D.append(D)

    if ncat > 1:
        from scipy.special import logsumexp
        logw = np.log(wts)[:, None] + cat_ll
        logZ = logsumexp(logw, axis=0)
        cat_post = np.exp(logw - logZ[None, :])
    else:
        logZ = cat_ll[0]
        cat_post = np.ones((1, ncol))

    posterior_sum = {n: np.zeros((ncol, 20)) for n in internal_nodes}
    for ci, r in enumerate(rates):
        D = per_cat_D[ci]
        U = _up_pass(work, model_obj, D, r, ncol)
        weight_col = cat_post[ci][:, None]
        for node in internal_nodes:
            post = D[node] * U[node]
            s = post.sum(axis=1, keepdims=True)
            s = np.where(s > 0, s, 1.0)
            posterior_sum[node] += (post / s) * weight_col

    sequences: Dict[str, str] = {
        leaf_name: aln.seqs[aln_idx[leaf_name]].upper() for leaf_name in leaves
    }
    posterior_out: Dict[str, "np.ndarray"] = {}
    confidence: Dict[str, "np.ndarray"] = {}
    for node in internal_nodes:
        post = posterior_sum[node]
        best = post.argmax(axis=1)
        conf = post.max(axis=1)
        chars = [AA_STATES[b] for b in best]
        if mask_below is not None:
            chars = [c if cv >= mask_below else "X" for c, cv in zip(chars, conf)]
        seq_str = "".join(chars)
        sequences[node.name] = seq_str
        posterior_out[node.name] = post
        confidence[node.name] = conf
        node.data["asr_sequence"] = seq_str
        node.data["asr_mean_confidence"] = float(conf.mean())

    return {
        "tree": work,
        "sequences": sequences,
        "posterior": posterior_out,
        "confidence": confidence,
        "mean_confidence": {n.name: float(confidence[n.name].mean())
                            for n in internal_nodes},
        "model": name,
        "gamma_shape": float(model_obj.shape) if model_obj.ncat > 1 else None,
        "n_gamma_categories": ncat,
        "logLik": float(logZ.sum()),
        "fitted_logLik": fitted_loglik,
        "all_gap_columns": all_gap_columns,
        "fit_model": fit_model,
    }


def ancestral_alignment(result: Dict[str, object],
                        nodes: Optional[Iterable[str]] = None,
                        include_tips: bool = False) -> Alignment:
    """Package reconstructed (and optionally observed) sequences from
    :func:`reconstruct_ancestral_sequences` as an :class:`Alignment`, so the
    existing :meth:`Alignment.to_fasta` writes them out -- rather than a
    second, parallel FASTA writer.

    ``nodes`` selects which internal nodes to include (by name, matching
    ``result["tree"]``'s node names); by default, all of them. Set
    ``include_tips=True`` to add the observed sequences alongside the
    reconstructed ones, e.g. for a combined tips+ancestors FASTA to plot
    against the tree.
    """
    tree = result["tree"]
    internal_names = [n.name for n in tree.traverse() if not n.is_leaf]
    wanted = list(nodes) if nodes is not None else internal_names
    unknown = sorted(set(wanted) - set(internal_names))
    if unknown:
        raise ValueError(
            f"ancestral_alignment: {unknown} are not internal node names in "
            f"this result; have {sorted(internal_names)[:10]}"
            f"{' ...' if len(internal_names) > 10 else ''}"
        )
    names = list(wanted)
    if include_tips:
        names = names + tree.leaf_names()
    seqs = [result["sequences"][n] for n in names]
    return Alignment(names, seqs)
