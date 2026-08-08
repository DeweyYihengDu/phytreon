"""One-call sequences -> tree pipeline, every stage configurable.

    align  ->  trim ("cut")  ->  infer  ->  bootstrap

Each stage is opt-in/opt-out and forwards a kwargs dict, so the same entry
point covers "quick NJ from raw sequences" and "trim hard, ML, 1000 boots".
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Union

from ..core.tree import Tree
from .align import Alignment, align, align_external, read_fasta

Records = List[Tuple[str, str]]
SeqInput = Union[str, Records, Alignment]


def _to_records(seqs: SeqInput) -> Union[Records, Alignment]:
    if isinstance(seqs, Alignment):
        return seqs
    if isinstance(seqs, str):
        return read_fasta(seqs)
    return list(seqs)


_METHODS = ("nj", "upgma", "ml", "parsimony", "mp")


def build_tree(sequences: SeqInput, *,
               aligner: str = "builtin",
               align_kw: Optional[Dict] = None,
               trim_kw: Optional[Dict] = None,
               method: str = "nj",
               dist_model: str = "jc69",
               nj_engine: str = "builtin",
               constraint: Optional[Union[Dict[str, object], Tree, str]] = None,
               root: Union[str, List[str]] = "none",
               bootstrap: int = 0,
               ml_engine: str = "native",
               ml_model: str = "HKY85",
               ml_gamma: int = 0,
               ml_search: bool = True,
               parsimony_model: str = "fitch",
               seed: Optional[int] = None,
               return_alignment: bool = False
               ) -> Union[Tree, Tuple[Tree, Alignment]]:
    """Build a tree from sequences.

    Parameters
    ----------
    sequences   FASTA path/string, list of ``(name, seq)``, or an
                pre-built :class:`Alignment`.
    aligner     ``"builtin"`` (pure Python), ``"mafft"``/``"muscle"``
                (external), or ``"none"`` (input already aligned).
    align_kw    forwarded to the aligner (``seqtype``, ``match``, ``gap`` ...).
    trim_kw     ``None`` to skip trimming, else forwarded to
                :func:`phytreon.infer.trim.trim` (``max_gap``,
                ``min_occupancy``, ``min_conservation`` ...).
    method      ``"nj"`` | ``"upgma"`` | ``"ml"`` | ``"parsimony"``/``"mp"``.
    nj_engine   ``"builtin"`` (Biopython's textbook NJ, O(n^3)) or
                ``"rapidnj"`` (:func:`phytreon.infer.ml.infer_rapidnj`,
                external) for ``method="nj"`` on alignments too large for
                the builtin engine -- somewhere in the low thousands of
                tips, exactly the range a large 16S ASV table reaches.
                RapidNJ computes its own distance matrix (a smaller model
                menu than ``dist_model``'s) and is incompatible with
                ``constraint`` (:func:`~phytreon.infer.distance.
                constrained_nj` needs phytreon's own matrix directly).
    root        ``"none"``, ``"midpoint"``, or a tip name / list of tip names
                to root on as the outgroup (:func:`phytreon.treeops.
                outgroup_root`) -- the more defensible choice whenever one is
                known, since midpoint rooting assumes every lineage drifts at
                about the same rate and an outgroup does not.
    bootstrap   number of bootstrap replicates (0 = none). Works under
                ``constraint`` too (every replicate is rebuilt under the same
                constraint, not without it) and for ``ml_engine="iqtree"``/
                ``"raxml-ng"`` (each engine's own built-in bootstrap, not a
                rebuild-from-scratch loop). No effect with
                ``ml_engine="fasttree"``, which reports its own per-branch
                support automatically instead of taking a replicate count.
    constraint  Force a taxonomy grouping (e.g. genus) to come out
                monophyletic, as a ``{tip_name: label}`` mapping (a tip
                missing or mapped to ``None`` is left free). Two different
                strengths, depending on ``method``:

                * ``method="nj"`` runs :func:`phytreon.infer.distance.
                  constrained_nj`: NJ inside each group, then NJ again to
                  place the groups (and any free tip) relative to each
                  other. This *forces* every group monophyletic -- there is
                  no way for the result to disagree, even where the
                  sequence data would.
                * ``method="ml"`` with ``ml_engine="iqtree"`` or
                  ``"raxml-ng"`` instead runs a *constrained* search: the
                  mapping is turned into polytomies (see :func:`phytreon.
                  infer.constraint.constraint_tree`) and passed to IQ-TREE's
                  ``-g`` / RAxML-NG's ``--tree-constraint``, which still
                  resolve everything -- within a group, between groups,
                  where an unlisted tip goes -- by likelihood, and only
                  refuse to break the groups named. Pass a :class:`~phytreon.
                  core.tree.Tree` or a constraint-file path here directly to
                  skip that conversion. Native ML, FastTree, RapidNJ, and the
                  distance/parsimony methods other than NJ have no
                  constrained search to hook into and reject ``constraint``.

    Protein sequences work the same way -- pass ``ml_model="JTT"``/``"WAG"``/
    ``"LG"`` for ``method="ml"`` (:func:`phytreon.infer.ml_native.ml_tree`
    validates the model matches the data's alphabet); ``dist_model="poisson"``
    for ``method="nj"``/``"upgma"``.

    Single-cell CRISPR lineage-tracing character matrices (see
    :func:`phytreon.infer.lineage.read_allele_table`) work through
    ``method="parsimony"`` too -- pass ``parsimony_model="camin_sokal"`` for
    the irreversible model (:func:`phytreon.infer.lineage.lineage_tree`)
    instead of the default reversible Fitch parsimony
    (:func:`phytreon.infer.parsimony.parsimony_tree`).
    """
    if method not in _METHODS:
        raise ValueError(f"unknown method {method!r}; choose one of {_METHODS}")
    data = _to_records(sequences)

    # 1. alignment ------------------------------------------------------
    if isinstance(data, Alignment) or aligner == "none":
        aln = data if isinstance(data, Alignment) else Alignment(
            [n for n, _ in data], [s for _, s in data])
    elif aligner == "builtin":
        aln = align(data, **(align_kw or {}))
    elif aligner in ("mafft", "muscle"):
        aln = align_external(data, tool=aligner, **(align_kw or {}))
    else:
        raise ValueError(f"unknown aligner {aligner!r}")

    # 2. trim ("cut") ---------------------------------------------------
    if trim_kw is not None:
        from .trim import trim
        aln = trim(aln, **trim_kw)

    # 3. inference ------------------------------------------------------
    if method == "ml":
        if ml_engine == "native":
            if constraint is not None:
                raise ValueError(
                    "the native ML engine has no constrained search to hook "
                    "into -- use ml_engine='iqtree' for a constrained search, "
                    "or method='nj' to force the grouping outright"
                )
            from .ml_native import ml_tree as _native_ml
            tree = _native_ml(aln, model=ml_model, gamma=ml_gamma, search=ml_search)
        else:                                   # external engine
            from .ml import infer_ml
            _constrainable = ("iqtree", "iqtree2", "raxml-ng", "raxmlng", "raxml")
            ml_kw = {}
            # Each engine's own bootstrap, not phytreon's generic resampling
            # loop in step 4 below: that loop rebuilds the tree from scratch
            # per replicate, which for an external ML engine means one full
            # subprocess search per replicate -- fine for the native engine,
            # a non-starter for a few hundred external searches. FastTree has
            # no equivalent flag (it reports its own per-branch support
            # automatically instead), so it gets neither.
            if bootstrap and ml_engine in _constrainable:
                ml_kw["bootstrap"] = bootstrap
            if constraint is not None:
                if ml_engine not in _constrainable:
                    raise ValueError(
                        f"ml_engine={ml_engine!r} has no constrained search to "
                        "hook into -- use ml_engine='iqtree' or 'raxml-ng', or "
                        "method='nj' to force the grouping outright"
                    )
                if isinstance(constraint, dict):
                    from .constraint import constraint_tree
                    constraint = constraint_tree(constraint)
                ml_kw["constraint"] = constraint
            tree = infer_ml(aln, tool=ml_engine, **ml_kw)
    elif method in ("parsimony", "mp"):
        if constraint is not None:
            raise ValueError(
                "parsimony search has no constrained mode here -- use "
                "method='nj' to force the grouping outright, or "
                "method='ml' with ml_engine='iqtree' for a constrained search"
            )
        if parsimony_model == "camin_sokal":
            from .lineage import lineage_tree
            tree = lineage_tree(aln, search=ml_search)
        else:
            from .parsimony import parsimony_tree
            tree = parsimony_tree(aln, search=ml_search)
    elif method == "nj" and nj_engine != "builtin":
        if constraint is not None:
            raise ValueError(
                "constrained_nj needs phytreon's own distance matrix "
                "directly -- use nj_engine='builtin' (the default) to force "
                "the grouping outright, or method='ml' with an external "
                "ml_engine for a constrained search"
            )
        if nj_engine != "rapidnj":
            raise ValueError(f"unknown nj_engine {nj_engine!r}; use "
                             "'builtin' or 'rapidnj'")
        from .ml import infer_rapidnj
        tree = infer_rapidnj(aln)
    else:                                    # method in ("nj", "upgma")
        from .bootstrap import distance_matrix_model
        from .distance import neighbor_joining, upgma
        names, D = distance_matrix_model(aln, dist_model)
        if constraint is None:
            tree = neighbor_joining(names, D) if method == "nj" else upgma(names, D)
        elif method == "upgma":
            raise ValueError(
                "constrained UPGMA is not implemented -- use method='nj' to "
                "force the grouping outright, or method='ml' with "
                "ml_engine='iqtree' for a constrained search"
            )
        elif not isinstance(constraint, dict):
            raise TypeError(
                "constrained_nj needs a {tip_name: label} mapping, not a "
                "constraint tree/file -- build the grouping dict yourself, "
                "or use method='ml' with ml_engine='iqtree' to search under "
                "an existing constraint tree"
            )
        else:
            from .distance import constrained_nj
            tree = constrained_nj(names, D, constraint)

    # 3b. rooting -------------------------------------------------------
    if root == "midpoint":
        from ..treeops import midpoint_root
        tree = midpoint_root(tree)
    elif root != "none":
        if isinstance(root, str):
            raise ValueError(
                f"unknown root mode {root!r}; use 'none', 'midpoint', or a "
                "tip name/list of tip names to root on that outgroup"
            )
        from ..treeops import outgroup_root
        tree = outgroup_root(tree, root)

    # 4. bootstrap (works for distance, parsimony and native-ML methods) -
    if bootstrap:
        from .bootstrap import bootstrap_support, nj_builder, upgma_builder
        if method == "nj" and constraint is not None:
            # a replicate rebuilt *without* the constraint would score how
            # well the constrained split does against data that never had to
            # honour it in the first place -- read as "genuine" support, that
            # is not a fair test of the tree actually being reported. Every
            # replicate goes through constrained_nj too, so a low number here
            # means what it should: the data resample-to-resample disagree
            # with the grouping the constraint enforced, not that this one
            # resample happened to draw the split differently on its own.
            from .bootstrap import distance_matrix_model as _dmm
            from .distance import constrained_nj
            builder = lambda a: constrained_nj(*_dmm(a, dist_model), constraint)  # noqa: E731
        elif method == "nj" and nj_engine != "builtin":
            # rebuilt per replicate rather than routed through RapidNJ's own
            # -b: undocumented exactly how it reports support, and unlike an
            # external ML search a fast approximate NJ rebuild is cheap
            # enough that phytreon's own resampling loop is not a bottleneck.
            from .ml import infer_rapidnj
            builder = lambda a: infer_rapidnj(a)                       # noqa: E731
        elif method == "nj":
            builder = lambda a: nj_builder(a, dist_model)            # noqa: E731
        elif method == "upgma":
            builder = lambda a: upgma_builder(a, dist_model)         # noqa: E731
        elif method in ("parsimony", "mp") and parsimony_model == "camin_sokal":
            from .lineage import lineage_tree
            builder = lambda a: lineage_tree(a, search=ml_search)  # noqa: E731
        elif method in ("parsimony", "mp"):
            from .parsimony import parsimony_tree
            builder = lambda a: parsimony_tree(a, search=ml_search)  # noqa: E731
        elif method == "ml" and ml_engine == "native":
            from .ml_native import ml_tree as _nml
            # replicates skip NNI for tractability (branch+model opt only)
            builder = lambda a: _nml(a, model=ml_model, gamma=ml_gamma, search=False)  # noqa: E731
        else:
            builder = None
        if builder is not None:
            tree, _ = bootstrap_support(aln, builder=builder, n=bootstrap,
                                        seed=seed, reference=tree)

    return (tree, aln) if return_alignment else tree
