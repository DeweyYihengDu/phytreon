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
               root: str = "none",
               bootstrap: int = 0,
               ml_engine: str = "native",
               ml_model: str = "HKY85",
               ml_gamma: int = 0,
               ml_search: bool = True,
               ml_tool: str = "iqtree",
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
    bootstrap   number of bootstrap replicates (0 = none; distance methods).
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
            from .ml_native import ml_tree as _native_ml
            tree = _native_ml(aln, model=ml_model, gamma=ml_gamma, search=ml_search)
        else:                                   # external engine (iqtree/fasttree)
            from .ml import infer_ml
            tree = infer_ml(aln, tool=ml_engine)
    elif method in ("parsimony", "mp"):
        from .parsimony import parsimony_tree
        tree = parsimony_tree(aln, search=ml_search)
    else:                                    # method in ("nj", "upgma")
        from .bootstrap import distance_matrix_model
        from .distance import neighbor_joining, upgma
        names, D = distance_matrix_model(aln, dist_model)
        tree = neighbor_joining(names, D) if method == "nj" else upgma(names, D)

    # 3b. rooting -------------------------------------------------------
    if root == "midpoint":
        from ..treeops import midpoint_root
        tree = midpoint_root(tree)
    elif root != "none":
        raise ValueError(f"unknown root mode {root!r}; use 'none' or 'midpoint'")

    # 4. bootstrap (works for distance, parsimony and native-ML methods) -
    if bootstrap:
        from .bootstrap import bootstrap_support, nj_builder, upgma_builder
        if method == "nj":
            builder = lambda a: nj_builder(a, dist_model)            # noqa: E731
        elif method == "upgma":
            builder = lambda a: upgma_builder(a, dist_model)         # noqa: E731
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
