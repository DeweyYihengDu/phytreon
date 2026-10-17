"""Phylogenetic inference via external engines: ML, and fast approximate NJ.

Pure-Python ML is impractical, and pure-Python NJ is the textbook O(n^3)
algorithm, impractical past a few thousand tips (a large 16S ASV table gets
there) -- so phytreon shells out to the standard tools for both. Each
wrapper is graceful: if the program is not on PATH it raises a clear,
actionable error rather than failing obscurely.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple

from ..core.tree import Tree
from ..core.io import parse_newick, to_newick
from .align import Alignment


def _require(tool: str, path: Optional[str]) -> str:
    exe = path or shutil.which(tool)
    if exe is None:
        raise RuntimeError(
            f"{tool!r} not found on PATH. Install it (e.g. conda install -c "
            f"bioconda {tool}) or pass path=. ML inference needs an external "
            f"engine; distance methods (nj/upgma) are built in."
        )
    return exe


def _run(cmd: List[str]) -> subprocess.CompletedProcess:
    """``subprocess.run`` that fails with the tool's own message attached.

    ``check=True`` alone raises ``CalledProcessError``, whose default message
    is just "Command [...] returned non-zero exit status 1" -- the *reason*
    (a model name the engine did not recognise, a constraint with too few
    taxa, an alignment it could not parse) is sitting right there in the
    stderr this already captures, and was being silently dropped on the floor
    every time a run actually failed. Contradicted the module's own claim to
    be graceful about exactly this.
    """
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(
            f"{cmd[0]} exited with status {proc.returncode}"
            + (f":\n{detail}" if detail else " (no output captured)")
        )
    return proc


# Characters every engine accepts in a sequence name. IQ-TREE rewrites
# anything else to '_' (so the tree it hands back no longer matches the input
# names or the metadata keyed on them), and Newick punctuation such as ':' or
# ',' in a name breaks parsing the tree that comes back.
_SAFE_NAME = re.compile(r"[A-Za-z0-9_.\-]+")


def _safe_names(aln: Alignment) -> Tuple[Alignment, Optional[Dict[str, str]]]:
    """Return ``aln`` unchanged if every name is engine-safe; otherwise a copy
    with placeholder names plus the ``{placeholder: original}`` map that
    :func:`_restore_names` uses to put the real names back on the tree."""
    if all(_SAFE_NAME.fullmatch(name) for name in aln.names):
        return aln, None
    placeholders = [f"phytreon_t{i}" for i in range(len(aln.names))]
    return (Alignment(placeholders, list(aln.seqs)),
            dict(zip(placeholders, aln.names)))


def _restore_names(tree: Tree, mapping: Optional[Dict[str, str]]) -> Tree:
    if mapping:
        for leaf in tree.leaves():
            leaf.name = mapping.get(leaf.name, leaf.name)
    return tree


def _constraint_file(constraint, tmp: str,
                     mapping: Optional[Dict[str, str]]) -> str:
    """Path of a constraint tree for the engine, renamed to match the
    alignment's placeholders when :func:`_safe_names` had to rename it."""
    if mapping is None and isinstance(constraint, str):
        return constraint
    gfile = os.path.join(tmp, "constraint.tre")
    if mapping is None:
        constraint.write(gfile)
        return gfile
    if isinstance(constraint, str):
        with open(constraint) as f:
            tree = parse_newick(f.read())
    else:
        tree = parse_newick(to_newick(constraint))  # a copy: leave the caller's
    to_placeholder = {orig: ph for ph, orig in mapping.items()}
    for leaf in tree.leaves():
        leaf.name = to_placeholder.get(leaf.name, leaf.name)
    with open(gfile, "w") as f:
        f.write(to_newick(tree))
    return gfile


def infer_iqtree(aln: Alignment, model: str = "MFP", bootstrap: int = 0,
                 seed: Optional[int] = None, path: Optional[str] = None,
                 extra_args: Optional[List[str]] = None,
                 constraint=None) -> Tree:
    """ML tree with IQ-TREE. ``model='MFP'`` runs ModelFinder Plus.

    ``constraint`` asks IQ-TREE to only search trees that respect a set of
    clades (``-g``) -- typically the polytomies :func:`~phytreon.infer.
    constraint.constraint_tree` builds from a taxonomy column, one per group,
    each left internally unresolved so IQ-TREE still picks the topology
    *within* every group and *between* groups by likelihood; only "no other
    tip may land inside this group" is fixed. Pass the :class:`~phytreon.
    core.tree.Tree` it returns directly, or a path to an existing constraint
    Newick file. ``seed`` makes a stochastic search (and UFBoot) reproducible
    (``-seed``); left unset, IQ-TREE picks its own, and two calls differ.
    """
    exe = _require("iqtree2", path) if (path or shutil.which("iqtree2")) else _require("iqtree", path)
    aln, names = _safe_names(aln)
    with tempfile.TemporaryDirectory() as tmp:
        infile = os.path.join(tmp, "aln.fasta")
        aln.to_fasta(infile)
        cmd = [exe, "-s", infile, "-m", model, "-redo", "-quiet"]
        if bootstrap:
            cmd += ["-bb", str(max(bootstrap, 1000))]
        if seed is not None:
            cmd += ["-seed", str(seed)]
        if constraint is not None:
            cmd += ["-g", _constraint_file(constraint, tmp, names)]
        cmd += extra_args or []
        _run(cmd)
        with open(infile + ".treefile") as f:
            return _restore_names(parse_newick(f.read()), names)


def infer_fasttree(aln: Alignment, nucleotide: Optional[bool] = None,
                   path: Optional[str] = None) -> Tree:
    """Approximate-ML tree with FastTree."""
    exe = _require("FastTree", path) if (path or shutil.which("FastTree")) else _require("fasttree", path)
    if nucleotide is None:
        from .align import guess_type
        nucleotide = guess_type(aln.seqs) == "nucleotide"
    aln, names = _safe_names(aln)
    with tempfile.TemporaryDirectory() as tmp:
        infile = os.path.join(tmp, "aln.fasta")
        aln.to_fasta(infile)
        cmd = [exe] + (["-nt"] if nucleotide else []) + [infile]
        out = _run(cmd).stdout
        return _restore_names(parse_newick(out), names)


def infer_raxmlng(aln: Alignment, model: str = "GTR+G", bootstrap: int = 0,
                  seed: Optional[int] = None, threads: int = 1,
                  path: Optional[str] = None,
                  extra_args: Optional[List[str]] = None,
                  constraint=None) -> Tree:
    """ML tree with RAxML-NG, IQ-TREE's other standard alternative.

    Unlike IQ-TREE's ``-m MFP``, RAxML-NG has no built-in "find the model for
    me" default, so ``model`` has to name one explicitly (``"GTR+G"`` for
    nucleotides, ``"LG+G"``/``"WAG+G"`` etc. for protein). ``constraint``
    works the same way as :func:`infer_iqtree`'s -- a :class:`~phytreon.core.
    tree.Tree` of polytomies (see :func:`~phytreon.infer.constraint.
    constraint_tree`) or a path to one, passed through to RAxML-NG's
    ``--tree-constraint`` -- and needs at least 4 taxa in it, a RAxML-NG
    restriction phytreon does not relax. ``bootstrap`` runs search and
    bootstrapping together (``--all``) rather than phytreon's rebuild-from-
    scratch loop, for the same reason :func:`infer_iqtree` reaches for its
    own ``-bb`` instead: one external search per replicate does not scale.
    ``seed`` makes the search reproducible (``--seed``); left unset it still
    passes a fixed one (``1``, not RAxML-NG's own random default) so that two
    calls with the same arguments agree, since a tree that changes between
    identical-looking runs is its own kind of surprise. ``threads`` is passed
    through explicitly (``--threads``) since older RAxML-NG releases error
    out without it; newer ones default to auto-detection, but there's no
    version probe here to tell which is running.
    """
    exe = _require("raxml-ng", path)
    aln, names = _safe_names(aln)
    with tempfile.TemporaryDirectory() as tmp:
        infile = os.path.join(tmp, "aln.fasta")
        aln.to_fasta(infile)
        prefix = os.path.join(tmp, "run")
        cmd = [exe, "--msa", infile, "--model", model, "--prefix", prefix,
              "--threads", str(threads), "--seed", str(1 if seed is None else seed)]
        cmd.append("--all" if bootstrap else "--search")
        if bootstrap:
            cmd += ["--bs-trees", str(bootstrap)]
        if constraint is not None:
            cmd += ["--tree-constraint", _constraint_file(constraint, tmp, names)]
        cmd += extra_args or []
        _run(cmd)
        # --all additionally writes <prefix>.raxml.support (the best tree
        # with bootstrap values mapped onto it); --search alone never
        # produces that file, so bestTree is the one that always exists.
        result = prefix + (".raxml.support" if bootstrap else ".raxml.bestTree")
        with open(result) as f:
            return _restore_names(parse_newick(f.read()), names)


def infer_rapidnj(aln: Alignment, evolution_model: str = "kim",
                  bootstrap: int = 0, nonneg: bool = True,
                  path: Optional[str] = None,
                  extra_args: Optional[List[str]] = None) -> Tree:
    """Fast approximate NJ with RapidNJ, for alignments too large for the
    built-in neighbour-joining (:func:`~phytreon.infer.distance.
    neighbor_joining` wraps Biopython's textbook O(n^3) implementation,
    impractical somewhere in the low thousands of tips -- exactly the range
    a large 16S ASV table reaches).

    ``evolution_model`` is RapidNJ's own, more limited menu (``"jc"`` /
    ``"kim"``, its own default) -- it computes its distance matrix internally
    from the alignment, rather than accepting phytreon's (:func:`~phytreon.
    infer.bootstrap.distance_matrix_model`, ``"k2p"``/``"poisson"``/...) the
    way :func:`~phytreon.infer.distance.constrained_nj` does; pass a
    precomputed phylip distance matrix via ``extra_args=["-i", "pd"]``
    instead (with ``aln`` then unused) if that model menu matters more than
    RapidNJ's speed advantage over building the matrix in Python first.
    ``bootstrap`` uses RapidNJ's own resampling (``-b``), in the same run
    rather than through phytreon's rebuild-from-scratch loop, for the same
    reason :func:`infer_iqtree` reaches for its own ``-bb`` instead.
    """
    exe = _require("rapidnj", path)
    aln, names = _safe_names(aln)
    with tempfile.TemporaryDirectory() as tmp:
        infile = os.path.join(tmp, "aln.fasta")
        # rapidnj's FASTA reader wants one line per sequence; to_fasta()
        # wraps at 60 columns by default, which it cannot parse.
        aln.to_fasta(infile, width=10 ** 9)
        cmd = [exe, infile, "-i", "fa", "-a", evolution_model]
        if nonneg:
            cmd.append("-n")
        if bootstrap:
            cmd += ["-b", str(bootstrap)]
        cmd += extra_args or []
        out = _run(cmd).stdout
        return _restore_names(parse_newick(out), names)


def infer_ml(aln: Alignment, tool: str = "iqtree", **kw) -> Tree:
    if tool in ("iqtree", "iqtree2"):
        return infer_iqtree(aln, **kw)
    if tool in ("fasttree", "FastTree"):
        return infer_fasttree(aln, **kw)
    if tool in ("raxml-ng", "raxmlng", "raxml"):
        return infer_raxmlng(aln, **kw)
    raise ValueError(f"unknown ML tool {tool!r}; use 'iqtree', 'fasttree', "
                     "or 'raxml-ng'")
