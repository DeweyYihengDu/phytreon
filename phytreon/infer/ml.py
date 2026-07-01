"""Maximum-likelihood / parsimony inference via external engines.

Pure-Python ML is impractical, so phytreon shells out to the standard tools
when the user wants ML.  Each wrapper is graceful: if the program is not on
PATH it raises a clear, actionable error rather than failing obscurely.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from typing import List, Optional

from ..core.tree import Tree
from ..core.io import parse_newick
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


def infer_iqtree(aln: Alignment, model: str = "MFP", bootstrap: int = 0,
                 path: Optional[str] = None, extra_args: Optional[List[str]] = None
                 ) -> Tree:
    """ML tree with IQ-TREE. ``model='MFP'`` runs ModelFinder Plus."""
    exe = _require("iqtree2", path) if (path or shutil.which("iqtree2")) else _require("iqtree", path)
    with tempfile.TemporaryDirectory() as tmp:
        infile = os.path.join(tmp, "aln.fasta")
        aln.to_fasta(infile)
        cmd = [exe, "-s", infile, "-m", model, "-redo", "-quiet"]
        if bootstrap:
            cmd += ["-bb", str(max(bootstrap, 1000))]
        cmd += extra_args or []
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        with open(infile + ".treefile") as f:
            return parse_newick(f.read())


def infer_fasttree(aln: Alignment, nucleotide: Optional[bool] = None,
                   path: Optional[str] = None) -> Tree:
    """Approximate-ML tree with FastTree."""
    exe = _require("FastTree", path) if (path or shutil.which("FastTree")) else _require("fasttree", path)
    if nucleotide is None:
        from .align import guess_type
        nucleotide = guess_type(aln.seqs) == "nucleotide"
    with tempfile.TemporaryDirectory() as tmp:
        infile = os.path.join(tmp, "aln.fasta")
        aln.to_fasta(infile)
        cmd = [exe] + (["-nt"] if nucleotide else []) + [infile]
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
        return parse_newick(out)


def infer_ml(aln: Alignment, tool: str = "iqtree", **kw) -> Tree:
    if tool in ("iqtree", "iqtree2"):
        return infer_iqtree(aln, **kw)
    if tool in ("fasttree", "FastTree"):
        return infer_fasttree(aln, **kw)
    raise ValueError(f"unknown ML tool {tool!r}; use 'iqtree' or 'fasttree'")
