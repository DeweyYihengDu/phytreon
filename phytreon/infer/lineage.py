"""Single-cell CRISPR lineage-tracing tree reconstruction.

Cas9-based lineage tracing (GESTALT, Cassiopeia-style "allele table"
experiments, etc.) engineers cells to accumulate small indel "scars" at a
handful of genomic target sites as they divide; cells that share a scar
inherited it from a common ancestor, so the scar pattern is itself a
phylogenetic character matrix -- just one with an **irreversible**
evolutionary model, unlike DNA/protein sequences: a site, once cut, can
never revert to unedited, and can't spontaneously become a *different* edit
either (Cas9 cannot recut a site whose target sequence the first edit
already destroyed).

This is a parallel addition alongside the existing (reversible) Fitch
parsimony in :mod:`phytreon.infer.parsimony`, which is not touched by
anything here:

* :func:`read_allele_table` turns a raw Cassiopeia-style allele table into
  an :class:`~phytreon.infer.align.Alignment`, reusing the existing
  :func:`phytreon.infer.matrix.read_character_matrix` for the actual
  state-recoding rather than reimplementing it.
* :func:`sankoff_score` is a general Sankoff (cost-matrix) parsimony engine;
  :func:`camin_sokal_score` is the irreversible ("Camin-Sokal") preset
  appropriate for lineage-barcode data -- multiple *independent* origins of
  the same derived state are allowed (the identical indel outcome can recur
  by chance), reversion is not.
* :func:`lineage_tree` is the NNI hill-climbing search minimizing
  :func:`camin_sokal_score`, mirroring
  :func:`phytreon.infer.parsimony.parsimony_tree`'s structure and reusing
  the same shared NNI machinery (:mod:`phytreon.infer._search`).
"""
from __future__ import annotations

import re
from typing import Dict, Optional, Sequence, Union

from ..core.tree import Tree
from .align import Alignment
from .matrix import read_character_matrix

_MISSING = "?"                  # read_character_matrix's own missing-code
_ANCESTRAL_LABEL = ""           # sorts before every real edit label, always
_ANCESTRAL_ROW = "__ancestral__"
_EDIT_RE = re.compile(r"\[([^\]]*)\]")


# --------------------------------------------------------------------------
# raw allele table -> Alignment
# --------------------------------------------------------------------------
def _classify(raw: str, missing_markers: Sequence[str]) -> Optional[str]:
    """One r1/r2/r3 cell value -> None (missing/dropout), ``""`` (ancestral/
    uncut), or the specific edit's label."""
    if any(m in raw for m in missing_markers):
        return None
    if "[None]" in raw:
        return _ANCESTRAL_LABEL
    m = _EDIT_RE.search(raw)
    return m.group(1) if m else raw


def read_allele_table(source: Union[str, "object"], *, cell_col: str = "cellBC",
                      intbc_col: str = "intBC",
                      site_cols: Sequence[str] = ("r1", "r2", "r3"),
                      missing_markers: Sequence[str] = ("NC",)) -> Alignment:
    """Build an :class:`Alignment` from a Cassiopeia-style CRISPR allele
    table (tab-separated path, or an existing :class:`pandas.DataFrame`).

    One character per ``(intBC, site)`` combination. Each cell's raw edit
    string (e.g. ``"CGCCG[111:3D]AATGG"``) is reduced to its bracketed edit
    label (``"111:3D"``), the ancestral/uncut state (``"[None]"``
    substring) is normalized to a dedicated ancestral label, and any value
    containing one of ``missing_markers`` (Cassiopeia uses ``"NC"``) is
    treated as missing/dropout. Combinations with no row at all (allele
    dropout at the row level, not just a missing-marker value) are also
    encoded as missing, not silently omitted.

    The character actually assigned to the ancestral state is guaranteed to
    be code ``"0"`` at *every* site -- including a site that happens to be
    edited in 100% of profiled cells (this occurs in real data; without the
    fix, :func:`phytreon.infer.matrix.read_character_matrix` would recode
    whichever real edit happens to sort first as if it were ancestral).
    """
    import pandas as pd

    df = source.copy() if isinstance(source, pd.DataFrame) else pd.read_csv(source, sep="\t")

    key = [cell_col, intbc_col]
    if df.duplicated(key).any():
        for col in site_cols:
            if (df.groupby(key)[col].nunique() > 1).any():
                raise ValueError(
                    f"duplicate {tuple(key)} rows disagree on column {col!r}; "
                    "resolve duplicates before calling read_allele_table")
        df = df.groupby(key, as_index=False).first()

    cells = sorted(df[cell_col].astype(str).unique())
    if _ANCESTRAL_ROW in cells:
        raise ValueError(f"reserved cell name {_ANCESTRAL_ROW!r} collides with real data; "
                         "rename that cell before calling read_allele_table")
    intbcs = sorted(df[intbc_col].astype(str).unique())
    all_sites = [f"{intbc}:{site}" for intbc in intbcs for site in site_cols]

    long_rows = []
    for _, row in df.iterrows():
        cell = str(row[cell_col])
        intbc = str(row[intbc_col])
        for site in site_cols:
            long_rows.append((cell, f"{intbc}:{site}", _classify(str(row[site]), missing_markers)))
    long_df = pd.DataFrame(long_rows, columns=["cell", "site", "state"])
    wide = long_df.pivot(index="cell", columns="site", values="state")
    wide = wide.reindex(index=cells, columns=all_sites)   # dropout -> NaN, not omitted

    # phantom all-ancestral row: guarantees "" (-> code "0") is present in
    # every column even where no real cell happens to carry the ancestral
    # state (see docstring)
    wide.loc[_ANCESTRAL_ROW] = _ANCESTRAL_LABEL

    aln = read_character_matrix(wide)
    keep = [i for i, n in enumerate(aln.names) if n != _ANCESTRAL_ROW]
    return Alignment([aln.names[i] for i in keep], [aln.seqs[i] for i in keep])


# --------------------------------------------------------------------------
# encoding: compressed site patterns, ancestral state always slot 0
# --------------------------------------------------------------------------
def _encode_lineage(aln: Alignment):
    """Compressed site patterns as integer state slots -- slot 0 is always
    the ancestral state (character ``"0"``, per :func:`read_allele_table`'s
    encoding convention), slots ``1..k-1`` are whichever distinct derived
    states that particular site happens to have (assignment doesn't need to
    be consistent site-to-site: :func:`camin_sokal_score`'s cost matrix
    treats every non-ancestral state symmetrically). Missing (``"?"``)
    encoded as ``-1``."""
    import numpy as np
    names = list(aln.names)
    ncol = aln.ncol
    patterns: Dict[tuple, int] = {}
    order = []
    for j in range(ncol):
        col = tuple(s[j] for s in aln.seqs)
        if col not in patterns:
            patterns[col] = 0
            order.append(col)
        patterns[col] += 1
    npat = len(order)

    max_k = 1
    site_maps = []
    for col in order:
        distinct = sorted({c for c in col if c not in (_MISSING, "0")})
        site_maps.append({"0": 0, **{c: i + 1 for i, c in enumerate(distinct)}})
        max_k = max(max_k, len(distinct) + 1)

    states = np.full((len(names), npat), -1, dtype=np.int16)
    for p, (col, smap) in enumerate(zip(order, site_maps)):
        for i, ch in enumerate(col):
            if ch != _MISSING:
                states[i, p] = smap[ch]
    weights = np.array([patterns[c] for c in order], dtype=float)
    return names, states, weights, max_k


def _camin_sokal_cost(k: int):
    """(k, k) cost matrix: same-state free, ancestral(0)->derived costs 1,
    a derived state can never revert to ancestral or convert to a different
    derived state (both forbidden -- ``inf``)."""
    import numpy as np
    cost = np.full((k, k), np.inf)
    np.fill_diagonal(cost, 0.0)
    cost[0, 1:] = 1.0
    return cost


# --------------------------------------------------------------------------
# Sankoff parsimony (general engine) + the Camin-Sokal preset
# --------------------------------------------------------------------------
def sankoff_score(tree: Tree, aln: Alignment, cost_matrix, *,
                  root_state: Optional[int] = 0, data=None) -> float:
    """Sankoff parsimony cost of ``tree`` for ``aln`` under an arbitrary
    ``(k, k)`` ``cost_matrix`` (state slots as encoded by
    :func:`_encode_lineage` -- slot 0 is the ancestral state).

    ``root_state`` fixes the tree's root to a specific state slot rather
    than optimizing over every possible root state, which is the right
    choice whenever the true ancestral state is known a priori -- as it is
    here: the clonal progenitor cell predates all CRISPR editing, so the
    root is always ancestral (slot 0, the default). Pass ``root_state=None``
    to instead optimize freely over the root state (standard unrooted
    Sankoff parsimony).
    """
    import numpy as np
    if data is None:
        data = _encode_lineage(aln)
    names, states, weights, k = data
    if cost_matrix.shape != (k, k):
        raise ValueError(f"cost_matrix must be {(k, k)} to match the encoded "
                         f"data (got {cost_matrix.shape})")
    idx = {n: i for i, n in enumerate(names)}
    npat = states.shape[1]
    state_range = np.arange(k)
    cache: Dict[int, "np.ndarray"] = {}

    for node in tree.traverse("postorder"):
        if node.is_leaf:
            s = states[idx[node.name]]
            dp = np.where(s[:, None] == state_range[None, :], 0.0, np.inf)
            dp[s < 0] = 0.0                      # missing: free at every state
            cache[id(node)] = dp
        else:
            dp = np.zeros((npat, k))
            for c in node.children:
                child_dp = cache[id(c)]
                combined = cost_matrix[None, :, :] + child_dp[:, None, :]
                dp = dp + combined.min(axis=2)
            cache[id(node)] = dp

    root_dp = cache[id(tree.root)]
    col = root_dp.min(axis=1) if root_state is None else root_dp[:, root_state]
    return float((weights * col).sum())


def camin_sokal_score(tree: Tree, aln: Alignment, data=None) -> float:
    """Irreversible ("Camin-Sokal") parsimony cost: the minimum number of
    independent scar-acquisition events needed to explain the observed
    lineage-barcode states, given that state ``"0"`` (ancestral/uncut) can
    only transition to a derived state -- never revert, never interconvert
    directly with a different derived state -- and the tree's root is fixed
    ancestral at every site (see :func:`sankoff_score`'s ``root_state``).

    Assumes ``aln`` encodes its ancestral state as character ``"0"`` at
    every site, which :func:`read_allele_table` guarantees.
    """
    if data is None:
        data = _encode_lineage(aln)
    _, _, _, k = data
    return sankoff_score(tree, aln, _camin_sokal_cost(k), root_state=0, data=data)


def _min_possible_score(data) -> float:
    """Lower bound: every distinct non-ancestral state needs >=1 origin,
    regardless of topology (assumes best case -- each occurrence forms one
    monophyletic clade)."""
    _, states, weights, _ = data
    total = 0.0
    for p in range(states.shape[1]):
        col = states[:, p]
        total += len({s for s in col.tolist() if s > 0}) * weights[p]
    return total


# --------------------------------------------------------------------------
# tree search
# --------------------------------------------------------------------------
def lineage_tree(aln: Alignment, start: Optional[Tree] = None,
                 search: bool = True, max_sweeps: int = 20) -> Tree:
    """Minimum-Camin-Sokal-cost lineage tree by NNI hill-climbing (from an
    NJ start on raw Hamming distance).

    Result carries ``tree.data['camin_sokal_score']`` (the reconstruction's
    cost), ``['min_possible_score']`` (topology-independent lower bound) and
    ``['excess_origins']`` (their difference -- homoplasy the tree couldn't
    avoid). Camin-Sokal's bounds don't correspond to Fitch parsimony's
    ``ci``/``ri`` (those formulas assume a reversible model), so those
    fields are intentionally not set here.
    """
    from .bootstrap import nj_builder
    from ..treeops import midpoint_root
    from ._search import internal_edges, nni_neighbors

    data = _encode_lineage(aln)
    _, _, _, k = data
    cost = _camin_sokal_cost(k)
    tree = start or midpoint_root(nj_builder(aln, model="raw"))
    best = sankoff_score(tree, aln, cost, root_state=0, data=data)

    if search:
        for _ in range(max_sweeps):
            improved = False
            for node in internal_edges(tree):
                for swap in list(nni_neighbors(node)):
                    swap()
                    sc = sankoff_score(tree, aln, cost, root_state=0, data=data)
                    if sc < best - 1e-9:
                        best = sc
                        improved = True
                        break
                    swap()                              # involution -> undo
            if not improved:
                break

    min_possible = _min_possible_score(data)
    tree.data["camin_sokal_score"] = best
    tree.data["min_possible_score"] = min_possible
    tree.data["excess_origins"] = best - min_possible
    tree.data["mp_search"] = "NNI-local"
    return tree
