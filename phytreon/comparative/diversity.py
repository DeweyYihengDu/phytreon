"""Phylogenetic diversity: Faith's PD and UniFrac.

Both answer a question ancestral-state reconstruction does not: not "what
did the ancestors look like", but "how much evolutionary history does this
set of taxa -- the ones actually observed in one sample, or in two samples
being compared -- collectively represent". The standard alpha- and
beta-diversity metrics for a community sitting on a tree (16S ASVs, say),
and the natural next step once that tree is built.
"""
from __future__ import annotations

from typing import Dict, Iterable, Set, Union

from ..core.tree import Node, Tree


def _leaves_by_name(tree: Tree) -> Dict[str, Node]:
    return {lf.name: lf for lf in tree.leaves()}


def _check_taxa(taxa: Set[str], leaves: Dict[str, Node]) -> None:
    missing = taxa - set(leaves)
    if missing:
        raise ValueError(f"taxa not found in tree: {sorted(missing)}")


def _check_table_columns(table, leaves: Dict[str, Node], caller: str) -> None:
    """Every column of a samples x taxa table has to be a tip of the tree.

    Checked up front, over all columns, rather than per row as a side effect of
    whichever taxa happened to be present in that row: an ASV table normally
    has *more* ASVs than the tree does tips, because tree building drops
    sequences (too short, failed alignment, chimeras filtered after the table
    was counted). Left unchecked that is not a loud failure but a quiet one --
    the extra taxa's abundance still lands in each sample's total, so every real
    taxon's fraction comes out too small and the distances are wrong rather
    than absent.
    """
    unknown = [str(c) for c in table.columns if c not in leaves]
    if unknown:
        shown = sorted(unknown)[:10]
        raise ValueError(
            f"{caller}: {len(unknown)} of {len(table.columns)} table columns "
            f"are not tips of the tree: {shown}"
            f"{' ...' if len(unknown) > len(shown) else ''}. Subset the table "
            f"to the tree's tips first, e.g. "
            f"table[[c for c in table.columns if c in set(tree.leaf_names())]]"
        )


def _edges_to_root(leaves: Dict[str, Node], taxa: Iterable[str]) -> Set[Node]:
    """Nodes whose edge-to-parent lies on some tip-to-root path in ``taxa``.

    One node per edge (a node's ``length`` is the branch above it), so this
    set *is* "the minimal subtree connecting taxa to the root", represented
    the cheapest way phytreon can: as the edges themselves rather than a
    materialised sub-tree. Shared once a further-out node is already in the
    set, every path above it is too, so each leaf's walk stops the moment it
    rejoins one already taken -- the total work is bounded by the number of
    edges in the tree, not by (number of taxa) x (tree depth).
    """
    used: Set[Node] = set()
    for name in taxa:
        node = leaves[name]
        while node.parent is not None and node not in used:
            used.add(node)
            node = node.parent
    return used


def faiths_pd(tree: Tree, taxa: Union[str, Iterable[str]]) -> float:
    """Faith's phylogenetic diversity: total branch length of the smallest
    subtree connecting ``taxa`` to the root.

    The standard phylogenetic alpha-diversity metric -- how much
    evolutionary history a sample's observed taxa collectively represent,
    as opposed to a plain species count, which treats a genus-full of near-
    identical ASVs the same as one that spans the whole tree. Measured from
    the root (not the sample's own MRCA), the usual convention when
    comparing several samples embedded in one shared reference tree, since
    it makes every sample's number a fraction of the same total.
    """
    taxa = {taxa} if isinstance(taxa, str) else set(taxa)
    leaves = _leaves_by_name(tree)
    _check_taxa(taxa, leaves)
    if not taxa:
        return 0.0
    return sum(n.length or 0.0 for n in _edges_to_root(leaves, taxa))


def faiths_pd_table(tree: Tree, table) -> "pd.Series":  # noqa: F821
    """Faith's PD for every sample (row) of a samples x taxa table at once.

    ``table`` is a :class:`pandas.DataFrame`, presence/absence or abundance
    -- only whether each entry is nonzero matters, since Faith's PD does not
    weight by how much of a taxon is present (see :func:`weighted_unifrac`
    for the metric that does).
    """
    import pandas as pd
    leaves = _leaves_by_name(tree)
    _check_table_columns(table, leaves, "faiths_pd_table")
    out = {}
    for sample, row in table.iterrows():
        present = row.index[row.to_numpy() > 0]
        out[sample] = faiths_pd(tree, list(present))
    return pd.Series(out, name="faiths_pd")


def unweighted_unifrac(tree: Tree, taxa_a: Iterable[str],
                       taxa_b: Iterable[str]) -> float:
    """Unweighted UniFrac: the fraction of branch length spanning two
    samples' taxa that belongs to only one of them.

    0 when the two samples' observed taxa span exactly the same branches
    (however different the two sets of taxa look), 1 when they share no
    branches at all. Presence/absence only -- see :func:`weighted_unifrac`
    to also account for how abundant each taxon is.
    """
    a, b = set(taxa_a), set(taxa_b)
    leaves = _leaves_by_name(tree)
    _check_taxa(a | b, leaves)
    if not a and not b:
        return 0.0
    in_a = _edges_to_root(leaves, a)
    in_b = _edges_to_root(leaves, b)
    total = sum(n.length or 0.0 for n in (in_a | in_b))
    if total <= 0.0:
        return 0.0
    unshared = sum(n.length or 0.0 for n in (in_a ^ in_b))
    return unshared / total


def _subtree_fractions(tree: Tree, abundance: Dict[str, float]) -> Dict[Node, float]:
    """For every non-root node, the fraction of ``abundance``'s total
    carried by tips in its subtree -- one post-order pass, reused for every
    branch's contribution to weighted UniFrac.
    """
    total = sum(abundance.values())
    per_node: Dict[Node, float] = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf:
            per_node[node] = abundance.get(node.name, 0.0)
        else:
            per_node[node] = sum(per_node[c] for c in node.children)
    if total <= 0.0:
        return {n: 0.0 for n in per_node if n.parent is not None}
    return {n: v / total for n, v in per_node.items() if n.parent is not None}


def weighted_unifrac(tree: Tree, abundance_a: Dict[str, float],
                     abundance_b: Dict[str, float],
                     normalized: bool = True) -> float:
    """Weighted UniFrac: like :func:`unweighted_unifrac`, but a branch's
    contribution is scaled by how differently abundant the two samples are
    on either side of it, not just whether they differ at all.

    ``abundance_a``/``abundance_b`` map tip name -> a count or relative
    abundance (need not already sum to 1; each is normalised to its own
    total internally, so the two samples' totals do not have to match).
    ``normalized=True`` (default) divides by the maximum possible value for
    this pair of samples, bounding the result to ``[0, 1]`` and making it
    comparable across sample pairs and trees the way the raw metric is not;
    ``normalized=False`` returns the original, unbounded branch-length units.
    """
    leaves = _leaves_by_name(tree)
    _check_taxa(set(abundance_a) | set(abundance_b), leaves)
    frac_a = _subtree_fractions(tree, abundance_a)
    frac_b = _subtree_fractions(tree, abundance_b)
    nodes = set(frac_a) | set(frac_b)
    if not nodes:
        return 0.0
    num = sum((n.length or 0.0) * abs(frac_a.get(n, 0.0) - frac_b.get(n, 0.0))
             for n in nodes)
    if not normalized:
        return num
    denom = sum((n.length or 0.0) * (frac_a.get(n, 0.0) + frac_b.get(n, 0.0))
               for n in nodes)
    return num / denom if denom > 0.0 else 0.0


def unifrac_matrix(tree: Tree, table, weighted: bool = False,
                   normalized: bool = True) -> "pd.DataFrame":  # noqa: F821
    """Pairwise UniFrac distances for every sample (row) of a samples x taxa
    table -- the usual next step being an ordination (PCoA) or PERMANOVA on
    the result, both outside phytreon's own scope.

    Each sample's per-branch data is computed once, into a row of a
    samples x edges array, so every pair is one numpy expression rather than a
    Python loop over a dict -- 6-15x faster than the latter depending on shape
    (967 samples on a 500-tip tree: about 10 seconds, from about 5 minutes).
    That is a constant factor, not a change of complexity: the work is still
    proportional to ``samples^2 x edges``, since a pair's branch-by-branch
    comparison genuinely has to look at every branch.
    """
    import numpy as np
    import pandas as pd
    leaves = _leaves_by_name(tree)
    _check_table_columns(table, leaves, "unifrac_matrix")
    samples = list(table.index)
    n = len(samples)
    mat = np.zeros((n, n))

    # one fixed edge ordering for the whole call, so each sample's per-branch
    # data becomes a row of an array and every pair is a numpy expression rather
    # than a Python loop over a dict
    edges = [nd for nd in tree.traverse("postorder") if nd.parent is not None]
    index = {nd: i for i, nd in enumerate(edges)}
    lengths = np.array([nd.length or 0.0 for nd in edges])

    if weighted:
        F = np.zeros((n, len(edges)))
        for si, s in enumerate(samples):
            for nd, frac in _subtree_fractions(tree, table.loc[s].to_dict()).items():
                F[si, index[nd]] = frac
        # the denominator separates -- sum(L * (F_i + F_j)) is just r_i + r_j --
        # so it costs one matrix-vector product for all pairs, not a pass per pair
        r = F @ lengths
        for i in range(n):
            num = np.abs(F[i] - F[i + 1:]) @ lengths
            if normalized:
                denom = r[i] + r[i + 1:]
                row = np.divide(num, denom, out=np.zeros_like(num), where=denom > 0.0)
            else:
                row = num
            mat[i, i + 1:] = mat[i + 1:, i] = row
    else:
        A = np.zeros((n, len(edges)), dtype=bool)
        for si, s in enumerate(samples):
            present = list(table.columns[table.loc[s].to_numpy() > 0])
            for nd in _edges_to_root(leaves, present):
                A[si, index[nd]] = True
        for i in range(n):
            total = (A[i] | A[i + 1:]) @ lengths
            unshared = (A[i] ^ A[i + 1:]) @ lengths
            row = np.divide(unshared, total, out=np.zeros_like(total),
                            where=total > 0.0)
            mat[i, i + 1:] = mat[i + 1:, i] = row
    return pd.DataFrame(mat, index=samples, columns=samples)
