"""Per-domain trees for multidomain proteins, and what their disagreement means.

One tree built from a whole multidomain protein compresses several different
histories into one topology. Domains are recombined, fused and lost
independently: a protein whose N-terminal domain is related to one family and
whose C-terminal domain belongs to an unrelated superfamily has two ancestries,
and averaging them produces a tree that is nobody's history.

The cyanobacterial orange carotenoid protein is the textbook case -- an
all-helical N-terminal domain whose relatives are the HCPs, fused to an
NTF2-like beta-barrel C-terminal domain whose relatives are the CTDHs. A single
OCP tree is a category error; two domain trees, and the places they disagree,
are the fusion event.

So the disagreement is the result, not a nuisance. :func:`compare_domain_trees`
reports it two ways, because they mean different things:

* a **global** distance between two domain trees (Robinson-Foulds), which says
  how much they differ overall;
* a **leave-one-out** test (:func:`rogue_taxon`) asking whether dropping a single
  lineage would explain the disagreement -- which is what separates a real
  recombination signal in one lineage from a domain simply too short to resolve
  anything, where the conflict is spread over the whole tree and no single
  removal helps.
"""
from __future__ import annotations

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from ..core.tree import Tree
from .matrix import Alignment


def residue_to_column(aln: Alignment, reference: str,
                      ranges: Mapping[str, Tuple[int, int]]
                     ) -> Dict[str, Tuple[int, int]]:
    """Translate per-residue domain boundaries into alignment column ranges.

    The bridge that is easy to get wrong by one. HMMER and Pfam report a domain
    as residue positions *on one sequence*, counting only real residues; an
    alignment is indexed by columns, which include gaps. This walks the
    reference sequence and maps one to the other.

    ``ranges`` maps a domain name to ``(start, end)`` residue positions on
    ``reference``, **1-based and inclusive**, which is how HMMER prints them.
    Returns column ranges that are **0-based and half-open**, ready for
    :meth:`Alignment.select_columns` and Python slicing -- the two conventions
    are deliberately different so that mixing them up fails loudly rather than
    shifting a boundary by one.
    """
    try:
        row = aln.seqs[aln.names.index(reference)]
    except ValueError:
        raise ValueError(
            f"reference {reference!r} is not in the alignment; have "
            f"{aln.names[:8]}{' ...' if len(aln.names) > 8 else ''}"
        ) from None
    # column index of each successive residue of the reference
    residue_col: List[int] = [i for i, ch in enumerate(row) if ch not in "-."]
    n_res = len(residue_col)

    out: Dict[str, Tuple[int, int]] = {}
    for name, (start, end) in ranges.items():
        if start < 1 or end < start:
            raise ValueError(
                f"domain {name!r}: ({start}, {end}) is not a 1-based inclusive "
                f"range with start <= end"
            )
        if end > n_res:
            raise ValueError(
                f"domain {name!r} ends at residue {end} but {reference!r} has "
                f"only {n_res} residues in this alignment"
            )
        out[name] = (residue_col[start - 1], residue_col[end - 1] + 1)
    return out


def split_domains(aln: Alignment, domains: Mapping[str, Tuple[int, int]],
                  reference: Optional[str] = None,
                  min_columns: int = 1) -> Dict[str, Alignment]:
    """Cut an alignment into one sub-alignment per domain.

    ``domains`` maps a domain name to a range. With ``reference`` given, the
    ranges are residue positions on that sequence (1-based inclusive, as HMMER
    reports) and are converted by :func:`residue_to_column`; without it they are
    taken as alignment columns directly (0-based, half-open).

    Sequences are not filtered: a taxon whose domain is missing comes through as
    an all-gap row, which is visible and can then be dropped deliberately rather
    than silently. :func:`domain_trees` does drop them, and says how many.
    """
    ranges = (residue_to_column(aln, reference, domains) if reference is not None
              else {k: (int(a), int(b)) for k, (a, b) in domains.items()})
    if not ranges:
        raise ValueError("no domains given")
    out: Dict[str, Alignment] = {}
    for name, (start, end) in ranges.items():
        if not (0 <= start < end <= aln.ncol):
            raise ValueError(
                f"domain {name!r} maps to columns [{start}, {end}) which is "
                f"outside the alignment's {aln.ncol} columns"
            )
        if end - start < min_columns:
            raise ValueError(
                f"domain {name!r} is {end - start} columns, below "
                f"min_columns={min_columns}"
            )
        out[name] = aln.select_columns(list(range(start, end)))
    return out


def _drop_empty_rows(aln: Alignment, min_residues: int) -> Tuple[Alignment, List[str]]:
    keep, dropped = [], []
    for i, (name, seq) in enumerate(aln.records()):
        if sum(1 for ch in seq if ch not in "-.") >= min_residues:
            keep.append(i)
        else:
            dropped.append(name)
    kept = Alignment([aln.names[i] for i in keep], [aln.seqs[i] for i in keep])
    return kept, dropped


def domain_trees(aln: Alignment, domains: Mapping[str, Tuple[int, int]],
                 reference: Optional[str] = None, min_residues: int = 1,
                 **build_kw) -> Dict[str, object]:
    """Build one tree per domain from a single multidomain alignment.

    ``build_kw`` is passed through to :func:`~phytreon.build_tree`, so the domain
    trees are built exactly as any other tree would be -- including
    ``method="ml"`` with an external ``ml_engine``, which is what you want here,
    since the models that matter for divergent protein families are the
    site-heterogeneous ones only the external engines implement.

    Taxa carrying fewer than ``min_residues`` residues in a domain are dropped
    from that domain's tree rather than entered as an all-gap row, and are
    listed in the result: a domain present in only some of the taxa is itself a
    finding, and one that quietly became a zero-length branch would not be.

    Returns ``{"trees": {domain: Tree}, "alignments": {domain: Alignment},
    "dropped": {domain: [names]}, "columns": {domain: (start, end)}}``.
    """
    ranges = (residue_to_column(aln, reference, domains) if reference is not None
              else {k: (int(a), int(b)) for k, (a, b) in domains.items()})
    parts = split_domains(aln, ranges)
    trees, alns, dropped = {}, {}, {}
    for name, sub in parts.items():
        kept, gone = _drop_empty_rows(sub, min_residues)
        if kept.nseq < 3:
            raise ValueError(
                f"domain {name!r} has only {kept.nseq} taxa with at least "
                f"{min_residues} residue(s); too few to build a tree"
            )
        alns[name] = kept
        dropped[name] = gone
        trees[name] = _build(kept, **build_kw)
    return {"trees": trees, "alignments": alns, "dropped": dropped,
            "columns": ranges}


def _build(aln: Alignment, **kw):
    from .pipeline import build_tree
    return build_tree(aln, **kw)


def _topological_distances(tree: Tree, taxa: Sequence[str]):
    """Tip-to-tip distances counting *edges*, not branch lengths.

    Branch lengths are the wrong currency for comparing topologies: a gene
    evolving twice as fast has every branch twice as long and the same shape, and
    that must not register as a taxon having moved. Counting edges makes the
    comparison depend on the tree's shape alone.
    """
    from ..core.tree import Tree as _Tree
    # round-tripped through newick so the caller's tree keeps its branch lengths
    unit = _Tree.from_newick(tree.write())
    for node in unit.traverse("postorder"):
        if node.parent is not None:
            node.length = 1.0
    from ..comparative.community import patristic_distances
    return patristic_distances(unit, list(taxa))


def _neighbour_sets(tree: Tree, taxa: Sequence[str], k: int,
                    topological: bool) -> Dict[str, set]:
    """Each taxon's nearest tips, ties included.

    "The k nearest" is not well defined when several tips are equidistant, which
    on a topological distance is the normal case rather than an edge case -- every
    member of a balanced clade is the same number of edges away. Taking everything
    within the k-th smallest distance keeps the set deterministic; slicing an
    argsort would break ties by array order and manufacture displacement out of
    nothing.
    """
    import numpy as np
    if topological:
        names, D = _topological_distances(tree, taxa)
    else:
        from ..comparative.community import patristic_distances
        names, D = patristic_distances(tree, list(taxa))
    out = {}
    for i, name in enumerate(names):
        others = np.array([j for j in range(len(names)) if j != i])
        dists = D[i, others]
        cutoff = np.sort(dists)[min(k, len(dists)) - 1]
        out[name] = {names[j] for j, d in zip(others, dists) if d <= cutoff}
    return out


def taxon_displacement(tree_a: Tree, tree_b: Tree, k: int = 3,
                       topological: bool = True) -> "pd.Series":  # noqa: F821
    """How differently each taxon is placed in two trees, taxon by taxon.

    For every taxon shared by both trees, compares the set of its nearest tips
    and returns ``1 - Jaccard`` of the two sets: 0 means the same close relatives
    in both trees, 1 means none in common.

    Neighbour sets rather than a split-by-split comparison because they answer
    the question actually being asked -- "did this taxon change who it sits
    with". One displaced taxon perturbs many splits at once, so a split-based
    measure spreads the blame across the whole tree instead of naming the taxon
    responsible.

    ``topological=True`` (the default) counts edges between tips rather than
    summing branch lengths. That matters: genes evolve at different rates, so two
    trees of identical shape routinely have quite different branch lengths, and
    a distance-based comparison reports that as displacement. Measured on a gene
    tree that matched its species tree exactly -- Robinson-Foulds distance of 0 --
    the distance-based version still assigned displacements up to 0.5. Pass
    ``topological=False`` only if branch lengths are genuinely part of the
    question.
    """
    import pandas as pd
    shared = sorted(set(tree_a.leaf_names()) & set(tree_b.leaf_names()))
    if len(shared) < k + 1:
        raise ValueError(
            f"need more than k={k} shared taxa to compare neighbourhoods, "
            f"found {len(shared)}"
        )
    na = _neighbour_sets(tree_a, shared, k, topological)
    nb = _neighbour_sets(tree_b, shared, k, topological)
    out = {}
    for name in shared:
        a, b = na[name], nb[name]
        union = a | b
        out[name] = 1.0 - (len(a & b) / len(union) if union else 1.0)
    return pd.Series(out, name="displacement").sort_values(ascending=False)


def rogue_taxon(tree_a: Tree, tree_b: Tree) -> Dict[str, object]:
    """Ask whether one taxon accounts for two trees' disagreement.

    For every shared taxon, drops it from both trees and recomputes the
    normalised Robinson-Foulds distance. If one taxon's removal collapses the
    conflict, the two trees agree about everything except where that lineage
    goes -- which is what a horizontal transfer, a hidden paralogue, or a
    contaminated sequence looks like. If removing any single taxon barely helps,
    the disagreement is spread over the whole tree and the gene has no signal to
    offer rather than a story to tell.

    Returns ``rf`` (the full distance), ``rf_without`` (per taxon, the distance
    after dropping it), ``worst_taxon`` (whose removal helps most), and
    ``explained_by_one`` -- the fraction of the conflict that one taxon accounts
    for, 1.0 when removing it leaves the trees identical and 0.0 when it changes
    nothing.

    This replaced a ratio of maximum to mean neighbourhood displacement, which
    sounded reasonable and did not work: moving one taxon between clades
    perturbs the neighbourhoods of *both* clades, so on a modest number of taxa a
    single genuine transfer is not numerically "concentrated" at all. On an
    eight-taxon test the ratio ranked a saturated, signal-free gene above a
    simulated transfer. Leave-one-out asks the question directly instead of
    hoping a summary statistic stands in for it.
    """
    from ..treeops import prune_to_taxa, robinson_foulds
    import pandas as pd
    shared = sorted(set(tree_a.leaf_names()) & set(tree_b.leaf_names()))
    if len(shared) < 5:
        raise ValueError(
            f"rogue_taxon needs at least 5 shared taxa (dropping one must still "
            f"leave a tree with a resolvable split), found {len(shared)}"
        )
    a = prune_to_taxa(tree_a, shared) if len(shared) < len(tree_a.leaf_names()) else tree_a
    b = prune_to_taxa(tree_b, shared) if len(shared) < len(tree_b.leaf_names()) else tree_b
    full = float(robinson_foulds(a, b, normalized=True))
    without = {}
    for name in shared:
        keep = [t for t in shared if t != name]
        without[name] = float(robinson_foulds(
            prune_to_taxa(a, keep), prune_to_taxa(b, keep), normalized=True))
    series = pd.Series(without, name="rf_without").sort_values()
    best = float(series.iloc[0])
    return {
        "rf": full,
        "rf_without": series,
        "worst_taxon": str(series.index[0]) if full > 0 else "",
        "explained_by_one": (full - best) / full if full > 0 else float("nan"),
    }


def compare_domain_trees(trees: Mapping[str, Tree], k: int = 3
                        ) -> Dict[str, object]:
    """Compare per-domain trees pairwise, globally and taxon by taxon.

    Returns

    ``rf``
        a square DataFrame of normalised Robinson-Foulds distances between every
        pair of domain trees -- how much two domains disagree overall.
    ``displacement``
        a DataFrame of per-taxon displacement (see :func:`taxon_displacement`),
        one column per domain pair, sorted by the largest.
    ``explained_by_one``
        per pair, the fraction of the conflict that removing a single taxon
        accounts for (:func:`rogue_taxon`). **This is what separates a signal
        from a shrug.** A domain recombined or horizontally acquired in one
        lineage disagrees about that lineage and nowhere else, so dropping it
        collapses the conflict and this approaches 1; a domain too short to
        resolve anything disagrees everywhere and no single removal helps. Both
        can show the same overall RF distance and they mean opposite things.
    ``worst_taxon``
        per pair, the taxon whose removal helps most -- the recombination
        candidate when ``explained_by_one`` is high.
    """
    import itertools
    import pandas as pd
    names = list(trees)
    if len(names) < 2:
        raise ValueError(f"need at least 2 domain trees to compare, got {names}")
    rf = pd.DataFrame(0.0, index=names, columns=names)
    disp, conc, worst = {}, {}, {}
    for a, b in itertools.combinations(names, 2):
        shared = sorted(set(trees[a].leaf_names()) & set(trees[b].leaf_names()))
        from ..treeops import prune_to_taxa
        ta = prune_to_taxa(trees[a], shared) if len(shared) < len(trees[a].leaf_names()) else trees[a]
        tb = prune_to_taxa(trees[b], shared) if len(shared) < len(trees[b].leaf_names()) else trees[b]
        from ..treeops import robinson_foulds
        rf.loc[a, b] = rf.loc[b, a] = float(robinson_foulds(ta, tb, normalized=True))
        disp[f"{a}|{b}"] = taxon_displacement(ta, tb, k=k)
        rogue = rogue_taxon(ta, tb)
        conc[f"{a}|{b}"] = rogue["explained_by_one"]
        worst[f"{a}|{b}"] = rogue["worst_taxon"]
    return {
        "rf": rf,
        "displacement": pd.DataFrame(disp),
        "explained_by_one": pd.Series(conc, name="explained_by_one"),
        "worst_taxon": pd.Series(worst, name="worst_taxon"),
    }
