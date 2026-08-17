"""Species trees from many marker genes, and what each gene's disagreement means.

The standard route to a species tree from genomes or metagenome-assembled
genomes: take single-copy marker proteins, concatenate them into a supermatrix,
and infer one tree from all of it. Then the interesting part -- comparing each
gene's own tree with that species tree, because **the conflicts are the biology**.
A gene whose tree puts one lineage somewhere the rest of the genome does not is a
horizontal-transfer candidate.

What matters, and what :func:`gene_tree_conflict` is built around, is that two
completely different situations produce the same overall conflict score:

* **concentrated** conflict -- one lineage badly misplaced, everything else in
  agreement. That is what horizontal transfer, or a hidden paralogue, looks like.
* **diffuse** conflict -- every taxon slightly off. That is a gene too short, too
  saturated, or too badly aligned to resolve anything, and it carries no signal
  at all.

Ranking genes by total conflict mixes the two together and puts the useless genes
at the top alongside the interesting ones. The separating question is asked
directly instead: **would dropping one lineage explain the disagreement?** If it
would, that lineage is the transfer candidate; if no single removal helps, the
gene has no signal to offer. Both that and the total are reported.
"""
from __future__ import annotations

from typing import Dict, Mapping, Optional, Sequence, Tuple

from ..core.tree import Tree
from .matrix import Alignment


def concatenate(alignments: Mapping[str, Alignment],
                taxa: Optional[Sequence[str]] = None,
                gap: str = "-") -> Dict[str, object]:
    """Concatenate per-marker alignments into one supermatrix.

    ``alignments`` maps marker name to an :class:`Alignment` whose sequence names
    are taxa. Taxa absent from a marker are filled with gaps for that marker's
    columns, because requiring every marker in every taxon is not an option on
    real data -- metagenome-assembled genomes are incomplete by construction.

    ``taxa`` fixes the row order and set; by default the union of all markers'
    taxa is used. Returns

    ``alignment``
        the supermatrix.
    ``partitions``
        ``{marker: (start, end)}`` column ranges, 0-based half-open -- what a
        partitioned model needs, and what tells you which columns came from where.
    ``occupancy``
        per taxon, the fraction of markers it was actually present in.
        **Read this before the tree.** A taxon at 0.1 occupancy is placed by a
        tenth of the data with the rest gap-filled, and gaps carry no
        phylogenetic signal -- such a taxon can land almost anywhere and will
        still get high bootstrap support, because every replicate resamples the
        same missing data.
    ``markers_per_taxon`` / ``taxa_per_marker``
        the same information counted both ways.
    """
    import pandas as pd
    if not alignments:
        raise ValueError("no alignments given")
    # No per-marker checks for ragged lengths or repeated taxon names here:
    # Alignment's own constructor already rejects both, so a marker that reached
    # this function has passed them. Repeating the checks would be unreachable
    # code that reads as though it were doing something.
    order = list(alignments)
    all_taxa = (list(taxa) if taxa is not None
                else sorted({n for a in alignments.values() for n in a.names}))
    missing_requested = [t for t in all_taxa
                         if not any(t in a.names for a in alignments.values())]
    if missing_requested:
        raise ValueError(
            f"taxa absent from every marker: {missing_requested[:10]}"
            f"{' ...' if len(missing_requested) > 10 else ''}"
        )

    rows = {t: [] for t in all_taxa}
    partitions: Dict[str, Tuple[int, int]] = {}
    present: Dict[str, set] = {}
    at = 0
    for name in order:
        aln = alignments[name]
        width = aln.ncol
        by_taxon = dict(aln.records())
        present[name] = set(by_taxon)
        for t in all_taxa:
            rows[t].append(by_taxon.get(t, gap * width))
        partitions[name] = (at, at + width)
        at += width

    supermatrix = Alignment(all_taxa, ["".join(rows[t]) for t in all_taxa])
    occupancy = pd.Series(
        {t: sum(t in present[m] for m in order) / len(order) for t in all_taxa},
        name="occupancy").sort_values()
    return {
        "alignment": supermatrix,
        "partitions": partitions,
        "occupancy": occupancy,
        "markers_per_taxon": (occupancy * len(order)).astype(int),
        "taxa_per_marker": pd.Series(
            {m: len(present[m]) for m in order}, name="taxa_per_marker"),
    }


def species_tree(alignments: Mapping[str, Alignment],
                 taxa: Optional[Sequence[str]] = None,
                 min_occupancy: float = 0.0, **build_kw) -> Dict[str, object]:
    """Concatenate markers and infer one tree from the supermatrix.

    ``min_occupancy`` drops taxa present in fewer than that fraction of markers
    before inferring, which is worth setting above 0: a taxon carried by a
    handful of markers and gap-filled for the rest is positioned by almost no
    data, and bootstrap will not tell you so, because every replicate resamples
    the same absence.

    ``build_kw`` goes to :func:`~phytreon.build_tree`. For proteins at any real
    divergence use ``method="ml"`` with an external ``ml_engine``: the
    site-heterogeneous models that matter here (profile mixtures such as LG+C60,
    or LG+PMSF on large matrices) are implemented by those engines and not by
    phytreon's own likelihood, which offers only the site-homogeneous LG, WAG and
    JTT. A site-homogeneous model on a deep protein matrix does not merely lose
    resolution -- it produces long-branch attraction with high support, which is
    worse than an unresolved answer.

    Returns the concatenation result plus ``"tree"`` and ``"excluded"``.
    """
    packed = concatenate(alignments, taxa)
    occ = packed["occupancy"]
    excluded = sorted(occ.index[occ < min_occupancy]) if min_occupancy > 0 else []
    if excluded:
        keep = [t for t in packed["alignment"].names if t not in set(excluded)]
        if len(keep) < 3:
            raise ValueError(
                f"min_occupancy={min_occupancy} leaves only {len(keep)} taxa"
            )
        packed = concatenate(alignments, keep)
    from .pipeline import build_tree
    packed["tree"] = build_tree(packed["alignment"], **build_kw)
    packed["excluded"] = excluded
    return packed


def gene_trees(alignments: Mapping[str, Alignment], min_taxa: int = 4,
               **build_kw) -> Dict[str, object]:
    """Build one tree per marker, skipping markers with too few taxa.

    Returns ``{"trees": {marker: Tree}, "skipped": {marker: reason}}``. Markers
    are skipped rather than raising, because on real marker sets a few always are
    too sparse and losing the whole run to one of them is not useful.
    """
    from .pipeline import build_tree
    trees, skipped = {}, {}
    for name, aln in alignments.items():
        if aln.nseq < min_taxa:
            skipped[name] = f"only {aln.nseq} taxa (min_taxa={min_taxa})"
            continue
        try:
            trees[name] = build_tree(aln, **build_kw)
        except Exception as exc:                       # noqa: BLE001
            skipped[name] = f"{type(exc).__name__}: {exc}"
    if not trees:
        raise ValueError(f"no marker produced a tree; reasons: {skipped}")
    return {"trees": trees, "skipped": skipped}


def astrid_tree(gene_trees: Mapping[str, Tree]) -> Dict[str, object]:
    """A coalescent-aware species tree from many gene trees (ASTRID/NJst:
    Liu & Yu 2011; Vachaspati & Warnow 2015, *BMC Genomics* 16(Suppl 10):S3).

    :func:`species_tree` concatenates markers and infers one tree from the
    supermatrix, which implicitly assumes every gene shares one true tree --
    a real assumption, not a formality, and one that incomplete lineage
    sorting (ILS) breaks: under the multispecies coalescent, individual gene
    trees routinely disagree with the species tree and with each other even
    with no error or recombination involved, purely from how lineages happen
    to coalesce. Concatenation under strong ILS is not merely less accurate;
    it can be a statistically *inconsistent* estimator (Degnan & Rosenberg
    2006's "anomaly zone") -- more data does not fix it, because it is
    converging on the wrong answer.

    This instead builds, for every pair of taxa, the *topological* distance
    between them (edges on the path, not branch length -- see
    :func:`~phytreon.comparative.taxon_displacement` for the same reasoning
    applied elsewhere) on every gene tree that contains both, and averages
    those distances over all such gene trees. Liu & Yu showed this "average
    internode distance" matrix converges to one that is additive for the true
    species tree as the number of genes grows, so running neighbour-joining
    on it is a statistically consistent species-tree estimator under the
    coalescent -- unlike concatenation, and unlike trusting any single gene
    tree, or even a majority vote among them.

    This is *not* full ASTRAL (Mirarab et al. 2014), which solves a
    constrained quartet-optimisation problem by dynamic programming; NJst/
    ASTRID is a simpler, separately published, still actively used method
    with the same core consistency guarantee, reached by neighbour-joining on
    a well-chosen distance matrix instead.

    Returns ``{"tree": Tree, "distance_matrix": DataFrame, "gene_trees_per_pair":
    DataFrame}`` -- the last one says how many gene trees each pairwise
    distance was actually averaged over, since a pair covered by only one or
    two genes is a much weaker estimate than one averaged over hundreds.
    """
    import numpy as np
    import pandas as pd
    from .distance import neighbor_joining
    from ..comparative.community import patristic_distances

    if not gene_trees:
        raise ValueError("astrid_tree: no gene trees given")
    all_taxa = sorted({t for tree in gene_trees.values() for t in tree.leaf_names()})
    n = len(all_taxa)
    idx = {t: i for i, t in enumerate(all_taxa)}
    dist_sum = np.zeros((n, n))
    count = np.zeros((n, n))

    for tree in gene_trees.values():
        taxa = tree.leaf_names()
        if len(taxa) < 2:
            continue
        unit = Tree.from_newick(tree.write())
        for node in unit.traverse("postorder"):
            if node.parent is not None:
                node.length = 1.0
        names, D = patristic_distances(unit, taxa)
        gidx = np.array([idx[t] for t in names])
        dist_sum[np.ix_(gidx, gidx)] += D
        count[np.ix_(gidx, gidx)] += 1.0

    never = [(all_taxa[i], all_taxa[j]) for i in range(n) for j in range(i + 1, n)
            if count[i, j] == 0]
    if never:
        raise ValueError(
            f"astrid_tree: these taxon pairs never co-occur in any gene tree, "
            f"so no distance between them can be estimated: {never[:10]}"
            f"{' ...' if len(never) > 10 else ''}"
        )
    avg = np.divide(dist_sum, count, out=np.zeros_like(dist_sum), where=count > 0)
    tree = neighbor_joining(all_taxa, avg)
    return {
        "tree": tree,
        "distance_matrix": pd.DataFrame(avg, index=all_taxa, columns=all_taxa),
        "gene_trees_per_pair": pd.DataFrame(count.astype(int), index=all_taxa,
                                            columns=all_taxa),
    }


def gene_tree_conflict(trees: Mapping[str, Tree], reference: Tree,
                       k: int = 3) -> "pd.DataFrame":  # noqa: F821
    """Rank each gene tree's disagreement with a species tree, and say what kind.

    One row per gene, with

    ``rf``
        normalised Robinson-Foulds distance to ``reference`` over their shared
        taxa -- total disagreement.
    ``max_displacement`` / ``mean_displacement``
        the worst-placed taxon and the average, from
        :func:`~phytreon.infer.domains.taxon_displacement`.
    ``explained_by_one``
        the fraction of the conflict that dropping a single taxon accounts for
        (:func:`~phytreon.infer.domains.rogue_taxon`). **The column to sort by
        when looking for transfer.** A horizontally acquired gene disagrees about
        the recipient lineage and nowhere else, so removing it collapses the
        conflict and this approaches 1; a gene with too little signal disagrees
        everywhere and no single removal helps. Both can produce the same ``rf``
        -- in testing the signal-free gene produced a *larger* one -- and only
        the first is a biological result.
    ``worst_taxon``
        the taxon whose removal helps most: the transfer candidate.
    ``max_displacement`` / ``mean_displacement``
        neighbourhood change per taxon, for describing *how* it moved once
        ``explained_by_one`` has identified that something did.
    ``n_shared``
        how many taxa the comparison used, which bounds how much any of it means.

    Sorted by ``explained_by_one``, then by ``rf``.
    """
    import pandas as pd
    from ..treeops import prune_to_taxa
    from .domains import taxon_displacement

    ref_taxa = set(reference.leaf_names())
    rows = {}
    for name, tree in trees.items():
        shared = sorted(set(tree.leaf_names()) & ref_taxa)
        if len(shared) < max(k + 2, 5):
            continue
        gene = (prune_to_taxa(tree, shared)
                if len(shared) < len(tree.leaf_names()) else tree)
        ref = (prune_to_taxa(reference, shared)
               if len(shared) < len(ref_taxa) else reference)
        from .domains import rogue_taxon
        disp = taxon_displacement(gene, ref, k=k)
        rogue = rogue_taxon(gene, ref)
        rows[name] = {
            "rf": rogue["rf"],
            "explained_by_one": rogue["explained_by_one"],
            "worst_taxon": rogue["worst_taxon"],
            "max_displacement": float(disp.max()),
            "mean_displacement": float(disp.mean()),
            "n_shared": len(shared),
        }
    if not rows:
        raise ValueError(
            f"no gene tree shared more than k+1={k + 1} taxa with the reference"
        )
    return pd.DataFrame.from_dict(rows, orient="index").sort_values(
        ["explained_by_one", "rf"], ascending=[False, False])
