"""ASTRID/NJst: a coalescent-aware species tree from many gene trees.

Two regimes, both necessary to make the case for this over
:func:`~phytreon.species_tree`'s concatenation:

* **No discordance** -- gene trees vary only in branch length, never in
  topology. Any sane method must recover the true tree exactly here; this is
  the floor, not the interesting case.
* **Genuine incomplete lineage sorting**, simulated with a real (if
  simplified) Kingman multispecies-coalescent process on a species tree with
  short internal branches: most *individual* gene trees disagree with the
  species tree, and the single most common gene-tree topology is itself
  wrong. ASTRID's whole reason to exist is recovering the true species tree
  anyway, by averaging internode distances rather than trusting or voting on
  any one gene tree -- so that is exactly what gets tested.
"""
import matplotlib
matplotlib.use("Agg")

from collections import Counter

import numpy as np
import pytest

import phytreon as pt
from phytreon.core.tree import Node, Tree


def _simulate_coalescent_gene_tree(species_tree: Tree, rng) -> Tree:
    """One gene tree under the multispecies coalescent on ``species_tree``,
    whose branch lengths are taken directly as coalescent time units (one
    sampled lineage per tip). Standard Kingman construction: within each
    species-tree branch, lineages entering it coalesce at
    Exp(rate=C(k,2))-distributed waiting times; whatever has not coalesced by
    the top of the branch carries forward into the parent branch; the root's
    own branch is treated as unbounded so every remaining lineage eventually
    coalesces into one.
    """
    pending = {}
    for node in species_tree.traverse("postorder"):
        if node.is_leaf:
            lineages = [(Node(name=node.name, length=0.0), 0.0)]
        else:
            lineages = []
            for c in node.children:
                lineages += pending.pop(c)
        branch_len = node.length if (node.length and not node.is_root) else np.inf
        t = 0.0
        while len(lineages) > 1:
            k = len(lineages)
            wait = rng.exponential(1.0 / (k * (k - 1) / 2))
            if t + wait > branch_len:
                break
            t += wait
            i, j = rng.choice(k, size=2, replace=False)
            (a, ta), (b, tb) = lineages[i], lineages[j]
            a.length = (a.length or 0.0) + (t - ta)
            b.length = (b.length or 0.0) + (t - tb)
            parent = Node(length=0.0)
            parent.add_child(a)
            parent.add_child(b)
            lineages = [x for idx, x in enumerate(lineages) if idx not in (i, j)]
            lineages.append((parent, t))
        out = []
        for g, ts in lineages:
            g.length = (g.length or 0.0) + ((branch_len if branch_len != np.inf else t) - ts)
            out.append((g, 0.0))
        pending[node] = out
    [(root_node, _)] = pending[species_tree.root]
    return Tree(root=root_node)


# A 5-taxon caterpillar with short internal branches (coalescent units) -- the
# classic recipe for substantial ILS (Degnan & Rosenberg 2006).
HIGH_ILS_SPECIES_TREE = (
    "(((((A:0.5,B:0.5):0.15,C:0.65):0.15,D:0.8):0.15,E:0.95):0.3);"
)


def _topo_signature(tree: Tree):
    return tuple(sorted(tuple(sorted(n.leaf_names())) for n in tree.traverse()
                        if not n.is_leaf))


def test_astrid_recovers_the_tree_exactly_when_gene_trees_only_vary_in_length():
    tr = pt.datasets.random_tree(15, seed=3)
    gene_trees = {}
    for i in range(30):
        unit = Tree.from_newick(tr.write())
        r = np.random.default_rng(i)
        for n in unit.traverse("postorder"):
            if n.parent is not None:
                n.length = float(r.uniform(0.5, 2.0))
        gene_trees[f"g{i}"] = unit
    res = pt.astrid_tree(gene_trees)
    assert pt.robinson_foulds(res["tree"], tr, normalized=True) == 0.0


def test_astrid_recovers_the_species_tree_under_genuine_ils():
    sp = Tree.from_newick(HIGH_ILS_SPECIES_TREE)
    rng = np.random.default_rng(0)
    gene_trees = {f"g{i}": _simulate_coalescent_gene_tree(sp, rng) for i in range(400)}

    # confirm the regime really has substantial ILS, not just claim it: most
    # individual gene trees must disagree with the true species tree
    matching = sum(pt.robinson_foulds(t, sp, normalized=True) == 0.0
                   for t in gene_trees.values())
    assert matching / len(gene_trees) < 0.3, (
        f"only {matching}/{len(gene_trees)} gene trees matched -- expected "
        f"discordance is what this test is supposed to be checking against"
    )

    result = pt.astrid_tree(gene_trees)
    assert pt.robinson_foulds(result["tree"], sp, normalized=True) == 0.0


def test_the_naive_majority_gene_tree_topology_is_wrong_on_this_data():
    # the point of the method, made concrete: on the same discordant data
    # above, simply trusting the single most common gene-tree topology --
    # the obvious naive alternative to averaging internode distances --
    # gives the WRONG answer. ASTRID succeeding here is not "does about as
    # well as the obvious baseline", it is "succeeds where the obvious
    # baseline fails".
    sp = Tree.from_newick(HIGH_ILS_SPECIES_TREE)
    rng = np.random.default_rng(0)
    gene_trees = {f"g{i}": _simulate_coalescent_gene_tree(sp, rng) for i in range(400)}
    counts = Counter(_topo_signature(t) for t in gene_trees.values())
    most_common_topology, _freq = counts.most_common(1)[0]
    assert most_common_topology != _topo_signature(sp)


def test_astrid_tree_reports_the_distance_matrix_and_gene_tree_coverage():
    sp = Tree.from_newick(HIGH_ILS_SPECIES_TREE)
    rng = np.random.default_rng(1)
    gene_trees = {f"g{i}": _simulate_coalescent_gene_tree(sp, rng) for i in range(50)}
    result = pt.astrid_tree(gene_trees)
    taxa = sp.leaf_names()
    assert list(result["distance_matrix"].index) == sorted(taxa)
    assert result["distance_matrix"].shape == (5, 5)
    assert (result["distance_matrix"].to_numpy() >= 0).all()
    assert np.allclose(np.diag(result["distance_matrix"].to_numpy()), 0.0)
    # every gene tree here has all 5 taxa, so every pair is covered by all 50
    coverage = result["gene_trees_per_pair"].to_numpy()
    off_diag = coverage[~np.eye(5, dtype=bool)]
    assert (off_diag == 50).all()


def test_astrid_handles_gene_trees_with_partial_taxon_overlap():
    # a pair covered by only some genes should still get a distance -- the
    # average is just taken over fewer gene trees, which gene_trees_per_pair
    # reports rather than hides
    full = pt.Tree.from_newick("((A:1,B:1):1,(C:1,D:1):1);")
    partial = pt.Tree.from_newick("(A:1,C:1);")   # missing B, D
    result = pt.astrid_tree({"g1": full, "g2": partial})
    assert set(result["distance_matrix"].index) == {"A", "B", "C", "D"}
    cov = result["gene_trees_per_pair"]
    assert cov.loc["A", "C"] == 2      # in both gene trees
    assert cov.loc["A", "B"] == 1      # only in the first


def test_astrid_rejects_no_input_and_taxon_pairs_that_never_co_occur():
    with pytest.raises(ValueError, match="no gene trees given"):
        pt.astrid_tree({})
    disjoint_a = pt.Tree.from_newick("(A:1,B:1);")
    disjoint_b = pt.Tree.from_newick("(C:1,D:1);")
    with pytest.raises(ValueError, match="never co-occur"):
        pt.astrid_tree({"g1": disjoint_a, "g2": disjoint_b})
