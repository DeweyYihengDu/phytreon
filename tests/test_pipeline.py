"""Tests for the sequences->tree pipeline and tree operations."""
import matplotlib
matplotlib.use("Agg")

import phytreon as pt
from phytreon.infer import align, trim, Alignment

SEQS = [
    ("A1", "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"),
    ("A2", "ATGGCCATTGTTATGGGCCGCTGAAAGGGTGCCCGATAG"),
    ("A3", "ATGGCCATTGTAATGGGCCGCTGTAAGGGTGCCGATAG"),
    ("B1", "ATGTCGATTCTAATGAACCGCTGAAAGCGTGACCTTTAG"),
    ("B2", "ATGTCGATTCTAATGAACCGCTGTAAGCGTGACCTTTAG"),
    ("B3", "ATGTCGATTCTAATGAACCGGCTGAAAGCGTGACCTTTAG"),
]


def test_builtin_align_rectangular():
    aln = align([(n, s) for n, s in SEQS], seqtype="nucleotide")
    assert aln.nseq == 6
    # all rows aligned to equal length
    assert len({len(s) for s in aln.seqs}) == 1
    assert aln.ncol >= max(len(s) for _, s in SEQS)


def test_trim_removes_gappy_columns():
    aln = align([(n, s) for n, s in SEQS], seqtype="nucleotide")
    trimmed = trim(aln, max_gap=0.3)
    assert trimmed.ncol <= aln.ncol
    # every kept column has gap fraction <= 0.3
    for j in range(trimmed.ncol):
        assert pt.infer.column_gap_fraction(trimmed, j) <= 0.3 + 1e-9


def test_build_tree_recovers_groups():
    tree = pt.build_tree(SEQS, method="nj",
                         root="midpoint")
    kids = [frozenset(c.leaf_names()) for c in tree.root.children]
    assert frozenset({"A1", "A2", "A3"}) in kids
    assert frozenset({"B1", "B2", "B3"}) in kids


def test_bootstrap_support_range_and_signal():
    tree = pt.build_tree(SEQS, method="nj",
                         root="midpoint", bootstrap=100, seed=0)
    sups = [n.support for n in tree.traverse() if not n.is_leaf and n.support is not None]
    assert sups and all(0 <= s <= 100 for s in sups)
    # the A clade should be strongly supported
    aclade = tree.get_mrca(["A1", "A2", "A3"])
    assert aclade.support >= 80


def test_cut_tree_k_on_primates():
    tr = pt.datasets.primates()
    c = pt.cut_tree(tr, k=2)
    groups = {}
    for name, cid in c.items():
        groups.setdefault(cid, set()).add(name)
    assert {"Macaque", "Baboon"} in groups.values()


def test_rotate_reverses_order():
    tr = pt.datasets.primates()
    before = tr.leaf_names()
    pt.rotate(tr, tr.root)
    after = tr.leaf_names()
    assert before != after
    assert set(before) == set(after)


def test_collapse_low_support():
    tr = pt.datasets.primates()      # internal supports 85..100
    n_before = sum(1 for n in tr.traverse() if not n.is_leaf)
    pt.collapse_low_support(tr, threshold=92)
    n_after = sum(1 for n in tr.traverse() if not n.is_leaf)
    assert n_after < n_before        # 85 and 90 nodes collapsed


def test_midpoint_root_balances():
    tree = pt.build_tree(SEQS, method="nj")
    rooted = pt.midpoint_root(tree)
    assert len(rooted.root.children) == 2
    assert set(rooted.leaf_names()) == {n for n, _ in SEQS}


def test_native_ml_recovers_groups():
    aln = pt.align(SEQS, seqtype="nucleotide")
    tree = pt.ml_tree(aln, model="HKY85", search=True)
    # the A|B split must be present somewhere (rooting-independent check)
    sides = [frozenset(n.leaf_names()) for n in tree.traverse() if not n.is_leaf]
    assert (frozenset({"A1", "A2", "A3"}) in sides
            or frozenset({"B1", "B2", "B3"}) in sides)
    assert "logL" in tree.data


def test_nj_recovers_additive_tree():
    # NJ must return the original tree from its own (additive) patristic distances
    true_tree = pt.datasets.random_tree(10, seed=3)
    depth = {n: n.depth(use_lengths=True) for n in true_tree.traverse()}
    leaves = true_tree.leaves()
    names = [l.name for l in leaves]
    D = [[0.0] * len(leaves) for _ in leaves]
    for i, a in enumerate(leaves):
        for j in range(i + 1, len(leaves)):
            b = leaves[j]
            m = true_tree.get_mrca([a.name, b.name])
            D[i][j] = D[j][i] = depth[a] + depth[b] - 2 * depth[m]
    nj = pt.neighbor_joining(names, D)
    assert pt.robinson_foulds(true_tree, nj) == 0


def test_corrected_distance_inflates():
    from phytreon.infer import distance_matrix_model
    aln = pt.align(SEQS, seqtype="nucleotide")
    _, raw = distance_matrix_model(aln, "raw")
    _, jc = distance_matrix_model(aln, "jc69")
    assert jc[0][3] >= raw[0][3] - 1e-12       # JC correction >= raw distance


def test_nj_no_negative_branches():
    tree = pt.build_tree(SEQS, method="nj")
    assert all((n.length or 0.0) >= 0.0 for n in tree.traverse())


def test_ml_reports_aic():
    aln = pt.align(SEQS, seqtype="nucleotide")
    t = pt.ml_tree(aln, model="JC69", search=False)
    assert {"AIC", "BIC", "free_params"} <= set(t.data)


def test_robinson_foulds():
    a = pt.Tree.from_newick("((A,B),C,(D,E));")
    b = pt.Tree.from_newick("((A,C),B,(D,E));")
    assert pt.robinson_foulds(a, a) == 0.0
    assert pt.robinson_foulds(a, b) > 0.0


def test_ml_bootstrap_support():
    tree = pt.build_tree(SEQS, method="ml", ml_model="JC69", ml_search=False,
                         bootstrap=10, seed=0)
    sups = [n.support for n in tree.traverse()
            if not n.is_leaf and n.support is not None]
    assert sups and all(0 <= s <= 100 for s in sups)


def test_parsimony_inference():
    aln = pt.align(SEQS, seqtype="nucleotide")
    tr = pt.parsimony_tree(aln, search=True)
    assert "parsimony_score" in tr.data
    sides = [frozenset(n.leaf_names()) for n in tr.traverse() if not n.is_leaf]
    assert (frozenset({"A1", "A2", "A3"}) in sides
            or frozenset({"B1", "B2", "B3"}) in sides)
    good = pt.Tree.from_newick("((A1,A2),A3,((B1,B2),B3));")
    bad = pt.Tree.from_newick("((A1,B1),A2,((A3,B2),B3));")
    assert pt.parsimony_score(good, aln) < pt.parsimony_score(bad, aln)


def test_gamma_rate_heterogeneity():
    aln = pt.align(SEQS, seqtype="nucleotide")
    tr = pt.midpoint_root(pt.infer.nj_builder(aln))
    for n in tr.traverse():
        if not n.is_root:
            n.length = 0.1
    ll0 = pt.log_likelihood(tr, aln, "HKY85", gamma=0)
    llg = pt.log_likelihood(tr, aln, "HKY85", gamma=4, shape=0.5)
    assert abs(llg - ll0) > 1e-6                 # gamma changes the likelihood
    trg = pt.ml_tree(aln, model="JC69", gamma=4, search=False)
    assert trg.data["gamma_shape"] is not None and "+G" in trg.data["model"]


def test_native_ml_likelihood_prefers_true_topology():
    aln = pt.align(SEQS, seqtype="nucleotide")
    good = pt.Tree.from_newick("((A1,A2),A3,((B1,B2),B3));")
    bad = pt.Tree.from_newick("((A1,B1),A2,((A3,B2),B3));")
    for t in (good, bad):
        for n in t.traverse():
            if not n.is_root:
                n.length = 0.1
    assert pt.log_likelihood(good, aln, "HKY85") > pt.log_likelihood(bad, aln, "HKY85")
