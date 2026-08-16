"""Per-domain trees and marker-gene species trees, with conflict analysis.

Both workflows exist for the same reason: a tree's *disagreement* with another
tree is often the result, not the noise. So the tests are built around simulated
data whose answer is known by construction -- sequences evolved down one tree for
one domain or marker and down a deliberately different tree for another, so the
lineage that "moved" is known in advance and can be checked for by name.

The important test is not that conflict is detected but that the two kinds of
conflict are told apart:

* one lineage genuinely misplaced (recombination, horizontal transfer)
* nothing resolvable at all (a saturated or too-short region)

They are distinguishable and they mean opposite things, and on this data the
signal-free case has the *larger* Robinson-Foulds distance -- so any ranking by
total conflict puts the useless region above the real finding.
"""
import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest

import phytreon as pt

AA = "ACDEFGHIKLMNPQRSTVWY"

# Two trees over the same eight taxa, differing only in where X sits: in TREE_A
# it is sister to C inside the first clade, in TREE_B it has moved into the
# second. Anything evolved down both should place X differently and nothing else.
TREE_A = "(((A:.1,B:.1):.3,(C:.1,X:.1):.3):.4,((D:.1,E:.1):.3,(F:.1,G:.1):.3):.4);"
TREE_B = "(((A:.1,B:.1):.3,C:.4):.4,((D:.1,E:.1):.3,((F:.1,G:.1):.2,X:.3):.1):.4);"


def _evolve(newick, length, rate, seed):
    """Evolve a protein down a tree: independent sites, uniform replacement."""
    tree = pt.Tree.from_newick(newick)
    rng = np.random.default_rng(seed)
    root = "".join(rng.choice(list(AA), length))
    seqs = {}

    def walk(node, parent):
        seq = list(parent)
        p = 1.0 - np.exp(-rate * (node.length or 0.0))
        for i in np.flatnonzero(rng.random(length) < p):
            seq[i] = rng.choice(list(AA))
        seq = "".join(seq)
        if node.is_leaf:
            seqs[node.name] = seq
        for child in node.children:
            walk(child, seq)

    for child in tree.root.children:
        walk(child, root)
    return seqs


@pytest.fixture(scope="module")
def evolved():
    vertical = _evolve(TREE_A, 300, 1.0, 0)      # follows TREE_A
    moved = _evolve(TREE_B, 300, 1.0, 1)         # follows TREE_B -- X relocated
    saturated = _evolve(TREE_A, 300, 25.0, 2)    # right tree, no signal left
    return {"vertical": vertical, "moved": moved, "saturated": saturated,
            "taxa": sorted(vertical)}


# --------------------------------------------------------------------------
# residue -> column, the boundary arithmetic
# --------------------------------------------------------------------------
def test_residue_ranges_map_to_the_right_alignment_columns():
    # HMMER counts residues on a sequence; an alignment is indexed by columns
    # that include gaps. Hand-checked: "AB--CDE-FG" has residues A B C D E F G at
    # columns 0 1 4 5 6 8 9.
    aln = pt.Alignment(["ref", "other"], ["AB--CDE-FG", "ABXXCDEXFG"])
    got = pt.residue_to_column(aln, "ref", {
        "first_two": (1, 2), "middle": (3, 5), "last_two": (6, 7), "all": (1, 7)})
    assert got == {"first_two": (0, 2), "middle": (4, 7),
                   "last_two": (8, 10), "all": (0, 10)}
    # and the columns it returns really do cut those residues out
    parts = pt.split_domains(aln, {"middle": (3, 5)}, reference="ref")
    assert parts["middle"].seqs[0] == "CDE"


def test_residue_ranges_are_one_based_inclusive_and_columns_are_half_open():
    # deliberately different conventions, so confusing them fails loudly rather
    # than shifting a domain boundary by one residue
    aln = pt.Alignment(["ref"], ["ABCDE"])
    assert pt.residue_to_column(aln, "ref", {"d": (1, 5)}) == {"d": (0, 5)}
    assert pt.residue_to_column(aln, "ref", {"d": (2, 2)}) == {"d": (1, 2)}


def test_residue_to_column_rejects_ranges_it_cannot_honour():
    aln = pt.Alignment(["ref", "other"], ["AB--CDE-FG", "ABXXCDEXFG"])
    with pytest.raises(ValueError, match="only 7 residues"):
        pt.residue_to_column(aln, "ref", {"d": (1, 99)})
    with pytest.raises(ValueError, match="1-based inclusive"):
        pt.residue_to_column(aln, "ref", {"d": (5, 2)})
    with pytest.raises(ValueError, match="1-based inclusive"):
        pt.residue_to_column(aln, "ref", {"d": (0, 3)})
    with pytest.raises(ValueError, match="not in the alignment"):
        pt.residue_to_column(aln, "absent", {"d": (1, 2)})


def test_split_domains_rejects_columns_outside_the_alignment():
    aln = pt.Alignment(["a", "b"], ["ACGT", "ACGT"])
    with pytest.raises(ValueError, match="outside the alignment"):
        pt.split_domains(aln, {"d": (0, 99)})
    with pytest.raises(ValueError, match="below min_columns"):
        pt.split_domains(aln, {"d": (0, 2)}, min_columns=3)
    with pytest.raises(ValueError, match="no domains given"):
        pt.split_domains(aln, {})


def test_domain_trees_drops_and_reports_taxa_missing_a_domain():
    # a domain present in only some taxa is a finding; entering an all-gap row
    # would bury it as a zero-length branch
    taxa = [f"t{i}" for i in range(6)]
    left = ["ACDEFGHIKL"] * 6
    right = ["MNPQRSTVWY"] * 5 + ["----------"]      # t5 lacks the second domain
    aln = pt.Alignment(taxa, [a + b for a, b in zip(left, right)])
    res = pt.domain_trees(aln, {"one": (0, 10), "two": (10, 20)}, method="nj")
    assert res["dropped"]["one"] == []
    assert res["dropped"]["two"] == ["t5"]
    assert set(res["trees"]["two"].leaf_names()) == set(taxa[:5])
    assert res["columns"] == {"one": (0, 10), "two": (10, 20)}


# --------------------------------------------------------------------------
# taxon_displacement: topological, not branch-length based
# --------------------------------------------------------------------------
def test_displacement_is_zero_between_a_tree_and_itself():
    tr = pt.datasets.random_tree(15, seed=1)
    assert pt.taxon_displacement(tr, tr).max() == pytest.approx(0.0)


def test_displacement_ignores_branch_lengths_by_default():
    # Genes evolve at different rates, so two trees of identical shape routinely
    # differ in branch length. A distance-based comparison calls that
    # displacement; the topological default must not. This was a real bug -- a
    # gene tree matching its species tree exactly (RF 0) still scored
    # displacements up to 0.5 before the measure was made topological.
    tr = pt.Tree.from_newick(
        "(((A:.1,B:.1):.3,(C:.1,D:.1):.3):.4,((E:.1,F:.1):.3,(G:.1,H:.1):.3):.4);")
    stretched = pt.Tree.from_newick(
        "(((A:.9,B:.05):.7,(C:.4,D:.02):.1):.9,((E:.7,F:.3):.2,(G:.1,H:.8):.6):.3);")
    assert pt.robinson_foulds(tr, stretched, normalized=True) == 0.0
    assert pt.taxon_displacement(tr, stretched).max() == pytest.approx(0.0)
    # with the topological guard off, the same pair does register a difference
    by_distance = pt.taxon_displacement(tr, stretched, topological=False)
    assert by_distance.max() > 0.0


def test_displacement_names_the_taxon_that_moved(evolved):
    a = pt.build_tree(pt.Alignment(evolved["taxa"],
                                   [evolved["vertical"][t] for t in evolved["taxa"]]),
                      method="nj")
    b = pt.build_tree(pt.Alignment(evolved["taxa"],
                                   [evolved["moved"][t] for t in evolved["taxa"]]),
                      method="nj")
    disp = pt.taxon_displacement(a, b)
    assert disp.index[0] == "X"
    assert disp["X"] > 0.0


def test_displacement_needs_more_taxa_than_k():
    tr = pt.Tree.from_newick("(A:1,(B:1,C:1):1);")
    with pytest.raises(ValueError, match="shared taxa"):
        pt.taxon_displacement(tr, tr, k=3)


# --------------------------------------------------------------------------
# rogue_taxon: does dropping one lineage explain the conflict?
# --------------------------------------------------------------------------
def test_rogue_taxon_identifies_a_single_relocated_lineage(evolved):
    taxa = evolved["taxa"]
    a = pt.build_tree(pt.Alignment(taxa, [evolved["vertical"][t] for t in taxa]),
                      method="nj")
    b = pt.build_tree(pt.Alignment(taxa, [evolved["moved"][t] for t in taxa]),
                      method="nj")
    res = pt.rogue_taxon(a, b)
    assert res["worst_taxon"] == "X"
    # removing X leaves the two trees identical -- the whole conflict was X
    assert res["explained_by_one"] == pytest.approx(1.0)
    assert res["rf_without"]["X"] == pytest.approx(0.0)
    # and no other single removal achieves that
    others = res["rf_without"].drop("X")
    assert others.min() > 0.0


def test_rogue_taxon_finds_no_single_culprit_for_diffuse_conflict(evolved):
    taxa = evolved["taxa"]
    a = pt.build_tree(pt.Alignment(taxa, [evolved["vertical"][t] for t in taxa]),
                      method="nj")
    noise = pt.build_tree(pt.Alignment(taxa, [evolved["saturated"][t] for t in taxa]),
                          method="nj")
    res = pt.rogue_taxon(a, noise)
    # the point of the measure: this pair's total conflict is LARGER than the
    # real single-lineage move above, yet no one taxon accounts for it
    assert res["rf"] > 0.5
    assert res["explained_by_one"] < 0.6


def test_rogue_taxon_reports_no_conflict_as_nan_not_as_a_number():
    tr = pt.datasets.random_tree(12, seed=2)
    res = pt.rogue_taxon(tr, tr)
    assert res["rf"] == 0.0
    assert np.isnan(res["explained_by_one"])
    assert res["worst_taxon"] == ""


def test_rogue_taxon_needs_enough_taxa_to_drop_one():
    tr = pt.Tree.from_newick("((A:1,B:1):1,(C:1,D:1):1);")
    with pytest.raises(ValueError, match="at least 5 shared taxa"):
        pt.rogue_taxon(tr, tr)


# --------------------------------------------------------------------------
# The domain workflow end to end
# --------------------------------------------------------------------------
def test_domain_trees_recover_a_recombined_lineage(evolved):
    # the OCP-shaped case: one protein, two domains, different histories for one
    # lineage. The disagreement is the finding.
    taxa = evolved["taxa"]
    fused = pt.Alignment(taxa, [evolved["vertical"][t] + evolved["moved"][t]
                                for t in taxa])
    res = pt.domain_trees(fused, {"NTD": (0, 300), "CTD": (300, 600)}, method="nj")
    assert set(res["trees"]) == {"NTD", "CTD"}
    cmp = pt.compare_domain_trees(res["trees"])
    assert cmp["rf"].loc["NTD", "CTD"] > 0.0
    assert cmp["rf"].loc["NTD", "CTD"] == cmp["rf"].loc["CTD", "NTD"]
    assert cmp["worst_taxon"]["NTD|CTD"] == "X"
    assert cmp["explained_by_one"]["NTD|CTD"] == pytest.approx(1.0)
    assert cmp["displacement"]["NTD|CTD"].index[0] == "X"


def test_a_domain_with_no_signal_is_not_mistaken_for_recombination(evolved):
    taxa = evolved["taxa"]
    aln = pt.Alignment(taxa, [evolved["vertical"][t] + evolved["saturated"][t]
                              for t in taxa])
    res = pt.domain_trees(aln, {"real": (0, 300), "noise": (300, 600)}, method="nj")
    cmp = pt.compare_domain_trees(res["trees"])
    # larger total conflict than the genuine recombination above ...
    assert cmp["rf"].loc["real", "noise"] > 0.5
    # ... and yet no single lineage explains it
    assert cmp["explained_by_one"]["real|noise"] < 0.6


def test_compare_domain_trees_needs_two_trees():
    tr = pt.datasets.random_tree(8, seed=1)
    with pytest.raises(ValueError, match="at least 2 domain trees"):
        pt.compare_domain_trees({"only": tr})


# --------------------------------------------------------------------------
# The marker workflow end to end
# --------------------------------------------------------------------------
def _marker_set(evolved, n_vertical=7):
    taxa = evolved["taxa"]
    markers = {f"rpl{i}": pt.Alignment(
        taxa, [_evolve(TREE_A, 250, 1.0, 10 + i)[t] for t in taxa])
        for i in range(n_vertical)}
    markers["hgt_gene"] = pt.Alignment(
        taxa, [_evolve(TREE_B, 250, 1.0, 99)[t] for t in taxa])
    markers["junk_gene"] = pt.Alignment(
        taxa, [_evolve(TREE_A, 250, 25.0, 7)[t] for t in taxa])
    partial = [t for t in taxa if t not in ("F", "G")]
    markers["patchy"] = pt.Alignment(
        partial, [_evolve(TREE_A, 250, 1.0, 3)[t] for t in partial])
    return markers


def test_concatenate_builds_a_supermatrix_with_partitions_and_occupancy(evolved):
    markers = _marker_set(evolved, n_vertical=2)
    packed = pt.concatenate(markers)
    aln = packed["alignment"]
    assert aln.nseq == len(evolved["taxa"])
    assert aln.ncol == sum(m.ncol for m in markers.values())
    # partitions tile the supermatrix end to end with no gap or overlap
    spans = sorted(packed["partitions"].values())
    assert spans[0][0] == 0
    assert spans[-1][1] == aln.ncol
    for (_, end), (start, _) in zip(spans, spans[1:]):
        assert end == start
    # the taxa missing from one marker show a lower occupancy, and the rest 1.0
    occ = packed["occupancy"]
    assert occ["F"] < 1.0 and occ["G"] < 1.0
    assert occ.drop(["F", "G"]).min() == pytest.approx(1.0)
    assert packed["taxa_per_marker"]["patchy"] == len(evolved["taxa"]) - 2


def test_concatenate_gap_fills_a_missing_taxon_for_that_markers_columns(evolved):
    taxa = ["a", "b", "c"]
    m1 = pt.Alignment(taxa, ["ACDE", "ACDE", "ACDE"])
    m2 = pt.Alignment(["a", "b"], ["FGHIK", "FGHIK"])       # c absent
    packed = pt.concatenate({"m1": m1, "m2": m2})
    rows = dict(packed["alignment"].records())
    assert rows["a"] == "ACDEFGHIK"
    assert rows["c"] == "ACDE" + "-" * 5
    assert packed["occupancy"]["c"] == pytest.approx(0.5)


def test_concatenate_rejects_input_it_cannot_use():
    with pytest.raises(ValueError, match="no alignments given"):
        pt.concatenate({})
    with pytest.raises(ValueError, match="absent from every marker"):
        pt.concatenate({"m": pt.Alignment(["a", "b"], ["AC", "AC"])},
                       taxa=["a", "b", "ghost"])


def test_alignment_itself_rejects_ragged_and_duplicated_input():
    # concatenate deliberately does not re-check these: Alignment's constructor
    # already does, so a marker that reaches concatenate has passed them, and
    # duplicating the checks there would be unreachable code. Asserted here so
    # that stays true -- if Alignment ever stops enforcing it, this fails rather
    # than concatenate silently accepting a ragged marker.
    with pytest.raises(ValueError, match="same length"):
        pt.Alignment(["a", "b"], ["ACGT", "ACG"])
    with pytest.raises(ValueError, match="unique"):
        pt.Alignment(["a", "a"], ["ACGT", "ACGT"])


def test_species_tree_recovers_the_tree_the_markers_evolved_down(evolved):
    markers = _marker_set(evolved)
    res = pt.species_tree(markers, method="nj")
    # seven vertical markers outvote one transferred and one saturated gene
    assert pt.robinson_foulds(res["tree"], pt.Tree.from_newick(TREE_A),
                              normalized=True) == 0.0
    assert res["excluded"] == []


def test_species_tree_min_occupancy_excludes_thinly_covered_taxa(evolved):
    markers = _marker_set(evolved)
    res = pt.species_tree(markers, min_occupancy=0.95, method="nj")
    # F and G are missing from the patchy marker, so they fall below 0.95
    assert set(res["excluded"]) == {"F", "G"}
    assert "F" not in res["tree"].leaf_names()
    assert res["alignment"].nseq == len(evolved["taxa"]) - 2


def test_gene_tree_conflict_ranks_the_transferred_gene_first(evolved):
    markers = _marker_set(evolved)
    reference = pt.species_tree(markers, method="nj")["tree"]
    built = pt.gene_trees(markers, method="nj")
    assert built["skipped"] == {}
    conflict = pt.gene_tree_conflict(built["trees"], reference)

    assert conflict.index[0] == "hgt_gene"
    assert conflict.loc["hgt_gene", "worst_taxon"] == "X"
    assert conflict.loc["hgt_gene", "explained_by_one"] == pytest.approx(1.0)

    # the vertical markers and the patchy one show no conflict at all
    for name in [c for c in conflict.index if c.startswith("rpl")] + ["patchy"]:
        assert conflict.loc[name, "rf"] == 0.0
        assert np.isnan(conflict.loc[name, "explained_by_one"])

    # THE point of the measure: the signal-free gene has a LARGER total conflict
    # than the real transfer, so ranking by rf would bury the biology under it
    assert conflict.loc["junk_gene", "rf"] > conflict.loc["hgt_gene", "rf"]
    assert conflict.loc["junk_gene", "explained_by_one"] < 0.6


def test_gene_trees_skips_rather_than_fails_on_a_marker_it_cannot_use(evolved):
    markers = _marker_set(evolved, n_vertical=1)
    markers["too_few"] = pt.Alignment(["a", "b"], ["ACDE", "ACDF"])
    built = pt.gene_trees(markers, method="nj")
    assert "too_few" in built["skipped"]
    assert "too_few" not in built["trees"]
    assert len(built["trees"]) >= 3


def test_gene_tree_conflict_needs_a_shared_taxon_set_worth_comparing(evolved):
    tiny = pt.build_tree(pt.Alignment(["a", "b", "c", "d"],
                                      ["ACDE", "ACDF", "ACDG", "ACDH"]),
                         method="nj")
    with pytest.raises(ValueError, match="no gene tree shared"):
        pt.gene_tree_conflict({"t": tiny}, pt.datasets.random_tree(9, seed=1))
