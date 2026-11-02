"""Rough edges a first-time user hits: indexed metadata, missing files,
alignment types, ambiguity codes, copying a tree before editing it."""
import math

import pytest

import phytreon as pt
from phytreon.core.io import to_newick
from phytreon.infer.bootstrap import distance_matrix_model


def test_join_data_accepts_a_table_indexed_by_tip_name():
    meta = pt.datasets.primates_metadata()          # indexed by "name"
    tree = pt.datasets.primates().join_data(meta)
    assert all("habitat" in leaf.data for leaf in tree.leaves())


def test_join_data_warns_when_nothing_matches():
    meta = pt.datasets.primates_metadata().rename(index=str.upper)
    with pytest.warns(UserWarning, match="none of the 7 rows matched"):
        pt.datasets.primates().join_data(meta)


def test_join_data_names_the_missing_key():
    with pytest.raises(KeyError, match="no column 'species'"):
        pt.datasets.primates().join_data(pt.datasets.primates_metadata(),
                                         on="species")


def test_missing_fasta_file_says_so():
    with pytest.raises(FileNotFoundError, match="nope.fasta"):
        pt.build_tree("nope.fasta")


def test_tree_copy_is_independent():
    tree = pt.datasets.primates().join_data(pt.datasets.primates_metadata())
    before = to_newick(tree)
    copy = tree.copy()
    assert to_newick(copy) == before
    pt.collapse_low_support(copy, 101)
    copy.leaves()[0].data["habitat"] = "changed"
    assert to_newick(tree) == before
    assert tree.leaves()[0].data["habitat"] != "changed"


def test_distance_matrix_takes_a_phytreon_alignment_and_model_names():
    aln = pt.Alignment(["x", "y"], ["ACGTACGTAC", "ACGTACGTAA"])
    names, d = pt.distance_matrix(aln, "jc69")
    assert names == ["x", "y"]
    assert d[0][1] == pytest.approx(-0.75 * math.log(1 - 4 / 3 * 0.1))
    assert pt.distance_matrix(aln)[1][0][1] == pytest.approx(0.1)  # identity


def test_an_ambiguity_code_does_not_switch_off_the_jc69_correction():
    with_r = pt.Alignment(["x", "y"], ["ACGTACGTRC", "ACGTACGTAA"])
    with_n = pt.Alignment(["x", "y"], ["ACGTACGTNC", "ACGTACGTAA"])
    assert distance_matrix_model(with_r, "jc69")[1][0][1] == \
        pytest.approx(distance_matrix_model(with_n, "jc69")[1][0][1])
    assert distance_matrix_model(with_r, "jc69")[1][0][1] > 1 / 9


@pytest.mark.parametrize("model", ["jc69", "k2p", "poisson"])
def test_corrected_distances_never_shrink_as_divergence_grows(model):
    ref = "A" * 100
    d = [distance_matrix_model(pt.Alignment(["a", "b"],
                                            [ref, "C" * k + "A" * (100 - k)]),
                               model)[1][0][1]
         for k in range(101)]
    assert all(x <= y + 1e-12 for x, y in zip(d, d[1:]))
