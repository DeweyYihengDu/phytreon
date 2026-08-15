"""Community phylogenetics: patristic distances, MPD/MNTD, NRI/NTI, betaNTI,
PERMANOVA and Mantel.

Two kinds of check. The raw metrics have closed forms, so they are held against
arithmetic done by hand on a four-tip tree small enough to verify by reading.
The standardised indices do not -- they are a comparison with a null model -- so
they are checked on constructed communities whose answer is known by
construction (a real clade must come out clustered), and by the property any
correctly specified standardised effect size has to have: on *random* samples it
centres on zero with unit standard deviation.

The sign conventions get their own tests. NRI and NTI carry a factor of -1
(Webb 2000) so positive means clustered; betaNTI does not (Stegen et al. 2012)
so positive means more turnover than expected. Getting one of them backwards
inverts the biological conclusion while leaving every magnitude plausible, which
no amount of testing the arithmetic would catch.
"""
import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest

import phytreon as pt


# ((A:1,B:3):1,(C:1,D:2):2);  -- distances by hand:
#   d(A,B) = 1+3 = 4        d(C,D) = 1+2 = 3
#   d(A,C) = 1+1+2+1 = 5    d(A,D) = 1+1+2+2 = 6
#   d(B,C) = 3+1+2+1 = 7    d(B,D) = 3+1+2+2 = 8
HAND_TREE = "((A:1,B:3):1,(C:1,D:2):2);"
HAND_D = [[0, 4, 5, 6],
          [4, 0, 7, 8],
          [5, 7, 0, 3],
          [6, 8, 3, 0]]


def test_patristic_distances_match_hand_computation():
    tr = pt.Tree.from_newick(HAND_TREE)
    names, D = pt.patristic_distances(tr, ["A", "B", "C", "D"])
    assert names == ["A", "B", "C", "D"]
    assert np.allclose(D, np.array(HAND_D, dtype=float))


def test_patristic_distances_agree_with_the_covariance_identity():
    # d_ij = V_ii + V_jj - 2 V_ij; phylo_vcv is built by a separate walk, so
    # this is two independent routes to the same numbers
    tr = pt.datasets.random_tree(25, seed=4)
    names, D = pt.patristic_distances(tr)
    _, V = pt.phylo_vcv(tr, names)
    d = np.diag(V)
    assert np.allclose(D, d[:, None] + d[None, :] - 2.0 * V)


def test_patristic_distances_do_not_depend_on_where_the_tree_is_rooted():
    # the path between two tips runs through their MRCA whatever the root, so
    # rerooting moves depths but not differences of them
    tr = pt.datasets.random_tree(20, seed=9)
    names, before = pt.patristic_distances(tr)
    pt.outgroup_root(tr, [names[0]])
    _, after = pt.patristic_distances(tr, names)
    assert np.allclose(before, after)


def test_patristic_distances_reject_an_unknown_taxon():
    tr = pt.Tree.from_newick(HAND_TREE)
    with pytest.raises(ValueError, match="not found in tree"):
        pt.patristic_distances(tr, ["A", "NotReal"])


# --------------------------------------------------------------------------
# Raw metrics, against hand arithmetic on the tree above
# --------------------------------------------------------------------------
def test_mpd_and_mntd_match_hand_computation():
    tr = pt.Tree.from_newick(HAND_TREE)
    every = ["A", "B", "C", "D"]
    # MPD = mean of the six pairwise distances = (4+5+6+7+8+3)/6
    assert pt.mpd(tr, every) == pytest.approx(33 / 6)
    # MNTD: nearest kin are A-B (4), B-A (4), C-D (3), D-C (3) -> 14/4
    assert pt.mntd(tr, every) == pytest.approx(3.5)
    # a two-taxon sample: both metrics are just that one distance
    assert pt.mpd(tr, ["A", "B"]) == pytest.approx(4.0)
    assert pt.mntd(tr, ["A", "B"]) == pytest.approx(4.0)
    # fewer than two taxa has no pair to average
    assert np.isnan(pt.mpd(tr, ["A"]))
    assert np.isnan(pt.mntd(tr, ["A"]))


def test_mpd_abundance_weighting_reduces_to_the_unweighted_case():
    # equal abundances weight every pair equally, so the weighted mean is the
    # plain mean -- a property that fails if the diagonal is not excluded
    tr = pt.Tree.from_newick(HAND_TREE)
    flat = {t: 1.0 for t in ["A", "B", "C", "D"]}
    assert pt.mpd(tr, flat, abundance_weighted=True) == pytest.approx(
        pt.mpd(tr, list(flat)))
    assert pt.mntd(tr, flat, abundance_weighted=True) == pytest.approx(
        pt.mntd(tr, list(flat)))


def test_mpd_weighting_is_invariant_to_the_scale_of_the_abundances():
    tr = pt.Tree.from_newick(HAND_TREE)
    a = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
    b = {k: 1000.0 * v for k, v in a.items()}      # same profile, bigger numbers
    assert pt.mpd(tr, a, abundance_weighted=True) == pytest.approx(
        pt.mpd(tr, b, abundance_weighted=True))


def test_beta_mntd_matches_hand_computation():
    tr = pt.Tree.from_newick(HAND_TREE)
    # one taxon each side: the single cross distance, both directions
    assert pt.beta_mntd(tr, {"A": 1.0}, {"C": 1.0}) == pytest.approx(5.0)
    # {A,B} vs {C,D}, unweighted: A->min(5,6)=5, B->min(7,8)=7, mean 6;
    # C->min(5,7)=5, D->min(6,8)=6, mean 5.5; averaged -> 5.75
    assert pt.beta_mntd(tr, {"A": 1.0, "B": 1.0}, {"C": 1.0, "D": 1.0},
                        abundance_weighted=False) == pytest.approx(5.75)


def test_beta_mntd_of_a_sample_against_itself_is_zero_and_it_is_symmetric():
    tr = pt.Tree.from_newick(HAND_TREE)
    a, b = {"A": 1.0, "B": 2.0}, {"C": 3.0, "D": 1.0}
    assert pt.beta_mntd(tr, a, a) == pytest.approx(0.0)
    assert pt.beta_mntd(tr, a, b) == pytest.approx(pt.beta_mntd(tr, b, a))


def test_the_metrics_accept_a_precomputed_distance_matrix():
    # mpd/mntd/beta_mntd take a tree *or* an already-computed matrix, so a caller
    # working through hundreds of samples need not repeat patristic_distances'
    # O(n^2) walk each time. Documented and worth pinning: it is the one place in
    # the comparative API whose first argument is not a tree, so the two routes
    # have to be kept agreeing.
    tr = pt.datasets.random_tree(30, seed=1)
    names, D = pt.patristic_distances(tr)
    sample = names[:8]
    assert pt.mpd(D, sample, names=names) == pytest.approx(pt.mpd(tr, sample))
    assert pt.mntd(D, sample, names=names) == pytest.approx(pt.mntd(tr, sample))
    a, b = {names[0]: 1.0}, {names[9]: 1.0}
    assert pt.beta_mntd(D, a, b, names=names) == pytest.approx(
        pt.beta_mntd(tr, a, b))
    # a matrix built for a subset works too, given its own names
    sub_names, sub_D = pt.patristic_distances(tr, sample)
    assert pt.mpd(sub_D, sample, names=sub_names) == pytest.approx(
        pt.mpd(tr, sample))


def test_a_precomputed_matrix_without_names_says_so():
    tr = pt.datasets.random_tree(12, seed=1)
    names, D = pt.patristic_distances(tr)
    with pytest.raises(ValueError, match="pass names="):
        pt.mpd(D, names[:5])


# --------------------------------------------------------------------------
# Standardised indices, and the sign conventions
# --------------------------------------------------------------------------
def _clade_and_random_table(n_tips=60, seed=3):
    import pandas as pd
    tr = pt.datasets.random_tree(n_tips, seed=seed)
    tips = tr.leaf_names()
    clade = next(nd.leaf_names() for nd in tr.traverse("postorder")
                 if not nd.is_leaf and 6 <= len(nd.leaf_names()) <= 10)
    rng = np.random.default_rng(0)
    rows = [{t: 1.0 if t in set(clade) else 0.0 for t in tips}]
    index = ["clade"]
    for k in range(30):
        pick = set(rng.choice(tips, len(clade), replace=False))
        rows.append({t: 1.0 if t in pick else 0.0 for t in tips})
        index.append(f"rand{k}")
    return tr, pd.DataFrame(rows, index=index)


def test_nri_and_nti_are_positive_for_a_real_clade():
    # the -1 in Webb's definition means POSITIVE is clustered. A clade's tips are
    # as clustered as a sample can be, so both indices must come out well above
    # zero -- if either convention were inverted this is the test that catches it
    tr, table = _clade_and_random_table()
    nri = pt.ses_mpd(tr, table, n_null=299, seed=1)
    nti = pt.ses_mntd(tr, table, n_null=299, seed=1)
    assert nri.loc["clade", "NRI"] > 2.0
    assert nti.loc["clade", "NTI"] > 2.0
    assert nri.loc["clade", "p"] < 0.05
    assert nti.loc["clade", "p"] < 0.05
    # the observed MPD must be *below* the null mean for a clustered sample --
    # the raw direction, before the sign flip that turns it into NRI
    assert nri.loc["clade", "mpd"] < nri.loc["clade", "null_mean"]


def test_random_samples_centre_the_indices_on_zero_with_unit_spread():
    # what any correctly specified standardised effect size has to do, and the
    # check that the null model itself is right rather than merely present: draw
    # samples with no phylogenetic structure and the index should have mean 0 and
    # standard deviation 1. A wrong null shows up here as a shifted centre or a
    # spread away from 1, even while the clustered case above still passes.
    tr, table = _clade_and_random_table()
    nri = pt.ses_mpd(tr, table, n_null=299, seed=1)["NRI"][1:]
    nti = pt.ses_mntd(tr, table, n_null=299, seed=1)["NTI"][1:]
    for label, values in (("NRI", nri), ("NTI", nti)):
        assert abs(values.mean()) < 0.6, f"{label} mean {values.mean():.3f}"
        assert 0.5 < values.std() < 1.8, f"{label} sd {values.std():.3f}"


def test_beta_nti_is_negative_within_a_clade_and_positive_across_clades():
    # betaNTI has NO sign flip (Stegen et al. 2012), the opposite convention to
    # NRI/NTI above: positive means MORE phylogenetic turnover than chance.
    # Two halves of one clade share close relatives, so turnover is less than
    # chance -> negative; two different clades -> positive.
    import pandas as pd
    tr = pt.datasets.random_tree(60, seed=3)
    tips = tr.leaf_names()
    clades = [nd.leaf_names() for nd in tr.traverse("postorder")
              if not nd.is_leaf and 8 <= len(nd.leaf_names()) <= 14]
    first, last = clades[0], clades[-1]
    first = [t for t in first if t not in set(last)]
    last = [t for t in last if t not in set(first)]

    def row(taxa):
        return {t: 1.0 if t in set(taxa) else 0.0 for t in tips}

    table = pd.DataFrame(
        [row(first[:len(first) // 2]), row(first[len(first) // 2:]),
         row(last[:len(last) // 2])],
        index=["one_a", "one_b", "other"])
    bn = pt.beta_nti(tr, table, n_null=299, seed=0)
    assert bn.loc["one_a", "one_b"] < 0.0
    assert bn.loc["one_a", "other"] > 2.0
    assert np.allclose(np.diag(bn), 0.0)
    assert np.allclose(bn.to_numpy(), bn.to_numpy().T)


def test_beta_nti_null_applies_one_relabelling_to_both_samples():
    # The null relabels tips once per draw and applies that *one* relabelling to
    # both samples, so a taxon present in both stays present in both. Two
    # identical samples are the sharp test of it: relabelling them together
    # leaves them identical, so every null draw gives betaMNTD exactly 0, the
    # null has zero variance, and betaNTI is undefined -- NaN is the correct
    # answer here and is itself the evidence.
    #
    # Permuting the two independently would instead give the null a spread
    # (measured: sd 0.070 around a mean of 0.302 on this tree) and so return a
    # confident-looking finite value of about -4.3. A number here would be the
    # bug; NaN is the fix.
    import pandas as pd
    tr = pt.datasets.random_tree(40, seed=5)
    tips = tr.leaf_names()
    same = {t: 1.0 if t in set(tips[:10]) else 0.0 for t in tips}
    table = pd.DataFrame([same, dict(same)], index=["x", "x_copy"])
    bn = pt.beta_nti(tr, table, n_null=199, seed=0)
    assert np.isnan(bn.loc["x", "x_copy"])

    # and the underlying reason, checked directly rather than inferred
    from phytreon.comparative.community import (
        _beta_mntd_at, _null_indices, _resolve_sample, patristic_distances)
    names, D = patristic_distances(tr)
    ia, wa = _resolve_sample({t: 1.0 for t in tips[:10]}, names)
    rng = np.random.default_rng(0)
    nulls = [_beta_mntd_at(D, p[ia], wa, p[ia], wa, True)
             for p in _null_indices(len(names), 50, rng)]
    assert all(v == 0.0 for v in nulls)


def test_ses_functions_reject_table_columns_that_are_not_tips():
    import pandas as pd
    tr = pt.datasets.random_tree(12, seed=1)
    tips = tr.leaf_names()
    table = pd.DataFrame([{**{t: 1.0 for t in tips}, "Ghost": 1.0}], index=["s1"])
    for call in (lambda: pt.ses_mpd(tr, table, n_null=9),
                 lambda: pt.ses_mntd(tr, table, n_null=9),
                 lambda: pt.beta_nti(tr, table, n_null=9)):
        with pytest.raises(ValueError, match="not tips of the tree"):
            call()


# --------------------------------------------------------------------------
# PERMANOVA and Mantel
# --------------------------------------------------------------------------
def _blocked_distances(n, sep, rng):
    import pandas as pd
    co = np.vstack([rng.normal(0, 1, (n // 2, 3)), rng.normal(sep, 1, (n // 2, 3))])
    D = np.linalg.norm(co[:, None] - co[None, :], axis=-1)
    return pd.DataFrame(D), ["g1"] * (n // 2) + ["g2"] * (n // 2)


def test_permanova_detects_groups_that_really_differ():
    rng = np.random.default_rng(0)
    D, labels = _blocked_distances(24, sep=4.0, rng=rng)
    res = pt.permanova(D, labels, n_perm=499, seed=0)
    assert res["p"] < 0.01
    assert res["R2"] > 0.3
    assert res["pseudo_F"] > 1.0
    assert res["groups"] == {"g1": 12, "g2": 12}


def test_permanova_r2_is_near_zero_when_the_groups_are_arbitrary():
    import pandas as pd
    rng = np.random.default_rng(1)
    co = rng.normal(0, 1, (24, 3))
    D = pd.DataFrame(np.linalg.norm(co[:, None] - co[None, :], axis=-1))
    res = pt.permanova(D, ["g1"] * 12 + ["g2"] * 12, n_perm=499, seed=0)
    assert res["R2"] < 0.15
    assert 0.0 < res["p"] <= 1.0


def test_permanova_needs_two_groups_with_two_samples_each():
    import pandas as pd
    rng = np.random.default_rng(2)
    co = rng.normal(0, 1, (6, 3))
    D = pd.DataFrame(np.linalg.norm(co[:, None] - co[None, :], axis=-1))
    with pytest.raises(ValueError, match="at least 2 groups"):
        pt.permanova(D, ["g"] * 6, n_perm=9)
    with pytest.raises(ValueError, match="2 samples per group"):
        pt.permanova(D, ["a"] * 5 + ["b"], n_perm=9)


def test_permanova_takes_groups_as_a_series_indexed_by_sample():
    import pandas as pd
    rng = np.random.default_rng(3)
    D, labels = _blocked_distances(20, sep=5.0, rng=rng)
    D.index = D.columns = [f"s{i}" for i in range(20)]
    # deliberately out of order -- a Series must be aligned by label, not position
    series = pd.Series(labels, index=D.index).sample(frac=1.0, random_state=0)
    res = pt.permanova(D, series, n_perm=299, seed=0)
    assert res["p"] < 0.01
    assert res["groups"] == {"g1": 10, "g2": 10}


def test_mantel_of_a_matrix_against_itself_is_exactly_one():
    rng = np.random.default_rng(0)
    D, _ = _blocked_distances(20, sep=3.0, rng=rng)
    res = pt.mantel(D, D, n_perm=299, seed=0)
    assert res["r"] == pytest.approx(1.0)
    assert res["p"] < 0.05
    assert res["n"] == 20


def test_mantel_finds_little_between_independent_matrices():
    import pandas as pd
    rng = np.random.default_rng(1)
    out = []
    for _ in range(3):
        a = rng.normal(0, 1, (20, 3))
        b = rng.normal(0, 1, (20, 3))
        Da = pd.DataFrame(np.linalg.norm(a[:, None] - a[None, :], axis=-1))
        Db = pd.DataFrame(np.linalg.norm(b[:, None] - b[None, :], axis=-1))
        out.append(abs(pt.mantel(Da, Db, n_perm=299, seed=0)["r"]))
    assert max(out) < 0.5


def test_mantel_spearman_sees_a_monotone_but_curved_relationship():
    import pandas as pd
    rng = np.random.default_rng(2)
    co = rng.normal(0, 1, (20, 3))
    D = np.linalg.norm(co[:, None] - co[None, :], axis=-1)
    curved = pd.DataFrame(np.exp(D))        # monotone in D, far from linear
    sp = pt.mantel(pd.DataFrame(D), curved, method="spearman", n_perm=299, seed=0)
    # a rank correlation is exactly 1 for any monotone transform
    assert sp["r"] == pytest.approx(1.0)
    assert sp["method"] == "spearman"


def test_mantel_rejects_an_unknown_method_and_mismatched_sizes():
    import pandas as pd
    rng = np.random.default_rng(3)
    co = rng.normal(0, 1, (10, 3))
    D = pd.DataFrame(np.linalg.norm(co[:, None] - co[None, :], axis=-1))
    with pytest.raises(ValueError, match="pearson.*spearman"):
        pt.mantel(D, D, method="kendall", n_perm=9)
    with pytest.raises(ValueError, match="different sizes"):
        pt.mantel(np.asarray(D), np.asarray(D)[:5, :5], n_perm=9)
