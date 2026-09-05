"""GC content, CpG islands, tandem repeats -- validated on their own terms,
not just "runs without error".

**A real calibration subtlety, found while writing these tests and worth
recording.** The first attempt at a CpG-island negative control used plain
i.i.d. random ACGT sequence and got a 30/30 false-positive rate -- which
looked like a broken implementation but was a wrong test instead. Uniform
random ACGT has ~50% GC and an observed/expected CpG ratio close to 1.0 by
construction (nothing is suppressing CpG), so it clears both
``cpg_islands`` thresholds by design, not by bug. Real genomic bulk
sequence is nothing like that: cytosine methylation drives CpG well below
random expectation over evolutionary time (typical obs/exp ~0.2-0.3), and
a CpG island is meaningful precisely because it is a region where that
suppression is absent. The tests below build a CpG-*suppressed* background
(reject any generated "CG" dinucleotide) as the actual negative control,
which is what makes the Type-I check below mean anything.

**Tandem repeats, similarly**: short repeats (period 1-2, a handful of
copies) occur constantly by chance in any random sequence -- real genomes
are full of them too. The default settings are not expected to have a low
false-positive rate at that scale, and the tests don't claim otherwise;
the meaningful negative control here is at a stricter `min_copies`, where
a long chance repeat becomes genuinely rare.
"""
import numpy as np
import pytest

import phytreon as pt

BASES = np.array(list("ACGT"))


def _random_seq(n, seed):
    return "".join(np.random.default_rng(seed).choice(BASES, size=n))


def _cpg_suppressed_seq(n, seed):
    # reject any 'CG' dinucleotide -- crude but real: mimics the
    # methylation-driven CpG suppression that makes actual genomic
    # background nothing like i.i.d. random sequence (see module docstring)
    rng = np.random.default_rng(seed)
    seq = list(rng.choice(BASES, size=n))
    for i in range(n - 1):
        while seq[i] == "C" and seq[i + 1] == "G":
            seq[i + 1] = rng.choice(BASES)
    return "".join(seq)


# --------------------------------------------------------------------------
# sequence_lengths
# --------------------------------------------------------------------------
def test_sequence_lengths_counts_ungapped_bases():
    aln = pt.Alignment(["a", "b"], ["AC-GT", "ACGT-"])
    lengths = pt.sequence_lengths(aln)
    assert lengths["a"] == 4
    assert lengths["b"] == 4


# --------------------------------------------------------------------------
# gc_content
# --------------------------------------------------------------------------
def test_gc_content_matches_hand_computation():
    # GGCC (4/4 GC) then AATT (0/4 GC), non-overlapping window=4
    aln = pt.Alignment(["s"], ["GGCCAATT"])
    res = pt.gc_content(aln, window=4)
    assert list(res["start"]) == [0, 4]
    assert list(res["gc"]) == pytest.approx([1.0, 0.0])


def test_gc_content_overlapping_windows_via_step():
    aln = pt.Alignment(["s"], ["GCGCAAAA"])   # GC then AT-only
    res = pt.gc_content(aln, window=4, step=2)
    assert list(res["start"]) == [0, 2, 4]
    # window [0,4)=GCGC->1.0, [2,6)=GCAA->0.5, [4,8)=AAAA->0.0
    assert list(res["gc"]) == pytest.approx([1.0, 0.5, 0.0])


def test_gc_content_skips_tips_shorter_than_window():
    # Alignment requires equal-length rows -- pad the short one with gaps
    # so it strips down to a genuinely shorter raw sequence
    aln = pt.Alignment(["short", "long"], ["ACG---------", "ACGTACGTACGT"])
    res = pt.gc_content(aln, window=6)
    assert set(res["tip"]) == {"long"}


def test_gc_content_rejects_bad_window_or_step():
    aln = pt.Alignment(["s"], ["ACGTACGT"])
    with pytest.raises(ValueError, match="window must be"):
        pt.gc_content(aln, window=0)
    with pytest.raises(ValueError, match="step must be"):
        pt.gc_content(aln, window=4, step=0)


# --------------------------------------------------------------------------
# cpg_islands
# --------------------------------------------------------------------------
def test_cpg_islands_matches_a_hand_constructed_region():
    # a pure CG-repeat block, flanked by AT: windows straddling the
    # boundary can still clear the thresholds (diluted GC, still >= 0.5),
    # so the reported island legitimately extends a bit past the exact
    # planted block -- independently recompute gc/obs_exp_cpg over the
    # function's OWN reported span (plain Python, not the module's numpy
    # path) rather than hand-guessing a number for a span whose exact
    # extent depends on the windowing itself.
    seq = "AT" * 60 + "CG" * 60 + "AT" * 60
    aln = pt.Alignment(["s"], [seq])
    res = pt.cpg_islands(aln, window=40, gc_min=0.5, oe_min=0.6, min_length=40)
    assert len(res) == 1
    row = res.iloc[0]
    start, end = int(row["start"]), int(row["end"])
    # island must fully cover the planted CG block [120, 240)
    assert start <= 120 and end >= 240
    span = seq[start:end]
    g, c = span.count("G"), span.count("C")
    cpg = sum(1 for i in range(len(span) - 1) if span[i:i + 2] == "CG")
    expected_gc = (g + c) / len(span)
    expected_oe = cpg / (g * c / len(span)) if g * c > 0 else 0.0
    assert row["gc"] == pytest.approx(expected_gc)
    assert row["obs_exp_cpg"] == pytest.approx(expected_oe)
    assert row["obs_exp_cpg"] > 1.0   # genuinely CpG-enriched, not just GC-rich


def test_cpg_islands_type_i_error_near_nominal_on_suppressed_background():
    # the meaningful negative control -- see module docstring for why
    # plain i.i.d. random sequence is NOT a valid one for this function
    n_reps, hits = 30, 0
    for rep in range(n_reps):
        aln = pt.Alignment(["s"], [_cpg_suppressed_seq(2000, seed=1000 + rep)])
        hits += len(pt.cpg_islands(aln)) > 0
    assert hits / n_reps < 0.1, f"false-positive rate {hits}/{n_reps} too high"


def test_cpg_islands_detects_a_planted_island_in_suppressed_background():
    n_reps, hits = 15, 0
    for rep in range(n_reps):
        seed = 2000 + rep * 3
        before = _cpg_suppressed_seq(400, seed)
        island = "".join(np.random.default_rng(seed + 1).choice(["C", "G"], size=300))
        after = _cpg_suppressed_seq(400, seed + 2)
        aln = pt.Alignment(["s"], [before + island + after])
        res = pt.cpg_islands(aln)
        hits += any((r.start <= 400 and r.end >= 700) for r in res.itertuples())
    assert hits / n_reps > 0.8, f"only {hits}/{n_reps} detected"


def test_cpg_islands_rejects_bad_params():
    aln = pt.Alignment(["s"], ["ACGT" * 20])
    with pytest.raises(ValueError, match="window must be"):
        pt.cpg_islands(aln, window=1)
    with pytest.raises(ValueError, match="min_length must be"):
        pt.cpg_islands(aln, min_length=0)


# --------------------------------------------------------------------------
# tandem_repeats
# --------------------------------------------------------------------------
def test_tandem_repeats_matches_a_hand_constructed_example():
    # 'AT' x5 = 10bp period-2 repeat, flanked by non-repeating context
    aln = pt.Alignment(["s"], ["GCGC" + "AT" * 5 + "GCGC"])
    res = pt.tandem_repeats(aln, min_period=2, max_period=2, min_copies=3)
    assert len(res) == 1
    row = res.iloc[0]
    assert (row["start"], row["end"]) == (4, 14)
    assert row["period"] == 2
    assert row["unit"] == "AT"
    assert row["copies"] == pytest.approx(5.0)


def test_tandem_repeats_resolves_overlapping_periods_by_keeping_longest():
    # a run of identical bases is trivially "periodic" at every period;
    # the greedy resolution must keep the period-1 call (more copies for
    # the same span), not report the same region multiple times
    aln = pt.Alignment(["s"], ["GC" + "A" * 12 + "GC"])
    res = pt.tandem_repeats(aln, min_period=1, max_period=4, min_copies=3)
    assert len(res) == 1
    assert res.iloc[0]["period"] == 1


def test_tandem_repeats_mismatch_tolerance():
    # one mismatch inside an otherwise-perfect period-4 repeat: invisible
    # to exact matching (breaks the run in two, each too short to call at
    # min_copies=4), recoverable once a small mismatch budget is allowed
    unit_run = "ACGT" * 3 + "ACGA" + "ACGT" * 3   # one mismatch at position 12
    aln = pt.Alignment(["s"], ["TTTT" + unit_run + "TTTT"])
    exact = pt.tandem_repeats(aln, min_period=4, max_period=4, min_copies=4,
                              max_mismatch_frac=0.0)
    tolerant = pt.tandem_repeats(aln, min_period=4, max_period=4, min_copies=4,
                                 max_mismatch_frac=0.15)
    assert not any((r.start <= 4 and r.end >= 32) for r in exact.itertuples())
    assert any((r.start <= 4 and r.end >= 32) for r in tolerant.itertuples())


def test_tandem_repeats_long_chance_repeats_are_rare_at_stricter_settings():
    # NOT a claim that default settings have few hits: short repeats are
    # common by chance (see module docstring), and even min_copies=6
    # across 6 candidate periods and ~1000 starting positions turned out
    # empirically NOT rare either (13/20 replicates, a multiple-comparisons
    # effect -- many candidate (position, period) pairs are tried). At
    # min_copies=9 a chance repeat that long is genuinely rare (0/20,
    # verified directly, not assumed), which is the meaningful check.
    n_reps, hits = 20, 0
    for rep in range(n_reps):
        aln = pt.Alignment(["s"], [_random_seq(1000, seed=3000 + rep)])
        hits += len(pt.tandem_repeats(aln, min_copies=9)) > 0
    assert hits / n_reps < 0.25, f"false-positive rate {hits}/{n_reps} too high"


def test_tandem_repeats_rejects_bad_params():
    aln = pt.Alignment(["s"], ["ACGT" * 10])
    with pytest.raises(ValueError, match="min_period must be"):
        pt.tandem_repeats(aln, min_period=0)
    with pytest.raises(ValueError, match="max_period.*must be"):
        pt.tandem_repeats(aln, min_period=5, max_period=2)
    with pytest.raises(ValueError, match="min_copies must be"):
        pt.tandem_repeats(aln, min_copies=1)


# --------------------------------------------------------------------------
# shared input validation
# --------------------------------------------------------------------------
def test_rejects_non_nucleotide_input():
    aln = pt.Alignment(["s"], ["MKVLA" * 10])   # protein-looking
    with pytest.raises(ValueError, match="nucleotide"):
        pt.gc_content(aln, window=5)
