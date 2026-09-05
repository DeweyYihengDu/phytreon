"""Per-tip DNA sequence features: GC content, CpG islands, tandem repeats.

Every function here reads each tip's own raw (ungapped) sequence out of an
alignment and scans it for a feature that lives *within* one sequence, not
between sequences the way the rest of :mod:`phytreon.infer` compares trees
or alignments to each other.

**What this is not: no TSS prediction.** Transcription-start-site calling
from raw sequence alone, with no other evidence (CAGE-seq, RNA-seq, an
annotated gene model), is not reliable enough to ship as a phytreon
algorithm -- promoter signals like the TATA box are absent from a large
fraction of real promoters, and a purely compositional heuristic would be
guessing, not measuring. This mirrors :mod:`phytreon.infer.recombination`'s
own choice to decline the published PHI test's "refined incompatibility"
statistic after failing to verify it with confidence, rather than ship
something that looks like an established method but isn't: nothing here
predicts a TSS. If you have TSS coordinates from other evidence, hand them
straight to :meth:`~phytreon.plot.figure.TreeFigure.sequence_features` as
zero-length intervals -- the plotting side draws whatever you give it.

**tandem_repeats' scope, stated plainly.** This finds exact (or, with
``max_mismatch_frac``, near-exact) periodic runs by direct sequence
comparison. It is not a reimplementation of Tandem Repeats Finder's
statistical scoring model (indel-aware alignment, a fitted probability of
each kind of change) -- just the simpler, honestly-scoped version of "does
this region repeat".

``cpg_islands`` uses the standard Gardiner-Garden & Frommer (1987,
*J. Mol. Biol.* 196(2):261-282) criteria, with the same defaults EMBOSS's
``newcpgreport`` uses.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

from .align import guess_type
from .matrix import Alignment


def _ungapped_nucleotide_seqs(aln: Alignment) -> Dict[str, str]:
    """Each tip's own raw sequence, gaps stripped, upper-cased.

    ``Alignment`` has no ungapped accessor of its own (only ``names``,
    ``seqs``, ``ncol``, ``nseq``, ``column``, ``select_columns``,
    ``records``), so every function in this module strips gaps here, once.
    """
    if guess_type(aln.seqs) != "nucleotide":
        raise ValueError(
            "seqfeatures: expected nucleotide sequences, but this alignment "
            "does not look like DNA/RNA"
        )
    return {name: seq.upper().replace("-", "").replace(".", "")
            for name, seq in zip(aln.names, aln.seqs)}


def sequence_lengths(aln: Alignment) -> "pd.Series":  # noqa: F821
    """Each tip's own ungapped sequence length.

    Exists so callers of :func:`cpg_islands`/:func:`tandem_repeats` --
    whose output can be *empty* for a tip with no called features -- have a
    correct length to pass to :meth:`~phytreon.plot.figure.TreeFigure.
    sequence_features`'s ``lengths=``. Inferring length from "the largest
    ``end`` seen in the results" would silently draw a zero-length track
    for a tip with nothing detected, which is wrong: a blank strip of the
    *right* length is the correct picture there, not an omission.
    """
    import pandas as pd
    seqs = _ungapped_nucleotide_seqs(aln)
    return pd.Series({name: len(seq) for name, seq in seqs.items()},
                     name="length")


def gc_content(aln: Alignment, window: int, step: Optional[int] = None
              ) -> "pd.DataFrame":  # noqa: F821
    """Sliding-window GC fraction along each tip's own sequence.

    ``step`` (default: ``window``, i.e. non-overlapping tiling windows)
    controls how far the window moves each step; smaller than ``window``
    gives overlapping windows. A tip shorter than ``window`` contributes no
    rows (there is no full window to report) rather than raising, since
    real alignments routinely mix sequence lengths.

    Returns a DataFrame with columns ``tip, start, end, gc`` (``start``/
    ``end`` are 0-based half-open coordinates into that tip's own ungapped
    sequence, ``gc`` in ``[0, 1]``).
    """
    import numpy as np
    import pandas as pd
    if window < 1:
        raise ValueError(f"gc_content: window must be >= 1, got {window}")
    step = window if step is None else step
    if step < 1:
        raise ValueError(f"gc_content: step must be >= 1, got {step}")

    rows = []
    for tip, seq in _ungapped_nucleotide_seqs(aln).items():
        n = len(seq)
        if n < window:
            continue
        is_gc = np.fromiter((ch in "GC" for ch in seq), dtype=np.int64, count=n)
        cum = np.concatenate(([0], np.cumsum(is_gc)))
        for start in range(0, n - window + 1, step):
            end = start + window
            rows.append({"tip": tip, "start": start, "end": end,
                        "gc": float(cum[end] - cum[start]) / window})
    return pd.DataFrame(rows, columns=["tip", "start", "end", "gc"])


def cpg_islands(aln: Alignment, window: int = 200, gc_min: float = 0.5,
                oe_min: float = 0.6, min_length: int = 200
               ) -> "pd.DataFrame":  # noqa: F821
    """CpG islands by the Gardiner-Garden & Frommer (1987) criteria.

    A window of length ``window`` is flagged when its GC fraction is
    ``>= gc_min`` *and* its observed/expected CpG ratio -- ``n_CpG /
    (n_C * n_G / window)``, the standard formula -- is ``>= oe_min``.
    Flagged windows are scanned at 1-bp resolution and merged wherever they
    overlap (this is a correctness requirement for getting island
    boundaries right, not a tunable parameter, so it is not exposed).
    Merged regions shorter than ``min_length`` are dropped. ``gc``/
    ``obs_exp_cpg`` in the result are recomputed over each island's final,
    merged span, not just the window that first flagged it.

    **What "expected" means here, stated plainly.** ``obs_exp_cpg``
    compares observed CpG count against the count expected under pure
    statistical independence of G and C -- not against real genomic
    background. In actual vertebrate genomes CpG is heavily suppressed by
    methylation-driven mutation (bulk sequence typically sits around
    0.2-0.3 obs/exp), so a real CpG island's ~1.0 stands out sharply
    against it; on already CpG-neutral input -- uniform-random test
    sequence being the clearest case, verified directly: i.i.d. random
    ACGT hits both thresholds essentially everywhere, a 30/30
    false-positive rate on one such check -- there is no suppressed
    baseline to stand out against, and large stretches can be flagged.
    This function reports the same ratio real CpG-island callers use; it
    does not (and cannot, without a wider reference) correct for whether a
    given input already has realistic methylation-driven suppression.

    Returns a DataFrame with columns ``tip, start, end, gc, obs_exp_cpg``.
    """
    import numpy as np
    import pandas as pd
    if window < 2:
        raise ValueError(f"cpg_islands: window must be >= 2, got {window}")
    if min_length < 1:
        raise ValueError(f"cpg_islands: min_length must be >= 1, got {min_length}")

    rows = []
    for tip, seq in _ungapped_nucleotide_seqs(aln).items():
        n = len(seq)
        if n < window:
            continue
        g = np.fromiter((ch == "G" for ch in seq), dtype=np.int64, count=n)
        c = np.fromiter((ch == "C" for ch in seq), dtype=np.int64, count=n)
        cpg = np.fromiter((seq[i] == "C" and seq[i + 1] == "G" for i in range(n - 1)),
                          dtype=np.int64, count=n - 1)
        cum_g = np.concatenate(([0], np.cumsum(g)))
        cum_c = np.concatenate(([0], np.cumsum(c)))
        cum_cpg = np.concatenate(([0], np.cumsum(cpg)))

        def window_stats(start: int, end: int) -> Tuple[float, float]:
            length = end - start
            g_n = float(cum_g[end] - cum_g[start])
            c_n = float(cum_c[end] - cum_c[start])
            cpg_n = float(cum_cpg[end - 1] - cum_cpg[start])
            gc = (g_n + c_n) / length
            denom = g_n * c_n / length
            oe = cpg_n / denom if denom > 0 else 0.0
            return gc, oe

        n_windows = n - window + 1
        flagged = np.zeros(n_windows, dtype=bool)
        for start in range(n_windows):
            gc, oe = window_stats(start, start + window)
            flagged[start] = gc >= gc_min and oe >= oe_min

        starts = np.flatnonzero(flagged)
        if starts.size == 0:
            continue
        merged: List[Tuple[int, int]] = []
        cur_s, cur_e = int(starts[0]), int(starts[0]) + window
        for s in starts[1:]:
            s, e = int(s), int(s) + window
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                merged.append((cur_s, cur_e))
                cur_s, cur_e = s, e
        merged.append((cur_s, cur_e))

        for s, e in merged:
            if e - s < min_length:
                continue
            gc, oe = window_stats(s, e)
            rows.append({"tip": tip, "start": s, "end": e,
                        "gc": gc, "obs_exp_cpg": oe})
    return pd.DataFrame(rows, columns=["tip", "start", "end", "gc", "obs_exp_cpg"])


def _periodic_runs(match: "np.ndarray", max_mismatch_frac: float  # noqa: F821
                   ) -> List[Tuple[int, int]]:
    """Maximal runs of ``match`` (bool array) tolerating up to
    ``max_mismatch_frac`` of the run being ``False``.

    Greedy left-to-right: skip non-matches, then extend a candidate run as
    far as possible while its cumulative mismatch fraction stays within
    budget, close it, and continue from right after. At
    ``max_mismatch_frac=0`` this reduces to maximal runs of all-``True``
    values, since a single mismatch immediately exceeds a zero budget.
    """
    runs = []
    i, n = 0, len(match)
    while i < n:
        if not match[i]:
            i += 1
            continue
        j, mismatches, last_good = i, 0, i
        while j < n:
            if not match[j]:
                mismatches += 1
            if mismatches / (j - i + 1) > max_mismatch_frac:
                break
            last_good = j
            j += 1
        runs.append((i, last_good))
        i = last_good + 1
    return runs


def tandem_repeats(aln: Alignment, min_period: int = 1, max_period: int = 6,
                   min_copies: int = 3, max_mismatch_frac: float = 0.0
                  ) -> "pd.DataFrame":  # noqa: F821
    """Tandem repeats by direct periodic-match scanning, one tip at a time.

    For each candidate period ``p`` in ``[min_period, max_period]``, finds
    maximal runs where position ``i`` matches position ``i + p`` (see
    :func:`_periodic_runs` for how ``max_mismatch_frac`` tolerance works),
    keeps runs spanning at least ``min_copies`` full repeats of the unit,
    then resolves calls that overlap across different periods by keeping
    the longest span (ties broken by more copies) and discarding anything
    that overlaps an already-kept call for that tip -- the greedy
    interval-scheduling rule that stops e.g. a run of ``AAAA...`` (a
    trivial period-1 repeat) from also being reported redundantly as
    period 2, 3, ....

    Short repeats at the default settings are common by chance, not a
    calibration problem -- real genomes are full of them too. Even a
    fairly demanding ``min_copies`` does not make chance hits rare on its
    own: scanning every period in ``[min_period, max_period]`` at every
    position is a multiple-comparisons search, so a false positive at
    *some* period *somewhere* stays common well past where a single-period
    check alone would predict (checked directly: at ``min_copies=6``,
    ``max_period=6``, 1000bp random sequence still hits 13/20 replicates;
    ``min_copies=9`` is needed for that to drop below 1/20). Treat a single
    call as suggestive, not significant on its own, unless ``min_copies``
    is set well above what a single period's chance rate alone would
    suggest.

    Returns a DataFrame with columns ``tip, start, end, period, unit, copies``.
    """
    import numpy as np
    import pandas as pd
    if min_period < 1:
        raise ValueError(f"tandem_repeats: min_period must be >= 1, got {min_period}")
    if max_period < min_period:
        raise ValueError(
            f"tandem_repeats: max_period ({max_period}) must be >= "
            f"min_period ({min_period})"
        )
    if min_copies < 2:
        raise ValueError(f"tandem_repeats: min_copies must be >= 2, got {min_copies}")

    rows = []
    for tip, seq in _ungapped_nucleotide_seqs(aln).items():
        n = len(seq)
        candidates = []
        for p in range(min_period, max_period + 1):
            if n - p < 1:
                continue
            match = np.fromiter((seq[i] == seq[i + p] for i in range(n - p)),
                                dtype=bool, count=n - p)
            for i, last in _periodic_runs(match, max_mismatch_frac):
                start, end = i, last + p + 1
                span = end - start
                copies = span / p
                if copies >= min_copies:
                    candidates.append((start, end, p, seq[start:start + p], copies))
        # greedy: longest span first, then more copies, keep non-overlapping
        candidates.sort(key=lambda r: (-(r[1] - r[0]), -r[4]))
        accepted: List[Tuple[int, int]] = []
        for start, end, p, unit, copies in candidates:
            if any(start < a_end and end > a_start for a_start, a_end in accepted):
                continue
            accepted.append((start, end))
            rows.append({"tip": tip, "start": start, "end": end, "period": p,
                        "unit": unit, "copies": round(copies, 2)})
    return pd.DataFrame(rows, columns=["tip", "start", "end", "period", "unit", "copies"])
