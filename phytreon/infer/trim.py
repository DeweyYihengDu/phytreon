"""Alignment trimming -- the "cut" step before tree building.

Removes poorly aligned / gap-rich columns before tree inference.
Every criterion is an independent, tunable knob; a column is kept only if it
passes *all* the active ones:

* ``max_gap``         -- drop columns with gap fraction above this.
* ``min_occupancy``   -- keep columns with at least this non-gap fraction.
* ``min_conservation``-- keep columns where the commonest residue reaches
  this fraction (of non-gap residues).

Returns the trimmed :class:`~phytreon.infer.align.Alignment`; pass
``return_mask=True`` to also get the boolean keep-mask for traceability.
"""
from __future__ import annotations

from collections import Counter
from typing import List, Optional, Tuple, Union

from .align import Alignment


def column_gap_fraction(aln: Alignment, j: int) -> float:
    col = aln.column(j)
    return sum(c == "-" for c in col) / len(col)


def column_conservation(aln: Alignment, j: int) -> float:
    col = [c for c in aln.column(j) if c != "-"]
    if not col:
        return 0.0
    return Counter(col).most_common(1)[0][1] / len(col)


def trim(aln: Alignment, *, max_gap: float = 0.5,
         min_occupancy: Optional[float] = None,
         min_conservation: Optional[float] = None,
         min_length: int = 1,
         return_mask: bool = False
         ) -> Union[Alignment, Tuple[Alignment, List[bool]]]:
    """Trim alignment columns by the active criteria. See module docstring."""
    keep_idx: List[int] = []
    mask: List[bool] = []
    for j in range(aln.ncol):
        gapf = column_gap_fraction(aln, j)
        ok = gapf <= max_gap
        if ok and min_occupancy is not None:
            ok = (1.0 - gapf) >= min_occupancy
        if ok and min_conservation is not None:
            ok = column_conservation(aln, j) >= min_conservation
        mask.append(ok)
        if ok:
            keep_idx.append(j)

    if len(keep_idx) < min_length:
        raise ValueError(
            f"trimming kept only {len(keep_idx)} columns (< min_length="
            f"{min_length}); relax the thresholds."
        )

    trimmed = aln.select_columns(keep_idx)
    return (trimmed, mask) if return_mask else trimmed


def trim_terminal_gaps(aln: Alignment, max_gap: float = 0.9) -> Alignment:
    """Trim only the leading/trailing gap-rich columns, keep the interior."""
    keep = [j for j in range(aln.ncol) if column_gap_fraction(aln, j) <= max_gap]
    if not keep:
        return aln.select_columns([])
    lo, hi = keep[0], keep[-1]
    return aln.select_columns(range(lo, hi + 1))
