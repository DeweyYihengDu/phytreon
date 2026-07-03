"""Discrete character/trait matrices -> :class:`~phytreon.infer.align.Alignment`.

A character matrix is a table with one row per taxon and one column per
discrete character -- e.g. a 0/1 gene presence/absence matrix, a coded
morphological character matrix, or any small per-column state space. This
recodes each column's states to single printable characters and packs them
into an :class:`Alignment`, so the result plugs directly into
:func:`phytreon.infer.parsimony.parsimony_tree` (or
``build_tree(..., method="parsimony")``).
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional, Union

from .align import Alignment

if TYPE_CHECKING:
    import pandas as pd

# codes assigned to states in order of first appearance: digits, then
# uppercase, then lowercase -- 62 distinct states per column before erroring.
_CODES = [chr(c) for c in range(48, 58)] + \
    [chr(c) for c in range(65, 91)] + \
    [chr(c) for c in range(97, 123)]


def read_character_matrix(source: Union[str, "pd.DataFrame"], *,
                          taxa_col: Optional[str] = None,
                          missing: Optional[object] = None) -> Alignment:
    """Build an :class:`Alignment` from a discrete character/trait matrix.

    ``source`` is a path to a CSV/TSV file, or an existing
    :class:`pandas.DataFrame` (taxa as rows, one column per character).
    Missing values (``NaN``, or an explicit ``missing`` sentinel) are
    encoded as ambiguous states, so they never force a false character
    change under parsimony scoring.

    ``taxa_col`` names the column holding taxon labels; if omitted, the
    DataFrame's existing index (or the file's first column) is used.

    Each column may use any small set of hashable states (numbers,
    strings, booleans, ...); states are recoded internally to single
    printable characters (up to 62 distinct states per column).
    """
    import pandas as pd

    if isinstance(source, pd.DataFrame):
        df = source.copy()
    else:
        sep = "\t" if str(source).lower().endswith((".tsv", ".tab")) else ","
        df = pd.read_csv(source, sep=sep, index_col=0 if taxa_col is None else None)
    if taxa_col is not None:
        df = df.set_index(taxa_col)
    names = [str(n) for n in df.index]

    rows = [[] for _ in names]
    for col in df.columns:
        values = df[col]
        is_missing = values.isna()
        if missing is not None:
            is_missing = is_missing | (values == missing)
        states = sorted({v for v, m in zip(values, is_missing) if not m}, key=str)
        if len(states) > len(_CODES):
            raise ValueError(
                f"column {col!r} has {len(states)} distinct states; "
                f"read_character_matrix supports at most {len(_CODES)} per column")
        code = {state: _CODES[i] for i, state in enumerate(states)}
        for row, value, m in zip(rows, values, is_missing):
            row.append("?" if m else code[value])

    return Alignment(names, ["".join(row) for row in rows])
