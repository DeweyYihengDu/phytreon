"""Pure-Python multiple sequence alignment.

A self-contained **progressive aligner** (no external program required):

1. cheap k-mer distances between the raw sequences,
2. a UPGMA *guide tree* from those distances,
3. profile-profile Needleman-Wunsch merging up the guide tree.

Scoring is fully configurable (match/mismatch/gap for nucleotides, or any
Biopython substitution matrix for proteins).  This is meant for small/medium
inputs; for large alignments plug in MAFFT/MUSCLE via :func:`align_external`.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

Record = Tuple[str, str]


# --------------------------------------------------------------------------
# alignment container
# --------------------------------------------------------------------------
@dataclass
class Alignment:
    names: List[str]
    seqs: List[str]                       # aligned, all equal length

    @property
    def nseq(self) -> int:
        return len(self.names)

    @property
    def ncol(self) -> int:
        return len(self.seqs[0]) if self.seqs else 0

    def column(self, j: int) -> List[str]:
        return [s[j] for s in self.seqs]

    def select_columns(self, idx: Sequence[int]) -> "Alignment":
        keep = list(idx)
        return Alignment(list(self.names),
                         ["".join(s[i] for i in keep) for s in self.seqs])

    def records(self) -> List[Record]:
        return list(zip(self.names, self.seqs))

    def to_fasta(self, path: Optional[str] = None, width: int = 60) -> Optional[str]:
        out = []
        for n, s in zip(self.names, self.seqs):
            out.append(f">{n}")
            out += [s[i:i + width] for i in range(0, len(s), width)] or [""]
        text = "\n".join(out) + "\n"
        if path:
            with open(path, "w") as f:
                f.write(text)
            return None
        return text

    @classmethod
    def from_fasta(cls, source: str) -> "Alignment":
        recs = read_fasta(source)
        return cls([n for n, _ in recs], [s for _, s in recs])

    def __repr__(self) -> str:
        return f"<Alignment nseq={self.nseq} ncol={self.ncol}>"


def read_fasta(source: str) -> List[Record]:
    """Read FASTA from a path or a raw string."""
    import os
    text = open(source).read() if os.path.exists(source) else source
    recs, name, buf = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(buf)))
            name, buf = line[1:].strip().split()[0] if line[1:].strip() else line[1:], []
        elif line.strip():
            buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


# --------------------------------------------------------------------------
# scoring
# --------------------------------------------------------------------------
class Scoring:
    def __init__(self, match=1.0, mismatch=-1.0, gap=-2.0, matrix=None):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self._matrix = matrix              # dict-like (a,b)->score, or None

    def sub(self, a: str, b: str) -> float:
        if self._matrix is not None:
            try:
                return float(self._matrix[a, b])
            except (KeyError, IndexError):
                return self.mismatch
        return self.match if a == b else self.mismatch

    @classmethod
    def for_type(cls, seqtype: str, **kw):
        if seqtype == "protein":
            from Bio.Align import substitution_matrices
            m = substitution_matrices.load("BLOSUM62")
            kw.setdefault("gap", -4.0)
            return cls(matrix=m, **kw)
        return cls(**kw)


def guess_type(seqs: Sequence[str]) -> str:
    sample = "".join(seqs)[:2000].upper()
    nuc = sum(c in "ACGTUN-" for c in sample)
    return "nucleotide" if sample and nuc / len(sample) > 0.9 else "protein"


# --------------------------------------------------------------------------
# progressive aligner
# --------------------------------------------------------------------------
class _Profile:
    """Aligned block: ordered names + columns (each column a list of chars)."""
    __slots__ = ("names", "cols")

    def __init__(self, names: List[str], cols: List[List[str]]):
        self.names = names
        self.cols = cols

    @classmethod
    def single(cls, name: str, seq: str) -> "_Profile":
        return cls([name], [[c] for c in seq])

    def to_records(self) -> List[Record]:
        seqs = ["".join(col[i] for col in self.cols) for i in range(len(self.names))]
        return list(zip(self.names, seqs))


def _col_score(ca: List[str], cb: List[str], sc: Scoring) -> float:
    tot, n = 0.0, 0
    for a in ca:
        if a == "-":
            continue
        for b in cb:
            if b == "-":
                continue
            tot += sc.sub(a, b)
            n += 1
    return tot / n if n else 0.0


def _align_profiles(P: _Profile, Q: _Profile, sc: Scoring) -> _Profile:
    m, n = len(P.cols), len(Q.cols)
    g = sc.gap
    # column scores are cached by content: conserved alignments have many
    # identical columns, so this collapses most of the O(m*n) score work.
    Pk = [tuple(c) for c in P.cols]
    Qk = [tuple(c) for c in Q.cols]
    cache: dict = {}

    def cs(a, b):
        v = cache.get((a, b))
        if v is None:
            v = _col_score(a, b, sc)
            cache[(a, b)] = v
        return v

    F = [[0.0] * (n + 1) for _ in range(m + 1)]
    T = [[0] * (n + 1) for _ in range(m + 1)]      # 0 diag, 1 up, 2 left
    for i in range(1, m + 1):
        F[i][0] = i * g
        T[i][0] = 1
    for j in range(1, n + 1):
        F[0][j] = j * g
        T[0][j] = 2
    for i in range(1, m + 1):
        Fi, Fim, Ti, Pki = F[i], F[i - 1], T[i], Pk[i - 1]
        for j in range(1, n + 1):
            d = Fim[j - 1] + cs(Pki, Qk[j - 1])
            u = Fim[j] + g
            l = Fi[j - 1] + g
            if d >= u and d >= l:
                Fi[j] = d
            elif u >= l:
                Fi[j] = u; Ti[j] = 1
            else:
                Fi[j] = l; Ti[j] = 2

    gapsP = ["-"] * len(P.names)
    gapsQ = ["-"] * len(Q.names)
    cols: List[List[str]] = []
    i, j = m, n
    while i > 0 or j > 0:
        t = T[i][j]
        if t == 0:
            cols.append(P.cols[i - 1] + Q.cols[j - 1]); i -= 1; j -= 1
        elif t == 1:
            cols.append(P.cols[i - 1] + gapsQ); i -= 1
        else:
            cols.append(gapsP + Q.cols[j - 1]); j -= 1
    cols.reverse()
    return _Profile(P.names + Q.names, cols)


def _kmer_counts(seq: str, k: int) -> Counter:
    s = seq.replace("-", "").upper()
    return Counter(s[i:i + k] for i in range(len(s) - k + 1))


def kmer_distance_matrix(names: Sequence[str], seqs: Sequence[str], k: int = 3):
    counts = [_kmer_counts(s, k) for s in seqs]
    n = len(seqs)
    D = [[0.0] * n for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            ca, cb = counts[a], counts[b]
            shared = sum((ca & cb).values())
            denom = min(sum(ca.values()), sum(cb.values())) or 1
            d = 1.0 - shared / denom
            D[a][b] = D[b][a] = d
    return D


def align(records: List[Record], seqtype: str = "auto", k: int = 3,
          **scoring_kw) -> Alignment:
    """Progressively align unaligned ``records`` -> :class:`Alignment`.

    ``seqtype`` is ``"auto"`` | ``"nucleotide"`` | ``"protein"``.  Extra
    keywords (``match``, ``mismatch``, ``gap``, ``matrix``) tune scoring.
    """
    records = [(n, s.replace("-", "")) for n, s in records]
    names = [n for n, _ in records]
    seqs = [s for _, s in records]
    if len(records) == 0:
        return Alignment([], [])
    if len(records) == 1:
        return Alignment(names, seqs)

    if seqtype == "auto":
        seqtype = guess_type(seqs)
    sc = Scoring.for_type(seqtype, **scoring_kw)
    kk = max(1, min(k, min(len(s) for s in seqs)))

    # guide tree via UPGMA on k-mer distances
    from .distance import upgma
    D = kmer_distance_matrix(names, seqs, k=kk)
    guide = upgma(names, D)

    seqmap = dict(records)
    prof: Dict[int, _Profile] = {}
    for node in guide.traverse("postorder"):
        if node.is_leaf:
            prof[id(node)] = _Profile.single(node.name, seqmap[node.name])
        else:
            acc = prof[id(node.children[0])]
            for child in node.children[1:]:
                acc = _align_profiles(acc, prof[id(child)], sc)
            prof[id(node)] = acc

    merged = prof[id(guide.root)]
    recs = dict(merged.to_records())
    # restore the caller's original order
    return Alignment(names, [recs[n] for n in names])


# --------------------------------------------------------------------------
# external aligner adapter (optional, for large inputs)
# --------------------------------------------------------------------------
def align_external(records: List[Record], tool: str = "mafft",
                   path: Optional[str] = None, extra_args: Optional[List[str]] = None
                   ) -> Alignment:
    """Align via MAFFT or MUSCLE if installed.  Raises if the tool is absent."""
    import shutil
    import subprocess
    import tempfile
    import os

    exe = path or shutil.which(tool)
    if exe is None:
        raise RuntimeError(
            f"{tool!r} not found on PATH. Install it, pass path=, or use the "
            f"built-in align()."
        )
    with tempfile.TemporaryDirectory() as tmp:
        infile = os.path.join(tmp, "in.fasta")
        Alignment([n for n, _ in records],
                  [s.replace("-", "") for _, s in records]).to_fasta(infile)
        if tool == "mafft":
            cmd = [exe, *(extra_args or ["--auto"]), infile]
            out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
            return Alignment.from_fasta(out)
        if tool == "muscle":
            outfile = os.path.join(tmp, "out.fasta")
            subprocess.run([exe, "-align", infile, "-output", outfile],
                           capture_output=True, text=True, check=True)
            return Alignment.from_fasta(outfile)
    raise ValueError(f"unsupported tool {tool!r}")
