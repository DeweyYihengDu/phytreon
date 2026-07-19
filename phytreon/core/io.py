"""Tree input/output.

We reuse Biopython's robust parsers/writers for the standard formats
(Newick, Nexus, PhyloXML, NeXML) and convert to/from our :class:`Tree`.
A tiny self-contained Newick parser is also provided so that the package
has a zero-config fast path for the common case (and so the data model
can be used even if Biopython is unavailable).
"""
from __future__ import annotations

import io as _io
import re
from typing import Optional

from .tree import Node, Tree

_BIO_FORMATS = {"newick", "nexus", "phyloxml", "nexml", "cdao"}


# --------------------------------------------------------------------------
# Biopython bridge
# --------------------------------------------------------------------------
def from_biopython(bp_tree) -> Tree:
    """Convert a ``Bio.Phylo`` tree into a phytreon :class:`Tree`."""
    def convert(clade) -> Node:
        node = Node(
            name=clade.name,
            length=clade.branch_length,
            support=getattr(clade, "confidence", None),
            comment=getattr(clade, "comment", None),
        )
        for child in clade.clades:
            node.add_child(convert(child))
        return node

    root = convert(bp_tree.root)
    return Tree(root=root, name=getattr(bp_tree, "name", None))


def to_biopython(tree: Tree):
    from Bio.Phylo.BaseTree import Clade, Tree as BioTree

    def convert(node: Node) -> Clade:
        clade = Clade(
            branch_length=node.length,
            name=node.name,
            confidence=node.support,
        )
        clade.clades = [convert(c) for c in node.children]
        return clade

    return BioTree(root=convert(tree.root), name=tree.name)


def read(path: str, fmt: str = "newick") -> Tree:
    # annotated NEXUS keeps the per-node estimates a plain reader discards
    if fmt.lower() in ("beast", "mrbayes", "nexus-annotated"):
        from .nexus import read_annotated_nexus
        return read_annotated_nexus(path)
    if fmt not in _BIO_FORMATS:
        raise ValueError(f"unsupported format {fmt!r}; choose from {_BIO_FORMATS}")
    from Bio import Phylo
    return from_biopython(Phylo.read(path, fmt))


def write(tree: Tree, path: Optional[str] = None, fmt: str = "newick") -> Optional[str]:
    """Write to ``path``; if ``path`` is ``None`` return the string."""
    if fmt == "newick" and path is None:
        return to_newick(tree)
    from Bio import Phylo
    bp = to_biopython(tree)
    if path is None:
        buf = _io.StringIO()
        Phylo.write(bp, buf, fmt)
        return buf.getvalue()
    Phylo.write(bp, path, fmt)
    return None


# --------------------------------------------------------------------------
# Lightweight self-contained Newick parser / writer
# --------------------------------------------------------------------------
_TOKEN = re.compile(r"\s*([(),;])\s*|\s*([^(),;]+)\s*")
_NEEDS_QUOTE = re.compile(r"[()\[\]{}/\\,;:=*'\s]")


def _quote_label(label: str) -> str:
    """Newick-quote ``label`` if it contains reserved punctuation or
    whitespace (``()[]{}/\\,;:=*'`` or any space); an embedded single quote
    is escaped by doubling it, per the standard Newick quoting convention."""
    if not _NEEDS_QUOTE.search(label):
        return label
    return "'" + label.replace("'", "''") + "'"


def parse_newick(newick: str) -> Tree:
    """Parse a Newick string into a :class:`Tree` (no external deps).

    Supports branch lengths (``name:0.1``), internal labels / support
    values, and quoted names.  Square-bracket comments (e.g. NHX) are
    captured verbatim onto ``node.comment``.
    """
    s = newick.strip()
    if not s.endswith(";"):
        s += ";"

    pos = 0
    n = len(s)

    def parse_clade() -> Node:
        nonlocal pos
        node = Node()
        if s[pos] == "(":
            pos += 1  # consume '('
            while True:
                node.add_child(parse_clade())
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
        _parse_label(node)
        return node

    def _parse_label(node: Node) -> None:
        nonlocal pos
        # optional name/support
        if pos < n and s[pos] == "'":
            # quoted label: reserved characters inside are literal; an
            # embedded '' is an escaped literal single quote
            pos += 1
            buf = []
            while pos < n:
                if s[pos] == "'":
                    if pos + 1 < n and s[pos + 1] == "'":
                        buf.append("'")
                        pos += 2
                        continue
                    pos += 1
                    break
                buf.append(s[pos])
                pos += 1
            label = "".join(buf)
            while pos < n and s[pos] not in "():,;[":  # trailing whitespace
                pos += 1
        else:
            start = pos
            while pos < n and s[pos] not in "():,;[":
                pos += 1
            label = s[start:pos].strip().strip("'\"")
        if label:
            if node.is_leaf:
                node.name = label
            else:
                # internal label is usually a support value
                try:
                    node.support = float(label)
                except ValueError:
                    node.name = label
        _maybe_comment(node)
        # optional branch length
        if pos < n and s[pos] == ":":
            pos += 1
            start = pos
            while pos < n and s[pos] not in "():,;[":
                pos += 1
            try:
                node.length = float(s[start:pos])
            except ValueError:
                pass
        # a second comment may follow the branch length -- BEAST writes
        # per-branch rates there. Left unconsumed it derails the whole parse.
        _maybe_comment(node)

    def _maybe_comment(node: Node) -> None:
        """Consume a ``[...]`` block, appending it to ``node.comment``."""
        nonlocal pos
        if pos >= n or s[pos] != "[":
            return
        depth, cstart = 0, pos + 1
        while pos < n:
            if s[pos] == "[":
                depth += 1
            elif s[pos] == "]":
                depth -= 1
                if depth == 0:
                    text = s[cstart:pos]
                    node.comment = (f"{node.comment},{text}"
                                    if node.comment else text)
                    pos += 1
                    return
            pos += 1

    root = parse_clade()
    return Tree(root=root)


def to_newick(tree: Tree, with_support: bool = True) -> str:
    """Serialise a :class:`Tree` back to a Newick string."""
    def fmt(node: Node) -> str:
        if node.is_leaf:
            label = _quote_label(node.name) if node.name else ""
        else:
            inner = ",".join(fmt(c) for c in node.children)
            sup = ""
            if with_support and node.support is not None:
                sup = _num(node.support)
            elif node.name:
                sup = _quote_label(node.name)
            label = f"({inner}){sup}"
        if node.length is not None:
            label += f":{_num(node.length)}"
        return label

    return fmt(tree.root) + ";"


def _num(x: float) -> str:
    return f"{x:g}"
