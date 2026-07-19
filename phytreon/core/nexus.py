"""Annotated NEXUS trees -- BEAST, MrBayes and friends.

Bayesian phylogenetics programs write their per-node estimates into NEXUS
comments::

    TREE t1 = [&R] ((1:1.0,2:1.0)[&height=1.0,height_95%_HPD={0.8,1.4},
                                  posterior=0.98]:2.0,3:3.0);

A plain NEXUS reader keeps the topology and throws all of that away, which is
exactly the part worth plotting -- node ages and their credible intervals,
posterior clade probabilities, per-branch rates.  This module keeps it, on
``node.data``::

    tree = pt.Tree.read("beast.tre", fmt="beast")
    tree.root.data["posterior"]            # 1.0
    tree.root.data["height_95_lower"]      # 2.1

``{lo, hi}`` intervals also get flattened to ``<name>_lower`` / ``<name>_upper``
so a 95% HPD lands on the keys
:meth:`~phytreon.plot.figure.TreeFigure.node_bars` reads by default.
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional

from .tree import Tree


# --------------------------------------------------------------------------
# the "&key=value,..." comment payload
# --------------------------------------------------------------------------
def _split_top_level(text: str) -> List[str]:
    """Split on commas that are not inside ``{}`` or quotes."""
    parts, buf, depth, quote = [], [], 0, ""
    for ch in text:
        if quote:
            buf.append(ch)
            if ch == quote:
                quote = ""
            continue
        if ch in "\"'":
            quote = ch
            buf.append(ch)
        elif ch in "{[":
            depth += 1
            buf.append(ch)
        elif ch in "}]":
            depth -= 1
            buf.append(ch)
        elif ch == "," and depth == 0:
            parts.append("".join(buf))
            buf = []
        else:
            buf.append(ch)
    if buf:
        parts.append("".join(buf))
    return [p.strip() for p in parts if p.strip()]


def _coerce(raw: str):
    """Turn one annotation value into a number, string or list."""
    raw = raw.strip()
    if raw.startswith("{") and raw.endswith("}"):
        return [_coerce(v) for v in _split_top_level(raw[1:-1])]
    if len(raw) > 1 and raw[0] == raw[-1] and raw[0] in "\"'":
        return raw[1:-1]
    try:
        return float(raw)
    except ValueError:
        return raw


#: ``height_95%_HPD`` -> ``height_95``, so the flattened bounds read
#: ``height_95_lower`` / ``height_95_upper``
_HPD_SUFFIX = re.compile(r"_?(HPD|hpd)$")


def _interval_base(key: str) -> str:
    return _HPD_SUFFIX.sub("", key.replace("%", "")).rstrip("_")


def parse_annotation(comment: Optional[str]) -> Dict[str, object]:
    """Parse a ``&key=value,...`` NEXUS node comment into a dict.

    Two-element numeric values (BEAST's ``{lower, upper}`` intervals) are also
    flattened to ``<base>_lower`` / ``<base>_upper``; the raw list is kept
    under the original key.
    """
    if not comment:
        return {}
    text = comment.strip()
    if text.startswith("&"):
        text = text[1:]
    out: Dict[str, object] = {}
    for item in _split_top_level(text):
        if "=" not in item:
            continue                      # a bare flag such as [&R]
        key, raw = item.split("=", 1)
        # "!" prefixes BEAST's styling keys; "&" opens each comment block, and
        # a node carrying two blocks (one before the branch length, one after)
        # has them joined, so the second block's "&" lands mid-string
        key = key.strip().lstrip("!&").strip()
        value = _coerce(raw)
        out[key] = value
        if (isinstance(value, list) and len(value) == 2
                and all(isinstance(v, float) for v in value)):
            base = _interval_base(key)
            lo, hi = sorted(value)
            out[f"{base}_lower"] = lo
            out[f"{base}_upper"] = hi
    return out


def annotate_from_comments(tree: Tree) -> Tree:
    """Parse every node's ``comment`` into ``node.data`` (in place)."""
    for node in tree.traverse():
        for key, value in parse_annotation(node.comment).items():
            node.data.setdefault(key, value)
    return tree


# --------------------------------------------------------------------------
# the NEXUS wrapper
# --------------------------------------------------------------------------
_TREE_LINE = re.compile(
    r"^\s*TREE\s+.*?=\s*(?:\[[^\]]*\]\s*)?(?P<newick>.+;)\s*$",
    re.IGNORECASE | re.MULTILINE)


def _translate_table(block: str) -> Dict[str, str]:
    m = re.search(r"TRANSLATE(?P<body>.*?);", block,
                  re.IGNORECASE | re.DOTALL)
    if not m:
        return {}
    table: Dict[str, str] = {}
    for entry in m.group("body").split(","):
        parts = entry.split()
        if len(parts) >= 2:
            table[parts[0]] = " ".join(parts[1:]).strip("'\"")
    return table


def read_annotated_nexus(path: str, *, tree_index: int = 0) -> Tree:
    """Read a NEXUS tree written by BEAST / MrBayes, keeping node annotations.

    ``tree_index`` picks which tree to take from a file holding several (a
    posterior sample); the default is the first, which for a TreeAnnotator
    summary file is the summary tree.
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    if "#NEXUS" not in text.upper():
        raise ValueError(f"{path!r} does not look like a NEXUS file")

    matches = _TREE_LINE.findall(text)
    if not matches:
        raise ValueError(f"no TREE statement found in {path!r}")
    if tree_index >= len(matches):
        raise IndexError(
            f"tree_index {tree_index} out of range: {path!r} holds "
            f"{len(matches)} tree(s)")

    from .io import parse_newick
    tree = annotate_from_comments(parse_newick(matches[tree_index]))

    table = _translate_table(text)
    if table:
        for leaf in tree.leaves():
            if leaf.name in table:
                leaf.name = table[leaf.name]
    return tree


def read_beast(path: str, **kwargs) -> Tree:
    """Alias for :func:`read_annotated_nexus` -- BEAST and MrBayes write the
    same annotated-NEXUS dialect."""
    return read_annotated_nexus(path, **kwargs)
