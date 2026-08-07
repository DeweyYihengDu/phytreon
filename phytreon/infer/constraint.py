"""Taxonomy-derived constraint trees, for a constrained ML search.

A constraint tree tells an external ML engine (IQ-TREE's ``-g``, RAxML-NG's
``--tree-constraint``) which clades it may not break -- here, the clades a
taxonomy column claims. The engine still searches everything else freely: it
resolves relationships between groups, within a group, and where to put any
tip the constraint says nothing about, under the actual likelihood. This is
the constraint asking the data to agree wherever it is free to; for the
stronger, structural version that never lets the data disagree, see
:func:`~phytreon.infer.distance.constrained_nj`.
"""
from __future__ import annotations

from typing import Dict, List, Optional

from ..core.tree import Node, Tree


def _isnan(v) -> bool:
    try:
        return v != v          # only NaN is unequal to itself
    except Exception:
        return False


def constraint_tree(groups, column: Optional[str] = None) -> Tree:
    """Build a multifurcating constraint tree from a taxonomy grouping.

    ``groups`` is a ``{tip_name: label}`` mapping (a genus column, say), or a
    table with a ``name`` column (or index) plus ``column=`` naming which
    column to read the label from.

    Every label with two or more tips becomes one polytomy: unresolved
    internally, so nothing about *how* those tips relate is dictated, only
    that a constrained search may not place another tip between them. A
    label with a single tip needs no polytomy of its own. A tip mapped to
    ``None``/NaN, or missing from ``groups`` altogether, is left out of the
    file entirely -- which is how IQ-TREE and RAxML-NG already read a
    constraint that does not name every taxon: that tip is placed freely
    rather than forced into whichever group happens to be nearby.

    The result carries no branch lengths (a constraint is topology only) and
    is written with :meth:`~phytreon.core.tree.Tree.write`, same as any other
    tree here.
    """
    if hasattr(groups, "columns"):
        if column is None:
            raise TypeError("column= is required when groups is a table, "
                            "not a {tip_name: label} mapping")
        groups = (dict(zip(groups["name"], groups[column]))
                 if "name" in groups.columns else groups[column].to_dict())

    by_label: Dict[object, List[str]] = {}
    for name, label in groups.items():
        if label is None or _isnan(label):
            continue
        by_label.setdefault(label, []).append(str(name))
    if not by_label:
        raise ValueError("no tip has a non-missing group label -- nothing "
                         "to constrain")

    root = Node()
    for label in sorted(by_label, key=str):
        tips = by_label[label]
        if len(tips) == 1:
            root.add_child(Node(name=tips[0]))
        else:
            clade = root.add_child(Node())
            for name in tips:
                clade.add_child(Node(name=name))
    return Tree(root=root)
