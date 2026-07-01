"""Rectangular ("phylogram") layout -- the default."""
from __future__ import annotations

from typing import List, Tuple

from ..core.tree import Node
from .base import Layout

XY = Tuple[float, float]


class RectangularLayout(Layout):
    """x = distance from root, y = leaf row.

    Branches are drawn as elbows: a horizontal segment from the parent's
    x to the node's x at the node's row, plus a vertical connector at the
    parent spanning its children's rows.
    """

    is_polar = False

    def compute(self) -> None:
        for node in self.tree.traverse():
            node.x = node.data["_rootd"]
            node.y = node.data["_row"]

    def branch_path(self, node: Node) -> List[XY]:
        if node.is_root:
            return []
        p = node.parent
        return [(p.x, node.y), (node.x, node.y)]

    def child_connector(self, node: Node) -> List[XY]:
        if node.is_leaf:
            return []
        ys = [c.y for c in node.children]
        return [(node.x, min(ys)), (node.x, max(ys))]


class SlantedLayout(RectangularLayout):
    """Same node coordinates as rectangular, but edges are drawn as direct
    diagonal segments parent->child (``layout="slanted"``)."""

    def branch_path(self, node: Node) -> List[XY]:
        if node.is_root:
            return []
        return [(node.parent.x, node.parent.y), (node.x, node.y)]

    def child_connector(self, node: Node) -> List[XY]:
        return []


class DendrogramLayout(Layout):
    """Top-down dendrogram (``layout="dendrogram"``): root at the top,
    tips along the bottom, height increasing downward -- the hclust look.

    It is a rectangular tree with the axes transposed: leaf order runs along
    x, distance-from-root runs down y.
    """

    kind = "dendrogram"
    invert_y = False                # coordinates already put the root on top

    def compute(self) -> None:
        for node in self.tree.traverse():
            node.x = node.data["_row"]
            node.y = -node.data["_rootd"]       # root at y=0 (top), tips below

    def branch_path(self, node: Node) -> List[XY]:
        if node.is_root:
            return []
        return [(node.x, node.parent.y), (node.x, node.y)]   # vertical drop

    def child_connector(self, node: Node) -> List[XY]:
        if node.is_leaf:
            return []
        xs = [c.x for c in node.children]
        return [(min(xs), node.y), (max(xs), node.y)]        # horizontal bar
