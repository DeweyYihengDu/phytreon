"""Circular / fan layout (``layout="circular"`` / ``"fan"``)."""
from __future__ import annotations

import math
from typing import List, Tuple

from ..core.tree import Node
from .base import Layout

XY = Tuple[float, float]


class CircularLayout(Layout):
    """Radius = distance from root, angle = leaf row mapped onto an arc.

    ``start`` / ``extent`` (degrees) control the fan opening, so a fan
    plot is just a circular plot with ``extent < 360``.  Branch elbows
    become a radial segment plus an arc connector, both sampled to
    polylines so the backends stay dumb.
    """

    is_polar = True
    kind = "polar"
    equal_aspect = True
    inward = False

    def __init__(self, use_branch_lengths: bool = True,
                 start: float = 0.0, extent: float = 350.0,
                 inner_radius: float = 0.0):
        super().__init__(use_branch_lengths)
        self.start = math.radians(start)
        self.extent = math.radians(extent)
        self.inner_radius = inner_radius
        self.center: XY = (0.0, 0.0)

    # -- polar helpers ---------------------------------------------------
    def _polar_to_xy(self, r: float, a: float) -> XY:
        return (r * math.cos(a), r * math.sin(a))

    def _angle_of_row(self, row: float) -> float:
        denom = max(self.n_leaves - 1, 1)
        return self.start + (row / denom) * self.extent

    def _radius_of(self, node: Node) -> float:
        return self.inner_radius + node.data["_rootd"]

    def compute(self) -> None:
        for node in self.tree.traverse():
            r = self._radius_of(node)
            a = self._angle_of_row(node.data["_row"])
            node._r = r
            node._angle = a
            node.x, node.y = self._polar_to_xy(r, a)

    def branch_path(self, node: Node) -> List[XY]:
        if node.is_root:
            return []
        # radial segment along this node's angle, from parent radius to node radius
        a = node._angle
        return [self._polar_to_xy(node.parent._r, a),
                self._polar_to_xy(node._r, a)]

    def child_connector(self, node: Node) -> List[XY]:
        if node.is_leaf:
            return []
        angles = [c._angle for c in node.children]
        return self._arc(node._r, min(angles), max(angles))

    def _arc(self, r: float, a0: float, a1: float, step: float = math.radians(2)) -> List[XY]:
        n = max(2, int(abs(a1 - a0) / step) + 1)
        return [self._polar_to_xy(r, a0 + (a1 - a0) * i / (n - 1)) for i in range(n)]


class InwardCircularLayout(CircularLayout):
    """Circular tree drawn inward (``layout="inward_circular"``): the root
    sits on the outer rim and the tips point toward the centre."""

    inward = True

    def _radius_of(self, node: Node) -> float:
        # reverse the radius: root at the outer edge, tips near the centre
        return self.inner_radius + (self.max_x - node.data["_rootd"])
