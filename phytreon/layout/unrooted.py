"""Unrooted layouts.

:class:`EqualAngleLayout` implements the equal-angle unrooted layout: every
subtree is allotted an angular wedge
proportional to its number of tips, and nodes are placed outward from the
root along the bisector of their wedge.  Branches are straight lines, so the
result is the classic unrooted "star/dendrite" tree.

``layout="unrooted"`` maps to the equal-*daylight* algorithm, an iterative
refinement of equal-angle that evens out the blank space between subtrees
(Felsenstein, *Inferring Phylogenies*, pp. 582-584).  Daylight is a documented
extension point on top of this class (see :meth:`compute`).
"""
from __future__ import annotations

import math
from typing import Dict, List, Tuple

from ..core.tree import Node
from .base import Layout

XY = Tuple[float, float]


class EqualAngleLayout(Layout):
    """Equal-angle unrooted layout (angles in half-rotation units)."""

    is_polar = False
    kind = "radial"
    equal_aspect = True

    def compute(self) -> None:
        tree = self.tree
        use_len = self.use_branch_lengths and tree.has_branch_lengths

        # tips-per-subtree (postorder so children are counted first)
        ntip: Dict[Node, int] = {}
        for node in tree.traverse("postorder"):
            ntip[node] = 1 if node.is_leaf else sum(ntip[c] for c in node.children)

        # root at the origin, owning the full circle [0, 2) half-rotations
        root = tree.root
        root.x = root.y = 0.0
        root._angle = 0.0
        wedge: Dict[Node, Tuple[float, float]] = {root: (0.0, 2.0)}

        # preorder: parent is always placed before its children
        for node in tree.traverse("preorder"):
            start, end = wedge[node]
            total = end - start
            acc = start
            for child in node.children:
                alpha = total * ntip[child] / ntip[node]
                beta = acc + alpha / 2.0           # bisector of this child's wedge
                length = (child.length if use_len else 1.0) or 1.0
                child.x = node.x + math.cos(math.pi * beta) * length
                child.y = node.y + math.sin(math.pi * beta) * length
                child._angle = math.pi * beta      # outward direction (for labels)
                wedge[child] = (acc, acc + alpha)
                acc += alpha

    def branch_path(self, node: Node) -> List[XY]:
        if node.is_root:
            return []
        return [(node.parent.x, node.parent.y), (node.x, node.y)]

    def child_connector(self, node: Node) -> List[XY]:
        return []                                   # straight edges, no connectors


class DaylightLayout(EqualAngleLayout):
    """Equal-*daylight* unrooted layout (Felsenstein, *Inferring Phylogenies*).

    Starts from equal-angle, then iteratively rotates the subtrees around each
    internal node so the blank wedges ("daylight") between them are equal --
    this spreads the tree out and removes the lopsidedness of equal-angle.
    """

    def __init__(self, use_branch_lengths: bool = True, iterations: int = 5,
                 min_change: float = 0.05):
        super().__init__(use_branch_lengths)
        self.iterations = iterations
        self.min_change = min_change          # radians; stop when avg change below

    def compute(self) -> None:
        super().compute()                      # equal-angle initialisation
        nodes = self.tree.nodes()
        internals = [n for n in nodes if not n.is_leaf]
        for _ in range(self.iterations):
            total = 0.0
            for v in internals:
                total += self._equalize(v, nodes)
            if internals and total / len(internals) <= self.min_change:
                break
        # refresh outward directions for tip-label rotation
        for n in self.tree.traverse():
            if not n.is_root:
                n._angle = math.atan2(n.y - n.parent.y, n.x - n.parent.x)

    def _equalize(self, v: Node, allnodes) -> float:
        sub_v = {id(n) for n in v.traverse()}
        groups = [list(c.traverse()) for c in v.children]
        if not v.is_root:
            groups.append([n for n in allnodes if id(n) not in sub_v])
        if len(groups) < 2:
            return 0.0
        vx, vy = v.x, v.y

        specs = []
        for g in groups:
            angs = sorted(math.atan2(n.y - vy, n.x - vx) for n in g)
            if len(angs) == 1:
                start, width = angs[0], 0.0
            else:
                gaps = [(angs[(i + 1) % len(angs)] - angs[i]) % (2 * math.pi)
                        for i in range(len(angs))]
                k = max(range(len(gaps)), key=lambda i: gaps[i])
                start = angs[(k + 1) % len(angs)]
                width = 2 * math.pi - gaps[k]
            specs.append([g, start, width])

        specs.sort(key=lambda s: s[1])
        daylight = 2 * math.pi - sum(s[2] for s in specs)
        if daylight < 0:
            return 0.0
        gap = daylight / len(specs)
        max_change = 0.0
        cursor = specs[0][1] + specs[0][2] + gap        # first subtree fixed
        for i in range(1, len(specs)):
            g, start, width = specs[i]
            delta = ((cursor - start + math.pi) % (2 * math.pi)) - math.pi
            self._rotate(g, vx, vy, delta)
            max_change = max(max_change, abs(delta))
            cursor = cursor + width + gap
        return max_change

    @staticmethod
    def _rotate(group, cx, cy, delta):
        ca, sa = math.cos(delta), math.sin(delta)
        for n in group:
            dx, dy = n.x - cx, n.y - cy
            n.x = cx + dx * ca - dy * sa
            n.y = cy + dx * sa + dy * ca
