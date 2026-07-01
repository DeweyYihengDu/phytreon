"""Layout base class.

A *layout* turns a tree topology into display coordinates.  To keep both
rendering backends identical, every layout writes final cartesian
coordinates onto ``node.x`` / ``node.y`` and (for polar layouts)
``node._r`` / ``node._angle``.  Geoms then read those coordinates and the
layout's small helper API (``is_polar``, ``branch_path`` ...) to emit
scene primitives.

Adding a new layout (slanted, unrooted/daylight, radial) means
subclassing :class:`Layout` and implementing :meth:`compute` -- nothing
in the plotting layer needs to change.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

from ..core.tree import Node, Tree

XY = Tuple[float, float]


class Layout:
    """Abstract base. Subclasses implement :meth:`compute`."""

    is_polar = False
    #: how geoms place labels: "rect" | "polar" | "radial"
    kind = "rect"
    #: keep a 1:1 data aspect (circular/unrooted) and do not flip the y axis
    equal_aspect = False
    #: rectangular renders with the y axis flipped (first leaf on top); a
    #: dendrogram already encodes orientation in its coordinates, so it opts out
    invert_y = True

    def __init__(self, use_branch_lengths: bool = True):
        self.use_branch_lengths = use_branch_lengths
        self.tree: Tree | None = None
        # populated by _assign_grid():
        self._order: List[Node] = []        # leaves top->bottom
        self.max_x: float = 1.0
        self.n_leaves: int = 0

    # -- public ----------------------------------------------------------
    def apply(self, tree: Tree) -> "Layout":
        self.tree = tree
        self._assign_grid()
        self.compute()
        return self

    def compute(self) -> None:                       # pragma: no cover
        raise NotImplementedError

    def branch_path(self, node: Node) -> List[XY]:   # pragma: no cover
        """Points of the branch connecting ``node`` to its parent."""
        raise NotImplementedError

    def child_connector(self, node: Node) -> List[XY]:
        """Points of the connector spanning ``node``'s children.

        For a rectangular tree this is the vertical bar; for a circular
        tree it is the arc.  Returns ``[]`` for leaves.
        """
        return []

    # -- shared helpers --------------------------------------------------
    def _assign_grid(self) -> None:
        """Compute the abstract grid: leaf order (rows) and root distance.

        ``node._row`` (stored on ``node.data``) is the leaf index for tips
        and the centred mid-row for internal nodes.  ``node._rootd`` is the
        distance from the root (branch-length sum, or edge count for a
        cladogram).  These are layout-independent; rectangular and circular
        both build on them.
        """
        assert self.tree is not None
        use_len = self.use_branch_lengths and self.tree.has_branch_lengths

        leaves = self.tree.leaves()
        self._order = leaves
        self.n_leaves = len(leaves)
        for i, leaf in enumerate(leaves):
            leaf.data["_row"] = float(i)

        # internal node row = midpoint of extreme children (centred)
        for node in self.tree.traverse("postorder"):
            if not node.is_leaf:
                rows = [c.data["_row"] for c in node.children]
                node.data["_row"] = 0.5 * (min(rows) + max(rows))

        # root distance (preorder so parent is set first)
        max_d = 0.0
        for node in self.tree.traverse("preorder"):
            if node.is_root:
                node.data["_rootd"] = 0.0
            else:
                step = (node.length or 0.0) if use_len else 1.0
                node.data["_rootd"] = node.parent.data["_rootd"] + step
            max_d = max(max_d, node.data["_rootd"])
        self.max_x = max_d or 1.0
