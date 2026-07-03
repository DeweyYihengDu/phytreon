"""Unified tree data model.

A small, dependency-light ``Node`` / ``Tree`` pair that is the single
representation used by every other sub-package (layout, inference,
comparative methods, plotting).  Parsing is delegated to Biopython
(see :mod:`phytreon.core.io`) but we immediately convert into our own
nodes so the rest of the library never touches Biopython internals.

Design notes
------------
* ``Node.data`` is a free-form ``dict`` used for everything that is not a
  first-class attribute: joined metadata columns, ancestral-state
  reconstructions, bootstrap values from external programs, etc.  This is
  the analogue of ``treeio``'s "treedata" -- arbitrary data travelling
  alongside the topology.
* Layout fills ``x`` / ``y`` (display coordinates) and, for polar
  layouts, ``_r`` / ``_angle``.  Keeping them on the node keeps geoms
  simple.
"""
from __future__ import annotations

from typing import Dict, Iterable, Iterator, List, Optional


class Node:
    """A single tree node (works for both tips and internal nodes)."""

    __slots__ = (
        "name", "length", "support", "comment",
        "parent", "children", "data",
        "x", "y", "_r", "_angle",
    )

    def __init__(
        self,
        name: Optional[str] = None,
        length: Optional[float] = None,
        support: Optional[float] = None,
        comment: Optional[str] = None,
    ):
        self.name = name
        self.length = length            # branch length to parent
        self.support = support          # node support (bootstrap/posterior)
        self.comment = comment
        self.parent: Optional["Node"] = None
        self.children: List["Node"] = []
        self.data: Dict[str, object] = {}
        # filled by layout:
        self.x: float = 0.0
        self.y: float = 0.0
        self._r: float = 0.0            # polar radius (circular layout)
        self._angle: float = 0.0        # polar angle in radians

    # -- structure -------------------------------------------------------
    @property
    def is_leaf(self) -> bool:
        return not self.children

    @property
    def is_root(self) -> bool:
        return self.parent is None

    def add_child(self, child: "Node") -> "Node":
        child.parent = self
        self.children.append(child)
        return child

    # -- traversal -------------------------------------------------------
    def traverse(self, order: str = "preorder") -> Iterator["Node"]:
        """Iterate over this node and all descendants.

        ``order`` is one of ``"preorder"``, ``"postorder"`` or
        ``"levelorder"``.
        """
        if order == "preorder":
            yield self
            for c in self.children:
                yield from c.traverse(order)
        elif order == "postorder":
            for c in self.children:
                yield from c.traverse(order)
            yield self
        elif order == "levelorder":
            queue = [self]
            while queue:
                node = queue.pop(0)
                yield node
                queue.extend(node.children)
        else:
            raise ValueError(f"unknown traversal order: {order!r}")

    def iter_leaves(self) -> Iterator["Node"]:
        for n in self.traverse("preorder"):
            if n.is_leaf:
                yield n

    def get_leaves(self) -> List["Node"]:
        return list(self.iter_leaves())

    def leaf_names(self) -> List[str]:
        return [n.name for n in self.iter_leaves()]

    # -- metrics ---------------------------------------------------------
    def depth(self, use_lengths: bool = True) -> float:
        """Distance from the root (sum of branch lengths or edge count)."""
        d, node = 0.0, self
        while node.parent is not None:
            d += (node.length or 0.0) if use_lengths else 1.0
            node = node.parent
        return d

    def __repr__(self) -> str:
        kind = "leaf" if self.is_leaf else f"node({len(self.children)} children)"
        return f"<Node {self.name!r} {kind}>"


class Tree:
    """A rooted tree: a thin wrapper around a root :class:`Node`."""

    def __init__(self, root: Optional[Node] = None, name: Optional[str] = None):
        self.root = root if root is not None else Node()
        self.name = name
        self.data: Dict[str, object] = {}      # tree-level metadata (logL, model, ...)

    # -- construction ----------------------------------------------------
    @classmethod
    def from_newick(cls, newick: str) -> "Tree":
        from .io import parse_newick
        return parse_newick(newick)

    @classmethod
    def read(cls, path: str, fmt: str = "newick") -> "Tree":
        from .io import read as _read
        return _read(path, fmt)

    @classmethod
    def from_biopython(cls, bp_tree) -> "Tree":
        from .io import from_biopython as _fb
        return _fb(bp_tree)

    def write(self, path: Optional[str] = None, fmt: str = "newick") -> Optional[str]:
        from .io import write as _write
        return _write(self, path, fmt)

    # -- traversal / access ---------------------------------------------
    def traverse(self, order: str = "preorder") -> Iterator[Node]:
        yield from self.root.traverse(order)

    def nodes(self, order: str = "preorder") -> List[Node]:
        return list(self.traverse(order))

    def leaves(self) -> List[Node]:
        return self.root.get_leaves()

    def leaf_names(self) -> List[str]:
        return self.root.leaf_names()

    @property
    def n_leaves(self) -> int:
        return sum(1 for _ in self.root.iter_leaves())

    def search_nodes(self, **kwargs) -> List[Node]:
        """Return nodes whose attributes / data match all kwargs."""
        out = []
        for n in self.traverse():
            ok = True
            for k, v in kwargs.items():
                cur = getattr(n, k, None) if hasattr(n, k) else n.data.get(k)
                if cur != v:
                    ok = False
                    break
            if ok:
                out.append(n)
        return out

    def get_node(self, name: str) -> Optional[Node]:
        for n in self.traverse():
            if n.name == name:
                return n
        return None

    def get_mrca(self, names: Iterable[str]) -> Optional[Node]:
        """Most recent common ancestor of the named tips."""
        targets = set(names)
        # path from each target tip up to the root
        paths = []
        for leaf in self.leaves():
            if leaf.name in targets:
                p, node = [], leaf
                while node is not None:
                    p.append(node)
                    node = node.parent
                paths.append(p)
        if not paths:
            return None
        common = set(paths[0])
        for p in paths[1:]:
            common &= set(p)
        # deepest common node = the one furthest from root
        return max(common, key=lambda n: n.depth(use_lengths=False)) if common else None

    # -- shape -----------------------------------------------------------
    def ladderize(self, ascending: bool = True) -> "Tree":
        """Order children by subtree size for a tidy 'laddered' look."""
        def size(node: Node) -> int:
            if node.is_leaf:
                return 1
            s = sum(size(c) for c in node.children)
            node.children.sort(key=size, reverse=not ascending)
            return s
        size(self.root)
        return self

    @property
    def has_branch_lengths(self) -> bool:
        return all(
            n.length is not None
            for n in self.traverse() if not n.is_root
        )

    def total_branch_length(self) -> float:
        return sum((n.length or 0.0) for n in self.traverse())

    # -- associated data (treeio-style) ---------------------------------
    def join_data(self, df, on: str = "name") -> "Tree":
        """Attach a :class:`pandas.DataFrame` of metadata to nodes.

        Each row is matched to a node by ``on`` (a column name; ``"name"``
        matches against ``node.name``).  Every other column becomes an
        entry in ``node.data``.  Returns ``self`` for chaining.
        """
        cols = [c for c in df.columns if c != on]
        index = {str(row[on]): row for _, row in df.iterrows()}
        for node in self.traverse():
            key = str(node.name) if on == "name" else str(node.data.get(on))
            row = index.get(key)
            if row is not None:
                for c in cols:
                    node.data[c] = row[c]
        return self

    def __repr__(self) -> str:
        return f"<Tree {self.name or ''!r} n_leaves={self.n_leaves}>"
