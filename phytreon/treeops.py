"""Tree manipulation: reshape topology and reading order without
re-inferring the tree.

These cover the common reshaping operations -- ``rotate`` / ``flip`` /
``collapse`` / ``scale_clade`` and ``cut_tree``.  Layout assigns tip rows by
depth-first leaf order, so reordering ``node.children`` is exactly what moves
branches up/down on the plot -- that is how "freely adjust branch positions"
works here.

All functions mutate the tree in place and return it (or the cluster map for
:func:`cut_tree`), so they chain.
"""
from __future__ import annotations

from typing import Dict, List, Optional

from .core.tree import Node, Tree


# --------------------------------------------------------------------------
# reordering branches (the "adjust positions" tools)
# --------------------------------------------------------------------------
def rotate(tree: Tree, node: Node) -> Tree:
    """Reverse the child order at ``node``.

    Flips the vertical arrangement of that clade's subtrees.
    """
    node.children.reverse()
    return tree


def swap_children(tree: Tree, node: Node, order: List[int]) -> Tree:
    """Set an arbitrary child order at ``node`` by index permutation."""
    if sorted(order) != list(range(len(node.children))):
        raise ValueError("order must be a permutation of child indices")
    node.children = [node.children[i] for i in order]
    return tree


def flip(tree: Tree, node_a: Node, node_b: Node) -> Tree:
    """Swap the vertical positions of two clades.

    ``node_a`` and ``node_b`` must be siblings *or* share a path; we swap them
    within their common parent's child list when they are siblings, else we
    swap the two ancestor branches descending from their MRCA.
    """
    pa, pb = node_a.parent, node_b.parent
    if pa is pb is not None:
        ia, ib = pa.children.index(node_a), pa.children.index(node_b)
        pa.children[ia], pa.children[ib] = pa.children[ib], pa.children[ia]
        return tree
    # general case: find the children of the MRCA leading to each node
    mrca = tree.get_mrca([*_leaf_names(node_a), *_leaf_names(node_b)])
    if mrca is None:
        raise ValueError("nodes are not in the same tree")
    ca = _child_towards(mrca, node_a)
    cb = _child_towards(mrca, node_b)
    ia, ib = mrca.children.index(ca), mrca.children.index(cb)
    mrca.children[ia], mrca.children[ib] = mrca.children[ib], mrca.children[ia]
    return tree


def ladderize(tree: Tree, ascending: bool = True) -> Tree:
    """Order every node's children by subtree size (delegates to Tree)."""
    return tree.ladderize(ascending=ascending)


# --------------------------------------------------------------------------
# collapse low-support edges -> polytomies
# --------------------------------------------------------------------------
def collapse_low_support(tree: Tree, threshold: float) -> Tree:
    """Contract internal edges whose ``support`` < ``threshold`` into
    polytomies.

    The collapsed node's children are re-parented to its parent and the
    collapsed branch length is added onto each child.
    """
    # postorder so we collapse deep nodes before their ancestors
    for node in list(tree.traverse("postorder")):
        if node.is_leaf or node.is_root or node.parent is None:
            continue
        if node.support is not None and node.support < threshold:
            _contract(node)
    return tree


def _contract(node: Node) -> None:
    parent = node.parent
    idx = parent.children.index(node)
    extra = node.length or 0.0
    for c in node.children:
        c.length = (c.length or 0.0) + extra
        c.parent = parent
    parent.children[idx:idx + 1] = node.children


# --------------------------------------------------------------------------
# scale a clade's branch lengths (display emphasis)
# --------------------------------------------------------------------------
def scale_clade(tree: Tree, node: Node, factor: float) -> Tree:
    """Multiply every branch length inside ``node``'s subtree by ``factor``."""
    for n in node.traverse("preorder"):
        if n is node:
            continue
        if n.length is not None:
            n.length *= factor
    return tree


# --------------------------------------------------------------------------
# grouping clades -> colour the tree by lineage
# --------------------------------------------------------------------------
def group_clade(tree: Tree, mapping: Dict[Node, str], key: str = "group",
                default: Optional[str] = None) -> Tree:
    """Label clades for colouring.

    ``mapping`` maps a node -> group label; every node in that node's subtree
    gets ``data[key] = label`` (later mappings win for nested clades).  Nodes
    not covered get ``default``.  Use with ``branches(color=key)`` /
    ``tip_labels(color=key)`` to colour by lineage.
    """
    for node in tree.traverse():
        node.data[key] = default
    # shallow clades first so deeper/nested labels override
    for anchor in sorted(mapping, key=lambda n: n.depth(use_lengths=False)):
        for d in anchor.traverse("preorder"):
            d.data[key] = mapping[anchor]
    return tree


def group_otu(tree: Tree, mapping: Dict[str, list], key: str = "group",
              default: Optional[str] = None) -> Tree:
    """Like :func:`group_clade` but keyed by tip-name sets.

    ``mapping`` maps a group label -> list of tip names; the label is applied
    to the smallest clade (MRCA subtree) containing those tips.
    """
    node_map: Dict[Node, str] = {}
    for label, tips in mapping.items():
        mrca = tree.get_mrca(tips)
        if mrca is not None:
            node_map[mrca] = label
    return group_clade(tree, node_map, key=key, default=default)


# --------------------------------------------------------------------------
# tree comparison
# --------------------------------------------------------------------------
def _bipartition_set(tree: Tree):
    leaves = frozenset(tree.leaf_names())
    anchor = min(leaves)
    n = len(leaves)
    s = set()
    for node in tree.traverse():
        if node.is_leaf or node.is_root:
            continue
        side = frozenset(node.leaf_names())
        if 2 <= len(side) <= n - 2:
            s.add(side if anchor not in side else (leaves - side))
    return s


def robinson_foulds(t1: Tree, t2: Tree, normalized: bool = False) -> float:
    """Robinson-Foulds (symmetric-difference) distance between two trees.

    Counts bipartitions present in one tree but not the other (rooting-
    independent).  ``normalized=True`` divides by the maximum possible (2n-6).
    """
    if frozenset(t1.leaf_names()) != frozenset(t2.leaf_names()):
        raise ValueError("trees must have the same taxon set")
    b1, b2 = _bipartition_set(t1), _bipartition_set(t2)
    rf = len(b1 ^ b2)
    if normalized:
        n = t1.n_leaves
        denom = 2 * (n - 3) or 1
        return rf / denom
    return float(rf)


# --------------------------------------------------------------------------
# rooting (essential for unrooted NJ trees)
# --------------------------------------------------------------------------
def _adjacency(tree: Tree):
    adj: Dict[Node, List] = {}
    for n in tree.traverse():
        if n.parent is not None:
            w = n.length or 0.0
            adj.setdefault(n, []).append((n.parent, w))
            adj.setdefault(n.parent, []).append((n, w))
    return adj


def _far_leaf(adj, src, leaves):
    dist = {src: 0.0}
    prev = {src: None}
    stack = [src]
    while stack:
        u = stack.pop()
        for v, w in adj.get(u, []):
            if v not in dist:
                dist[v] = dist[u] + w
                prev[v] = u
                stack.append(v)
    far = max(leaves, key=lambda l: dist.get(l, 0.0))
    return far, dist, prev


def _rebuild_rooted(tree, adj, u, v, du, dv) -> Tree:
    def build(cur, excl, length):
        node = Node(name=cur.name, length=length, support=cur.support)
        node.data = dict(cur.data)
        for nb, w in adj[cur]:
            if nb is excl:
                continue
            node.add_child(build(nb, cur, w))
        return node
    root = Node()
    root.add_child(build(u, v, du))
    root.add_child(build(v, u, dv))
    return Tree(root=root, name=tree.name)


def midpoint_root(tree: Tree) -> Tree:
    """Re-root on the midpoint of the longest leaf-to-leaf path.

    The standard way to give an unrooted (NJ) tree a sensible root; returns a
    new tree (the input is left unchanged).
    """
    leaves = tree.leaves()
    if len(leaves) < 2:
        return tree
    adj = _adjacency(tree)
    a, _, _ = _far_leaf(adj, leaves[0], leaves)
    b, dist, prev = _far_leaf(adj, a, leaves)
    diameter = dist[b]
    target = diameter / 2.0

    # path a..b, then walk accumulating until we cross the midpoint
    path = []
    node = b
    while node is not None:
        path.append(node)
        node = prev[node]
    path.reverse()                                  # a -> b
    acc = 0.0
    for i in range(len(path) - 1):
        u, v = path[i], path[i + 1]
        w = next(ww for nb, ww in adj[u] if nb is v)
        if acc + w >= target:
            return _rebuild_rooted(tree, adj, u, v, target - acc, acc + w - target)
        acc += w
    return tree


# --------------------------------------------------------------------------
# cut the tree into clusters (cutree)
# --------------------------------------------------------------------------
def cut_tree(tree: Tree, height: Optional[float] = None,
             k: Optional[int] = None) -> Dict[str, int]:
    """Cut the tree into clusters and return ``{tip_name: cluster_id}``.

    Provide exactly one of:
      * ``height`` -- cut at a fixed root-distance; each maximal subtree whose
        stem crosses the line becomes a cluster.
      * ``k`` -- choose the height that yields ``k`` clusters (by collapsing
        the deepest internal splits first, like hierarchical ``cutree(k=)``).
    """
    if (height is None) == (k is None):
        raise ValueError("provide exactly one of height= or k=")

    depth = {n: n.depth(use_lengths=True) for n in tree.traverse()}

    if k is not None:
        # grow clusters by repeatedly expanding the shallowest internal node
        # (= cutting the k-1 shallowest edges). Exact for binary trees; with
        # multifurcations the count may overshoot to the next achievable k.
        roots: List[Node] = [tree.root]
        while len(roots) < k:
            internal = [n for n in roots if not n.is_leaf]
            if not internal:
                break
            node = min(internal, key=lambda n: depth[n])
            roots.remove(node)
            roots.extend(node.children)
    else:
        # height cut: a cluster root is any node whose edge crosses `height`
        roots = []
        for n in tree.traverse():
            pd = depth[n.parent] if n.parent is not None else float("-inf")
            if pd < height <= depth[n]:
                roots.append(n)

    clusters: Dict[str, int] = {}
    for cid, root in enumerate(roots):
        for leaf in root.get_leaves():
            clusters.setdefault(leaf.name, cid)
    # leaves above the cut (none of the roots covers them) -> singletons
    nxt = len(roots)
    for leaf in tree.leaves():
        if leaf.name not in clusters:
            clusters[leaf.name] = nxt
            nxt += 1
    return clusters


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def _leaf_names(node: Node) -> List[str]:
    return [l.name for l in node.get_leaves()]


def _child_towards(ancestor: Node, descendant: Node) -> Node:
    """The child of ``ancestor`` that lies on the path to ``descendant``."""
    node = descendant
    while node.parent is not None and node.parent is not ancestor:
        node = node.parent
    return node
