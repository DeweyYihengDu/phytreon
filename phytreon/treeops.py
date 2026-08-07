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

from collections import Counter
from itertools import permutations
from typing import Callable, Dict, Iterable, List, Optional, Union

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


def sort_by(tree: Tree, key: Union[str, Dict[str, object], Callable[[Node], object]],
            order: Optional[List[str]] = None, rounds: int = 4) -> Tree:
    """Reorder sibling branches so tips sharing a category sit together.

    ``key`` reads the grouping label per tip: a column name looked up in
    ``leaf.data`` (e.g. after ``tree.join_data(meta, on="name")``), a
    ``{tip_name: label}`` mapping, or a function of the leaf node.

    Only ``node.children`` order changes -- the same move :func:`ladderize`
    and :func:`untangle` already make, no branch moves and no split is added
    or dropped -- which bounds what this can do. Rotating a fork brings two
    clades that are already siblings closer together; it cannot join clades
    the tree keeps apart. If a label is not monophyletic under this tree --
    common for genus-level 16S calls, which a single hypervariable region
    often cannot resolve -- its tips collect into as few separate runs as the
    topology allows rather than one, and no amount of reordering closes the
    rest of the gap. A tree that keeps splitting a genus in two is telling
    you something the label does not; forcing it into one run by editing the
    topology would tell the reader the opposite of what the data support.
    ``highlight(by=key)`` makes those remaining runs visible as separate
    bands instead of hiding the split.

    This is :func:`untangle`'s own search, aimed at a different scorecard:
    that one counts crossings against a reference tree, this one counts
    label changes walking the tips left to right, and a rotation is kept only
    when it lowers that count -- so, unlike sorting every node by a per-node
    summary in one pass (the first cut at this, which on a real 106-taxon 16S
    tree made the phylum grouping *worse* -- 33 runs against 31 plainly
    ladderized, because a summary taken in isolation at one node knows
    nothing about what its own neighbours need), the result is never worse
    than what the tree started with, at the cost of being a greedy hill climb
    rather than a provably optimal reordering: it finds a good arrangement,
    not the best one over every possible rotation. Polytomies of up to six
    children are searched exactly (6! = 720 arrangements); wider ones sort by
    the subtree's majority category, which is not guaranteed to help and is
    kept only if it measurably does.

    ``order`` breaks ties among arrangements that score equally on label
    changes, preferring the one whose categories run in this left-to-right
    sequence (unlisted ones follow, alphabetically) -- the same convention as
    ``highlight(order=)``/``ring(order=)``. Unlabelled tips (``key`` gives
    ``None``) count as their own group.
    """
    if callable(key):
        catfun = key
    elif isinstance(key, dict):
        catfun: Callable[[Node], object] = lambda leaf: key.get(leaf.name)
    elif isinstance(key, str):
        catfun = lambda leaf: leaf.data.get(key)
    else:
        raise TypeError("key must be a column name, a {tip_name: label} "
                        "mapping, or a function of the leaf node")

    cats: Dict[str, Optional[str]] = {}
    for leaf in tree.leaves():
        v = catfun(leaf)
        cats[leaf.name] = None if v is None else str(v)

    present = sorted({c for c in cats.values() if c is not None})
    if order:
        wanted = [str(c) for c in order]
        present = [c for c in wanted if c in present] + \
            [c for c in present if c not in wanted]
    rank = {c: i for i, c in enumerate(present)}
    none_rank = len(present)          # unlabelled tips trail, as one group

    # The majority category of a node's own subtree never changes across
    # anything this function does -- reordering children moves leaves left
    # and right, never into or out of a subtree -- so it is computed once,
    # bottom-up, and reused for every candidate tried in every round rather
    # than re-walked from scratch each time.
    votes: Dict[Node, Counter] = {}
    for node in tree.traverse("postorder"):
        votes[node] = (Counter([cats[node.name]]) if node.is_leaf
                       else Counter())
        if not node.is_leaf:
            for child in node.children:
                votes[node].update(votes[child])
    node_rank = {n: (none_rank if v.most_common(1)[0][0] is None
                     else rank[v.most_common(1)[0][0]])
                for n, v in votes.items()}

    def transitions() -> int:
        names = tree.leaf_names()
        t = 0
        prev = cats.get(names[0]) if names else None
        for nm in names[1:]:
            c = cats.get(nm)
            if c != prev:
                t += 1
            prev = c
        return t

    for _ in range(rounds):
        changed = False
        for node in tree.traverse("preorder"):
            k = len(node.children)
            if node.is_leaf or k < 2:
                continue
            current = node.children
            best_order = current
            best_key = (transitions(), tuple(node_rank[c] for c in current))
            if k <= 6:
                for perm in permutations(current):
                    perm = list(perm)
                    if perm == current:
                        continue
                    node.children = perm
                    cand = (transitions(), tuple(node_rank[c] for c in perm))
                    if cand < best_key:
                        best_key, best_order = cand, perm
            else:
                cand_order = sorted(current, key=lambda c: node_rank[c])
                if cand_order != current:
                    node.children = cand_order
                    cand = (transitions(),
                           tuple(node_rank[c] for c in cand_order))
                    if cand < best_key:
                        best_key, best_order = cand, cand_order
            node.children = best_order
            changed = changed or best_order != current
        if not changed:
            break
    return tree


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
        if n < 4:
            return 0.0          # no non-trivial bipartitions possible below 4 taxa
        return rf / (2 * (n - 3))
    return float(rf)


# --------------------------------------------------------------------------
# tanglegram support: crossings between two trees' tip orders
# --------------------------------------------------------------------------
def _inversions(seq: List[int]) -> int:
    """Number of inverted pairs in ``seq``, by merge sort -- O(n log n)."""
    def sort(a: List[int]):
        if len(a) < 2:
            return a, 0
        mid = len(a) // 2
        left, il = sort(a[:mid])
        right, ir = sort(a[mid:])
        merged: List[int] = []
        cross = 0
        i = j = 0
        while i < len(left) and j < len(right):
            if left[i] <= right[j]:
                merged.append(left[i])
                i += 1
            else:
                # right[j] jumps ahead of every left element still pending
                cross += len(left) - i
                merged.append(right[j])
                j += 1
        merged.extend(left[i:])
        merged.extend(right[j:])
        return merged, il + ir + cross
    return sort(list(seq))[1]


def _shared_rank_sequence(t1: Tree, t2: Tree) -> List[int]:
    """``t2``'s tip order expressed as ranks in ``t1``'s order (shared tips)."""
    shared = set(t1.leaf_names()) & set(t2.leaf_names())
    rank = {name: i for i, name in
            enumerate(n for n in t1.leaf_names() if n in shared)}
    return [rank[n] for n in t2.leaf_names() if n in shared]


def crossing_number(t1: Tree, t2: Tree) -> int:
    """Number of crossing tip-to-tip links when ``t1`` and ``t2`` face each
    other in a tanglegram.

    Two links cross exactly when their tips appear in opposite vertical order
    in the two trees, so this is the number of inverted pairs between the two
    tip orderings.  Tips present in only one tree are ignored.

    Unlike :func:`robinson_foulds`, this measures *display* discordance -- it
    depends on how each tree's children happen to be ordered, so minimise it
    with :func:`untangle` before reading anything into the value.

    Zero crossings does **not** mean the two trees are the same: it only means
    some tip order satisfies both, which conflicting splits can still allow.
    Always read :func:`robinson_foulds` alongside it for the actual
    topological difference.
    """
    return _inversions(_shared_rank_sequence(t1, t2))


def _untangle_against(target: Tree, reference: Tree, rounds: int) -> int:
    """Greedily reverse child order in ``target`` to reduce crossings."""
    best = crossing_number(reference, target)
    for _ in range(rounds):
        improved = False
        for node in target.traverse("preorder"):
            if node.is_leaf or len(node.children) < 2:
                continue
            node.children.reverse()
            score = crossing_number(reference, target)
            if score < best:
                best = score
                improved = True
            else:
                node.children.reverse()          # no gain -- put it back
        if not improved:
            break                                # local optimum reached
    return best


def untangle(t1: Tree, t2: Tree, *, fix: Optional[str] = "left",
             rounds: int = 3) -> int:
    """Rotate nodes so the two trees' shared tips line up, and return the
    resulting :func:`crossing_number`.

    Rotating a node reverses the order of its children, which slides that
    clade up or down the plot without changing the topology -- so untangling
    only changes how the trees *read*, never what they say.  The search is a
    greedy hill-climb over single rotations (the standard heuristic; it finds
    a good arrangement, not a provably optimal one).

    ``fix="left"`` holds ``t1``'s tip order and rotates only ``t2``;
    ``fix="right"`` does the reverse; ``fix=None`` alternates between the two,
    which usually untangles further but leaves neither tree in its original
    order.  Both trees are modified in place.
    """
    if fix not in ("left", "right", None):
        raise ValueError("fix must be 'left', 'right' or None")
    if fix == "left":
        return _untangle_against(t2, t1, rounds)
    if fix == "right":
        return _untangle_against(t1, t2, rounds)
    best = crossing_number(t1, t2)
    for _ in range(rounds):
        _untangle_against(t2, t1, 1)
        score = _untangle_against(t1, t2, 1)
        if score >= best:
            break                                # alternating stopped helping
        best = score
    return best


# --------------------------------------------------------------------------
# collapse a clade for display (drawn as a triangle)
# --------------------------------------------------------------------------
def collapse_clade(tree: Tree, node: Node, *, name: Optional[str] = None) -> Tree:
    """Collapse ``node``'s subtree into a single tip, drawn as a triangle.

    This is a *display* operation, the standard way to compress a large tree
    down to the clades you actually want to discuss.  The clade's children are
    dropped and the summary needed to draw it is recorded on
    ``node.data["_collapsed"]``:

    ``n`` (tips collapsed), ``near`` / ``far`` (branch length from ``node`` to
    its closest and farthest leaf) and ``leaves`` (the names).  A triangle
    whose two sides use ``near`` and ``far`` therefore shows how deep and how
    ragged the hidden clade is, the same convention iTOL uses.

    Collapsing a clade that already contains a collapsed one is accounted for:
    the outer summary counts the tips hidden inside the inner clade and reaches
    to their true depth, rather than treating the inner clade as a single tip
    sitting at its own node.

    The tree is modified in place, so work on a copy
    (``Tree.from_newick(tree.write())``) to keep the original.  Renders via
    :meth:`~phytreon.plot.figure.TreeFigure.collapsed_clades`.
    """
    if node.is_leaf:
        raise ValueError("cannot collapse a leaf")
    base = node.depth(use_lengths=True)
    base_e = node.depth(use_lengths=False)
    depths: List[float] = []
    edges: List[float] = []
    leaves: List[str] = []
    total = 0
    for leaf in node.get_leaves():
        offset = leaf.depth(use_lengths=True) - base
        offset_e = leaf.depth(use_lengths=False) - base_e
        inner = leaf.data.get("_collapsed")
        if inner:                       # an already-collapsed clade: unfold it
            depths += [offset + inner["near"], offset + inner["far"]]
            edges += [offset_e + inner["near_edges"], offset_e + inner["far_edges"]]
            leaves += list(inner["leaves"])
            total += inner["n"]
        else:
            depths.append(offset)
            edges.append(offset_e)
            leaves.append(leaf.name)
            total += 1
    node.data["_collapsed"] = {
        "n": total,
        "near": min(depths),
        "far": max(depths),
        # the same span counted in edges, for cladogram layouts where the
        # drawing axis is depth-in-edges rather than branch length
        "near_edges": min(edges),
        "far_edges": max(edges),
        "leaves": leaves,
    }
    if name is not None:
        node.name = name
    elif not node.name:
        node.name = f"{leaves[0]} +{total - 1}"
    node.children = []
    return tree


# --------------------------------------------------------------------------
# restrict to a leaf subset
# --------------------------------------------------------------------------
def prune_to_taxa(tree: Tree, taxa: Iterable[str], *, strict: bool = True) -> Tree:
    """Restrict ``tree`` to only the given leaf names, returning a new tree
    (the input is left unchanged).

    Subtrees with no kept leaves are dropped entirely; internal nodes left
    with a single surviving child are collapsed away (their branch length is
    added onto that child), so the result is a clean, minimal tree over
    exactly ``taxa`` rather than one padded with unary "pass-through" nodes.
    ``strict=True`` (default) raises ``ValueError`` listing any requested
    name that isn't a leaf of ``tree`` (mirrors :meth:`Tree.get_mrca`); pass
    ``strict=False`` to silently keep only whichever names are present.
    """
    keep = set(taxa)
    missing = keep - set(tree.leaf_names())
    if missing and strict:
        raise ValueError(f"taxa not found in tree: {sorted(missing)}")
    keep -= missing

    def build(node: Node) -> Optional[Node]:
        if node.is_leaf:
            return _copy_node(node) if node.name in keep else None
        kids = [k for k in (build(c) for c in node.children) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            kids[0].length = (kids[0].length or 0.0) + (node.length or 0.0)
            return kids[0]
        new = _copy_node(node)
        for k in kids:
            new.add_child(k)
        return new

    new_root = build(tree.root) or Node()
    return Tree(root=new_root, name=tree.name)


def _copy_node(node: Node) -> Node:
    new = Node(name=node.name, length=node.length, support=node.support, comment=node.comment)
    new.data = dict(node.data)
    return new


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
    far = max(leaves, key=lambda leaf: dist.get(leaf, 0.0))
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
    return [leaf.name for leaf in node.get_leaves()]


def _child_towards(ancestor: Node, descendant: Node) -> Node:
    """The child of ``ancestor`` that lies on the path to ``descendant``."""
    node = descendant
    while node.parent is not None and node.parent is not ancestor:
        node = node.parent
    return node
