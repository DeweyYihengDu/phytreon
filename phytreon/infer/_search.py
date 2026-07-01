"""Shared topology-search primitives (NNI), used by ML and parsimony search."""
from __future__ import annotations


def internal_edges(tree):
    """Internal nodes whose parent is also internal -- the NNI-able edges
    (both endpoints binary)."""
    return [n for n in tree.traverse()
            if not n.is_leaf and not n.is_root and not n.parent.is_root
            and len(n.children) == 2 and len(n.parent.children) == 2]


def nni_neighbors(node):
    """Yield the two NNI swap closures across the edge between ``node`` and its
    parent.  Each closure is its own inverse (an involution)."""
    parent = node.parent
    sib = parent.children[1] if parent.children[0] is node else parent.children[0]
    for ci in (0, 1):
        child = node.children[ci]

        def make(a=sib, b=child):
            def swap():
                pa, pb = a.parent, b.parent
                ia, ib = pa.children.index(a), pb.children.index(b)
                pa.children[ia], pb.children[ib] = b, a
                a.parent, b.parent = pb, pa
            return swap
        yield make()
