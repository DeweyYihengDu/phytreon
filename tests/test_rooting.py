"""Rerooting moves the root, never the data on the branches."""
import pytest

import phytreon as pt
from phytreon.core.io import parse_newick

TREES = [
    # unrooted-style: three children at the display root
    "(((A:1,B:1)90:1,(C:1,D:1)80:1)70:1,(E:1,F:1)60:1,G:2);",
    # rooted binary: the old root becomes a pass-through node when rerooted
    "(((A:1,B:1)90:1,(C:1,D:1)80:1)70:1,(E:1,(F:1,G:10)55:1)70:1);",
]


def _splits(tree):
    """{non-trivial bipartition: supports found on it}, root-independent."""
    every = frozenset(tree.leaf_names())
    out = {}
    for node in tree.traverse():
        if node.is_leaf or node.is_root:
            continue
        side = frozenset(node.leaf_names())
        key = min(side, every - side, key=lambda s: (len(s), sorted(s)))
        if 1 < len(key) < len(every) - 1:
            out.setdefault(key, set()).add(node.support)
    return out


@pytest.mark.parametrize("src", TREES)
@pytest.mark.parametrize("reroot", [
    pt.midpoint_root,
    lambda t: pt.outgroup_root(t, "E"),
    lambda t: pt.outgroup_root(t, ["C", "D"]),
    lambda t: pt.outgroup_root(t, "G"),
])
def test_support_stays_on_its_split_and_no_pass_through_node_is_left(src, reroot):
    before = parse_newick(src)
    after = reroot(parse_newick(src))
    assert _splits(after) == _splits(before)
    assert after.total_branch_length() == pytest.approx(before.total_branch_length())
    assert not [n for n in after.traverse() if len(n.children) == 1]


def test_deep_ladder_tree_can_be_rerooted():
    import sys
    text = "t0:1"
    for i in range(1, sys.getrecursionlimit() + 200):
        text = f"({text},t{i}:1):1"
    tree = parse_newick(text + ";")
    assert pt.midpoint_root(tree).n_leaves == tree.n_leaves
