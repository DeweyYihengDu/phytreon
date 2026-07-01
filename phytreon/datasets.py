"""Tiny bundled example data for demos and tests."""
from __future__ import annotations

from .core.tree import Tree

# A small primate tree with branch lengths and internal support values.
PRIMATES_NEWICK = (
    "(((((Human:0.06,Chimp:0.06)95:0.04,Gorilla:0.10)90:0.04,"
    "Orangutan:0.16)98:0.06,Gibbon:0.20)100:0.10,"
    "(Macaque:0.20,Baboon:0.18)85:0.12);"
)


def primates() -> Tree:
    """A rooted primate tree (7 tips) with support values."""
    return Tree.from_newick(PRIMATES_NEWICK)


def random_tree(n: int = 60, seed: int = 0) -> Tree:
    """A random bifurcating tree with branch lengths (for demos/large layouts)."""
    import random
    from .core.tree import Node
    rng = random.Random(seed)
    pool = [Node(name=f"t{i:02d}", length=rng.uniform(0.02, 0.12)) for i in range(n)]
    while len(pool) > 1:
        rng.shuffle(pool)
        a, b = pool.pop(), pool.pop()
        p = Node(length=rng.uniform(0.04, 0.18))
        p.add_child(a)
        p.add_child(b)
        pool.append(p)
    root = pool[0]
    root.length = None
    return Tree(root=root)


def primates_metadata():
    """Per-tip metadata as a DataFrame indexed by tip name."""
    import pandas as pd
    return pd.DataFrame(
        {
            "name": ["Human", "Chimp", "Gorilla", "Orangutan",
                     "Gibbon", "Macaque", "Baboon"],
            "habitat": ["urban", "forest", "forest", "forest",
                        "forest", "savanna", "savanna"],
            "body_mass_kg": [62, 45, 160, 75, 8, 11, 25],
        }
    ).set_index("name")
