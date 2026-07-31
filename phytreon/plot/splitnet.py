"""Split networks: a phylogeny drawn as a network rather than a tree.

A tree can only say that one grouping is right. Real data often supports two
groupings at once -- recombination, hybridisation, incomplete lineage sorting
all leave that signature -- and a tree resolves the conflict silently, by
picking a winner. A split network refuses to: every conflicting pair of splits
is drawn as a pair of parallel edges, so the conflict appears as a **box**. A
dataset with a clean tree signal draws as a tree; a reticulate one draws a
lattice, and the size of the boxes is the size of the disagreement.

    net = pt.SplitNetwork.from_distances(names, matrix)
    net.color_by(group).save("network.pdf")

This is the picture SplitsTree produces, and it is the standard companion to a
tree wherever recombination is on the table.

Implementation: splits are extracted by neighbour-joining over the distance
matrix (each internal edge of the NJ tree is one split, weighted by its branch
length), then drawn by the split-decomposition convention -- start from a
point and add each split as a displacement shared by every taxon on one side.
Conflicting splits therefore open into boxes. This is the drawing rule of
SplitsTree's ``EqualAngle``, not the full NeighborNet circular ordering; it
recovers the same boxes for the conflicts an NJ tree plus its residuals
expose, but it is not a reimplementation of NeighborNet and will not find
every split that program does.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

from ..scene import Label, Marker, Path, Scene
from .figure import RenderContext, _Renderable, build_color_scale

XY = Tuple[float, float]


class _NetLayout:
    is_polar = False
    equal_aspect = True
    invert_y = False
    kind = "rect"
    use_branch_lengths = False

    def __init__(self, max_x: float = 1.0):
        self.max_x = max_x

    @staticmethod
    def _collapsed_span(node, use_len: bool):
        return (0.0, 0.0)


def splits_from_tree(tree, names: Sequence[str]) -> List[Tuple[frozenset, float]]:
    """Every internal edge of ``tree`` as ``(taxon subset, edge length)``."""
    index = set(names)
    out = []
    for node in tree.traverse():
        if node.is_leaf or node.is_root or node.parent is None:
            continue
        side = frozenset(n for n in node.leaf_names() if n in index)
        if 1 < len(side) < len(index):
            out.append((side, float(node.length or 0.0)))
    return out


def conflicting(a: frozenset, b: frozenset, universe: frozenset) -> bool:
    """True when two splits cannot both sit on one tree.

    Two splits are compatible if one of the four intersections between their
    sides is empty; if all four are populated they conflict, and that is
    exactly what opens a box in the drawing.
    """
    a2, b2 = universe - a, universe - b
    return all((a & b, a & b2, a2 & b, a2 & b2))


class SplitNetwork(_Renderable):
    """A split network -- conflicting splits drawn as boxes.

    Build from a distance matrix (:meth:`from_distances`), an alignment
    (:meth:`from_alignment`), or a set of trees (:meth:`from_trees`), where
    every split's weight is how many of the trees contain it. That last one
    turns a bootstrap or posterior sample straight into a picture of which
    groupings the sample disagrees about.

    ``max_splits`` is the readability knob, and it matters more than it looks.
    A split network only stays legible while the conflict is modest; past that
    it degenerates into a mesh in which nothing can be read, and the median
    closure that draws it grows steeply too. Measured on a 60-replicate 16S
    bootstrap set (31 distinct splits, 82 conflicting pairs in total):

    ======  =========  =====  ======  ========
    cap     conflicts  boxes  verts   time
    ======  =========  =====  ======  ========
    16      1          1      18      instant
    20      11         15     34      0.1 s
    24      31         82     80      1.1 s
    31      82         474    284     80 s
    ======  =========  =====  ======  ========

    Hence the default of 20. Raise it to chase weaker conflicting signal, and
    expect both the picture and the wait to degrade; lower it for a cleaner
    figure that shows only the conflicts among the strongest splits.
    """

    def __init__(self, names: Sequence[str],
                 splits: Sequence[Tuple[frozenset, float]], *,
                 min_weight: float = 0.0, max_splits: int = 20,
                 color: str = "#37618e", width: float = 1.0,
                 tip_labels: bool = True, label_size: float = 8.0,
                 node_size: float = 5.0):
        self.names = list(names)
        universe = frozenset(self.names)
        kept = [(frozenset(s) & universe, float(w)) for s, w in splits]
        kept = [(s, w) for s, w in kept
                if 1 <= len(s) < len(universe) and w > min_weight]
        kept.sort(key=lambda sw: -sw[1])
        self.splits = self._select(kept, max_splits, universe)
        self.color = color
        self.width = width
        self.tip_labels = tip_labels
        self.label_size = label_size
        self.node_size = node_size
        self.title: Optional[str] = None
        self._groups: Optional[Dict[str, object]] = None
        self._group_title = "group"
        self._baseline = None
        self._pos: Optional[Dict[str, XY]] = None

    @staticmethod
    def _select(ranked, max_splits: int, universe: frozenset):
        """Choose which splits to draw, keeping the conflicts.

        Taking the strongest ``max_splits`` and stopping would be wrong here:
        the splits that appear in more than half the trees are exactly the
        majority-rule consensus, and a majority consensus is compatible *by
        construction*. A pure top-N cut therefore draws a tree no matter how
        reticulate the data is -- on a real 16S bootstrap set it discarded all
        82 conflicting pairs and produced a plain tree.

        So the budget is split: most of it goes to the strongest splits, and
        the rest is reserved for the strongest splits that actually conflict
        with something already kept. If the data really is tree-like nothing
        conflicts, the reserve goes unused, and the drawing is a tree because
        the data says so rather than because the selection said so.
        """
        if len(ranked) <= max_splits:
            return ranked
        core = max(1, int(max_splits * 0.7))
        kept = list(ranked[:core])
        chosen = {s for s, _ in kept}
        rest = [sw for sw in ranked[core:] if sw[0] not in chosen]

        # Sweep repeatedly rather than once: a split admitted on this pass can
        # be the very thing a later split conflicts with, and one pass would
        # miss those. Without the repeat the reserve found nothing at all on a
        # real 16S set whose full split system holds 82 conflicting pairs.
        added = True
        while added and len(kept) < max_splits:
            added = False
            for split, weight in rest:
                if len(kept) >= max_splits:
                    break
                if split in chosen:
                    continue
                if any(conflicting(split, other, universe) for other, _ in kept):
                    kept.append((split, weight))
                    chosen.add(split)
                    added = True

        # any budget left over goes back to the next-strongest splits
        for split, weight in rest:
            if len(kept) >= max_splits:
                break
            if split not in chosen:
                kept.append((split, weight))
                chosen.add(split)
        return kept

    # -- constructors ----------------------------------------------------
    @classmethod
    def from_distances(cls, names: Sequence[str], matrix, **kwargs) -> "SplitNetwork":
        from ..infer.distance import neighbor_joining
        tree = neighbor_joining(list(names), matrix)
        return cls(names, splits_from_tree(tree, names), **kwargs)

    @classmethod
    def from_alignment(cls, alignment, *, model: str = "identity",
                       **kwargs) -> "SplitNetwork":
        from ..infer.distance import distance_matrix
        from .network import _as_distance_matrix
        names, mat = _as_distance_matrix(alignment, model, distance_matrix)
        return cls.from_distances(names, mat, **kwargs)

    @classmethod
    def from_trees(cls, trees: Sequence[object], **kwargs) -> "SplitNetwork":
        """Weight each split by the fraction of ``trees`` that contain it.

        A bootstrap or posterior sample becomes a picture of exactly where the
        sample disagrees: a split found in every tree draws long, one found in
        half of them draws short and boxed against its rival.
        """
        trees = list(trees)
        if not trees:
            raise ValueError("from_trees() needs at least one tree")
        names = list(trees[0].leaf_names())
        universe = frozenset(names)
        counts: Dict[frozenset, int] = {}
        for tree in trees:
            # A split and its complement are one split, and a rooted tree
            # carries both sides as separate edges -- so canonicalise, then
            # deduplicate *within* the tree before counting. Without the
            # dedupe every balanced split is tallied twice and its weight can
            # exceed 1, which is meant to be a fraction of the trees.
            here = set()
            for side, _ in splits_from_tree(tree, names):
                here.add(min(side, universe - side,
                             key=lambda s: (len(s), sorted(s))))
            for key in here:
                counts[key] = counts.get(key, 0) + 1
        splits = [(s, c / len(trees)) for s, c in counts.items()]
        return cls(names, splits, **kwargs)

    # -- composition -----------------------------------------------------
    def color_by(self, groups, *, title: str = "group",
                 baseline=None) -> "SplitNetwork":
        self._groups = (dict(groups) if isinstance(groups, dict)
                        else dict(zip(self.names, groups)))
        self._group_title = title
        self._baseline = baseline
        return self

    def titled(self, title: str) -> "SplitNetwork":
        self.title = title
        return self

    def conflicts(self) -> List[Tuple[frozenset, frozenset]]:
        """Pairs of splits that cannot coexist on a tree -- the boxes."""
        universe = frozenset(self.names)
        out = []
        for i, (a, _) in enumerate(self.splits):
            for b, _ in self.splits[i + 1:]:
                if conflicting(a, b, universe):
                    out.append((a, b))
        return out

    # -- the network -----------------------------------------------------
    def _signatures(self) -> Dict[str, Tuple[int, ...]]:
        """Each taxon as a 0/1 vector: which side of each split it sits on."""
        return {nm: tuple(1 if nm in side else 0 for side, _ in self.splits)
                for nm in self.names}

    def _median_network(self):
        """Nodes and edges of the median network over the split system.

        Nodes are 0/1 vectors; two nodes are joined when they differ in
        exactly one coordinate, i.e. by crossing exactly one split. The vertex
        set is closed under the coordinate-wise majority of any three vectors
        -- those medians are the extra corners that turn a pair of conflicting
        splits into a box rather than a crossing.
        """
        verts = set(self._signatures().values())
        if not verts:
            return [], []
        # Closure under medians, capped. The closure of many conflicting
        # splits grows combinatorially -- and the triple loop that finds the
        # medians is cubic in the vertex count on top of that, so an
        # unbounded run on a heavily reticulate split system takes minutes to
        # produce a mesh nobody could read. The cap bounds both.
        limit = 1200
        frontier = list(verts)
        while frontier and len(verts) < limit:
            new = set()
            pool = list(verts)
            for i, a in enumerate(frontier):
                for j, b in enumerate(pool):
                    for c in pool[j + 1:]:
                        med = tuple(1 if (x + y + z) >= 2 else 0
                                    for x, y, z in zip(a, b, c))
                        if med not in verts:
                            new.add(med)
                            if len(verts) + len(new) >= limit:
                                break
                    if len(verts) + len(new) >= limit:
                        break
                if len(verts) + len(new) >= limit:
                    break
            if not new:
                break
            verts |= new
            frontier = list(new)

        verts = sorted(verts)
        index = {v: i for i, v in enumerate(verts)}
        edges = []
        for a in verts:
            for k in range(len(a)):
                b = list(a)
                b[k] = 1 - b[k]
                b = tuple(b)
                if b in index and index[a] < index[b]:
                    edges.append((index[a], index[b], k))
        return verts, edges

    @property
    def positions(self) -> Dict[str, XY]:
        """Taxon coordinates -- the split embedding.

        Every split is assigned a direction; a vertex's position is the sum of
        the displacements of the splits it lies on the far side of. Compatible
        splits then displace nested groups and the drawing stays tree-like,
        while two conflicting splits displace overlapping groups along
        different directions -- and the four corners that produces are the box.
        """
        if self._pos is None:
            verts, _ = self._median_network()
            coords = self._vertex_coords(verts)
            sig = self._signatures()
            index = {v: i for i, v in enumerate(verts)}
            # Taxa that no kept split separates share a network vertex, and
            # drawn there they would sit exactly on top of each other. Give
            # each a short pendant, fanned out around the shared vertex --
            # which is also what the vertex means: these taxa are together.
            span = self._span(coords)
            at_vertex: Dict[int, List[str]] = {}
            for nm in self.names:
                at_vertex.setdefault(index[sig[nm]], []).append(nm)
            pos = {}
            for vi, members in at_vertex.items():
                vx, vy = coords[vi]
                if len(members) == 1:
                    pos[members[0]] = (vx, vy)
                    continue
                r = 0.06 * span
                for k, nm in enumerate(sorted(members)):
                    a = 2 * math.pi * k / len(members)
                    pos[nm] = (vx + r * math.cos(a), vy + r * math.sin(a))
            self._pos = pos
        return self._pos

    @staticmethod
    def _span(coords) -> float:
        if not coords:
            return 1.0
        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        return max(max(xs) - min(xs), max(ys) - min(ys)) or 1.0

    def _split_angles(self) -> List[float]:
        """One direction per split, spread over the half-circle.

        Splits are already sorted by weight, so ordering the angles the same
        way keeps the strongest structure aligned rather than scattered.
        """
        n = len(self.splits)
        return [math.pi * i / max(n, 1) for i in range(n)]

    def _vertex_coords(self, verts) -> List[XY]:
        angles = self._split_angles()
        total = sum(w for _, w in self.splits) or 1.0
        out = []
        for v in verts:
            x = y = 0.0
            for k, bit in enumerate(v):
                if bit:
                    step = self.splits[k][1] / total
                    x += math.cos(angles[k]) * step
                    y += math.sin(angles[k]) * step
            out.append((x, y))
        return out

    def _build(self) -> RenderContext:
        verts, edges = self._median_network()
        coords = self._vertex_coords(verts)
        scene = Scene()

        for i, j, k in edges:
            weight = self.splits[k][1]
            scene.add(Path([coords[i], coords[j]], color=self.color,
                           width=self.width * (0.4 + 1.0 * min(weight, 1.0)),
                           opacity=0.7, zorder=0.6))

        pos = self.positions
        sig = self._signatures()
        index = {v: i for i, v in enumerate(verts)}
        span = self._span(coords)
        cfunc = self._node_colors(scene)
        for nm in self.names:
            x, y = pos[nm]
            vx, vy = coords[index[sig[nm]]]
            if (x, y) != (vx, vy):        # pendant edge back to its vertex
                scene.add(Path([(vx, vy), (x, y)], color=self.color,
                               width=self.width * 0.5, opacity=0.6, zorder=0.6))
            scene.add(Marker(x, y, size=self.node_size, color=cfunc(nm),
                             edgecolor=cfunc(nm), zorder=3, label=nm))
            if self.tip_labels:
                # anchor the text on the side facing its vertex so it grows
                # outward; a centred anchor puts half the name back over the
                # network, and taxa fanned around one vertex then collide
                dx, dy = x - vx, y - vy
                pad = 0.02 * span
                if abs(dx) > abs(dy) * 0.5 and (dx or dy):
                    ha = "left" if dx > 0 else "right"
                    scene.add(Label(x + (pad if dx > 0 else -pad), y, nm,
                                    size=self.label_size, color="#333333",
                                    ha=ha, va="center"))
                else:
                    va = "bottom" if dy >= 0 else "top"
                    scene.add(Label(x, y + (pad if dy >= 0 else -pad), nm,
                                    size=self.label_size, color="#333333",
                                    ha="center", va=va))

        ctx = RenderContext(_NamelessTree(self.names), _NetLayout(span))
        ctx.scene = scene
        return ctx

    def _node_colors(self, scene: Scene):
        if not self._groups:
            return lambda nm: "#37618e"
        scale = build_color_scale(self._group_title,
                                  [self._groups.get(n) for n in self.names],
                                  baseline=self._baseline)
        scene.add_legend(scale.title, scale.legend)
        scene.legend_swatch[scale.title] = scale.swatch
        return lambda nm: scale.color(self._groups.get(nm))

    def _default_figsize(self, ctx: RenderContext = None):
        return (6.0, 5.0)


class _NamelessTree:
    def __init__(self, names: Sequence[str]):
        self._names = list(names)

    @property
    def n_leaves(self) -> int:
        return len(self._names)

    def leaves(self):
        return []

    def traverse(self, order: str = "preorder"):
        return iter(())

    def nodes(self, order: str = "preorder"):
        return []
