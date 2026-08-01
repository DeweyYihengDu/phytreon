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

Implementation, in three steps, because the middle one is what makes the
picture readable rather than a tangle:

1. **A circular ordering.** The taxa are arranged around a circle so that as
   many splits as possible cut it as a single arc -- so each split becomes a
   *chord*. This is the idea Neighbor-Net turns on, and without it the split
   directions are arbitrary and the drawing crosses itself (measured on a real
   16S split system: 48 edges, 87 crossings). From a distance matrix the
   ordering is built by agglomeration on the distances themselves
   (:func:`neighbornet_ordering`); from a set of trees, by nesting the
   compatible splits and reordering their children.
2. **Splits.** From a set of trees, each split is weighted by the fraction of
   trees holding it. From a distance matrix, every split the ordering can draw
   is fitted by non-negative least squares against the distances -- the
   Neighbor-Net estimation step, without which a distance matrix yields no
   boxes at all, since splits read off one tree are compatible by
   construction.
3. **The arrangement.** Those chords cut the disc into cells; the network is
   the dual -- one node per cell, one edge per shared chord segment. Two
   chords that cross give a cell on all four of their sides, which is the box.
   Drawing each edge perpendicular to its own chord makes the result planar.

The last step is the one that has to be got right. Taking instead the median
closure of the taxon signatures -- the Buneman graph -- overshoots: for three
mutually conflicting splits it returns the whole 3-cube, eight nodes, which in
two dimensions can only be drawn as a wireframe cube with its hidden edges
crossing the visible ones. The chord arrangement returns seven cells, which is
the hexagon of three rhombi that SplitsTree draws.

Not every split system can be drawn this way, and the ones that cannot are
worth understanding rather than papering over. Four taxa admit three ways of
splitting two against two, and a circle can only show two of them -- the third
pair would have to sit opposite each other. A dataset that supports all three
resolutions of some quartet is therefore not circular under *any* ordering,
and the splits that lose out land in :attr:`SplitNetwork.dropped` rather than
being drawn wrong. On a 60-replicate 16S bootstrap set, 392 of the 3060
quartets carried all three resolutions; the seven splits that had to go held
0.8% of the weight between them, and a search from 31 independent starting
orderings found no ordering that saved any of them.
"""
from __future__ import annotations

import math
import warnings
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


def splits_from_tree(tree, names: Sequence[str], *,
                     trivial: bool = False) -> List[Tuple[frozenset, float]]:
    """Every edge of ``tree`` as ``(taxon subset, edge length)``.

    ``trivial`` includes the terminal edges -- the splits that separate one
    taxon from all the others. They say nothing about topology, since every
    tree has all of them, but they are what holds the taxa apart in the
    drawing: without them a taxon sits directly on an internal node and its
    name lands on its neighbours'.
    """
    index = set(names)
    out = []
    for node in tree.traverse():
        if node.is_root or node.parent is None:
            continue
        side = frozenset(n for n in node.leaf_names() if n in index)
        low = 0 if trivial else 1
        if low < len(side) < len(index):
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


def is_circular(side: frozenset, order: Sequence[str]) -> bool:
    """True when ``side`` is one contiguous arc of the circular ``order``."""
    n = len(order)
    marks = [nm in side for nm in order]
    # count runs around the circle; one run (of either state) means one arc
    runs = sum(1 for i in range(n) if marks[i] and not marks[i - 1])
    return runs <= 1


#: How far a chord end is nudged off the exact midpoint of its gap, as a
#: fraction of the spacing between two taxa. Chord ends otherwise sit on a
#: regular grid, which makes three chords meeting at a single point ordinary
#: rather than freakish -- three splits that each cover half the taxa are three
#: diameters through the centre -- and a cell pinched to zero area there is a
#: box that silently disappears from the drawing. The nudge only enters the
#: test for which cells exist, never the coordinates, so it cannot be seen.
_JITTER = 0.02


def _gap_angle(gap: int, n: int) -> float:
    """Angle of the gap that follows position ``gap``, in general position.

    Wrapped first: the gap before position 0 and the gap after position n-1
    are one gap, and two chords that end there must end at the *same* point or
    they cross each other for no reason.
    """
    gap %= n
    return 2 * math.pi * (gap + 0.5 + _JITTER * math.sin(3.0 * gap)) / n


def _clip(poly: Sequence[XY], p: XY, q: XY, keep: float) -> List[XY]:
    """The part of convex ``poly`` on one side of the line ``pq``."""
    dx, dy = q[0] - p[0], q[1] - p[1]

    def where(pt: XY) -> float:
        return keep * (dx * (pt[1] - p[1]) - dy * (pt[0] - p[0]))

    out: List[XY] = []
    for i, a in enumerate(poly):
        b = poly[(i + 1) % len(poly)]
        fa, fb = where(a), where(b)
        if fa >= 0:
            out.append(a)
        if (fa > 0) != (fb > 0) and fa != fb:
            t = fa / (fa - fb)
            out.append((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
    return out


def _area(poly: Sequence[XY]) -> float:
    total = 0.0
    for i, a in enumerate(poly):
        b = poly[(i + 1) % len(poly)]
        total += a[0] * b[1] - b[0] * a[1]
    return abs(total) / 2.0


def _nudge_apart(pos: Dict[str, XY], least: float,
                 rounds: int = 20) -> Dict[str, XY]:
    """Push taxa drawn nearer than ``least`` just far enough to both be seen.

    Two taxa a gene cannot separate belong at the same point, and the network
    is right to put them there -- but one marker then hides the other, and a
    reader counts one taxon where there are two. The nudge is a fraction of the
    figure's own width, on the order of the marker it is keeping visible, so it
    cannot be misread as a distance: nothing is legible at that scale anyway.
    """
    names = list(pos)
    if len(names) < 2 or least <= 0:
        return pos
    out = {nm: list(xy) for nm, xy in pos.items()}
    for _ in range(rounds):
        moved = False
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                ax, ay = out[a]
                bx, by = out[b]
                dx, dy = bx - ax, by - ay
                gap = math.hypot(dx, dy)
                if gap >= least:
                    continue
                if gap < 1e-12:            # exactly coincident: pick a way
                    ang = 2 * math.pi * (hash((a, b)) % 360) / 360.0
                    dx, dy, gap = math.cos(ang), math.sin(ang), 1.0
                push = (least - gap) / 2.0
                ux, uy = dx / gap, dy / gap
                out[a] = [ax - ux * push, ay - uy * push]
                out[b] = [bx + ux * push, by + uy * push]
                moved = True
        if not moved:
            break
    return {nm: (xy[0], xy[1]) for nm, xy in out.items()}


def _compatible_core(names, splits) -> List[frozenset]:
    """The strongest splits that all fit on one tree, heaviest first."""
    universe = frozenset(names)
    core: List[frozenset] = []
    for side, _ in sorted(splits, key=lambda sw: -sw[1]):
        side = min(side, universe - side, key=lambda s: (len(s), sorted(s)))
        if not 1 < len(side) < len(universe) or side in core:
            continue
        if all(not conflicting(side, other, universe) for other in core):
            core.append(side)
    return core


def _hierarchy(names, core):
    """Nest the compatible splits: ``{group: [children]}``, leaves are names."""
    universe = frozenset(names)
    children: Dict[frozenset, List] = {universe: []}
    # largest first, so every superset is already in place when a set arrives
    for side in sorted(core, key=lambda s: -len(s)):
        if side == universe or side in children:
            continue
        parent = min((p for p in children if side < p), key=len)
        children[parent].append(side)
        children[side] = []
    for name in names:
        deepest = min((p for p in children if name in p), key=len)
        children[deepest].append(name)
    return universe, children


def _flatten(node, children) -> List[str]:
    if isinstance(node, str):
        return [node]
    out: List[str] = []
    for child in children[node]:
        out.extend(_flatten(child, children))
    return out


def circular_ordering(names: Sequence[str],
                      splits: Sequence[Tuple[frozenset, float]]) -> List[str]:
    """An order of the taxa around a circle cut by as many splits as possible
    as a single arc.

    This is the step NeighborNet turns on, and it is what makes the drawing
    planar: once a split is one arc of the circle it is one chord, and chords
    only ever cross -- never tangle. Without an ordering the split directions
    are arbitrary and the picture crosses itself 87 times on a 48-edge network.

    Built in two stages. The compatible splits -- the strongest ones that all
    fit on a single tree -- are nested into a hierarchy, and *any* depth-first
    walk of a hierarchy already lays every one of them out as an arc, because
    each group stays contiguous whatever order its children take. That freedom
    is then spent on the rest: the children of each group are reordered to
    bring as much weight as possible from the *conflicting* splits into arcs
    too. Reordering children can never break the hierarchy's own splits, so
    the second stage can only improve on the first.
    """
    names = list(names)
    if len(names) <= 3 or not splits:
        return names
    root, children = _hierarchy(names, _compatible_core(names, splits))
    order = _flatten(root, children)

    def score(order_: Sequence[str]) -> float:
        return sum(w for side, w in splits if is_circular(side, order_))

    best = score(order)
    internal = [k for k in children if children[k]]
    for _ in range(4):
        improved = False
        for node in internal:
            kids = children[node]
            if len(kids) < 2:
                continue
            # reversal always, pair swaps while the node is small enough that
            # the quadratic sweep stays cheap
            moves = [list(reversed(kids))]
            if len(kids) <= 10:
                for i in range(len(kids)):
                    for j in range(i + 1, len(kids)):
                        trial = list(kids)
                        trial[i], trial[j] = trial[j], trial[i]
                        moves.append(trial)
            for trial in moves:
                children[node] = trial
                cand = _flatten(root, children)
                value = score(cand)
                if value > best + 1e-12:
                    best, order, kids = value, cand, trial
                    improved = True
                children[node] = kids
        if not improved:
            break
    return order


#: Above this many taxa the weight fit is skipped by default. The problem is
#: square in the number of taxon *pairs*, so it grows as the fourth power:
#: measured at 0.02 s for 30 taxa, 0.8 s for 60, 4.5 s for 80.
FIT_MAX_TAXA = 80

#: How a chain's end node divides its attention between the taxon at the end
#: and the one behind it, once the chain is too long to track node by node.
#: The end taxon is what the next link will attach to, so it carries the
#: larger share.
_END_WEIGHT = 2.0 / 3.0


class _Chain:
    """A run of taxa already placed next to each other, and its two ends.

    Only the ends can take another link, so only the ends need distances --
    ``nodes`` holds one bookkeeping node per end (one while the chain is a
    single taxon). Everything in the middle is settled and no longer consulted.
    """

    __slots__ = ("taxa", "nodes")

    def __init__(self, taxa, nodes):
        self.taxa = list(taxa)
        self.nodes = list(nodes)

    def reversed_(self) -> "_Chain":
        return _Chain(self.taxa[::-1], self.nodes[::-1])


def _mixture_distance(dist, left, right) -> float:
    """Distance between two nodes that each stand for a weighted mixture.

    A reduced node is a mixture over the nodes it replaced, so every distance
    it takes part in is the expected distance under that mixture. One rule
    covers both the node-to-node and node-to-outsider cases, which is why the
    reduction needs no special constants of its own.
    """
    return sum(wa * wb * dist[(a, b)] if a != b else 0.0
               for a, wa in left for b, wb in right)


def neighbornet_ordering(names: Sequence[str], matrix) -> List[str]:
    """The agglomerative circular ordering of Bryant and Moulton.

    Neighbour joining picks the two nodes to *merge*, and merging them away is
    exactly what costs it: from then on the pair is a single subtree and no
    later step can put anything between them. Neighbor-Net picks the two nodes
    to stand *next to each other* and merges nothing. Clusters are therefore
    chains of taxa rather than subtrees, they grow at both ends, and when one
    circle is left that circle is the ordering.

    The freedom that buys is the whole point of using it over a tree's leaf
    order. Any ordering read off a tree is constrained by the tree's splits;
    a chain can seat two taxa side by side that no single tree groups, which
    is precisely the situation a split network exists to draw.

    Selection is neighbour joining's own criterion, applied twice: once over
    the chains to decide which two to join, and again over their end nodes to
    decide which ends meet. A chain longer than two nodes is then reduced back
    to two, each new node standing for a mixture of the ones it replaces.
    """
    names = list(names)
    n = len(names)
    if n <= 3:
        return names

    dist: Dict[Tuple[int, int], float] = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                dist[(i, j)] = float(matrix[i][j])
    # what each bookkeeping node stands for, as (node, weight) pairs
    stands_for: Dict[int, List[Tuple[int, float]]] = {
        i: [(i, 1.0)] for i in range(n)}
    chains = [_Chain([nm], [i]) for i, nm in enumerate(names)]
    fresh = n

    def between(a: _Chain, b: _Chain) -> float:
        return (sum(dist[(x, y)] for x in a.nodes for y in b.nodes)
                / (len(a.nodes) * len(b.nodes)))

    # all the way down to one chain, not three: with three left their
    # orientations are still undecided, and concatenating them as they happen
    # to lie would seat three arbitrary pairs of taxa next to each other
    while len(chains) > 1:
        r = len(chains)
        gap = [[0.0] * r for _ in range(r)]
        for i in range(r):
            for j in range(i + 1, r):
                gap[i][j] = gap[j][i] = between(chains[i], chains[j])
        away = [sum(gap[i]) for i in range(r)]
        scale = max(r - 2, 1)
        pick = min((scale * gap[i][j] - away[i] - away[j], i, j)
                   for i in range(r) for j in range(i + 1, r))
        _, i, j = pick
        first, second = chains[i], chains[j]
        others = [c for k, c in enumerate(chains) if k not in (i, j)]

        # second round: the same criterion over the ends that could meet,
        # with the two chosen chains opened up into their own end nodes
        entities = ([[(x, 1.0)] for x in first.nodes]
                    + [[(y, 1.0)] for y in second.nodes]
                    + [[(x, 1.0 / len(c.nodes)) for x in c.nodes]
                       for c in others])
        m = len(entities)
        reach = [sum(_mixture_distance(dist, e, f)
                     for f in entities if f is not e) for e in entities]
        span = max(m - 2, 1)
        ends = min((span * dist[(x, y)] - reach[a] - reach[len(first.nodes) + b],
                    a, b)
                   for a, x in enumerate(first.nodes)
                   for b, y in enumerate(second.nodes))
        _, a, b = ends
        # face them so the two chosen ends touch, then join
        left = first if a == len(first.nodes) - 1 else first.reversed_()
        right = second if b == 0 else second.reversed_()
        merged = _Chain(left.taxa + right.taxa, left.nodes + right.nodes)

        if len(merged.nodes) > 2:
            # fold the inner nodes into the two ends: the end node keeps the
            # larger share, since it is what the next link attaches to
            head, tail = merged.nodes[0], merged.nodes[-1]
            inner_head, inner_tail = merged.nodes[1], merged.nodes[-2]
            new_head = [(head, _END_WEIGHT), (inner_head, 1 - _END_WEIGHT)]
            new_tail = [(tail, _END_WEIGHT), (inner_tail, 1 - _END_WEIGHT)]
            u, v = fresh, fresh + 1
            fresh += 2
            for node, parts in ((u, new_head), (v, new_tail)):
                stands_for[node] = parts
                for other in others:
                    for w in other.nodes:
                        value = _mixture_distance(dist, parts, [(w, 1.0)])
                        dist[(node, w)] = dist[(w, node)] = value
            value = _mixture_distance(dist, new_head, new_tail)
            dist[(u, v)] = dist[(v, u)] = value
            merged = _Chain(merged.taxa, [u, v])

        chains = others + [merged]

    return [nm for chain in chains for nm in chain.taxa]


def circular_splits(order: Sequence[str]) -> List[frozenset]:
    """Every split a given circular ordering can draw: all of its arcs.

    There are exactly ``n*(n-1)/2`` of them, the same as the number of pairs
    of taxa -- which is why fitting weights to them against a distance matrix
    is a square problem rather than an under-determined one.
    """
    order = list(order)
    n = len(order)
    universe = frozenset(order)
    seen, out = set(), []
    for size in range(1, n):
        for start in range(n):
            side = frozenset(order[(start + t) % n] for t in range(size))
            key = min(side, universe - side, key=lambda s: (len(s), sorted(s)))
            if key not in seen:
                seen.add(key)
                out.append(side)
    return out


def circular_split_weights(names: Sequence[str], matrix, order: Sequence[str],
                           *, tol: float = 1e-8) -> List[Tuple[frozenset, float]]:
    """Fit a weight to *every* split the ordering can draw.

    This is the estimation step of NeighborNet, and it is a different question
    from "which splits did a tree happen to contain". A tree hands over n-3
    internal splits and nothing else; here every one of the ``n*(n-1)/2``
    circular splits is a candidate, and each is given the weight that best
    explains the observed distances:

        minimise ||A w - d||   subject to   w >= 0

    where ``d`` is the pairwise distances and ``A[pair, split]`` says whether
    that split separates that pair. The non-negativity is what does the work:
    a split the data does not support is driven to exactly zero rather than to
    a small negative number, so the answer is sparse and every split that
    survives has earned its place.

    The point of running this is that it recovers conflict a tree cannot
    report. Splits taken from a tree are compatible by construction except for
    the disagreements *between* trees, so a single distance matrix yields no
    boxes at all; fitting all circular splits finds the conflicting signal
    that was in the matrix the whole time.
    """
    import numpy as np
    from scipy.optimize import nnls

    order = list(order)
    n = len(order)
    if n < 4:
        return []
    if n > FIT_MAX_TAXA:
        raise ValueError(
            "fitting all circular splits for %d taxa means a %d x %d "
            "non-negative least squares problem; measured cost is 0.8 s at 60 "
            "taxa, 4.5 s at 80, and it climbs from there. Pass estimate=False "
            "to read the splits off the tree instead."
            % (n, n * (n - 1) // 2, n * (n - 1) // 2))

    splits = circular_splits(order)
    index = {nm: i for i, nm in enumerate(names)}
    at = {nm: i for i, nm in enumerate(order)}

    member = np.zeros((len(splits), n), dtype=bool)
    for row, side in enumerate(splits):
        for nm in side:
            member[row, at[nm]] = True

    left, right, dist = [], [], []
    for a in range(n):
        for b in range(a + 1, n):
            left.append(a)
            right.append(b)
            dist.append(float(matrix[index[order[a]]][index[order[b]]]))
    design = (member[:, left] != member[:, right]).T.astype(float)

    weights, _ = nnls(design, np.asarray(dist, dtype=float))
    return [(side, float(w)) for side, w in zip(splits, weights) if w > tol]


class SplitNetwork(_Renderable):
    """A split network -- conflicting splits drawn as boxes.

    Build from a distance matrix (:meth:`from_distances`), an alignment
    (:meth:`from_alignment`), or a set of trees (:meth:`from_trees`), where
    every split's weight is how many of the trees contain it. That last one
    turns a bootstrap or posterior sample straight into a picture of which
    groupings the sample disagrees about.

    ``max_splits`` caps how many *informative* splits are drawn. Terminal
    splits -- one taxon against all the rest -- are exempt, because they can
    never conflict with anything and so can never be the point of the picture,
    while without them every taxon sits on an internal node and the names pile
    up. So this is a readability knob, not a speed one: each conflict opens a
    box, and past a couple of dozen boxes the drawing is a mesh. Measured on a
    60-replicate 16S bootstrap set (47 distinct splits, 7 of which no circular
    ordering can draw):

    ======  =========  =====  ======  ======
    cap     conflicts  boxes  nodes   time
    ======  =========  =====  ======  ======
    12      0          0      31      0.1 s
    16      2          2      37      0.2 s
    20      9          9      48      0.3 s
    28      13         13     54      0.6 s
    ======  =========  =====  ======  ======

    Hence the default of 20. Raise it to chase weaker conflicting signal,
    lower it for a cleaner figure showing only the strongest disagreements.
    Boxes and conflicts match exactly at every setting, which is the check
    that the drawing neither invents a conflict nor swallows one.

    A split that no circular ordering lays out as a single arc has no chord
    and cannot be drawn without crossings; those are listed in
    :attr:`dropped` rather than drawn wrong.
    """

    def __init__(self, names: Sequence[str],
                 splits: Sequence[Tuple[frozenset, float]], *,
                 min_weight: float = 0.0, max_splits: int = 20,
                 color: str = "#37618e", width: float = 1.0,
                 tip_labels: bool = True, label_size: float = 8.0,
                 node_size: float = 5.0, order: Optional[Sequence[str]] = None,
                 label_ring: Optional[bool] = None):
        self.names = list(names)
        universe = frozenset(self.names)
        cand = [(frozenset(s) & universe, float(w)) for s, w in splits]
        cand = [(s, w) for s, w in cand
                if 1 <= len(s) < len(universe) and w > min_weight]
        cand.sort(key=lambda sw: -sw[1])
        # A split and its complement are one split, and drawn as two they are
        # two chords on top of each other with no room between them -- which
        # loses the cells that should sit there. A rooted tree hands over both
        # sides of its root edge, so this is the ordinary case, not a rare one.
        seen, unique = set(), []
        for side, weight in cand:
            key = min(side, universe - side, key=lambda s: (len(s), sorted(s)))
            if key not in seen:
                seen.add(key)
                unique.append((side, weight))
        cand = unique
        self.order = list(order) if order else circular_ordering(self.names, cand)
        # Only a split that is one arc of the ordering is a chord, and only
        # chords can be drawn without crossings. The rest are set aside rather
        # than drawn wrong -- see :attr:`dropped`.
        arcs = [is_circular(s, self.order) for s, _ in cand]
        self.dropped = [s for (s, _), ok in zip(cand, arcs) if not ok]
        drawable = [sw for sw, ok in zip(cand, arcs) if ok]
        # Terminal splits -- one taxon against the rest -- are exempt from the
        # budget. They can never conflict with anything, so they can never be
        # the interesting part of the picture, and they cost one node each;
        # spending the budget on them would only crowd out real conflicts.
        terminal = [sw for sw in drawable if min(len(sw[0]),
                                                 len(universe) - len(sw[0])) < 2]
        rest = [sw for sw in drawable if sw not in terminal]
        self.splits = self._select(rest, max_splits, universe) + terminal
        self.color = color
        self.width = width
        self.tip_labels = tip_labels
        self.label_size = label_size
        self.label_ring = label_ring
        self.node_size = node_size
        self.title: Optional[str] = None
        self._groups: Optional[Dict[str, object]] = None
        self._group_title = "group"
        self._baseline = None
        self._pos: Optional[Dict[str, XY]] = None
        #: True when the split weights came from the circular fit rather than
        #: from a tree -- see :meth:`from_distances`.
        self.estimated = False

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
    def from_distances(cls, names: Sequence[str], matrix, *,
                       estimate="auto", ordering: str = "neighbornet",
                       **kwargs) -> "SplitNetwork":
        """Build from a distance matrix. By default, this is Neighbor-Net.

        Two steps, and both are needed. The taxa are put in a circular
        ordering by agglomeration on the distances themselves
        (:func:`neighbornet_ordering`), and then *every* split that ordering
        can draw is given the weight that best explains the distances, subject
        to being non-negative. Splits the data does not support come back at
        exactly zero and disappear.

        ``ordering="tree"`` takes the ordering from a neighbour-joining tree
        instead, which is cheaper and worse: measured on 40 distance matrices
        built from known circular split systems, the agglomeration made every
        split drawable in 40 cases out of 40, while the tree's leaf order
        managed 3 and left a fifth of the split weight undrawable. A tree's
        leaf order can only ever respect the tree's own splits.

        That second step is the whole point. Splits read off a single tree are
        compatible with one another by construction, so taking them and
        stopping produces no boxes at all -- the "network" is a tree however
        conflicted the data is. Measured on the 18-taxon 16S distance matrix:
        33 splits and **0 boxes** from the tree, against 40 splits and **20
        boxes** from the fit, reproducing the distances to 4.6%. The conflict
        was in the matrix the whole time; only the fit can report it.

        ``estimate`` is ``True``, ``False``, or ``"auto"`` (the default), which
        fits up to :data:`FIT_MAX_TAXA` taxa and reads the tree above that,
        warning as it does so -- a drawing with no boxes in it is a claim about
        the data, and it should not be made on the quiet because the fit was
        skipped. For more taxa than that, reach for :meth:`from_trees` and a
        bootstrap set, which carries conflict without needing the fit.
        Whichever route ran is recorded in :attr:`estimated`.
        """
        from ..infer.distance import neighbor_joining
        names = list(names)
        if ordering not in ("neighbornet", "tree"):
            raise ValueError("ordering must be 'neighbornet' or 'tree', "
                             "not %r" % (ordering,))
        splits = None
        if ordering == "tree" or estimate is False:
            tree = neighbor_joining(names, matrix)
            splits = splits_from_tree(tree, names, trivial=True)
        if estimate == "auto":
            estimate = len(names) <= FIT_MAX_TAXA
            if not estimate:
                warnings.warn(
                    "%d taxa is past the %d the circular split fit is worth "
                    "running for, so the splits come from the tree instead -- "
                    "and splits from one tree are compatible with each other, "
                    "so this drawing cannot show a box however conflicted the "
                    "data is. Use from_trees() with a bootstrap sample, or "
                    "pass estimate=True and wait."
                    % (len(names), FIT_MAX_TAXA), stacklevel=2)
        if not estimate:
            if splits is None:
                tree = neighbor_joining(names, matrix)
                splits = splits_from_tree(tree, names, trivial=True)
            net = cls(names, splits, **kwargs)
            net.estimated = False
            return net
        order = kwargs.pop("order", None)
        if order is None:
            order = (neighbornet_ordering(names, matrix) if ordering == "neighbornet"
                     else circular_ordering(names, splits))
        net = cls(names, circular_split_weights(names, matrix, order),
                  order=order, **kwargs)
        net.estimated = True
        return net

    @classmethod
    def neighbor_net(cls, names: Sequence[str], matrix,
                     **kwargs) -> "SplitNetwork":
        """Neighbor-Net, by name: the agglomerative ordering plus the fit."""
        kwargs.setdefault("estimate", True)
        kwargs.setdefault("ordering", "neighbornet")
        return cls.from_distances(names, matrix, **kwargs)

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
            for side, _ in splits_from_tree(tree, names, trivial=True):
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

    def _arcs(self) -> List[Tuple[int, int]]:
        """Each split as ``(first, last)`` positions of its arc in the order."""
        n = len(self.order)
        at = {nm: i for i, nm in enumerate(self.order)}
        out = []
        for side, _ in self.splits:
            idx = sorted(at[nm] for nm in side if nm in at)
            if not idx:
                out.append((0, 0))
                continue
            # the arc may straddle the end of the list; the largest step
            # between consecutive members is the gap *outside* it
            gaps = [(idx[(k + 1) % len(idx)] - idx[k]) % n
                    for k in range(len(idx))]
            g = gaps.index(max(gaps))
            out.append((idx[(g + 1) % len(idx)], idx[g]))
        return out

    def _chords(self):
        """Each split as a chord of the unit circle: two ends and a normal.

        The taxa sit at angles ``2*pi*i/n``; a split that covers positions
        ``first..last`` is cut off by the chord joining the two *gaps* that
        bound it, and the direction that separates its two sides is that
        chord's normal, pointing at the middle of the arc.
        """
        n = len(self.order)
        out = []
        for first, last in self._arcs():
            a = _gap_angle(first - 1, n)
            b = _gap_angle(last, n)
            span = (last - first) % n
            theta = 2 * math.pi * ((first + span / 2.0) % n) / n
            out.append(((math.cos(a), math.sin(a)),
                        (math.cos(b), math.sin(b)), theta))
        return out

    def _network(self):
        """Nodes and edges of the split network: the dual of the chords.

        The chords cut the disc into cells. Each cell lies on a definite side
        of every chord, so a cell *is* a 0/1 vector over the splits; two cells
        that differ in one coordinate share a piece of that chord and are
        joined by an edge. Two chords that cross leave a cell on each of their
        four sides, and those four cells are the box.

        Cells are found by walking outwards from the ones the taxa occupy,
        flipping one coordinate at a time and keeping the vector when the
        corresponding region of the disc is really non-empty. Each side of a
        chord is a convex piece of the disc, so a candidate cell is just an
        intersection of convex pieces and the test is a polygon clip.
        """
        m = len(self.splits)
        if not m:
            return [], []
        chords = self._chords()
        n = len(self.order)
        # A 2n-gon whose corners land exactly on the taxa *and* on the gaps,
        # so every chord runs corner to corner and nothing is clipped askew.
        disc = []
        for i in range(n):
            for angle in (2 * math.pi * i / n, _gap_angle(i, n)):
                disc.append((math.cos(angle), math.sin(angle)))
        sides = []
        for (p, q, theta) in chords:
            inside = (math.cos(theta), math.sin(theta))
            dx, dy = q[0] - p[0], q[1] - p[1]
            ref = dx * (inside[1] - p[1]) - dy * (inside[0] - p[0])
            sides.append((p, q, 1.0 if ref >= 0 else -1.0))

        def feasible(vec) -> bool:
            poly = disc
            for k, bit in enumerate(vec):
                p, q, sign = sides[k]
                poly = _clip(poly, p, q, sign if bit else -sign)
                if len(poly) < 3:
                    return False
            return _area(poly) > 1e-9

        verts = set(self._signatures().values())
        frontier = list(verts)
        while frontier:
            nxt = []
            for v in frontier:
                for k in range(m):
                    w = list(v)
                    w[k] = 1 - w[k]
                    w = tuple(w)
                    if w not in verts and feasible(w):
                        verts.add(w)
                        nxt.append(w)
            frontier = nxt

        verts = sorted(verts)
        index = {v: i for i, v in enumerate(verts)}
        edges = []
        for a in verts:
            for k in range(m):
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
            verts, _ = self._network()
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
            outward = self._outward()
            pos = {}
            for vi, members in at_vertex.items():
                vx, vy = coords[vi]
                if len(members) == 1:
                    pos[members[0]] = (vx, vy)
                    continue
                # Fan them out, but around the directions the circular
                # ordering already gives them, and in that order: taxa no
                # split separates are neighbours on the circle, so a fan
                # centred on where they belong keeps the picture's rotational
                # order intact instead of scrambling it alphabetically.
                rank = {nm: i for i, nm in enumerate(self.order)}
                members.sort(key=lambda nm: rank.get(nm, 0))
                mean = math.atan2(sum(outward[nm][1] for nm in members),
                                  sum(outward[nm][0] for nm in members))
                step = min(math.radians(55), 2 * math.pi * 0.85 / len(members))
                first = mean - step * (len(members) - 1) / 2.0
                r = 0.06 * span
                for k, nm in enumerate(members):
                    a = first + step * k
                    pos[nm] = (vx + r * math.cos(a), vy + r * math.sin(a))
            self._pos = _nudge_apart(pos, 0.025 * span)
        return self._pos

    def _ring_anchors(self, coords, pos) -> Optional[Dict[str, XY]]:
        """Where each name goes when the labels are moved out to a ring.

        Past a dozen or so taxa the names stop fitting beside their nodes: the
        drawing bunches taxa wherever the splits between them are short, and a
        name three centimetres long then lies across its neighbour's node. So
        the names move out to a circle and a hairline leader points back.

        The circle's angles come from the circular ordering, evenly spaced, so
        no two names can ever land on top of each other. The leaders cannot
        tangle either: the drawing keeps the taxa in that same rotational
        order around its outer face, so a node and its label are always in
        step with their neighbours.
        """
        want = self.label_ring
        if want is None:
            want = len(self.names) >= 12
        if not want or len(self.names) < 3:
            return None
        pts = list(coords) + list(pos.values())
        cx = (max(p[0] for p in pts) + min(p[0] for p in pts)) / 2.0
        cy = (max(p[1] for p in pts) + min(p[1] for p in pts)) / 2.0
        radius = max(math.hypot(x - cx, y - cy) for x, y in pos.values())
        radius = (radius or self._span(coords)) * 1.08
        out = {}
        for nm, (ux, uy) in self._outward().items():
            out[nm] = (cx + radius * ux, cy + radius * uy)
        return out

    def _outward(self) -> Dict[str, XY]:
        """Which way is *away from the network* for each taxon.

        The circular ordering already answers this: taxon ``i`` of ``n`` sits
        at angle ``2*pi*i/n`` on the circle the chords were cut from, and the
        drawing preserves that rotational order around its outer face. Using
        it means labels radiate outward and cannot pile up on one side --
        which is what happens if every solitary taxon simply puts its name
        above its node.
        """
        n = len(self.order)
        out = {}
        for i, nm in enumerate(self.order):
            a = 2 * math.pi * i / n
            out[nm] = (math.cos(a), math.sin(a))
        return out

    @staticmethod
    def _span(coords) -> float:
        if not coords:
            return 1.0
        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        return max(max(xs) - min(xs), max(ys) - min(ys)) or 1.0

    def _split_angles(self) -> List[float]:
        """One direction per split: the normal of the chord it cuts."""
        return [theta for _, _, theta in self._chords()]

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
        verts, edges = self._network()
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
        outward = self._outward()
        ring = self._ring_anchors(coords, pos)
        cfunc = self._node_colors(scene)
        for nm in self.names:
            x, y = pos[nm]
            vx, vy = coords[index[sig[nm]]]
            if (x, y) != (vx, vy):        # pendant edge back to its vertex
                scene.add(Path([(vx, vy), (x, y)], color=self.color,
                               width=self.width * 0.5, opacity=0.6, zorder=0.6))
            scene.add(Marker(x, y, size=self.node_size, color=cfunc(nm),
                             edgecolor=cfunc(nm), zorder=3, label=nm))
            if self.tip_labels and ring:
                scene.add(Path([(x, y), ring[nm]], color="#b8bec7",
                               width=self.width * 0.4, opacity=0.8, zorder=0.4))
                x, y = ring[nm]
                vx, vy = x - outward[nm][0], y - outward[nm][1]
            if self.tip_labels:
                # Anchor the text on the side that faces away from the
                # network, so it grows outward: a centred anchor puts half the
                # name back over the drawing, and neighbouring taxa collide.
                # A taxon fanned off a shared vertex already points outward;
                # a solitary one takes the direction the circular ordering
                # gives it, rather than defaulting to straight up.
                dx, dy = x - vx, y - vy
                if not (dx or dy):
                    dx, dy = outward[nm]
                norm = math.hypot(dx, dy) or 1.0
                ux, uy = dx / norm, dy / norm
                pad = 0.02 * span
                # Grow the name left or right on the faintest horizontal hint.
                # Two taxa side by side near the top of the drawing point
                # almost straight up; centre both names over their nodes and
                # they overlap, while sending one left and the other right
                # separates them by the whole width of both.
                ha = "center" if abs(ux) < 0.15 else ("left" if ux > 0
                                                      else "right")
                va = "center" if abs(uy) < 0.5 else ("bottom" if uy > 0
                                                     else "top")
                scene.add(Label(x + pad * ux, y + pad * uy, nm,
                                size=self.label_size, color="#333333",
                                ha=ha, va=va))

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


def neighbor_net(names: Sequence[str], matrix, **kwargs) -> SplitNetwork:
    """Neighbor-Net: a split network straight from a distance matrix.

    The taxa are placed in a circular ordering by agglomeration on the
    distances, then every split that ordering can draw is fitted to those
    distances by non-negative least squares. What comes back is planar, and
    its boxes are the conflict the matrix contains -- which no tree drawn from
    the same matrix can show, since a tree's splits are compatible with one
    another by construction.

        net = pt.neighbor_net(names, matrix)
        net.color_by(groups).save("network.pdf")
    """
    return SplitNetwork.neighbor_net(names, matrix, **kwargs)
