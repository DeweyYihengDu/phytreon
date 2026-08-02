"""Sequence-similarity networks -- the CLANS-style cluster map.

When a protein family is too divergent for a trustworthy alignment, a
phylogenetic tree stops being honest: branch order is then an artefact of
alignment error rather than a record of descent. The standard alternative is
to drop the tree and draw the *sequence space* itself -- one node per
sequence, an edge wherever a pairwise search finds significant similarity,
laid out by a force-directed algorithm so that groups of mutually similar
sequences fall into visible clusters::

    net = pt.SequenceNetwork.from_alignment(aln, cutoff=0.4)
    net.color_by(family_dict)
    net.save("clusters.pdf")

This is what CLANS does (Frickey & Lupas 2004), using the Fruchterman-Reingold
layout; the picture it produces is a familiar sight in comparative-genomics
papers where a family's internal structure is real but its deep branching
order is not recoverable.

Reading one: a tight ball of nodes is a group whose members all detect each
other; a thin bridge of edges between two balls means the groups are related
but only distantly; an isolated node found nothing above the cutoff. Distance
on the page is *not* an evolutionary distance -- it is only the layout's
compromise between many pairwise attractions, so read clusters and
connections, never a ruler.
"""
from __future__ import annotations

import math
import random
from typing import Dict, List, Optional, Sequence, Tuple

from ..scene import MIN_STROKE_PT, Label, Marker, Path, Scene
from .figure import RenderContext, _Renderable, build_color_scale

XY = Tuple[float, float]


class _NetworkLayout:
    """Layout shim so a network scene can drive the ordinary backends."""

    is_polar = False
    equal_aspect = True          # a force layout has no meaningful aspect
    invert_y = False
    kind = "rect"
    use_branch_lengths = False

    def __init__(self, max_x: float = 1.0):
        self.max_x = max_x

    @staticmethod
    def _collapsed_span(node, use_len: bool):
        return (0.0, 0.0)


def fruchterman_reingold(nodes: Sequence[str], edges: Sequence[Tuple[int, int, float]],
                         *, iterations: int = 300, seed: int = 0,
                         k: Optional[float] = None,
                         gravity: float = 0.8) -> List[XY]:
    """Force-directed layout: repulsion between all nodes, attraction along
    edges, with a temperature that cools each round (Fruchterman & Reingold
    1991 -- the algorithm CLANS uses).

    ``edges`` are ``(i, j, weight)``; heavier weights pull harder, so a
    stronger similarity ends up shorter. Returns one ``(x, y)`` per node.

    ``gravity`` pulls everything gently toward the centroid. Without it, a
    graph with disconnected pieces blows apart: nothing but repulsion acts
    between the pieces, so isolated nodes drift off, and since they then set
    the scale of the picture the connected core collapses into an unreadable
    dot. Sequence networks are disconnected by construction -- that is the
    whole point of a cutoff -- so this matters here more than in the textbook
    algorithm.

    It is a genuine trade-off, measured on a real 16S graph and a synthetic
    three-family one: raising gravity from 0.25 to 0.8 grows the main cluster
    from 16% to 27% of the frame while cluster separation only falls from
    4.0x to 3.6x. Past that, separation degrades faster than legibility
    improves. Raise it if isolated sequences push your clusters into a corner;
    lower it if distinct families are being squeezed together.
    """
    n = len(nodes)
    if n == 0:
        return []
    if n == 1:
        return [(0.0, 0.0)]

    rng = random.Random(seed)
    # start on a circle plus jitter: a pure random cloud often folds the graph
    # over itself and the layout then has to spend iterations untangling it
    pos = []
    for i in range(n):
        a = 2 * math.pi * i / n
        pos.append([math.cos(a) + rng.uniform(-0.01, 0.01),
                    math.sin(a) + rng.uniform(-0.01, 0.01)])

    area = 1.0
    k = k if k is not None else math.sqrt(area / n)   # ideal edge length
    temp = 0.1
    cool = temp / (iterations + 1)

    for _ in range(iterations):
        disp = [[0.0, 0.0] for _ in range(n)]

        # repulsion: every pair pushes apart, weaker with distance
        for i in range(n):
            xi, yi = pos[i]
            for j in range(i + 1, n):
                dx = xi - pos[j][0]
                dy = yi - pos[j][1]
                d2 = dx * dx + dy * dy
                if d2 < 1e-12:                 # coincident: nudge apart
                    dx, dy = rng.uniform(-1e-3, 1e-3), rng.uniform(-1e-3, 1e-3)
                    d2 = dx * dx + dy * dy
                d = math.sqrt(d2)
                force = k * k / d
                ux, uy = dx / d, dy / d
                disp[i][0] += ux * force
                disp[i][1] += uy * force
                disp[j][0] -= ux * force
                disp[j][1] -= uy * force

        # attraction: only along edges, scaled by the similarity weight
        for i, j, w in edges:
            dx = pos[i][0] - pos[j][0]
            dy = pos[i][1] - pos[j][1]
            d = math.hypot(dx, dy)
            if d < 1e-12:
                continue
            force = (d * d / k) * w
            ux, uy = dx / d, dy / d
            disp[i][0] -= ux * force
            disp[i][1] -= uy * force
            disp[j][0] += ux * force
            disp[j][1] += uy * force

        # gravity toward the centroid, so disconnected pieces stay in frame
        if gravity:
            gx = sum(p[0] for p in pos) / n
            gy = sum(p[1] for p in pos) / n
            for i in range(n):
                disp[i][0] -= (pos[i][0] - gx) * gravity * k
                disp[i][1] -= (pos[i][1] - gy) * gravity * k

        # step, capped by the current temperature so late rounds settle
        for i in range(n):
            dx, dy = disp[i]
            d = math.hypot(dx, dy)
            if d > 1e-12:
                step = min(d, temp)
                pos[i][0] += dx / d * step
                pos[i][1] += dy / d * step
        temp -= cool

    # Scale uniformly about the centroid so the result fills a unit frame.
    # The absolute scale carries no meaning (see the module docstring: page
    # distance is not an evolutionary distance), but it does decide whether a
    # dense component renders as a readable cloud or an unresolvable dot, and
    # that depends on node count in a way no fixed figure size can anticipate.
    # Uniform, not per-axis: stretching x and y independently would distort
    # the very cluster shapes the layout exists to show.
    cx = sum(p[0] for p in pos) / n
    cy = sum(p[1] for p in pos) / n
    reach = max(math.hypot(p[0] - cx, p[1] - cy) for p in pos)
    if reach > 1e-12:
        for p in pos:
            p[0] = (p[0] - cx) / reach
            p[1] = (p[1] - cy) / reach
    return [(p[0], p[1]) for p in pos]


def _components(n: int, edges) -> List[List[int]]:
    """Connected node indices, largest group first."""
    adj: Dict[int, List[int]] = {i: [] for i in range(n)}
    for i, j, _ in edges:
        adj[i].append(j)
        adj[j].append(i)
    seen = [False] * n
    out = []
    for start in range(n):
        if seen[start]:
            continue
        seen[start] = True
        stack, group = [start], [start]
        while stack:
            v = stack.pop()
            for w in adj[v]:
                if not seen[w]:
                    seen[w] = True
                    stack.append(w)
                    group.append(w)
        out.append(sorted(group))
    out.sort(key=len, reverse=True)
    return out


def _pack_discs(radii: Sequence[float]) -> List[XY]:
    """Centres for discs of the given radii: biggest in the middle, rest in
    rings around it. Sizes come in largest-first."""
    centres = [(0.0, 0.0)]
    edge = radii[0]
    i = 1
    while i < len(radii):
        widest = max(radii[i:])
        ring = edge + widest * 1.15
        angle = 0.0
        placed = 0
        while i < len(radii):
            step = 2 * math.asin(min(1.0, radii[i] * 1.15 / ring))
            if placed and angle + step > 2 * math.pi:
                break
            centres.append((ring * math.cos(angle + step / 2),
                            ring * math.sin(angle + step / 2)))
            angle += step
            i += 1
            placed += 1
        edge = ring + widest
    return centres


def layout_by_component(nodes: Sequence[str], edges, **kwargs) -> List[XY]:
    """Lay each connected piece out on its own, then pack the pieces.

    Running one force layout over the whole graph lets the pieces set each
    other's scale, and for a sequence network -- disconnected by construction,
    because that is what a cutoff does -- the result is decided by whichever
    isolated sequences drift furthest. Measured on a 106-sequence 16S graph at
    identity 0.78: eight components, and the one holding 90 of the 106
    sequences came out occupying 33% of the frame's width, with half of all
    nodes inside 13% of the radius. The picture was mostly empty paper with
    the data in one corner.

    Laid out separately and packed, each piece is sized by its own node count
    (area proportional to it, so node density is comparable between pieces)
    and the big one takes the room it deserves. A graph in one piece is
    untouched -- there is nothing to pack.
    """
    groups = _components(len(nodes), edges)
    if len(groups) <= 1:
        return fruchterman_reingold(nodes, edges, **kwargs)

    laid, radii = [], []
    for group in groups:
        local = {g: i for i, g in enumerate(group)}
        sub = [(local[i], local[j], w) for i, j, w in edges
               if i in local and j in local]
        pts = fruchterman_reingold([nodes[g] for g in group], sub, **kwargs)
        # area with the node count, so a big cluster is not drawn at the same
        # size as a lone sequence and every piece keeps a comparable density
        radius = math.sqrt(len(group))
        laid.append([(x * radius, y * radius) for x, y in pts])
        radii.append(radius)

    out: List[XY] = [(0.0, 0.0)] * len(nodes)
    for group, pts, (cx, cy) in zip(groups, laid, _pack_discs(radii)):
        for g, (x, y) in zip(group, pts):
            out[g] = (cx + x, cy + y)

    # Centre on the bounding box, not on the origin the packer happened to
    # start from. Packing puts the biggest piece at the origin and rings the
    # rest around it, so a graph with one big piece and a few stragglers ends
    # up lopsided: normalising by distance-from-origin then sizes the frame
    # for the far side and leaves the near side as margin.
    xs = [x for x, _ in out]
    ys = [y for _, y in out]
    cx = (max(xs) + min(xs)) / 2.0
    cy = (max(ys) + min(ys)) / 2.0
    reach = max(math.hypot(x - cx, y - cy) for x, y in out) or 1.0
    return [((x - cx) / reach, (y - cy) / reach) for x, y in out]


class SequenceNetwork(_Renderable):
    """A CLANS-style sequence-similarity cluster map.

    Build one from an alignment (:meth:`from_alignment`), a distance matrix
    (:meth:`from_distances`), or an explicit edge list (the constructor).
    ``cutoff`` is the similarity below which no edge is drawn -- the knob that
    decides how much of the sequence space you actually see, exactly like
    CLANS's E-value slider.
    """

    def __init__(self, names: Sequence[str],
                 edges: Sequence[Tuple[int, int, float]], *,
                 iterations: int = 300, seed: int = 0,
                 node_size: float = 3.2, node_color: str = "#37618e",
                 edge_color: str = "#9aa3ad", edge_width: float = 0.3,
                 edge_alpha: float = 0.45, weight_edges: bool = True):
        self.names = list(names)
        self.edges = [(int(i), int(j), float(w)) for i, j, w in edges]
        self.iterations = iterations
        self.seed = seed
        self.node_size = node_size
        self.node_color = node_color
        self.edge_color = edge_color
        self.edge_width = edge_width
        self.edge_alpha = edge_alpha
        #: draw stronger similarities as darker/thicker lines
        self.weight_edges = weight_edges
        self.title: Optional[str] = None
        self._groups: Optional[Dict[str, object]] = None
        self._group_title = "group"
        self._baseline = None
        self._order = None
        self._cluster_labels: List[Tuple[str, List[str]]] = []
        self._label_names: List[str] = []
        self._label_size = 8.0
        self._pos: Optional[List[XY]] = None

    # -- constructors ----------------------------------------------------
    @classmethod
    def from_distances(cls, names: Sequence[str], matrix, *,
                       cutoff: float = 0.5, **kwargs) -> "SequenceNetwork":
        """Build from a symmetric distance matrix.

        Distances are turned into similarities as ``1 - d`` and an edge is kept
        where that exceeds ``cutoff``.
        """
        n = len(names)
        edges = []
        for i in range(n):
            row = matrix[i]
            for j in range(i + 1, n):
                sim = 1.0 - float(row[j])
                if sim > cutoff:
                    edges.append((i, j, sim))
        return cls(names, edges, **kwargs)

    @classmethod
    def from_alignment(cls, alignment, *, cutoff: float = 0.4,
                       model: str = "identity", **kwargs) -> "SequenceNetwork":
        """Build from an alignment by all-against-all pairwise identity.

        This mirrors what CLANS gets from all-against-all BLAST, at the
        resolution an alignment can give: it needs the sequences to be
        alignable in the first place, so for a family too divergent to align
        globally, supply your own edges from a real search instead.
        """
        from ..infer.distance import distance_matrix
        names, mat = _as_distance_matrix(alignment, model, distance_matrix)
        return cls.from_distances(names, mat, cutoff=cutoff, **kwargs)

    @classmethod
    def from_pairs(cls, pairs, *, names: Optional[Sequence[str]] = None,
                   **kwargs) -> "SequenceNetwork":
        """Build from ``(name1, name2, similarity)`` triples -- e.g. parsed
        straight out of a BLAST tabular report."""
        rows = [tuple(p) for p in pairs]
        if names is None:
            seen: List[str] = []
            for a, b, *_ in rows:
                for nm in (str(a), str(b)):
                    if nm not in seen:
                        seen.append(nm)
            names = seen
        index = {nm: i for i, nm in enumerate(names)}
        missing = {str(nm) for a, b, *_ in rows for nm in (a, b)
                   if str(nm) not in index}
        if missing:
            raise ValueError(
                f"pairs reference names not in `names`: {sorted(missing)[:5]}")
        edges = [(index[str(a)], index[str(b)],
                  float(rest[0]) if rest else 1.0) for a, b, *rest in rows]
        return cls(names, edges, **kwargs)

    # -- composition -----------------------------------------------------
    def color_by(self, groups, *, title: str = "group", baseline=None,
                 order=None) -> "SequenceNetwork":
        """Colour nodes by a mapping ``{name: group}`` (or a list parallel to
        ``names``). ``baseline`` greys out the levels you want in the
        background; ``order`` fixes the legend order."""
        if isinstance(groups, dict):
            self._groups = dict(groups)
        else:
            self._groups = dict(zip(self.names, groups))
        self._group_title = title
        self._baseline = baseline
        self._order = order
        return self

    def label_clusters(self, labels: Dict[str, Sequence[str]]) -> "SequenceNetwork":
        """Write a name beside each cluster, as CLANS figures do.

        ``labels`` maps a display name to the sequences it covers; the text is
        placed just outside that group's centroid.
        """
        self._cluster_labels = [(k, list(v)) for k, v in labels.items()]
        return self

    def label_nodes(self, names: Optional[Sequence[str]] = None, *,
                    size: float = 8.0) -> "SequenceNetwork":
        """Label individual nodes (all of them, or just the ones named)."""
        self._label_names = list(names) if names is not None else list(self.names)
        self._label_size = size
        return self

    def titled(self, title: str) -> "SequenceNetwork":
        self.title = title
        return self

    # -- diagnostics -----------------------------------------------------
    @property
    def positions(self) -> List[XY]:
        """Node coordinates, computing the layout on first use."""
        if self._pos is None:
            self._pos = layout_by_component(
                self.names, self.edges,
                iterations=self.iterations, seed=self.seed)
        return self._pos

    def components(self) -> List[List[str]]:
        """Connected components, largest first -- the clusters an eye would
        pick out, computed rather than guessed at."""
        n = len(self.names)
        adj: Dict[int, List[int]] = {i: [] for i in range(n)}
        for i, j, _ in self.edges:
            adj[i].append(j)
            adj[j].append(i)
        seen = [False] * n
        out = []
        for start in range(n):
            if seen[start]:
                continue
            stack, group = [start], []
            seen[start] = True
            while stack:
                cur = stack.pop()
                group.append(self.names[cur])
                for nxt in adj[cur]:
                    if not seen[nxt]:
                        seen[nxt] = True
                        stack.append(nxt)
            out.append(group)
        out.sort(key=len, reverse=True)
        return out

    # -- building --------------------------------------------------------
    def _build(self) -> RenderContext:
        pos = self.positions
        scene = Scene()

        weights = [w for _, _, w in self.edges]
        wmin = min(weights) if weights else 0.0
        wmax = max(weights) if weights else 1.0
        wspan = (wmax - wmin) or 1.0

        # A stronger hit reads as a firmer line. The range it is drawn over
        # *starts* at the thinnest stroke that prints, rather than being scaled
        # freely and then clipped there: clipping would give every edge below
        # the floor the same width, so the weakest hits -- the ones a cutoff
        # exists to include -- would all come out looking identical. Mapping
        # into the range keeps them apart and printable at the same time, and
        # the fading is carried by opacity, which the press reproduces at any
        # value.
        lo = MIN_STROKE_PT
        hi = max(lo * 1.6, self.edge_width * 1.5)

        for i, j, w in self.edges:
            if self.weight_edges:
                frac = (w - wmin) / wspan
                width = lo + (hi - lo) * frac
                alpha = self.edge_alpha * (0.45 + 0.55 * frac)
            else:
                width, alpha = max(self.edge_width, lo), self.edge_alpha
            scene.add(Path([pos[i], pos[j]], color=self.edge_color,
                           width=width, opacity=alpha, zorder=0.5))

        cfunc = self._node_color_function(scene)
        for idx, name in enumerate(self.names):
            x, y = pos[idx]
            scene.add(Marker(x, y, size=self.node_size, color=cfunc(name),
                             edgecolor=cfunc(name), zorder=3, label=name))

        self._add_cluster_labels(scene, pos)
        for name in self._label_names:
            if name in self.names:
                x, y = pos[self.names.index(name)]
                scene.add(Label(x, y + 0.012, name, size=self._label_size,
                                color="#333333", ha="center", va="bottom"))

        xs = [p[0] for p in pos] or [0.0]
        span = (max(xs) - min(xs)) or 1.0
        ctx = RenderContext(_NamelessTree(self.names), _NetworkLayout(span))
        ctx.scene = scene
        return ctx

    def _node_color_function(self, scene: Scene):
        if not self._groups:
            return lambda name: self.node_color
        values = [self._groups.get(nm) for nm in self.names]
        scale = build_color_scale(self._group_title, values,
                                  baseline=self._baseline, order=self._order,
                                  swatch="point")
        scene.add_legend(scale.title, scale.legend)
        scene.legend_swatch[scale.title] = scale.swatch
        return lambda name: scale.color(self._groups.get(name))

    def _add_cluster_labels(self, scene: Scene, pos: Sequence[XY]) -> None:
        if not self._cluster_labels:
            return
        index = {nm: i for i, nm in enumerate(self.names)}
        cx = sum(p[0] for p in pos) / len(pos)
        cy = sum(p[1] for p in pos) / len(pos)
        xs = [p[0] for p in pos]
        ys = [p[1] for p in pos]
        span = max(max(xs) - min(xs), max(ys) - min(ys)) or 1.0
        for text, members in self._cluster_labels:
            pts = [pos[index[m]] for m in members if m in index]
            if not pts:
                continue
            mx = sum(p[0] for p in pts) / len(pts)
            my = sum(p[1] for p in pts) / len(pts)
            # push the label away from the figure's centre so it sits outside
            # its own cluster instead of on top of the nodes
            dx, dy = mx - cx, my - cy
            d = math.hypot(dx, dy)
            if d < 1e-9:            # cluster sits on the centroid: go up
                dx, dy, d = 0.0, 1.0, 1.0
            reach = max(math.hypot(p[0] - mx, p[1] - my) for p in pts)
            off = reach + 0.06 * span
            lx, ly = mx + dx / d * off, my + dy / d * off
            # anchor the text on the side facing the cluster, so it grows
            # outward instead of back over the nodes. A centred anchor would
            # put half the string inside the cluster, and how much that is
            # depends on the text length -- which is not knowable here.
            if abs(dx) > abs(dy) * 0.5:
                ha = "left" if dx > 0 else "right"
            else:
                ha = "center"
            va = "bottom" if (ha == "center" and dy > 0) else \
                 ("top" if ha == "center" else "center")
            scene.add(Label(lx, ly, text, size=10.0, color="#333333",
                            ha=ha, va=va))

    def _default_figsize(self, ctx: RenderContext = None):
        return (8.0, 7.0)


class _NamelessTree:
    """Minimal stand-in so ``RenderContext`` has the handful of tree
    attributes it touches; a network has no tree behind it."""

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


def _as_distance_matrix(alignment, model, distance_matrix):
    """Accept a phytreon ``Alignment`` or a Biopython alignment."""
    if hasattr(alignment, "names") and hasattr(alignment, "seqs"):
        names = list(alignment.names)
        seqs = list(alignment.seqs)
        n = len(names)
        mat = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                a, b = seqs[i], seqs[j]
                comparable = sum(1 for x, y in zip(a, b)
                                 if x not in "-." and y not in "-.")
                same = sum(1 for x, y in zip(a, b)
                           if x == y and x not in "-.")
                d = 1.0 - (same / comparable if comparable else 0.0)
                mat[i][j] = mat[j][i] = d
        return names, mat
    return distance_matrix(alignment, model)
