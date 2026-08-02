"""Visual elements -- the drawing vocabulary of a :class:`TreeFigure`.

Every element is a small :class:`~phytreon.plot.figure._Element` that reads
node coordinates from the layout and appends primitives to the scene.
Elements branch on ``ctx.layout.is_polar`` where a circular tree needs
different geometry (rotated labels, wedge highlights, arc clade bars).

Elements are not used directly; they are added through the fluent
``TreeFigure`` methods (``.tip_labels()``, ``.heatmap()`` …).
"""
from __future__ import annotations

import math
import warnings
from typing import Optional, Sequence

from ..core.tree import Node
from ..scene import Label, Marker, Path, Polygon, Raster
from .figure import _Element, RenderContext, build_color_scale, is_numeric

def readable_on(fill) -> str:
    """Black or white, whichever actually reads on ``fill``.

    A name printed inside its own coloured block was always white, which is
    right on a dark block and close to invisible on a pale one -- measured at
    2.2:1 against some of the palette, where 3:1 is the floor for large text
    and 4.5:1 for body text. Picking by the background's luminance keeps every
    block above 4.5:1 whatever colour it lands on.
    """
    from matplotlib.colors import to_rgba

    def channel(c):
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4

    r, g, b, _ = to_rgba(fill)
    lum = 0.2126 * channel(r) + 0.7152 * channel(g) + 0.0722 * channel(b)
    # Compare the two inks actually used, not white against an idealised black:
    # a near-black at #1a1a1a is itself luminous enough to lose a comparison it
    # appears to win, and the label then comes out at 3.9:1 on a mid blue.
    # Against true black the worst any fill can force is 4.58:1, at luminance
    # 0.179 where the two choices cross.
    return "#ffffff" if 1.05 / (lum + 0.05) >= (lum + 0.05) / 0.05 else "#000000"


#: How much of the fan opening a ring's name needs before it is worth drawing.
#: The name lies along the spoke, so what has to fit is its *height* in angle,
#: and that falls as the radius grows -- which is why this is measured rather
#: than derived. Counting collisions against the ring sectors on a 16S tree
#: while sweeping the tip count: the name lands on data up to about two degrees
#: of free wedge and is clear from roughly two and a half. Asking for *smaller*
#: names does not lower it, because a figure drawn with smaller names is itself
#: smaller by the same rule and the two cancel; asking for larger ones does
#: raise it. Measured at 6, 8, 10 and 12 pt.
_RING_NAME_MIN_GAP = math.radians(2.5)

#: above this many tips, per-tip cells are too thin to carry a separator
#: stroke -- see :class:`_Ring` / :class:`_Heatmap`
_RING_DENSE_TIPS = 150


# --------------------------------------------------------------------------
# tree skeleton
# --------------------------------------------------------------------------
class _Branches(_Element):
    def __init__(self, color="black", size: float = 1.0):
        self.color = color
        self.size = size

    def apply(self, ctx: RenderContext) -> None:
        nodes = ctx.tree.nodes()
        cfunc, scale = ctx.resolve_color(self.color, nodes, default="black")
        lay = ctx.layout
        for node in nodes:
            col = cfunc(node)
            bp = lay.branch_path(node)
            if bp:
                ctx.scene.add(Path(bp, color=col, width=self.size, zorder=1))
            conn = lay.child_connector(node)
            if conn:
                ctx.scene.add(Path(conn, color=cfunc(node), width=self.size, zorder=1))
        if scale is not None:
            ctx.add_scale(scale)


# --------------------------------------------------------------------------
# tip / node labels
# --------------------------------------------------------------------------
#: Past this many names an unrooted tree cannot keep them apart. Measured by
#: counting real glyph collisions on a 106-taxon 16S tree: none up to about 20
#: tips, 8 at 25, 21 at 40, 81 at 106 -- and growing the canvas does not fix
#: it, since the collisions are between taxa the layout puts at the same point.
_RADIAL_LABEL_LIMIT = 25


class _TipLabels(_Element):
    def __init__(self, color="black", size: float = 10.0,
                 offset: Optional[float] = None, italic: bool = False,
                 align: bool = False, max_labels: Optional[int] = None):
        self.color = color
        self.size = size
        self.offset = offset
        self.italic = italic
        self.align = align
        self.max_labels = max_labels       # show ~this many evenly-spaced labels

    def apply(self, ctx: RenderContext) -> None:
        lay = ctx.layout
        tips = ctx.tree.leaves()
        cfunc, scale = ctx.resolve_color(self.color, tips, default="black")
        kind = getattr(lay, "kind", "rect")
        # thin labels on large trees so they do not overlap
        step = 1
        if self.max_labels and len(tips) > self.max_labels:
            step = math.ceil(len(tips) / self.max_labels)
        # polar/radial labels sit at the tip radius, so they need a larger gap
        # to clear the tip marker (markers are sized in points, not data units)
        default_off = (0.06 if kind in ("polar", "radial") else 0.02) * lay.max_x
        off = self.offset if self.offset is not None else default_off
        # middle of the drawing, so an unrooted tree can point its names away
        # from it rather than along whichever way each branch happens to run
        cx = sum(t.x for t in tips) / len(tips) if tips else 0.0
        cy = sum(t.y for t in tips) / len(tips) if tips else 0.0
        drawn = 0
        for i, tip in enumerate(tips):
            text = tip.name or ""
            if not text or (step > 1 and i % step != 0):
                continue
            drawn += 1
            # a collapsed clade is drawn as a triangle reaching out to its
            # farthest hidden leaf; its label has to clear that, not sit on it.
            # Ask the layout so the span is in the units it draws in.
            far = lay._collapsed_span(
                tip, lay.use_branch_lengths and ctx.tree.has_branch_lengths)[1]
            if kind == "polar" and getattr(lay, "inward", False):
                # tips point toward the centre: label sits further inward
                a = tip._angle
                x, y = lay._polar_to_xy(tip._r - off - far, a)
                deg = math.degrees(a)
                if 90 < (deg % 360) < 270:
                    rot, ha = deg, "left"
                else:
                    rot, ha = deg + 180, "right"
                ctx.scene.add(Label(x, y, text, size=self.size, color=cfunc(tip),
                                    ha=ha, va="center", rotation=rot, italic=self.italic))
            elif kind == "polar":
                # sit outside any ring tracks (ctx.outer_radius) when present
                rings = ctx.outer_radius > ctx.ring_base
                r = (ctx.outer_radius if (self.align or rings)
                     else tip._r + far) + off
                a = tip._angle
                x, y = lay._polar_to_xy(r, a)
                deg = math.degrees(a)
                if 90 < (deg % 360) < 270:
                    rot, ha = deg + 180, "right"
                else:
                    rot, ha = deg, "left"
                ctx.scene.add(Label(x, y, text, size=self.size, color=cfunc(tip),
                                    ha=ha, va="center", rotation=rot, italic=self.italic))
            elif kind == "dendrogram":
                # tips along the bottom; labels drop below, rotated upright
                ctx.scene.add(Label(tip.x, tip.y - off, text, size=self.size,
                                    color=cfunc(tip), ha="right", va="center",
                                    rotation=90, italic=self.italic))
            elif kind == "radial":
                # Point the name away from the middle of the tree, not along
                # its own branch. Following the branch is the prettier idea and
                # it is what this used to do, but a terminal branch can point
                # anywhere -- including back across the drawing -- so names on
                # neighbouring tips ran parallel and straight through each
                # other. Measured on a 106-taxon 16S tree, counting real glyph
                # collisions: 164 with the branch direction, 84 pointing
                # outward, at the same canvas size.
                a = math.atan2(tip.y - cy, tip.x - cx) if (tip.x, tip.y) != (cx, cy) \
                    else tip._angle
                x = tip.x + (off + far) * math.cos(a)
                y = tip.y + (off + far) * math.sin(a)
                deg = math.degrees(a)
                if 90 < (deg % 360) < 270:
                    rot, ha = deg + 180, "right"
                else:
                    rot, ha = deg, "left"
                ctx.scene.add(Label(x, y, text, size=self.size, color=cfunc(tip),
                                    ha=ha, va="center", rotation=rot, italic=self.italic))
            else:
                x = (lay.max_x if self.align else tip.x + far) + off
                ctx.scene.add(Label(x, tip.y, text, size=self.size, color=cfunc(tip),
                                    ha="left", va="center", italic=self.italic,
                                    role="tiplab"))
        if drawn:
            ctx.tip_label_load = (drawn, self.size)
        if kind == "radial" and drawn > _RADIAL_LABEL_LIMIT:
            warnings.warn(
                "an unrooted tree cannot seat %d names without them running "
                "into each other -- it spaces its tips by branch geometry, so "
                "two taxa a gene cannot separate land on the same spot and so "
                "do their names, at any figure size. Counting real glyph "
                "collisions on a 106-taxon 16S tree: none up to about 20 tips, "
                "8 at 25, 21 at 40, 81 at 106. Thin them with "
                "tip_labels(max_labels=...), or use layout='circular', which "
                "seats every name on its own arc and stays clean."
                % drawn, stacklevel=2)
        if scale is not None:
            ctx.add_scale(scale)


class _NodeLabels(_Element):
    """Label internal nodes -- by default their support values.

    ``attr`` may name several keys, which are then printed as one combined
    label (``"88/95/1.0"``). A tree whose topology was checked by more than one
    method carries more than one support value, and reporting them together on
    the branch is how those papers present it -- drawing them as separate
    labels would just stack them on the same point.

    ``min_value`` hides weakly supported nodes, so a dense tree is not covered
    in numbers that say only "this branch is unreliable".

    ``stack=True`` writes the values on separate lines instead of joining them,
    and ``prefixes`` labels each line (``["p", "b", "n"]`` -> ``p:1.00`` /
    ``b:100`` / ``n:0.98``), which is the other common way these trees are
    annotated -- readable with four values where a slash-joined string is not.
    """

    def __init__(self, attr="support", size: float = 7.0,
                 color="#666666", offset: float = 0.0, fmt: str = "{:g}",
                 sep: str = "/", min_value: Optional[float] = None,
                 stack: bool = False, prefixes: Optional[Sequence[str]] = None):
        #: one key, or several to combine into a single label
        self.attrs = [attr] if isinstance(attr, str) else list(attr)
        self.size = size
        self.color = color
        self.offset = offset
        self.fmt = fmt
        self.sep = sep
        self.min_value = min_value
        self.stack = stack
        if prefixes is not None and len(prefixes) != len(self.attrs):
            raise ValueError(
                f"got {len(prefixes)} prefixes for {len(self.attrs)} keys")
        self.prefixes = list(prefixes) if prefixes else None

    @property
    def attr(self):
        """The first key -- kept for callers that set a single attribute."""
        return self.attrs[0]

    def _value(self, node, key):
        val = getattr(node, key, None)
        if val is None:
            val = node.data.get(key)
        return val

    def apply(self, ctx: RenderContext) -> None:
        lay = ctx.layout
        for node in ctx.tree.traverse():
            if node.is_leaf or node.is_root:
                continue
            values = [self._value(node, key) for key in self.attrs]
            if all(v is None for v in values):
                continue
            if self.min_value is not None:
                numeric = [v for v in values if is_numeric(v)]
                # keep the node only if something about it is worth reading
                if numeric and max(numeric) < self.min_value:
                    continue
            parts = []
            for i, val in enumerate(values):
                if val is None:
                    txt = "-"
                elif isinstance(val, (int, float)):
                    txt = self.fmt.format(val)
                else:
                    txt = str(val)
                if self.prefixes:
                    txt = f"{self.prefixes[i]}:{txt}"
                parts.append(txt)
            text = "\n".join(parts) if self.stack else self.sep.join(parts)
            if getattr(lay, "kind", "rect") != "rect":
                x, y, ha, va = node.x, node.y, "center", "bottom"
            else:
                # sit just above the branch leading into the node (no overlap
                # with the vertical connector or the node itself)
                x = 0.5 * (node.parent.x + node.x) + self.offset
                y = node.y - 0.3
                ha, va = "center", "center"
                if self.stack and len(parts) > 1:
                    # a stack grows downward from its anchor, so lift it clear
                    # of the branch instead of straddling it
                    va = "bottom"
                    y = node.y - 0.15
            ctx.scene.add(Label(x, y, text, size=self.size, color=self.color,
                                ha=ha, va=va))


# --------------------------------------------------------------------------
# points
# --------------------------------------------------------------------------
class _Points(_Element):
    def __init__(self, which: str = "tip", color="black", size=6.0,
                 marker: str = "o", shape=None, edgecolor: Optional[str] = None,
                 palette: str = "curated", cmap=None, baseline=None,
                 order=None):
        self.which = which                 # tip | node | all
        self.color = color
        self.size = size
        self.marker = marker
        self.shape = shape                 # categorical column -> marker shape
        self.edgecolor = edgecolor
        self.palette = palette
        self.cmap = cmap
        self.baseline = baseline
        self.order = order

    def _select(self, ctx):
        if self.which == "tip":
            return ctx.tree.leaves()
        if self.which == "node":
            return [n for n in ctx.tree.traverse() if not n.is_leaf]
        return ctx.tree.nodes()

    def apply(self, ctx: RenderContext) -> None:
        nodes = self._select(ctx)
        cfunc, cscale = ctx.resolve_color(self.color, nodes, default="black",
                                          palette=self.palette, cmap=self.cmap,
                                          baseline=self.baseline,
                                          order=self.order)
        sfunc, _ = _resolve_size(self.size, nodes)
        shfunc, shleg = ctx.resolve_shape(self.shape, nodes, default=self.marker)
        for n in nodes:
            hover = n.name or None
            ctx.scene.add(Marker(n.x, n.y, size=sfunc(n), color=cfunc(n),
                                 marker=shfunc(n), edgecolor=self.edgecolor,
                                 label=hover, zorder=3))
        if cscale is not None:
            ctx.add_scale(cscale)
        if shleg is not None:
            ctx.scene.add_legend(shleg[0], shleg[1])


def _resolve_size(spec, nodes):
    if isinstance(spec, str) and any(spec in n.data for n in nodes):
        vals = [n.data.get(spec) for n in nodes if is_numeric(n.data.get(spec))]
        lo, hi = (min(vals), max(vals)) if vals else (0, 1)
        rng = (hi - lo) or 1.0

        def f(n, _lo=lo, _rng=rng):
            v = n.data.get(spec)
            if not is_numeric(v):
                return 6.0
            return 4.0 + 12.0 * (v - _lo) / _rng
        return f, True
    return (lambda n: float(spec)), False


# --------------------------------------------------------------------------
# clade highlight
# --------------------------------------------------------------------------
class _Highlight(_Element):
    """Shade the band (or wedge) a clade occupies, behind the branches.

    Three ways to say which clade: ``node``, ``taxa`` (their common ancestor),
    or ``by`` -- the name of a joined column, which shades *every* group in it,
    each in its own colour, with a legend. One call for the whole figure rather
    than one per clade, and a reader can tell which shade means what.

    ``reach`` is how far out the colour goes, as a fraction of the distance
    from the band's inner edge to where the tip labels actually end: ``1.0``
    (the default) covers the whole name, ``0.7`` stops seven tenths of the way,
    ``0.5`` halfway. It is resolved while rendering, because how much room a
    name takes depends on the font and the figure size rather than on the tree
    -- so the same number gives the same *proportion* whatever the figure.

    A group whose taxa are not monophyletic is drawn as several bands, one per
    run of adjacent tips, rather than as one band over their common ancestor.
    The ancestor of a scattered group reaches down over other groups' taxa, so
    a single band there would colour in tips that do not belong to it -- it
    would look tidier and say something false. Several bands say what is
    actually there.
    """

    def __init__(self, node: Optional[Node] = None, taxa=None, by=None,
                 fill="#fdbf6f", alpha: float = 0.3, extend: float = 0.0,
                 palette: str = "curated", order=None, baseline=None,
                 span: str = "aligned", reach: float = 1.0):
        if span not in ("clade", "aligned", "full"):
            raise ValueError("span must be 'clade', 'aligned' or 'full', "
                             "not %r" % (span,))
        if not 0.0 < reach <= 1.5:
            raise ValueError("reach is a fraction of the way out to the end of "
                             "the tip labels; got %r" % (reach,))
        self.span = span
        self.reach = reach
        self.node = node
        self.taxa = taxa
        self.by = by
        self.fill = fill
        self.alpha = alpha
        self.extend = extend
        self.palette = palette
        self.order = order
        self.baseline = baseline

    def _target(self, ctx) -> Optional[Node]:
        if self.node is not None:
            return self.node
        if self.taxa is not None:
            return ctx.tree.get_mrca(self.taxa)
        return None

    @staticmethod
    def _runs(tips):
        """Split tips into groups adjacent in the drawn order."""
        tips = sorted(tips, key=lambda t: t.data["_row"])
        out, run = [], [tips[0]]
        for prev, cur in zip(tips, tips[1:]):
            if cur.data["_row"] - prev.data["_row"] == 1:
                run.append(cur)
            else:
                out.append(run)
                run = [cur]
        out.append(run)
        return out

    def apply(self, ctx: RenderContext) -> None:
        if self.by is not None:
            self._apply_by_column(ctx)
            return
        node = self._target(ctx)
        if node is not None:
            leaves = node.get_leaves()
            own = self._own_start(leaves, node)
            self._band(ctx, leaves, self.fill, node,
                       self._left_edge(own, None))

    def _apply_by_column(self, ctx: RenderContext) -> None:
        tips = ctx.tree.leaves()
        if not ctx.is_data_column(self.by, tips):
            raise ValueError(
                "highlight(by=%r): no tip carries that column -- join it on "
                "first with tree.join_data(df, on='name')" % (self.by,))
        scale = build_color_scale(self.by, [t.data.get(self.by) for t in tips],
                                 palette=self.palette, order=self.order,
                                 baseline=self.baseline, swatch="patch")
        groups: dict = {}
        for tip in tips:
            value = tip.data.get(self.by)
            if value is not None:
                groups.setdefault(value, []).append(tip)
        plan = []
        for value, members in groups.items():
            for run in self._runs(members):
                node = (ctx.tree.get_mrca([t.name for t in run])
                        if len(run) > 1 else run[0])
                plan.append((run, node, value, self._own_start(run, node)))
        flush = min((own for *_, own in plan), default=None)
        for run, node, value, own in plan:
            self._band(ctx, run, scale.color(value), node,
                       self._left_edge(own, flush))
        ctx.add_scale(scale)
        ctx.scene.legend_swatch[scale.title] = scale.swatch

    @staticmethod
    def _own_start(leaves, node) -> float:
        """Where this clade itself begins, measured from the root."""
        return node.x if len(leaves) > 1 else (
            node.parent.x if node.parent else node.x)

    def _left_edge(self, own: float, flush: Optional[float]) -> float:
        """The x a band starts at, under the chosen ``span``.

        ``"clade"`` starts at the group's own common ancestor, so each left
        edge is a fact about that clade -- and the edges come out ragged,
        because clades begin at different depths.

        ``"aligned"`` (the default) starts every band at the shallowest of
        those. Flush, and no band reaches further rootward than the drawing
        already did somewhere, so the deep backbone stays outside the colour
        and remains readable as a backbone.

        ``"full"`` runs from the root. Also flush, and honest for a different
        reason: a band then reads as *these rows*, which is exactly what it
        covers -- iTOL's coloured ranges are drawn this way. The cost is that
        the trunk shared with every other group is under colour too, so the
        deep structure stops standing out.

        Aligning needs to know every group at once, so it only applies to
        ``by=``; a lone ``node=``/``taxa=`` band has nothing to line up with
        and hugs its clade.
        """
        if self.span == "full":
            return 0.0
        if self.span == "aligned" and flush is not None:
            return flush
        return own

    def _band(self, ctx, leaves, fill, node, x0: float) -> None:
        lay = ctx.layout
        if lay.is_polar:
            angles = [lf._angle for lf in leaves]
            a0, a1 = min(angles), max(angles)
            da = (a1 - a0) / max(len(leaves) - 1, 1) / 2 + 1e-9
            if len(leaves) == 1:                      # no span to divide
                da = lay.extent / max(ctx.tree.n_leaves - 1, 1) / 2
            inner = lay.inner_radius + x0
            outer = lay.inner_radius + lay.max_x + self.extend
            pts = lay._arc(outer, a0 - da, a1 + da)
            pts += lay._arc(inner, a1 + da, a0 - da)
            ctx.scene.add(Polygon(pts, facecolor=fill, edgecolor=None,
                                  alpha=self.alpha, zorder=0,
                                  reach=self.reach))
        else:
            rows = [lf.y for lf in leaves]
            x1 = lay.max_x + self.extend
            y0, y1 = min(rows) - 0.45, max(rows) + 0.45
            pts = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
            ctx.scene.add(Polygon(pts, facecolor=fill, edgecolor=None,
                                  alpha=self.alpha, zorder=0, rounded=True,
                                  reach=self.reach))


# --------------------------------------------------------------------------
# clade label (bracket + text)
# --------------------------------------------------------------------------
class _CladeLabel(_Element):
    def __init__(self, label: str, node: Optional[Node] = None, taxa=None,
                 offset: float = 0.0, color="black", size: float = 11.0,
                 barsize: float = 2.0):
        self.label = label
        self.node = node
        self.taxa = taxa
        self.offset = offset
        self.color = color
        self.size = size
        self.barsize = barsize

    def apply(self, ctx: RenderContext) -> None:
        node = self.node or (ctx.tree.get_mrca(self.taxa) if self.taxa else None)
        if node is None:
            return
        lay = ctx.layout
        leaves = node.get_leaves()
        pad = self.offset + 0.04 * lay.max_x
        if lay.is_polar:
            angles = [lf._angle for lf in leaves]
            r = lay.inner_radius + lay.max_x + pad
            pts = lay._arc(r, min(angles), max(angles))
            ctx.scene.add(Path(pts, color=self.color, width=self.barsize, zorder=2))
            amid = 0.5 * (min(angles) + max(angles))
            x, y = lay._polar_to_xy(r + 0.03 * lay.max_x, amid)
            ctx.scene.add(Label(x, y, self.label, size=self.size, color=self.color,
                                ha="center", va="center",
                                rotation=math.degrees(amid)))
        else:
            rows = [lf.y for lf in leaves]
            gap = self.offset + 0.02 * lay.max_x
            x = lay.max_x + gap
            # align=True -> renderer shifts these to just past the tip labels
            ctx.scene.add(Path([(x, min(rows) + 0.1), (x, max(rows) - 0.1)],
                               color=self.color, width=self.barsize, zorder=2,
                               align=True))
            # horizontal label to the right of the bar
            ctx.scene.add(Label(x + 0.02 * lay.max_x, 0.5 * (min(rows) + max(rows)),
                                self.label, size=self.size, color=self.color,
                                ha="left", va="center", rotation=0, align=True))


# --------------------------------------------------------------------------
# heatmap alongside the tree
# --------------------------------------------------------------------------
class _Heatmap(_Element):
    """Draw a matrix of coloured cells aligned to the tips.

    ``data`` is a :class:`pandas.DataFrame` indexed by tip name (or having
    a column matching tip names).  Each column gets its own colour scale by
    default (``shared_scale=True`` for one scale across all columns).
    Rectangular layouts only.

    ``separators`` controls the hairline between cells: ``None`` (default)
    turns it off automatically past ~150 tips, where the stroke would be wider
    than the row itself and the block would read as stripes rather than
    colour. Force it with ``True``/``False``.
    """

    def __init__(self, data, offset: float = 0.0, width: float = 0.4,
                 cmap=None, palette: str = "curated", shared_scale: bool = False,
                 colnames: bool = True, colname_size: float = 9.0,
                 cell_gap: float = 0.05, separators: Optional[bool] = None,
                 baseline=None, order=None):
        self.data = _index_by_name(data)
        self.offset = offset
        self.width = width
        self.cmap = cmap
        self.palette = palette
        self.shared_scale = shared_scale   # one scale for all columns vs per-column
        self.colnames = colnames
        self.colname_size = colname_size
        self.cell_gap = cell_gap
        self.separators = separators   # None = auto (off once rows are thin)
        self.baseline = baseline
        self.order = order

    def apply(self, ctx: RenderContext) -> None:
        lay = ctx.layout
        if lay.is_polar:
            raise NotImplementedError(
                "heatmap() supports rectangular layouts; "
                "use ring() for circular trees."
            )
        df = self.data
        tips = {t.name: t for t in ctx.tree.leaves()}
        cols = list(df.columns)
        ncol = len(cols)
        total_w = self.width * lay.max_x
        cell_w = total_w / ncol
        # start at the running track cursor so multiple tracks stack rightward;
        # align=True -> renderer shifts the whole block past the tip labels
        x0 = ctx.track_cursor + (self.offset + 0.02) * lay.max_x

        # scales: one shared across all columns, or one per column (default).
        if self.shared_scale:
            flat = [df.iloc[i][c] for i in range(len(df)) for c in cols]
            shared = build_color_scale("value", flat, cmap=self.cmap,
                                       palette=self.palette,
                                       baseline=self.baseline,
                                       order=self.order, swatch="patch")
            scales = {c: shared for c in cols}
        else:
            scales = {c: build_color_scale(str(c), list(df[c]), cmap=self.cmap,
                                           palette=self.palette,
                                           baseline=self.baseline,
                                           order=self.order, swatch="patch")
                      for c in cols}

        # on a tall tree each row is thinner than its own separator stroke, so
        # the block reads as stripes rather than colour; drop the stroke and
        # let rows meet (same reasoning as the circular rings)
        dense = (len(tips) > _RING_DENSE_TIPS) if self.separators is None \
            else not self.separators
        for _, row in df.iterrows():
            name = str(row.name)
            tip = tips.get(name)
            if tip is None:
                continue
            y0, y1 = tip.y - 0.5, tip.y + 0.5
            for j, c in enumerate(cols):
                val = row[c]
                cx0 = x0 + j * cell_w
                cx1 = cx0 + cell_w - self.cell_gap * cell_w
                pts = [(cx0, y0), (cx1, y0), (cx1, y1), (cx0, y1)]
                fc = scales[c].color(val)
                ctx.scene.add(Polygon(pts, facecolor=fc,
                                      edgecolor=fc if dense else "white",
                                      width=0.4 if dense else 0.3, alpha=1.0,
                                      zorder=2, label=f"{name} | {c}: {val}",
                                      align=True))

        if self.colnames:
            ymax = max(t.y for t in ctx.tree.leaves())
            for j, c in enumerate(cols):
                cx = x0 + (j + 0.5) * cell_w
                ctx.scene.add(Label(cx, ymax + 0.8, str(c), size=self.colname_size,
                                    ha="right", va="top", rotation=45, align=True))
        for c in (["value"] if self.shared_scale else cols):
            ctx.add_scale(scales[cols[0]] if self.shared_scale else scales[c])
        ctx.track_cursor = x0 + total_w


# --------------------------------------------------------------------------
# rectangular bar track (aligned bars to the right of the tree)
# --------------------------------------------------------------------------
class _BarTrack(_Element):
    """Horizontal bar chart aligned to the tips, to the right of the tree.

    One numeric ``column`` -> one bar per tip (length encodes the value).
    Stacks after any earlier tracks; rectangular layouts only.
    """

    def __init__(self, data, column, width: float = 0.4, offset: float = 0.04,
                 fill: str = "#5b7897", bar_height: float = 0.8, colname=True,
                 colname_size: float = 9.0, axis: bool = True):
        self.data = _index_by_name(data)
        self.column = column
        self.width = width
        self.offset = offset
        self.fill = fill
        self.bar_height = bar_height
        self.colname = colname
        self.colname_size = colname_size
        self.axis = axis

    def apply(self, ctx: RenderContext) -> None:
        lay = ctx.layout
        if lay.is_polar:
            raise NotImplementedError("bar_track() is for rectangular layouts; "
                                      "use ring(geom='bar') for circular.")
        tips = ctx.tree.leaves()
        vals = [float(self.data.loc[t.name, self.column]) for t in tips
                if t.name in self.data.index]
        vmin, vmax = min(vals + [0.0]), max(vals + [0.0])
        rng = (vmax - vmin) or 1.0
        total_w = self.width * lay.max_x
        x0 = ctx.track_cursor + (self.offset + 0.02) * lay.max_x
        h = self.bar_height / 2.0
        for t in tips:
            if t.name not in self.data.index:
                continue
            v = float(self.data.loc[t.name, self.column])
            blen = (v - vmin) / rng * total_w
            y0, y1 = t.y - h, t.y + h
            pts = [(x0, y0), (x0 + blen, y0), (x0 + blen, y1), (x0, y1)]
            ctx.scene.add(Polygon(pts, facecolor=self.fill, edgecolor=None,
                                  alpha=1.0, zorder=2, align=True,
                                  label=f"{t.name} | {self.column}: {v:g}"))
        if self.axis:
            # Without a scale the bars are decoration: the reader can compare
            # two of them and read nothing off either. It matters most exactly
            # where the bars look alike -- 16S lengths run 1238 to 1584, and
            # against the zero baseline a bar chart owes them, every bar is
            # between 78% and 100% of full width. The axis is what says the
            # bars start at zero and where they end, so "nearly equal" reads
            # as the finding it is rather than as a drawing that failed.
            # the rectangular layout counts rows downward, so "above the
            # bars" is the *smaller* y and the labels anchor from their
            # bottom edge, or they hang back down over the rule they caption
            ylo = min(t.y for t in tips) - 0.6 - h
            ctx.scene.add(Path([(x0, ylo), (x0 + total_w, ylo)],
                               color="#555555", width=0.6, zorder=2,
                               align=True))
            for frac, value, side in ((0.0, vmin, "left"), (1.0, vmax, "right")):
                x = x0 + frac * total_w
                ctx.scene.add(Path([(x, ylo), (x, ylo + 0.22)], color="#555555",
                                   width=0.6, zorder=2, align=True))
                ctx.scene.add(Label(x, ylo - 0.12, ("%g" % value),
                                    size=self.colname_size - 2, color="#555555",
                                    ha=side, va="bottom", align=True))
        if self.colname:
            ymax = max(t.y for t in tips)
            ctx.scene.add(Label(x0 + total_w / 2, ymax + 0.8, str(self.column),
                                size=self.colname_size, ha="right", va="top",
                                rotation=45, align=True))
        ctx.track_cursor = x0 + total_w


# --------------------------------------------------------------------------
# multiple sequence alignment track
# --------------------------------------------------------------------------
NUC_COLORS = {"A": "#33a02c", "C": "#1f78b4", "G": "#ff7f00", "T": "#e31a1c",
              "U": "#e31a1c", "-": "#ffffff", "N": "#d9d9d9"}


def _residue_palette(seqs):
    chars = set("".join(seqs).upper())
    if chars <= set("ACGTUN-."):
        base = dict(NUC_COLORS)
    else:                                    # protein: distinct hues per residue
        from .palettes import hue_palette
        aa = sorted(c for c in chars if c not in "-.")
        cols = hue_palette(len(aa)) if aa else []
        base = {c: cols[i] for i, c in enumerate(aa)}
    base.setdefault("-", "#ffffff")
    base.setdefault(".", "#ffffff")
    return base


class _Alignment(_Element):
    """Render a multiple sequence alignment as a residue-coloured track.

    ``alignment`` is a :class:`~phytreon.infer.align.Alignment`, a FASTA
    path/string, or ``{name: aligned_seq}``.  Drawn as one raster (fast even
    for thousands of columns) aligned to the tip rows, right of the tree.
    """

    def __init__(self, alignment, width: float = 1.0, offset: float = 0.05,
                 colors=None, window=None):
        self.alignment = alignment
        self.width = width
        self.offset = offset
        self.colors = colors
        self.window = window                 # (start, end) column slice, optional

    def _seqmap(self):
        aln = self.alignment
        if isinstance(aln, dict):
            return dict(aln)
        if isinstance(aln, str):
            from .. infer.align import read_fasta
            return dict(read_fasta(aln))
        return dict(zip(aln.names, aln.seqs))   # Alignment

    def apply(self, ctx: RenderContext) -> None:
        import numpy as np
        lay = ctx.layout
        if lay.is_polar:
            raise NotImplementedError("alignment() is for rectangular layouts.")
        seqmap = self._seqmap()
        tips = ctx.tree.leaves()
        ncol = max((len(s) for s in seqmap.values()), default=0)
        if ncol == 0:
            return
        lo, hi = (self.window or (0, ncol))
        cols = hi - lo
        colormap = self.colors or _residue_palette(list(seqmap.values()))
        # ordered palette + char->code map (unknown/gap -> the '-' code = white)
        chars = sorted(colormap)
        code = {c: i for i, c in enumerate(chars)}
        palette = [colormap[c] for c in chars]
        gap_code = code.get("-", 0)

        codes = np.full((len(tips), cols), gap_code, dtype=np.int16)
        for i, t in enumerate(tips):
            s = seqmap.get(t.name, "")
            for j in range(lo, min(hi, len(s))):
                codes[i, j - lo] = code.get(s[j].upper(), gap_code)

        total_w = self.width * lay.max_x
        x0 = ctx.track_cursor + (self.offset + 0.02) * lay.max_x
        ys = [t.y for t in tips]
        ctx.scene.add(Raster(codes, palette, x0, x0 + total_w,
                             min(ys) - 0.5, max(ys) + 0.5, zorder=2, align=True))
        ctx.track_cursor = x0 + total_w


# --------------------------------------------------------------------------
# concentric metadata rings around a circular tree
# --------------------------------------------------------------------------
class _Ring(_Element):
    """Draw metadata as concentric coloured rings outside a circular tree.

    ``data`` is a DataFrame indexed by tip name (or carrying a ``name``
    column); each chosen column becomes one ring of per-tip sectors, coloured
    by that column's own scale (categorical palette or continuous cmap).
    Rings stack outward with a gap, every column gets its own legend, and tip
    labels are pushed outside all rings -- so nothing overlaps.

    Customisable: ``columns`` (which/what order), ``width``/``gap``/``offset``
    (radial geometry, fractions of the tree radius), ``pad_angle`` (gap between
    sectors, degrees), ``palette``/``cmap`` (per type), ``colnames``.

    Column names are written into the fan opening, along the spoke, so they
    never sit on top of the data they name. That opening has to be wide enough
    to hold them: the first and last sectors each hang half a sector into it,
    so on a tree with few tips -- where sectors are wide -- there may be
    nothing left. Widen it with the layout's ``extent`` (a 60-tip tree with two
    rings wants about 345 rather than the default 350), or turn the names off
    with ``colnames=False`` and rely on the legend.

    ``separators`` controls the hairline between neighbouring sectors: ``None``
    (default) turns it off automatically past ~150 tips, where the stroke would
    be wider than the sector itself and the ring would read as a comb of
    slivers instead of solid colour bands. Force it with ``True``/``False``.

    ``leaders=True`` draws a faint dotted guide from each tip out to the first
    ring. On a phylogram the tips sit at very different radii, so most stop
    well short of the rings and it stops being obvious which sector belongs to
    which tip; the alternative is to drop branch lengths entirely
    (``TreeFigure(tree, layout="circular", use_branch_lengths=False)``), which
    puts every tip on the same radius.
    """

    def __init__(self, data, columns=None, geom: str = "tile", width: float = 0.12,
                 gap: float = 0.02, offset: float = 0.04, pad_angle: float = 0.0,
                 cmap=None, palette: str = "curated", fill: str = "#5b7897",
                 bar_pad: float = 0.25, colnames: bool = True,
                 colname_size: float = 8.0, separators: Optional[bool] = None,
                 leaders: bool = False, leader_color: str = "#cccccc",
                 leader_width: float = 0.4, baseline=None, order=None):
        self.data = _index_by_name(data)
        self.columns = list(columns) if columns is not None else list(self.data.columns)
        self.geom = geom               # "tile" (heatmap ring) | "bar" (radial bars)
        self.width = width
        self.gap = gap
        self.offset = offset
        self.pad_angle = pad_angle
        self.cmap = cmap
        self.palette = palette
        self.fill = fill               # constant bar colour
        self.bar_pad = bar_pad         # fraction of the sector left blank around bars
        self.colnames = colnames
        self.colname_size = colname_size
        self.separators = separators   # None = auto (off once sectors are thin)
        self.leaders = leaders         # dotted tip -> ring guides
        self.leader_color = leader_color
        self.leader_width = leader_width
        self.baseline = baseline       # level(s) greyed out as the default state
        self.order = order             # explicit legend order

    def reserved_extent(self, layout) -> float:
        """Radial space (data units) this element claims, for the label pre-pass."""
        if not getattr(layout, "is_polar", False):
            return 0.0
        return (self.offset + len(self.columns) * (self.width + self.gap)) * layout.max_x

    def apply(self, ctx: RenderContext) -> None:
        lay = ctx.layout
        if not lay.is_polar:
            raise NotImplementedError(
                "ring() draws rings around a circular/fan tree; use a "
                "circular layout, or heatmap() for rectangular."
            )
        tips = ctx.tree.leaves()
        n = len(tips)
        step = lay.extent / max(n - 1, 1)
        # A hairline separator between sectors reads well until the sectors
        # get thin: past a few hundred tips the stroke is as wide as the
        # sector itself and the ring turns into a comb of slivers instead of
        # solid colour bands. Past that point drop the stroke and let
        # neighbours overlap very slightly, which also kills anti-alias seams.
        dense = (n > _RING_DENSE_TIPS) if self.separators is None \
            else not self.separators
        # clamp so a large pad_angle cannot eat the whole sector (or invert it)
        half = max(step / 2 - math.radians(self.pad_angle) / 2, step * 0.05)
        if dense and not self.pad_angle:
            half = step / 2                    # sectors meet edge to edge
        w = self.width * lay.max_x
        g = self.gap * lay.max_x
        r0 = ctx.ring_cursor + self.offset * lay.max_x

        # angle sitting in the fan opening (so column names never hit the rings)
        gap_angle = lay.start - (2 * math.pi - lay.extent) / 2

        # On a phylogram the tips sit at very different radii, so most of them
        # stop well short of the rings and the eye cannot tell which sector
        # belongs to which tip. A faint dotted guide from each tip out to the
        # first ring closes that gap (ggtree draws the same thing alongside
        # aligned tip labels).
        if self.leaders:
            for tip in tips:
                if tip.name not in self.data.index:
                    continue
                a = tip._angle
                ctx.scene.add(Path([lay._polar_to_xy(tip._r, a),
                                    lay._polar_to_xy(r0, a)],
                                   color=self.leader_color,
                                   width=self.leader_width, dash="dot",
                                   zorder=0.3))

        for ci, col in enumerate(self.columns):
            inner = r0 + ci * (w + g)
            outer = inner + w
            colvals = [self.data.loc[i, col] for i in self.data.index]

            if self.geom == "bar":
                # radial bars: length encodes the (numeric) value, baseline 0
                nums = [float(v) for v in colvals if v is not None]
                vmin, vmax = min(nums + [0.0]), max(nums + [0.0])
                rng = (vmax - vmin) or 1.0
                hbar = half * (1.0 - self.bar_pad)
                for tip in tips:
                    if tip.name not in self.data.index:
                        continue
                    val = float(self.data.loc[tip.name, col])
                    blen = (val - vmin) / rng * w
                    a = tip._angle
                    pts = lay._arc(inner + blen, a - hbar, a + hbar) + \
                        lay._arc(inner, a + hbar, a - hbar)
                    ctx.scene.add(Polygon(pts, facecolor=self.fill, edgecolor=None,
                                          alpha=1.0, zorder=2,
                                          label=f"{tip.name} | {col}: {val:g}"))
            else:
                scale = build_color_scale(str(col), colvals,
                                          cmap=self.cmap, palette=self.palette,
                                          baseline=self.baseline,
                                          order=self.order, swatch="patch")
                for tip in tips:
                    if tip.name not in self.data.index:
                        continue
                    val = self.data.loc[tip.name, col]
                    a = tip._angle
                    pts = lay._arc(outer, a - half, a + half) + \
                        lay._arc(inner, a + half, a - half)
                    # dense: stroke each sector in its own colour so abutting
                    # sectors have no anti-aliased hairline between them
                    fc = scale.color(val)
                    ctx.scene.add(Polygon(pts, facecolor=fc,
                                          edgecolor=fc if dense else "white",
                                          width=0.4 if dense else 0.3, alpha=1.0,
                                          zorder=2, label=f"{tip.name} | {col}: {val}"))
                ctx.add_scale(scale)
            if self.colnames:
                # Name each ring in the empty fan gap, running *along* the
                # spoke rather than across it. Across the spoke the name takes
                # up its own width in angle, and the gap is only the few
                # degrees the fan leaves open -- a name any longer than that
                # spilled out of the gap and printed itself over the ring's own
                # data. Along the spoke it takes up only its height, which
                # fits any gap; it grows back toward the centre from the ring's
                # outer edge, so consecutive rings' names do not run together.
                # Sit neighbouring rings' names on opposite sides of the gap as
                # well. Anchoring each to its own ring's outer edge separates
                # them only by the ring spacing, and a name longer than that
                # still reaches back into the one before it; putting
                # consecutive rings on either side of the gap separates them by
                # an amount that does not depend on how long the names are.
                # The wedge that is really free is narrower than the fan
                # opening: the first and last sectors each hang half a sector
                # past their tip, so that much of the opening is already
                # coloured in at both ends. On a tree with few tips the sectors
                # are wide enough to swallow the opening whole, leaving the
                # name nowhere to stand -- and rather than print it across the
                # ring it is naming, leave it out. The column already titles
                # its own legend, so nothing goes unlabelled.
                free = (2 * math.pi - lay.extent) / 2 - half
                if free <= _RING_NAME_MIN_GAP * max(1.0, self.colname_size / 8.0):
                    ctx.ring_slot += 1
                    continue
                spread = free * 0.5
                angle = gap_angle + (spread if ctx.ring_slot % 2 else -spread)
                x, y = lay._polar_to_xy(outer, angle)
                deg = math.degrees(angle)
                if 90 < (deg % 360) < 270:          # keep it the right way up
                    rot, ha = deg + 180, "left"
                else:
                    rot, ha = deg, "right"
                ctx.scene.add(Label(x, y, str(col), size=self.colname_size,
                                    color="#444444", ha=ha, va="center",
                                    rotation=rot))
            ctx.ring_slot += 1

        ctx.ring_cursor = r0 + len(self.columns) * (w + g)


def _index_by_name(data):
    """Return the DataFrame indexed by tip name (use a 'name' column if present)."""
    if "name" in getattr(data, "columns", []):
        return data.set_index("name")
    return data


# --------------------------------------------------------------------------
# painted branches -- colour each branch by stochastic-map state segments
# --------------------------------------------------------------------------
def _point_at(points, seglens, target):
    acc = 0.0
    for i, L in enumerate(seglens):
        if acc + L >= target:
            t = (target - acc) / L if L > 0 else 0.0
            x = points[i][0] + (points[i + 1][0] - points[i][0]) * t
            y = points[i][1] + (points[i + 1][1] - points[i][1]) * t
            return (x, y), i
        acc += L
    return points[-1], len(seglens) - 1


def _split_polyline(points, segs):
    """Cut a polyline into coloured pieces by fractional length (segs sum~1)."""
    if len(points) < 2:
        return [(segs[0][0], list(points))] if segs else []
    seglens = [math.dist(points[i], points[i + 1]) for i in range(len(points) - 1)]
    total = sum(seglens) or 1.0
    out = []
    start_pt, start_i, cum = points[0], 0, 0.0
    for st, frac in segs:
        cum += frac
        end_pt, end_i = _point_at(points, seglens, cum * total)
        sub = [start_pt] + [points[v] for v in range(start_i + 1, end_i + 1)] + [end_pt]
        out.append((st, sub))
        start_pt, start_i = end_pt, end_i
    return out


class _PaintedBranches(_Element):
    """Paint branches by stochastic-map state (run :func:`phytreon.stochastic_map` first).

    Each branch is split into segments proportional to the average time spent
    in each state; the child connector is drawn in the node's modal state.
    """

    def __init__(self, palette: str = "curated", size: float = 2.0):
        self.palette = palette
        self.size = size

    def apply(self, ctx: RenderContext) -> None:
        lay = ctx.layout
        nodes = ctx.tree.nodes()
        all_states = sorted({st for n in nodes
                             for st, _ in n.data.get("paint_segments", [])})
        if not all_states:
            raise ValueError("no painted-branch data on the tree; "
                             "call phytreon.stochastic_map() first")
        scale = build_color_scale("state", all_states, palette=self.palette)
        for node in nodes:
            bp = lay.branch_path(node)
            segs = node.data.get("paint_segments")
            if bp and segs:
                for st, sub in _split_polyline(bp, segs):
                    ctx.scene.add(Path(sub, color=scale.color(st), width=self.size,
                                       zorder=1))
            elif bp:
                ctx.scene.add(Path(bp, color="black", width=self.size, zorder=1))
            conn = lay.child_connector(node)
            if conn:
                modal = node.data.get("ace_state")
                ctx.scene.add(Path(conn, color=scale.color(modal) if modal else "black",
                                   width=self.size, zorder=1))
        ctx.add_scale(scale)


# --------------------------------------------------------------------------
# pie charts of ancestral-state probabilities at nodes
# --------------------------------------------------------------------------
class _NodePies(_Element):
    """Draw a pie chart at each internal node from a probability dict.

    Reads ``node.data[attr]`` (default ``'ace_probs'``, as written by
    :func:`phytreon.ace_ml` / :func:`phytreon.stochastic_map`) -- a
    ``{state: prob}`` mapping -- and draws a small pie wedge per state.
    Rectangular/slanted.
    """

    def __init__(self, attr: str = "ace_probs", radius: float = 0.4,
                 palette: str = "curated", tips: bool = False):
        self.attr = attr
        self.radius = radius            # in tip-row units
        self.palette = palette
        self.tips = tips

    def apply(self, ctx: RenderContext) -> None:
        lay = ctx.layout
        nodes = [n for n in ctx.tree.traverse()
                 if (self.tips or not n.is_leaf) and self.attr in n.data]
        states = sorted({s for n in nodes for s in n.data[self.attr]})
        if not states:
            return
        scale = build_color_scale("state", states, palette=self.palette)
        # keep pies visually circular: the rectangular layout has unequal x/y
        # data scales, so compensate using the default figure aspect (assumes
        # the default figsize from render_mpl; custom figsize may distort pies).
        n = max(ctx.tree.n_leaves, 1)
        fig_w, fig_h = 8.0, max(2.6, min(0.34 * n, 30))
        x_span = 1.21 * (lay.max_x or 1.0)
        ry = self.radius                              # rows (y data units)
        rxx = ry * (fig_h / fig_w) * (x_span / n)     # x data units
        for node in nodes:
            probs = node.data[self.attr]
            cx, cy = node.x, node.y
            a0 = 0.0
            for st in states:
                p = float(probs.get(st, 0.0))
                if p <= 0:
                    continue
                a1 = a0 + p * 2 * math.pi
                pts = [(cx, cy)]
                steps = max(2, int((a1 - a0) / math.radians(10)) + 1)
                for i in range(steps + 1):
                    a = a0 + (a1 - a0) * i / steps
                    pts.append((cx + rxx * math.cos(a), cy + ry * math.sin(a)))
                ctx.scene.add(Polygon(pts, facecolor=scale.color(st),
                                      edgecolor="white", width=0.3, zorder=3,
                                      label=f"{st}: {p:.2f}"))
                a0 = a1
        ctx.add_scale(scale)


# --------------------------------------------------------------------------
# time axis + geological time scale
# --------------------------------------------------------------------------
# ICS Phanerozoic periods: (name, young_Ma, old_Ma, colour)
GEO_PERIODS = [
    ("Quaternary", 0.0, 2.58, "#F9F97F"),
    ("Neogene", 2.58, 23.03, "#FFE619"),
    ("Paleogene", 23.03, 66.0, "#FD9A52"),
    ("Cretaceous", 66.0, 145.0, "#7FC64E"),
    ("Jurassic", 145.0, 201.4, "#34B2C9"),
    ("Triassic", 201.4, 251.9, "#812B92"),
    ("Permian", 251.9, 298.9, "#F04028"),
    ("Carboniferous", 298.9, 358.9, "#67A599"),
    ("Devonian", 358.9, 419.2, "#CB8C37"),
    ("Silurian", 419.2, 443.8, "#B3E1B6"),
    ("Ordovician", 443.8, 485.4, "#009270"),
    ("Cambrian", 485.4, 538.8, "#7FA056"),
]


def _nice_ticks(lo, hi, n):
    span = hi - lo
    if span <= 0:
        return [lo]
    raw = span / max(n - 1, 1)
    mag = 10 ** math.floor(math.log10(raw))
    step = next(m * mag for m in (1, 2, 2.5, 5, 10) if raw <= m * mag)
    start = math.ceil(lo / step) * step
    ticks, v = [], start
    while v <= hi + 1e-9:
        ticks.append(round(v, 10))
        v += step
    return ticks


class _TimeAxis(_Element):
    """A time axis below a (time-calibrated) rectangular tree.

    Branch lengths are assumed to be time; the most recent tip is at
    ``present`` (default 0) and time increases toward the root.  ``geo=True``
    shades the geological periods (Phanerozoic) behind the tree.
    """

    #: this element defines where the present sits on the x axis
    defines_present = True

    def __init__(self, geo: bool = False, n_ticks: int = 6, gridlines: bool = False,
                 present: float = 0.0, unit: str = "Mya", band_alpha: float = 0.3,
                 fontsize: float = 8.0):
        self.geo = geo
        self.n_ticks = n_ticks
        self.gridlines = gridlines
        self.present = present
        self.unit = unit
        self.band_alpha = band_alpha
        self.fontsize = fontsize

    def apply(self, ctx: RenderContext) -> None:
        lay = ctx.layout
        if lay.is_polar or getattr(lay, "kind", "rect") != "rect":
            raise NotImplementedError("time_axis() is for rectangular layouts.")
        maxx = lay.max_x
        n = ctx.tree.n_leaves
        ytop, ybot = -0.5, n - 0.5

        if self.geo:
            for name, y0, y1, col in GEO_PERIODS:
                xa = max(0.0, min(maxx, maxx - (y1 - self.present)))
                xb = max(0.0, min(maxx, maxx - (y0 - self.present)))
                if xb - xa <= 0:
                    continue
                ctx.scene.add(Polygon([(xa, ytop), (xb, ytop), (xb, ybot), (xa, ybot)],
                                      facecolor=col, edgecolor=None,
                                      alpha=self.band_alpha, zorder=0))
                if (xb - xa) > 0.03 * maxx:
                    ctx.scene.add(Label((xa + xb) / 2, ytop - 0.4, name,
                                        size=self.fontsize - 1, color="#555555",
                                        ha="center", va="bottom", rotation=90))

        ybase = ybot + 0.3
        ctx.scene.add(Path([(0, ybase), (maxx, ybase)], color="#333333", width=1.0,
                           zorder=4))
        for tbp in _nice_ticks(0, maxx, self.n_ticks):
            xt = maxx - tbp
            ctx.scene.add(Path([(xt, ybase), (xt, ybase + 0.25)], color="#333333",
                               width=1.0, zorder=4))
            ctx.scene.add(Label(xt, ybase + 0.45, f"{self.present + tbp:g}",
                                size=self.fontsize, color="#333333",
                                ha="center", va="top"))
            if self.gridlines:
                ctx.scene.add(Path([(xt, ytop), (xt, ybot)], color="#dddddd",
                                   width=0.6, dash="dot", zorder=0))
        if self.unit:
            ctx.scene.add(Label(maxx / 2, ybase + 1.2, self.unit, size=self.fontsize,
                                color="#333333", ha="center", va="top"))


# --------------------------------------------------------------------------
# collapsed clades (triangles)
# --------------------------------------------------------------------------
class _CollapsedClades(_Element):
    """Draw a triangle for every clade collapsed by
    :func:`phytreon.treeops.collapse_clade`.

    The two sides run to the collapsed clade's nearest and farthest leaf, so
    the wedge shows both how deep the hidden clade is and how uneven it is --
    the convention iTOL uses. Set ``scale_height=True`` to also let the width
    grow with the number of hidden tips, so a big clade reads as a big block.
    """

    def __init__(self, color="#8494a8", alpha: float = 1.0,
                 height: float = 0.8, scale_height: bool = False,
                 edgecolor: Optional[str] = None, min_extent: float = 0.04):
        self.color = color
        self.alpha = alpha
        self.height = height           # rows spanned (before any scaling)
        self.scale_height = scale_height
        self.edgecolor = edgecolor
        #: shortest triangle to draw, as a fraction of the tree depth. On a
        #: cladogram (or a clade of zero-length branches) the hidden leaves sit
        #: at the collapsed node itself, so an extent taken straight from the
        #: branch lengths would be a zero-size, invisible triangle.
        self.min_extent = min_extent

    def apply(self, ctx: RenderContext) -> None:
        lay = ctx.layout
        nodes = [n for n in ctx.tree.traverse() if "_collapsed" in n.data]
        if not nodes:
            return
        cfunc, scale = ctx.resolve_color(self.color, nodes, default="#8494a8")
        biggest = max(n.data["_collapsed"]["n"] for n in nodes)
        for node in nodes:
            info = node.data["_collapsed"]
            h = self.height / 2
            if self.scale_height:
                h *= 0.35 + 0.65 * (info["n"] / biggest)
            # ask the layout for the span in the units it actually draws in;
            # on a cladogram that is edges, not branch length
            near, far = lay._collapsed_span(node, lay.use_branch_lengths
                                            and ctx.tree.has_branch_lengths)
            floor = self.min_extent * lay.max_x
            if far < floor:            # cladogram / zero-length clade
                near, far = (near / far * floor if far else floor * 0.6), floor
            if lay.is_polar:
                a = node._angle
                da = h * (lay.extent / max(lay.n_leaves - 1, 1))
                pts = [lay._polar_to_xy(node._r, a),
                       lay._polar_to_xy(node._r + near, a - da),
                       lay._polar_to_xy(node._r + far, a + da)]
            else:
                pts = [(node.x, node.y),
                       (node.x + near, node.y - h),
                       (node.x + far, node.y + h)]
            ctx.scene.add(Polygon(pts, facecolor=cfunc(node),
                                  edgecolor=self.edgecolor, alpha=self.alpha,
                                  width=0.6 if self.edgecolor else 0.0,
                                  zorder=1.5,
                                  label=f"{node.name} ({info['n']} tips)"))
        if scale is not None:
            ctx.add_scale(scale)


# --------------------------------------------------------------------------
# node age / confidence bars (95% HPD)
# --------------------------------------------------------------------------
class _NodeBars(_Element):
    """Horizontal bars showing an interval at each internal node.

    The standard way to show divergence-time uncertainty -- a bar spanning the
    95% HPD of every node's age, as FigTree's "node bars" and ggtree's
    ``geom_range`` draw it. ``lower``/``upper`` name per-node data keys (e.g.
    written by BEAST/TreeAnnotator).

    Values are read as **ages** on the same scale as
    :meth:`~phytreon.plot.figure.TreeFigure.time_axis`: distance back from
    ``present``, increasing toward the root. Pass ``as_age=False`` if the two
    keys already hold plot x coordinates instead.
    """

    def __init__(self, lower: str = "height_95_lower",
                 upper: str = "height_95_upper", color: str = "#3a7ac1",
                 width: float = 3.0, alpha: float = 0.55,
                 present: Optional[float] = None, as_age: bool = True):
        self.lower = lower
        self.upper = upper
        self.color = color
        self.width = width
        self.alpha = alpha
        self.present = present
        self.as_age = as_age

    def apply(self, ctx: RenderContext) -> None:
        lay = ctx.layout
        if lay.is_polar or getattr(lay, "kind", "rect") != "rect":
            raise NotImplementedError(
                "node_bars() is for rectangular layouts (the bar runs along "
                "the time axis).")
        maxx = lay.max_x
        # follow the figure's time axis unless the caller pinned it; the two
        # used to default to 0 independently, so setting it on one silently
        # shifted the bars off the axis they are read against
        present = self.present if self.present is not None else ctx.present
        drawn = 0
        for node in ctx.tree.traverse():
            lo, hi = node.data.get(self.lower), node.data.get(self.upper)
            if not (is_numeric(lo) and is_numeric(hi)):
                continue
            if self.as_age:
                # age -> x: the present sits at max_x, ages run back to the root
                x0 = maxx - (float(hi) - present)
                x1 = maxx - (float(lo) - present)
            else:
                x0, x1 = float(lo), float(hi)
            ctx.scene.add(Path([(x0, node.y), (x1, node.y)], color=self.color,
                               width=self.width, opacity=self.alpha, zorder=2.5))
            drawn += 1
        if not drawn:
            raise ValueError(
                f"no node carries both {self.lower!r} and {self.upper!r}; "
                f"node_bars() needs an interval per node -- e.g. the "
                f"height_95%_HPD annotations on a BEAST summary tree")


# --------------------------------------------------------------------------
# connections between nodes (HGT, co-occurrence, host-symbiont)
# --------------------------------------------------------------------------
class _Connections(_Element):
    """Curved links drawn between arbitrary pairs of tips/nodes.

    iTOL's ``DATASET_CONNECTION``: how horizontal gene transfer, gene sharing,
    co-occurrence or host-symbiont pairings get shown on a tree. On a circular
    layout the curves bend toward the centre (iTOL's ``CENTER_CURVES``), which
    is what keeps a dense set of links readable; on a rectangular one they bow
    out past the tips so they clear the tree.

    ``pairs`` is an iterable of ``(name1, name2)``, optionally
    ``(name1, name2, value)``, or a DataFrame with those columns. Pass
    ``color="value"`` to colour each link by its third field.
    """

    def __init__(self, pairs, color="#c1553b", width: float = 0.9,
                 alpha: float = 0.55, curvature: float = 0.55,
                 dash: Optional[str] = None, cmap=None,
                 palette: str = "curated", label: str = "connection"):
        self.pairs = pairs
        self.color = color
        self.width = width
        self.alpha = alpha
        self.curvature = curvature     # 0 = straight, 1 = bends to the centre
        self.dash = dash
        self.cmap = cmap
        self.palette = palette
        self.label = label

    def _rows(self):
        data = self.pairs
        if hasattr(data, "itertuples"):           # DataFrame
            ncol = len(data.columns)
            for row in data.itertuples(index=False):
                vals = tuple(row)
                yield (str(vals[0]), str(vals[1]),
                       vals[2] if ncol > 2 else None)
            return
        for item in data:
            item = tuple(item)
            yield (str(item[0]), str(item[1]),
                   item[2] if len(item) > 2 else None)

    @staticmethod
    def _bezier(p0, p1, ctrl, n: int = 32):
        return [((1 - t) ** 2 * p0[0] + 2 * (1 - t) * t * ctrl[0] + t ** 2 * p1[0],
                 (1 - t) ** 2 * p0[1] + 2 * (1 - t) * t * ctrl[1] + t ** 2 * p1[1])
                for t in (i / (n - 1) for i in range(n))]

    def apply(self, ctx: RenderContext) -> None:
        lay = ctx.layout
        by_name = {n.name: n for n in ctx.tree.traverse() if n.name}
        rows = list(self._rows())
        missing = {nm for a, b, _ in rows for nm in (a, b) if nm not in by_name}
        if missing:
            shown = sorted(missing)[:5]
            raise ValueError(
                f"connections() got names that are not in the tree: {shown}"
                f"{' ...' if len(missing) > 5 else ''}")

        vals = [v for _, _, v in rows if v is not None]
        scale = None
        if vals and self.color == "value":
            scale = build_color_scale(self.label, vals, cmap=self.cmap,
                                      palette=self.palette)

        for a, b, val in rows:
            na, nb = by_name[a], by_name[b]
            p0, p1 = (na.x, na.y), (nb.x, nb.y)
            if lay.is_polar:
                # pull the control point toward the centre: the chord look
                ctrl = ((p0[0] + p1[0]) / 2 * (1 - self.curvature),
                        (p0[1] + p1[1]) / 2 * (1 - self.curvature))
            else:
                # bow sideways, past the deeper tip, so links clear the tree
                ctrl = (max(p0[0], p1[0]) + self.curvature * lay.max_x,
                        (p0[1] + p1[1]) / 2)
            col = scale.color(val) if scale is not None else self.color
            ctx.scene.add(Path(self._bezier(p0, p1, ctrl), color=col,
                               width=self.width, opacity=self.alpha,
                               dash=self.dash, zorder=0.8))
        if scale is not None:
            ctx.add_scale(scale)


# --------------------------------------------------------------------------
# compact scale bar
# --------------------------------------------------------------------------
class _ScaleBar(_Element):
    """A short bar giving the branch-length scale (ggtree's ``geom_treescale``).

    Unlike :class:`_TimeAxis` this assumes nothing about branch lengths being
    time and works on any layout, which is what a plain substitutions/site
    phylogram needs. The length defaults to a round number near a tenth of the
    tree's depth.
    """

    def __init__(self, length: Optional[float] = None, label: Optional[str] = None,
                 x: Optional[float] = None, y: Optional[float] = None,
                 color: str = "#333333", width: float = 1.4,
                 fontsize: float = 8.0):
        self.length = length
        self.label = label
        self.x = x
        self.y = y
        self.color = color
        self.width = width
        self.fontsize = fontsize

    @staticmethod
    def _nice(v: float) -> float:
        """Round down to 1, 2 or 5 times a power of ten."""
        if v <= 0:
            return 1.0
        mag = 10 ** math.floor(math.log10(v))
        for m in (5.0, 2.0, 1.0):
            if v >= m * mag:
                return m * mag
        return mag

    def apply(self, ctx: RenderContext) -> None:
        lay = ctx.layout
        length = self.length if self.length is not None \
            else self._nice(lay.max_x / 10.0)
        xmin, ymin, xmax, ymax = ctx.scene.bounds()
        span = (ymax - ymin) or 1.0
        x0 = self.x if self.x is not None else xmin
        y0 = self.y if self.y is not None else ymax + 0.05 * span
        ctx.scene.add(Path([(x0, y0), (x0 + length, y0)], color=self.color,
                           width=self.width, zorder=4))
        tick = 0.012 * span
        for xt in (x0, x0 + length):                 # end ticks
            ctx.scene.add(Path([(xt, y0 - tick), (xt, y0 + tick)],
                               color=self.color, width=self.width, zorder=4))
        text = self.label if self.label is not None else f"{length:g}"
        ctx.scene.add(Label(x0 + length / 2, y0 + 0.03 * span, text,
                            size=self.fontsize, color=self.color,
                            ha="center", va="top"))


# --------------------------------------------------------------------------
# domain architecture / gene neighbourhood beside the tips
# --------------------------------------------------------------------------
class _DomainTrack(_Element):
    """Draw each tip's domain architecture (or gene neighbourhood) beside it.

    A tree of protein-domain sequences is usually only half the story: the
    other half is what each protein is *built out of*, and putting the two
    side by side is what lets a reader see a domain being gained, lost or
    swapped along a clade. Papers in this area rarely show one without the
    other.

    ``data`` maps a tip name to its architecture, given either as plain names::

        {"protein_A": ["ParB", "DUF262", "HTH"], ...}

    or as ``(name, length)`` pairs when the relative sizes matter::

        {"protein_A": [("ParB", 120), ("DUF262", 200), ("HTH", 60)], ...}

    ``arrows=True`` draws each element as a block arrow instead of a
    rectangle, the convention for a gene neighbourhood; pass negative lengths
    for genes transcribed the other way and the arrow flips.
    """

    def __init__(self, data, width: float = 0.5, offset: float = 0.04,
                 height: float = 0.7, gap: float = 0.06, arrows: bool = False,
                 palette: str = "curated", to_scale: bool = True,
                 labels: bool = False, label_size: float = 6.0,
                 edgecolor: Optional[str] = "white"):
        self.data = {str(k): list(v) for k, v in dict(data).items()}
        self.width = width
        self.offset = offset
        self.height = height
        self.gap = gap                 # blank between adjacent elements
        self.arrows = arrows
        self.palette = palette
        #: honour the supplied lengths; False spaces the elements evenly
        self.to_scale = to_scale
        self.labels = labels
        self.label_size = label_size
        self.edgecolor = edgecolor

    def reserved_extent(self, layout) -> float:
        return 0.0                     # rectangular only; nothing radial

    @staticmethod
    def _parts(items):
        """Normalise an architecture to ``[(name, length), ...]``."""
        out = []
        for item in items:
            if isinstance(item, (tuple, list)):
                out.append((str(item[0]), float(item[1])))
            else:
                out.append((str(item), 1.0))
        return out

    def _block(self, x0, x1, y0, y1, forward=True):
        """A rectangle, or a block arrow pointing the way the gene reads."""
        if not self.arrows:
            return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        head = min(abs(x1 - x0) * 0.35, (y1 - y0) * 0.9)
        ym = 0.5 * (y0 + y1)
        if forward:
            return [(x0, y0), (x1 - head, y0), (x1, ym),
                    (x1 - head, y1), (x0, y1)]
        return [(x1, y0), (x0 + head, y0), (x0, ym), (x0 + head, y1), (x1, y1)]

    def apply(self, ctx: RenderContext) -> None:
        lay = ctx.layout
        if lay.is_polar:
            raise NotImplementedError(
                "domains() is for rectangular layouts (the architecture runs "
                "along x beside each tip).")
        tips = [t for t in ctx.tree.leaves() if t.name in self.data]
        if not tips:
            keys = sorted(self.data)[:3]
            raise ValueError(
                f"no tip name matches the architecture data (it has "
                f"{keys}{' ...' if len(self.data) > 3 else ''})")

        names = sorted({nm for items in self.data.values()
                        for nm, _ in self._parts(items)})
        scale = build_color_scale("domain", names, palette=self.palette,
                                  swatch="patch")

        total_w = self.width * lay.max_x
        x0 = ctx.track_cursor + (self.offset + 0.02) * lay.max_x
        longest = max(sum(abs(ln) for _, ln in self._parts(v))
                      for v in self.data.values()) or 1.0
        gap_w = self.gap * total_w / max(len(names), 1)
        half = self.height / 2

        for tip in tips:
            parts = self._parts(self.data[tip.name])
            span = sum(abs(ln) for _, ln in parts) or 1.0
            # to scale: every architecture shares one ruler, so a longer
            # protein really draws longer. Otherwise each fills the track,
            # which compares composition rather than size.
            scale_w = (total_w / longest) if self.to_scale else (total_w / span)
            cursor = x0
            for nm, ln in parts:
                w = abs(ln) * scale_w - gap_w
                if w <= 0:
                    w = abs(ln) * scale_w
                pts = self._block(cursor, cursor + w,
                                  tip.y - half, tip.y + half, forward=ln >= 0)
                ctx.scene.add(Polygon(pts, facecolor=scale.color(nm),
                                      edgecolor=self.edgecolor,
                                      width=0.4 if self.edgecolor else 0.0,
                                      alpha=1.0, zorder=2, align=True,
                                      label=f"{tip.name} | {nm}"))
                if self.labels and w > 0:
                    ctx.scene.add(Label(cursor + w / 2, tip.y, nm,
                                        size=self.label_size,
                                        color=readable_on(scale.color(nm)),
                                        ha="center", va="center", align=True))
                cursor += abs(ln) * scale_w
        ctx.track_cursor = x0 + total_w
        ctx.add_scale(scale)
