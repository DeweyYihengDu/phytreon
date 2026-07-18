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
from typing import Optional

from ..core.tree import Node
from ..scene import Label, Marker, Path, Polygon, Raster
from .figure import _Element, RenderContext, build_color_scale, is_numeric

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
        for i, tip in enumerate(tips):
            text = tip.name or ""
            if not text or (step > 1 and i % step != 0):
                continue
            # a collapsed clade is drawn as a triangle reaching out to its
            # farthest hidden leaf; its label has to clear that, not sit on it
            far = tip.data.get("_collapsed", {}).get("far", 0.0)
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
                # offset outward along the branch direction, rotate to match
                a = tip._angle
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
        if scale is not None:
            ctx.add_scale(scale)


class _NodeLabels(_Element):
    """Label internal nodes -- by default their support values."""

    def __init__(self, attr: str = "support", size: float = 7.0,
                 color="#666666", offset: float = 0.0, fmt: str = "{:g}"):
        self.attr = attr
        self.size = size
        self.color = color
        self.offset = offset
        self.fmt = fmt

    def apply(self, ctx: RenderContext) -> None:
        lay = ctx.layout
        for node in ctx.tree.traverse():
            if node.is_leaf or node.is_root:
                continue
            val = getattr(node, self.attr, None)
            if val is None:
                val = node.data.get(self.attr)
            if val is None:
                continue
            text = self.fmt.format(val) if isinstance(val, (int, float)) else str(val)
            if getattr(lay, "kind", "rect") != "rect":
                x, y, ha, va = node.x, node.y, "center", "bottom"
            else:
                # sit just above the branch leading into the node (no overlap
                # with the vertical connector or the node itself)
                x = 0.5 * (node.parent.x + node.x) + self.offset
                y = node.y - 0.3
                ha, va = "center", "center"
            ctx.scene.add(Label(x, y, text, size=self.size, color=self.color,
                                ha=ha, va=va))


# --------------------------------------------------------------------------
# points
# --------------------------------------------------------------------------
class _Points(_Element):
    def __init__(self, which: str = "tip", color="black", size=6.0,
                 marker: str = "o", shape=None, edgecolor: Optional[str] = None,
                 palette: str = "curated", cmap=None):
        self.which = which                 # tip | node | all
        self.color = color
        self.size = size
        self.marker = marker
        self.shape = shape                 # categorical column -> marker shape
        self.edgecolor = edgecolor
        self.palette = palette
        self.cmap = cmap

    def _select(self, ctx):
        if self.which == "tip":
            return ctx.tree.leaves()
        if self.which == "node":
            return [n for n in ctx.tree.traverse() if not n.is_leaf]
        return ctx.tree.nodes()

    def apply(self, ctx: RenderContext) -> None:
        nodes = self._select(ctx)
        cfunc, cscale = ctx.resolve_color(self.color, nodes, default="black",
                                          palette=self.palette, cmap=self.cmap)
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
    """Shade the rectangle / wedge occupied by a clade."""

    def __init__(self, node: Optional[Node] = None, taxa=None,
                 fill="#fdbf6f", alpha: float = 0.3, extend: float = 0.0):
        self.node = node
        self.taxa = taxa
        self.fill = fill
        self.alpha = alpha
        self.extend = extend

    def _target(self, ctx) -> Optional[Node]:
        if self.node is not None:
            return self.node
        if self.taxa is not None:
            return ctx.tree.get_mrca(self.taxa)
        return None

    def apply(self, ctx: RenderContext) -> None:
        node = self._target(ctx)
        if node is None:
            return
        lay = ctx.layout
        leaves = node.get_leaves()
        if lay.is_polar:
            angles = [lf._angle for lf in leaves]
            a0, a1 = min(angles), max(angles)
            da = (a1 - a0) / max(len(leaves) - 1, 1) / 2 + 1e-9
            inner = node.parent._r if node.parent else 0.0
            outer = lay.inner_radius + lay.max_x + self.extend
            pts = lay._arc(outer, a0 - da, a1 + da)
            pts += lay._arc(inner, a1 + da, a0 - da)
            ctx.scene.add(Polygon(pts, facecolor=self.fill, edgecolor=None,
                                  alpha=self.alpha, zorder=0))
        else:
            rows = [lf.y for lf in leaves]
            # hug the clade: start at the MRCA node, not the root-ward parent
            x0 = node.x
            x1 = lay.max_x + self.extend
            y0, y1 = min(rows) - 0.45, max(rows) + 0.45
            pts = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
            ctx.scene.add(Polygon(pts, facecolor=self.fill, edgecolor=None,
                                  alpha=self.alpha, zorder=0, rounded=True))


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
                 cell_gap: float = 0.05, separators: Optional[bool] = None):
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
                                       palette=self.palette)
            scales = {c: shared for c in cols}
        else:
            scales = {c: build_color_scale(str(c), list(df[c]), cmap=self.cmap,
                                           palette=self.palette) for c in cols}

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
                 colname_size: float = 9.0):
        self.data = _index_by_name(data)
        self.column = column
        self.width = width
        self.offset = offset
        self.fill = fill
        self.bar_height = bar_height
        self.colname = colname
        self.colname_size = colname_size

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

    ``separators`` controls the hairline between neighbouring sectors: ``None``
    (default) turns it off automatically past ~150 tips, where the stroke would
    be wider than the sector itself and the ring would read as a comb of
    slivers instead of solid colour bands. Force it with ``True``/``False``.
    """

    def __init__(self, data, columns=None, geom: str = "tile", width: float = 0.12,
                 gap: float = 0.02, offset: float = 0.04, pad_angle: float = 0.0,
                 cmap=None, palette: str = "curated", fill: str = "#5b7897",
                 bar_pad: float = 0.25, colnames: bool = True,
                 colname_size: float = 8.0, separators: Optional[bool] = None):
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
                                          cmap=self.cmap, palette=self.palette)
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
                # label each ring in the empty fan gap, oriented tangentially
                # (perpendicular to the spoke) so the names separate radially
                # instead of stacking on one another.
                rmid = (inner + outer) / 2
                x, y = lay._polar_to_xy(rmid, gap_angle)
                rot = (math.degrees(gap_angle) + 90) % 180 - 0  # keep upright
                ctx.scene.add(Label(x, y, str(col), size=self.colname_size,
                                    color="#444444", ha="center", va="center",
                                    rotation=rot))

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
                 edgecolor: Optional[str] = None):
        self.color = color
        self.alpha = alpha
        self.height = height           # rows spanned (before any scaling)
        self.scale_height = scale_height
        self.edgecolor = edgecolor

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
            near, far = info["near"], info["far"]
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
                 present: float = 0.0, as_age: bool = True):
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
        drawn = 0
        for node in ctx.tree.traverse():
            lo, hi = node.data.get(self.lower), node.data.get(self.upper)
            if not (is_numeric(lo) and is_numeric(hi)):
                continue
            if self.as_age:
                # age -> x: the present sits at max_x, ages run back to the root
                x0 = maxx - (float(hi) - self.present)
                x1 = maxx - (float(lo) - self.present)
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
