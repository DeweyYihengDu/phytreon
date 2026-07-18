"""Tanglegrams: two facing trees with their shared tips linked.

The standard way to show that two datasets disagree about a phylogeny --
e.g. a tree inferred from genomic data next to one inferred from
transcriptomic data, or a host tree beside its parasites'::

    fig = pt.TangleFigure(genome_tree, transcript_tree,
                          titles=("genome", "transcriptome"))
    fig.untangle()                       # rotate clades to line the tips up
    fig.connect(highlight_discordant=True)
    fig.save("tanglegram.pdf")

Each side is an ordinary :class:`~phytreon.plot.figure.TreeFigure`, reachable
as ``fig.left`` / ``fig.right``, so every element (tip points, clade
highlights, support labels, colour scales …) works on either tree::

    fig.left.tip_points(color="clade")
    fig.right.support_labels()

Implementation note: each side is built independently into its own scene, the
right-hand scene is then mirrored and shifted so the two trees face each
other, and the two are merged into one scene.  That way the tanglegram reuses
the entire element/layout/colour stack rather than reimplementing it.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Dict, List, Optional, Sequence, Tuple

from ..scene import Label, Path, Scene
from ..treeops import crossing_number, untangle as _untangle
from .figure import RenderContext, TreeFigure, _Renderable, build_color_scale

#: links whose tips sit in conflicting positions in the two trees
DISCORDANT_COLOR = "#c1553b"
CONCORDANT_COLOR = "#b8b8b8"


class _TangleLayout:
    """Layout shim describing the merged two-tree scene to a backend.

    The merged scene is already in final display coordinates, so the backend
    only needs the handful of hints it reads off a layout.  ``kind="rect"``
    keeps the ordinary (non-equal-aspect, y-flipped) framing.
    """

    is_polar = False
    equal_aspect = False
    invert_y = True
    kind = "rect"

    def __init__(self, max_x: float):
        self.max_x = max_x


def _mirror_scene(scene: Scene, x0: float, yscale: float) -> Scene:
    """Reflect a scene about ``x = x0/2`` and rescale y.

    ``x -> x0 - x`` turns the right-hand tree around so its tips face the
    middle; horizontal text anchors flip with it, otherwise every label would
    run back over its own tree.
    """
    flip_ha = {"left": "right", "right": "left", "center": "center"}
    out = Scene()
    for p in scene.paths:
        out.add(replace(p, points=[(x0 - x, y * yscale) for x, y in p.points]))
    for poly in scene.polygons:
        out.add(replace(poly, points=[(x0 - x, y * yscale) for x, y in poly.points]))
    for m in scene.markers:
        out.add(replace(m, x=x0 - m.x, y=m.y * yscale))
    for lb in scene.labels:
        out.add(replace(lb, x=x0 - lb.x, y=lb.y * yscale,
                        ha=flip_ha.get(lb.ha, lb.ha), rotation=-lb.rotation))
    for r in scene.rasters:
        out.add(replace(r, x0=x0 - r.x1, x1=x0 - r.x0,
                        y0=r.y0 * yscale, y1=r.y1 * yscale))
    return out


def _crossing_flags(order_l: List[str], order_r: List[str]) -> Dict[str, bool]:
    """Which shared tips have a link that crosses at least one other link.

    A pair of links crosses when the two tips appear in opposite order on the
    two sides, so this walks every pair -- O(n^2), which is fine at the sizes
    a tanglegram stays readable at (a few hundred tips).
    """
    rank_r = {name: i for i, name in enumerate(order_r)}
    flags = {name: False for name in order_l}
    for i, a in enumerate(order_l):
        for b in order_l[i + 1:]:
            # a is above b on the left; crossing iff a is below b on the right
            if rank_r[a] > rank_r[b]:
                flags[a] = True
                flags[b] = True
    return flags


class TangleFigure(_Renderable):
    """Two trees drawn face to face, their shared tips joined by links.

    ``left`` / ``right`` accept either a :class:`~phytreon.core.tree.Tree` (a
    default :class:`~phytreon.plot.figure.TreeFigure` is built for it) or a
    ready-made ``TreeFigure``, so the two sides can be styled independently.

    Only tips whose names appear in *both* trees are linked; unmatched tips
    are still drawn, just left unconnected.

    The middle band holds both the tip labels and the links.  By default its
    width is estimated from the longest label so the names fit; pass ``gap``
    (a fraction of the two trees' combined depth) to set it yourself.
    ``tip_labels`` selects which sides are named -- ``"both"`` (default),
    ``"left"``, ``"right"`` or ``False``.
    """

    def __init__(self, left, right, *, gap: Optional[float] = None,
                 titles: Optional[Sequence[str]] = None,
                 tip_labels: str = "both", layout: str = "rectangular",
                 **layout_kwargs):
        lab = self._label_sides(tip_labels)
        self.left = self._as_figure(left, layout, "left" in lab, **layout_kwargs)
        self.right = self._as_figure(right, layout, "right" in lab, **layout_kwargs)
        self.gap = gap
        self.titles = tuple(titles) if titles else None
        self.title: Optional[str] = None
        self._link_style = dict(color=CONCORDANT_COLOR, width=0.7, dash="dot",
                                opacity=1.0)
        self._link_color_by: Optional[str] = None
        self._highlight_discordant = False
        self._link_inset: Optional[Tuple[float, float]] = None

    @staticmethod
    def _label_sides(spec) -> set:
        """Normalise the ``tip_labels`` spec to a set of sides."""
        if spec is True or spec == "both":
            return {"left", "right"}
        if spec in (False, None, "none"):
            return set()
        if spec in ("left", "right"):
            return {spec}
        raise ValueError("tip_labels must be 'left', 'right', 'both' or False")

    @staticmethod
    def _as_figure(obj, layout: str, tip_labels: bool, **layout_kwargs):
        if isinstance(obj, TreeFigure):
            return obj                       # caller styled it themselves
        fig = TreeFigure(obj, layout=layout, **layout_kwargs)
        if tip_labels:
            fig.tip_labels()
        return fig

    # -- composition -----------------------------------------------------
    def untangle(self, *, fix: Optional[str] = "left",
                 rounds: int = 3) -> "TangleFigure":
        """Rotate clades so the two tip orders line up as closely as possible.

        Delegates to :func:`phytreon.treeops.untangle`; both trees are
        modified in place (rotation reorders children, so the topology and
        every branch length are untouched -- only the reading order moves).
        """
        _untangle(self.left.tree, self.right.tree, fix=fix, rounds=rounds)
        return self

    def connect(self, color: Optional[str] = None, width: float = 0.7,
                dash: Optional[str] = "dot", opacity: float = 1.0,
                highlight_discordant: bool = False,
                inset: Optional[Tuple[float, float]] = None) -> "TangleFigure":
        """Style the tip-to-tip links.

        ``color`` may be a literal colour or the name of a data column on the
        left tree's tips (which also emits a legend).  ``highlight_discordant``
        instead colours every link that crosses another one -- after
        :meth:`untangle` those are the taxa the two trees order differently.
        It flags *crossings*, so it can come up empty on trees that still
        conflict (see :func:`phytreon.treeops.crossing_number`); check
        :func:`phytreon.treeops.robinson_foulds` before concluding the two
        trees agree.

        ``inset`` is a ``(left, right)`` pair of fractions of the middle band
        to leave clear at each end, so the links start past the tip labels
        rather than running behind them.  It defaults to the space the labels
        are estimated to need; pass an explicit pair to tune it (or ``(0, 0)``
        to run the links tip to tip).
        """
        self._link_style = dict(color=color or CONCORDANT_COLOR, width=width,
                                dash=dash, opacity=opacity)
        self._link_color_by = color
        self._highlight_discordant = highlight_discordant
        self._link_inset = inset
        return self

    def titled(self, title: str) -> "TangleFigure":
        self.title = title
        return self

    # -- diagnostics -----------------------------------------------------
    def crossings(self) -> int:
        """Current number of crossing links (see
        :func:`phytreon.treeops.crossing_number`)."""
        return crossing_number(self.left.tree, self.right.tree)

    # -- building --------------------------------------------------------
    def _build(self) -> RenderContext:
        lctx = self.left._build()
        rctx = self.right._build()
        for side, ctx in (("left", lctx), ("right", rctx)):
            # mirroring a tree only reads as "facing" for layouts whose depth
            # runs along x; reflecting a circular tree just turns it inside out
            if ctx.layout.is_polar or getattr(ctx.layout, "kind", "rect") != "rect":
                raise ValueError(
                    f"tanglegrams need a layout whose depth runs along x; the "
                    f"{side} tree uses {type(ctx.layout).__name__}. Use "
                    f"'rectangular' or 'slanted'.")
        lmax = lctx.layout.max_x or 1.0
        rmax = rctx.layout.max_x or 1.0

        # put the right tree on the left tree's row grid, so the two tip
        # columns span the same height even with different taxon counts
        ln = max(self.left.tree.n_leaves - 1, 1)
        rn = max(self.right.tree.n_leaves - 1, 1)
        yscale = ln / rn

        # the middle band carries both the tip labels and the links
        band, res_l, res_r = self._band_and_reserves(lctx, rctx, lmax, rmax)
        x0 = lmax + band + rmax                  # right tree's root lands here

        scene = Scene()
        for prim in (list(lctx.scene.paths) + list(lctx.scene.polygons)
                     + list(lctx.scene.markers) + list(lctx.scene.labels)
                     + list(lctx.scene.rasters)):
            scene.add(prim)
        mirrored = _mirror_scene(rctx.scene, x0, yscale)
        for prim in (list(mirrored.paths) + list(mirrored.polygons)
                     + list(mirrored.markers) + list(mirrored.labels)
                     + list(mirrored.rasters)):
            scene.add(prim)

        self._add_links(scene, x0, yscale, lmax + res_l, x0 - rmax - res_r)

        # merge colour keys from both sides, dropping exact duplicates (both
        # trees coloured by the same column would otherwise repeat a legend)
        for src in (lctx.scene, rctx.scene):
            for title, entries in src.legends:
                if (title, entries) not in scene.legends:
                    scene.legends.append((title, entries))
            for cb in src.colorbars:
                if cb not in scene.colorbars:
                    scene.colorbars.append(cb)

        self._add_titles(scene, lmax, x0, ln)

        ctx = RenderContext(self.left.tree, _TangleLayout(x0))
        ctx.scene = scene
        return ctx

    @staticmethod
    def _has_tip_labels(ctx: RenderContext) -> bool:
        return any(lb.role == "tiplab" for lb in ctx.scene.labels)

    # Geometry of a rendered rectangular figure, measured once rather than
    # guessed: the axes spans ~77.5% of the figure width, and ``_set_limits``
    # pads the x range to ~1.21x the layout's max_x.  Together they convert a
    # text width in points into a fraction of the tanglegram's total width.
    _AXES_FRACTION = 0.775
    _XLIM_FACTOR = 1.21

    @staticmethod
    def _widest_label_points(ctx: RenderContext) -> float:
        """Width of this side's widest tip label, in points."""
        labs = [lb for lb in ctx.scene.labels if lb.role == "tiplab" and lb.text]
        if not labs:
            return 0.0
        # the longest string is not always the widest one (glyph widths vary),
        # so measure the few longest rather than trusting character count
        labs.sort(key=lambda lb: len(lb.text), reverse=True)
        try:
            from matplotlib.font_manager import FontProperties
            from matplotlib.textpath import TextPath
            return max(
                TextPath((0, 0), lb.text,
                         prop=FontProperties(size=lb.size)).get_extents().width
                for lb in labs[:4])
        except Exception:                      # no font machinery -- estimate
            return max(len(lb.text) * 0.55 * lb.size for lb in labs[:4])

    @classmethod
    def _label_width_fraction(cls, label_pts: float,
                              fig_width_in: float) -> float:
        """Fraction of the tanglegram's total width a label of ``label_pts``
        takes up once rendered at ``fig_width_in`` inches."""
        axes_pts = cls._AXES_FRACTION * 72.0 * fig_width_in
        return label_pts * cls._XLIM_FACTOR / axes_pts if axes_pts else 0.0

    def _band_and_reserves(self, lctx, rctx, lmax: float,
                           rmax: float) -> Tuple[float, float, float]:
        """``(band, left_reserve, right_reserve)`` in data units.

        Tip labels and links share the middle band. Solving
        ``band = (f_l + f_r) * total + link_min`` with
        ``total = lmax + rmax + band`` sizes the band so the labels actually
        fit, instead of a fixed fraction that overflows on long taxon names.
        """
        trees = lmax + rmax
        pts_l = self._widest_label_points(lctx)
        pts_r = self._widest_label_points(rctx)
        width_in = self._figure_width(pts_l + pts_r)
        self._fig_width = width_in                   # reused by _default_figsize
        if self.gap is not None:                     # explicit override
            band = self.gap * trees
            f_l = f_r = 0.0
        else:
            f_l = self._label_width_fraction(pts_l, width_in)
            f_r = self._label_width_fraction(pts_r, width_in)
            # a label takes ``f`` of the *total* width, and the total depends
            # on the band, so solve band = (f_l + f_r) * (trees + band) + links
            link_min = 0.34 * trees                  # keep crossings legible
            denom = max(1.0 - f_l - f_r, 0.3)        # guard absurd label widths
            band = ((f_l + f_r) * trees + link_min) / denom
        total = trees + band
        if self._link_inset is not None:
            il, ir = self._link_inset
            return band, band * il, band * ir
        if self.gap is not None:                     # proportional fallback
            lab_l = self._has_tip_labels(lctx)
            lab_r = self._has_tip_labels(rctx)
            share = 0.34 if (lab_l and lab_r) else 0.55
            return (band, band * (share if lab_l else 0.02),
                    band * (share if lab_r else 0.02))
        pad = 0.012 * trees
        return band, f_l * total + pad, f_r * total + pad

    def _add_links(self, scene: Scene, x0: float, yscale: float,
                   xl: float, xr: float) -> None:
        left_tips = {n.name: n for n in self.left.tree.leaves() if n.name}
        right_tips = {n.name: n for n in self.right.tree.leaves() if n.name}
        order_l = [n.name for n in self.left.tree.leaves()
                   if n.name in right_tips]
        if not order_l:
            return
        order_r = [n.name for n in self.right.tree.leaves()
                   if n.name in left_tips]

        flags = (_crossing_flags(order_l, order_r)
                 if self._highlight_discordant else {})

        color_fn = None
        spec = self._link_color_by
        if spec and any(spec in left_tips[n].data for n in order_l):
            scale = build_color_scale(
                spec, [left_tips[n].data.get(spec) for n in order_l])

            def color_fn(name):
                return scale.color(left_tips[name].data.get(spec))

            scene.add_legend(scale.title, scale.legend)

        style = self._link_style
        for name in order_l:
            lt, rt = left_tips[name], right_tips[name]
            if color_fn is not None:
                col = color_fn(name)
            elif self._highlight_discordant:
                col = DISCORDANT_COLOR if flags[name] else CONCORDANT_COLOR
            else:
                col = style["color"]
            scene.add(Path([(xl, lt.y), (xr, rt.y * yscale)],
                           color=col, width=style["width"], dash=style["dash"],
                           opacity=style["opacity"],
                           zorder=0.4))          # under branches and labels

    def _add_titles(self, scene: Scene, lmax: float, x0: float,
                    ln: float) -> None:
        if not self.titles:
            return
        y = -0.05 * max(ln, 1) - 0.6            # just above the first row
        left_t = self.titles[0] if len(self.titles) > 0 else ""
        right_t = self.titles[1] if len(self.titles) > 1 else ""
        if left_t:
            scene.add(Label(0.0, y, str(left_t), size=11, ha="left",
                            va="center", color="#333333"))
        if right_t:
            scene.add(Label(x0, y, str(right_t), size=11, ha="right",
                            va="center", color="#333333"))

    @classmethod
    def _figure_width(cls, label_pts: float) -> float:
        """Figure width in inches.

        Long taxon names would otherwise squeeze the two trees into slivers,
        so the figure widens until the labels take at most ~45% of the axes.
        """
        needed = label_pts / 0.45 / (cls._AXES_FRACTION * 72.0)
        return max(11.0, min(needed, 30.0))

    def _default_figsize(self, ctx: RenderContext = None):
        # ~0.26in per row keeps 10pt labels clear of each other without the
        # very tall, narrow page a larger per-row allowance gives on big trees
        n = max(self.left.tree.n_leaves, self.right.tree.n_leaves)
        height = max(3.0, min(0.26 * n, 40.0))
        return (getattr(self, "_fig_width", 11.0), height)
