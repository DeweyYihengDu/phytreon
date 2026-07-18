"""The :class:`TreeFigure` builder.

A figure is assembled by chaining methods on a tree::

    (TreeFigure(tr, layout="circular")
        .highlight(node=mrca, fill="#a6cee3")
        .tip_points(color="habitat", size=9)
        .tip_labels()
        .support_labels()
        .save("tree.pdf"))        # static  (matplotlib)
                                  # .save("tree.html") -> interactive (plotly)

Each method adds one visual element and returns the figure, so calls chain.
A figure is *lazy*: nothing is computed until a backend asks for the scene.
Layout computation and rendering are fully decoupled (see
:mod:`phytreon.scene`), which is what lets one figure drive both backends.
"""
from __future__ import annotations

from typing import Callable, List, Optional, Tuple

from ..core.tree import Tree
from ..layout import get_layout, Layout
from ..scene import Scene


# --------------------------------------------------------------------------
# colour scales
# --------------------------------------------------------------------------
class ColorScale:
    """Maps raw data values to colours and emits legend entries."""

    def __init__(self, title: str, mapping: Callable[[object], str],
                 legend: List[Tuple[str, str]], continuous: bool,
                 vmin: float = None, vmax: float = None):
        self.title = title
        self._map = mapping
        self.legend = legend
        self.continuous = continuous
        self.vmin = vmin
        self.vmax = vmax

    def color(self, value) -> str:
        return self._map(value)

    def gradient(self, n: int = 32) -> List[str]:
        """n hex colours sampled low->high across the continuous range."""
        if self.vmin is None or self.vmax is None:
            return [c for _, c in self.legend]
        span = self.vmax - self.vmin
        return [self.color(self.vmin + span * i / (n - 1)) for i in range(n)]


def _isnan(v) -> bool:
    try:
        return v != v          # only NaN is unequal to itself
    except Exception:
        return False


def is_numeric(v) -> bool:
    """True for real numbers, including numpy int64/float64 (which fail
    ``isinstance(v, (int, float))``); False for bool (a Number subclass)."""
    import numbers
    return isinstance(v, numbers.Number) and not isinstance(v, bool)


def build_color_scale(title: str, values, cmap=None,
                      palette: str = "curated") -> ColorScale:
    """Build a colour scale (categorical or continuous) from raw values.

    * ``palette`` -- name of the qualitative palette (default ``"curated"`` =
      the eight-hue colourblind-safe default; also ``"hue"`` for the raw HCL
      hue wheel, ``"set2"``, ``"dark2"``, ``"tab10"``).
    * ``cmap`` -- continuous colour spec: ``None`` for the default blue
      gradient, a ``(low, high)`` hex pair, or a matplotlib colormap name.
    """
    from .palettes import categorical_palette, continuous_mapper

    present = [v for v in values if v is not None and not _isnan(v)]
    numeric = len(present) > 0 and all(is_numeric(v) for v in present)
    if numeric:
        vmin, vmax = min(present), max(present)
        mapping = continuous_mapper(vmin, vmax, cmap)
        legend = []
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            val = vmin + frac * (vmax - vmin)
            legend.append((f"{val:g}", mapping(val)))
        return ColorScale(title, mapping, legend, continuous=True,
                          vmin=vmin, vmax=vmax)

    # categorical -- one colour per level
    cats = sorted({str(v) for v in present})
    colors = categorical_palette(len(cats), palette)
    lut = {c: colors[i] for i, c in enumerate(cats)}

    def mapping(v, _lut=lut):
        return _lut.get(str(v), "#cccccc") if v is not None else "#cccccc"

    legend = [(c, lut[c]) for c in cats]
    return ColorScale(title, mapping, legend, continuous=False)


# --------------------------------------------------------------------------
# render context
# --------------------------------------------------------------------------
class RenderContext:
    """Shared state handed to every element while a figure is built."""

    def __init__(self, tree: Tree, layout: Layout):
        self.tree = tree
        self.layout = layout
        self.scene = Scene()
        # radial bookkeeping for circular ring tracks so that rings stack
        # outward without overlap and tip labels sit outside them.
        base = getattr(layout, "inner_radius", 0.0) + layout.max_x
        self.ring_base = base          # where the first ring starts
        self.ring_cursor = base        # running outer edge while drawing rings
        self.outer_radius = base       # final outer edge (labels go beyond this)
        # rectangular right-side tracks (heatmap / bars / alignment) stack along x
        self.track_cursor = layout.max_x

    def add_scale(self, scale) -> None:
        """Register a ColorScale: continuous -> colorbar, categorical -> legend."""
        if scale is None:
            return
        if scale.continuous and scale.vmin is not None:
            self.scene.add_colorbar(scale.title, scale.vmin, scale.vmax,
                                    scale.gradient(32))
        else:
            self.scene.add_legend(scale.title, scale.legend)

    # -- aesthetic resolution -------------------------------------------
    def is_data_column(self, spec, nodes) -> bool:
        return isinstance(spec, str) and any(spec in n.data for n in nodes)

    def resolve_shape(self, spec, nodes, default: str = "o"):
        """Map a categorical column to marker shapes.

        Returns ``(func, legend)`` where ``func(node) -> marker`` and
        ``legend`` is a list of ``(label, color, marker)`` triples (or None
        for a constant shape).
        """
        markers = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">"]
        if spec is None:
            return (lambda n: default), None
        if isinstance(spec, str) and any(spec in n.data for n in nodes):
            cats = sorted({str(n.data.get(spec)) for n in nodes
                           if n.data.get(spec) is not None})
            lut = {c: markers[i % len(markers)] for i, c in enumerate(cats)}
            legend = [(c, "#444444", lut[c]) for c in cats]
            return (lambda n: lut.get(str(n.data.get(spec)), default)), (spec, legend)
        return (lambda n: spec), None        # literal marker

    def resolve_color(self, spec, nodes, default: str = "black",
                      palette: str = "curated", cmap=None):
        """Return ``(func, scale)``.

        ``func(node) -> color``.  ``scale`` is a :class:`ColorScale` if the
        spec mapped a data column (so the caller can register a legend),
        else ``None``.  ``palette`` / ``cmap`` choose the categorical /
        continuous colours when ``spec`` names a data column.
        """
        if spec is None:
            return (lambda n: default), None
        if callable(spec):
            return (lambda n: spec(n)), None
        if self.is_data_column(spec, nodes):
            scale = build_color_scale(spec, [n.data.get(spec) for n in nodes],
                                      cmap=cmap, palette=palette)
            return (lambda n: scale.color(n.data.get(spec))), scale
        # literal colour
        return (lambda n: spec), None


# --------------------------------------------------------------------------
# elements
# --------------------------------------------------------------------------
class _Element:
    """Base class for every visual element added to a figure."""

    def apply(self, ctx: RenderContext) -> None:    # pragma: no cover
        raise NotImplementedError


# --------------------------------------------------------------------------
# shared render / export plumbing
# --------------------------------------------------------------------------
class _Renderable:
    """Draw/save/show for anything that can build a :class:`RenderContext`.

    Shared by :class:`TreeFigure` and
    :class:`~phytreon.plot.tangle.TangleFigure` so both get identical backend
    dispatch and export behaviour (notably the editable-text SVG handling).
    """

    title: Optional[str] = None

    def _build(self) -> RenderContext:              # pragma: no cover
        raise NotImplementedError

    def _default_figsize(self, ctx: RenderContext):
        """Figure size to use when the caller did not pass one (None = let
        the backend decide)."""
        return None

    def draw(self, backend: str = "mpl", **kwargs):
        ctx = self._build()
        if backend in ("mpl", "matplotlib", "static"):
            from .backends import render_mpl
            if kwargs.get("figsize") is None:
                figsize = self._default_figsize(ctx)
                if figsize is not None:
                    kwargs["figsize"] = figsize
            return render_mpl(ctx, title=self.title, **kwargs)
        if backend in ("plotly", "interactive", "html"):
            from .backends import render_plotly
            return render_plotly(ctx, title=self.title, **kwargs)
        raise ValueError(f"unknown backend {backend!r}")

    def save(self, path: str, dpi: int = 300, **kwargs) -> str:
        """Save to ``path``; backend chosen from the file extension.

        Supports ``.png`` / ``.jpg`` (raster), ``.pdf`` / ``.svg`` (vector),
        and ``.html`` (interactive plotly). ``dpi`` controls raster
        resolution; remaining kwargs (e.g. ``figsize``) are forwarded to the
        renderer.

        SVG output keeps every label as a real ``<text>`` element (rather than
        outlined vector paths), so the figure stays fully editable after
        importing it into PowerPoint (Insert → Picture, then Graphics Format →
        Convert to Shape), Illustrator, Inkscape, etc. -- you can recolour,
        move, and re-type the text.
        """
        ext = path.lower().rsplit(".", 1)[-1]
        if ext == "html":
            self.draw(backend="plotly", **kwargs).write_html(path)
        else:  # pdf/svg/png/jpg -> matplotlib
            import matplotlib as mpl
            fig = self.draw(backend="mpl", **kwargs)
            extra = getattr(fig, "_phytreon_extra_artists", None)
            # for SVG, emit editable <text> (not glyph outlines) so labels
            # remain real, re-typeable text in PowerPoint / vector editors
            rc = {"svg.fonttype": "none"} if ext == "svg" else {}
            with mpl.rc_context(rc):
                fig.savefig(path, bbox_inches="tight", dpi=dpi,
                            bbox_extra_artists=extra)
        return path

    def show(self, backend: str = "mpl", **kwargs):
        fig = self.draw(backend=backend, **kwargs)
        if backend in ("mpl", "matplotlib", "static"):
            import matplotlib.pyplot as plt
            plt.show()
        else:
            fig.show()
        return fig


# --------------------------------------------------------------------------
# the figure object
# --------------------------------------------------------------------------
class TreeFigure(_Renderable):
    """A composable tree figure.

    ``TreeFigure(tree, layout=...)`` starts with the branch skeleton already
    drawn; chain methods (``tip_labels``, ``tip_points``, ``heatmap`` …) to
    add elements, then ``save``/``show``.  Set ``skeleton=False`` to start
    empty (e.g. to draw branches yourself with a custom colour).
    """

    def __init__(self, tree: Tree, layout: str = "rectangular",
                 *, skeleton: bool = True, **layout_kwargs):
        self.tree = tree
        self.layout_name = layout
        self.layout_kwargs = layout_kwargs
        self._elements: List[_Element] = []
        self.title: Optional[str] = None
        if skeleton:
            self.branches()

    def add(self, element: _Element) -> "TreeFigure":
        """Append a custom :class:`_Element` (for extensions)."""
        if not isinstance(element, _Element):
            raise TypeError(f"expected an element, got {type(element)!r}")
        self._elements.append(element)
        return self

    # -- composition methods (each returns self so calls chain) ----------
    def branches(self, color="black", size: float = 1.0) -> "TreeFigure":
        """Draw (or restyle) the tree skeleton. ``size`` is a single width
        applied to every branch in the tree -- there is no per-branch/data
        variation. A figure only ever has one skeleton layer: calling this
        again (e.g. to override the ``skeleton=True`` default added by
        ``__init__``) replaces it in place rather than drawing a second,
        overlapping set of lines.
        """
        from .elements import _Branches
        self._elements = [e for e in self._elements if not isinstance(e, _Branches)]
        return self.add(_Branches(color=color, size=size))

    def tip_labels(self, **kwargs) -> "TreeFigure":
        from .elements import _TipLabels
        return self.add(_TipLabels(**kwargs))

    def node_labels(self, **kwargs) -> "TreeFigure":
        from .elements import _NodeLabels
        return self.add(_NodeLabels(**kwargs))

    def support_labels(self, **kwargs) -> "TreeFigure":
        from .elements import _NodeLabels
        kwargs.setdefault("attr", "support")
        return self.add(_NodeLabels(**kwargs))

    def tip_points(self, **kwargs) -> "TreeFigure":
        from .elements import _Points
        kwargs["which"] = "tip"
        return self.add(_Points(**kwargs))

    def node_points(self, **kwargs) -> "TreeFigure":
        from .elements import _Points
        kwargs["which"] = "node"
        return self.add(_Points(**kwargs))

    def points(self, **kwargs) -> "TreeFigure":
        from .elements import _Points
        kwargs.setdefault("which", "all")
        return self.add(_Points(**kwargs))

    def highlight(self, **kwargs) -> "TreeFigure":
        from .elements import _Highlight
        return self.add(_Highlight(**kwargs))

    def clade_label(self, label: str, **kwargs) -> "TreeFigure":
        from .elements import _CladeLabel
        return self.add(_CladeLabel(label=label, **kwargs))

    def heatmap(self, data, **kwargs) -> "TreeFigure":
        from .elements import _Heatmap
        return self.add(_Heatmap(data, **kwargs))

    def ring(self, data, **kwargs) -> "TreeFigure":
        from .elements import _Ring
        return self.add(_Ring(data, **kwargs))

    def bar_track(self, data, column, **kwargs) -> "TreeFigure":
        from .elements import _BarTrack
        return self.add(_BarTrack(data, column, **kwargs))

    def alignment(self, alignment, **kwargs) -> "TreeFigure":
        from .elements import _Alignment
        return self.add(_Alignment(alignment, **kwargs))

    def painted_branches(self, **kwargs) -> "TreeFigure":
        from .elements import _PaintedBranches
        return self.add(_PaintedBranches(**kwargs))

    def node_pies(self, **kwargs) -> "TreeFigure":
        from .elements import _NodePies
        return self.add(_NodePies(**kwargs))

    def time_axis(self, **kwargs) -> "TreeFigure":
        from .elements import _TimeAxis
        return self.add(_TimeAxis(**kwargs))

    def titled(self, title: str) -> "TreeFigure":
        self.title = title
        return self

    # -- building / rendering -------------------------------------------
    def _build(self) -> RenderContext:
        layout = get_layout(self.layout_name, **self.layout_kwargs)
        layout.apply(self.tree)
        ctx = RenderContext(self.tree, layout)
        # pre-pass: reserve radial space for every ring track so tip labels
        # (whatever order they were added) land outside all the rings.
        for el in self._elements:
            reserve = getattr(el, "reserved_extent", None)
            if reserve is not None:
                ctx.outer_radius += reserve(layout)
        for el in self._elements:
            el.apply(ctx)
        return ctx
