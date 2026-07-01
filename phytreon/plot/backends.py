"""Renderers: turn a :class:`~phytreon.scene.Scene` into a figure.

Two backends consume the *same* scene:

* :func:`render_mpl`    -> a ``matplotlib.figure.Figure`` (publication).
* :func:`render_plotly` -> a ``plotly.graph_objects.Figure`` (interactive).

Because layout already produced final cartesian coordinates and colours
are pre-resolved to hex, these functions are pure translation -- no
phylogenetic logic lives here.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Optional

from .figure import RenderContext

_DASH_MPL = {None: "-", "dash": "--", "dot": ":"}
_VA_MPL = {"top": "top", "center": "center", "bottom": "bottom"}


# --------------------------------------------------------------------------
# matplotlib
# --------------------------------------------------------------------------
def render_mpl(ctx: RenderContext, title: Optional[str] = None,
               figsize=None, ax=None):
    import matplotlib.pyplot as plt

    scene = ctx.scene
    equal = getattr(ctx.layout, "equal_aspect", ctx.layout.is_polar)
    max_x = ctx.layout.max_x
    n = ctx.tree.n_leaves

    if ax is None:
        if figsize is None:
            figsize = (8, 8) if equal else (8, max(2.6, min(0.34 * n, 30)))
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # split "aligned" tracks (clade bars, heatmaps) from the base plot so we
    # can place them just past the tip labels after measuring label widths.
    base_polys = [p for p in scene.polygons if not p.align]
    al_polys = [p for p in scene.polygons if p.align]
    base_paths = [p for p in scene.paths if not p.align]
    al_paths = [p for p in scene.paths if p.align]
    base_labels = [l for l in scene.labels if not l.align]
    al_labels = [l for l in scene.labels if l.align]
    base_rasters = [r for r in scene.rasters if not r.align]
    al_rasters = [r for r in scene.rasters if r.align]

    # -- base layer ------------------------------------------------------
    for r in base_rasters:
        _draw_raster(ax, r)
    for poly in sorted(base_polys, key=lambda p: p.zorder):
        _draw_polygon(ax, poly)
    for p in base_paths:
        _draw_path(ax, p)
    for m in scene.markers:
        ax.plot([m.x], [m.y], linestyle="None", marker=m.marker,
                markersize=m.size, markerfacecolor=m.color,
                markeredgecolor=m.edgecolor or m.color, markeredgewidth=0.6,
                zorder=m.zorder, alpha=m.opacity)
    base_text = [_draw_label(ax, lb) for lb in base_labels]
    tiplabs = [t for t, lb in zip(base_text, base_labels) if lb.role == "tiplab"]

    # frame
    ax.set_axis_off()
    if equal:
        ax.set_aspect("equal")
    elif getattr(ctx.layout, "invert_y", True):
        ax.invert_yaxis()      # first leaf on top (rectangular); dendrogram opts out

    # right-side aligned tracks only make sense for the standard rectangular
    # layout (x = depth); other layouts use plain bounds-based limits.
    kind = getattr(ctx.layout, "kind", "rect")
    has_tracks = (kind == "rect" and not equal
                  and bool(al_polys or al_paths or al_labels or al_rasters))
    track_w = _aligned_extent(al_polys, al_paths, al_labels, al_rasters, max_x) \
        if has_tracks else 0.0

    # limits
    delta = 0.0
    if equal or kind != "rect":
        _set_equal_limits(ax, scene, pad=0.06)
    else:
        # provisional right margin so labels have room to render before measuring
        _set_limits(ax, scene, max_x, right_pad=(track_w + 0.8 * max_x) if has_tracks
                    else 0.18 * max_x)

    # -- measure tip labels, then drop aligned tracks just past them ------
    if has_tracks and tiplabs:
        right = _measure_right(fig, ax, tiplabs)          # pass 1
        if right is not None:
            _set_limits(ax, scene, max_x, right_pad=(right - max_x) + track_w + 0.06 * max_x)
            right = _measure_right(fig, ax, tiplabs) or right   # pass 2 (near-final scale)
            delta = max(0.0, right - max_x) + 0.015 * max_x

    al_text = []
    for poly in sorted(al_polys, key=lambda p: p.zorder):
        _draw_polygon(ax, _shift_poly(poly, delta))
    for p in al_paths:
        _draw_path(ax, _shift_path(p, delta))
    for lb in al_labels:
        al_text.append(_draw_label(ax, replace(lb, x=lb.x + delta)))
    for r in al_rasters:
        _draw_raster(ax, replace(r, x0=r.x0 + delta, x1=r.x1 + delta))

    if has_tracks:
        _set_limits(ax, scene, max_x, right_pad=delta + track_w + 0.04 * max_x)

    if title:
        ax.set_title(title, fontsize=13, pad=10)

    # -- legends (stacked top-right, recorded for the tight bbox) ---------
    from matplotlib.legend import Legend
    from matplotlib.lines import Line2D
    # place legends to the right of the tip labels so they never overlap them
    legend_x = 1.02
    if scene.legends and tiplabs:
        r = _measure_right(fig, ax, tiplabs)
        x0, x1 = ax.get_xlim()
        if r is not None and x1 != x0:
            legend_x = max(1.02, (r - x0) / (x1 - x0) + 0.06)
    extra = []
    y = 1.0
    for lt, entries in scene.legends:
        handles, labels = [], []
        for e in entries:
            lab, col = e[0], e[1]
            mk = e[2] if len(e) > 2 else "o"
            handles.append(Line2D([0], [0], marker=mk, linestyle="None",
                                  markerfacecolor=col, markeredgecolor=col))
            labels.append(str(lab))
        leg = Legend(ax, handles, labels, title=lt,
                     loc="upper left", bbox_to_anchor=(legend_x, y),
                     frameon=False, fontsize=8, title_fontsize=9)
        leg._legend_box.align = "left"
        ax.add_artist(leg)
        extra.append(leg)
        y -= 0.065 * (len(entries) + 2)

    # -- continuous colorbars (stacked below the legends) ----------------
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.cm import ScalarMappable
    for title, vmin, vmax, stops in scene.colorbars:
        cmap = LinearSegmentedColormap.from_list(title or "cb", list(stops))
        sm = ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap=cmap)
        cax = ax.inset_axes([legend_x + 0.01, max(y - 0.22, 0.0), 0.02, 0.2],
                            transform=ax.transAxes)
        cb = fig.colorbar(sm, cax=cax)
        cb.outline.set_visible(False)
        cb.ax.tick_params(labelsize=7, length=2)
        cb.set_label(title, fontsize=9)
        extra.append(cax)
        y -= 0.30

    fig._phytreon_extra_artists = extra + base_text + al_text
    return fig


# -- primitive drawing helpers ---------------------------------------------
def _draw_polygon(ax, poly):
    from matplotlib.patches import Polygon as MplPolygon
    ax.add_patch(MplPolygon(
        list(poly.points), closed=True,
        facecolor=poly.facecolor if poly.facecolor else "none",
        edgecolor=poly.edgecolor if poly.edgecolor else "none",
        alpha=poly.alpha, linewidth=poly.width, zorder=poly.zorder,
        joinstyle="round",
    ))


def _draw_path(ax, p):
    xs = [pt[0] for pt in p.points]
    ys = [pt[1] for pt in p.points]
    ax.plot(xs, ys, color=p.color, linewidth=p.width,
            linestyle=_DASH_MPL.get(p.dash, "-"), solid_capstyle="round",
            solid_joinstyle="round", zorder=p.zorder, alpha=p.opacity)


def _draw_raster(ax, r):
    # build an RGB image from the categorical codes + palette, then imshow.
    # row 0 is the first tip (smallest y); with the y-axis inverted,
    # origin="upper" puts it at the top, matching the tree.
    import numpy as np
    pal = np.array([[int(h.lstrip("#")[i:i + 2], 16) for i in (0, 2, 4)]
                    for h in r.palette], dtype=np.uint8)
    codes = np.clip(np.asarray(r.codes), 0, len(pal) - 1)
    rgb = pal[codes]
    ax.imshow(rgb, extent=(r.x0, r.x1, r.y1, r.y0), origin="upper",
              aspect="auto", interpolation="nearest", zorder=r.zorder)


def _draw_label(ax, lb):
    return ax.text(lb.x, lb.y, lb.text, fontsize=lb.size, color=lb.color,
                   ha=lb.ha, va=_VA_MPL.get(lb.va, "center"),
                   rotation=lb.rotation, rotation_mode="anchor",
                   style="italic" if lb.italic else "normal", zorder=lb.zorder)


# -- alignment helpers ------------------------------------------------------
def _shift_poly(poly, dx):
    return replace(poly, points=[(x + dx, y) for x, y in poly.points])


def _shift_path(p, dx):
    return replace(p, points=[(x + dx, y) for x, y in p.points])


def _aligned_extent(polys, paths, labels, rasters, max_x):
    xs = []
    for p in polys:
        xs += [pt[0] for pt in p.points]
    for p in paths:
        xs += [pt[0] for pt in p.points]
    for l in labels:
        xs.append(l.x)
    for r in rasters:
        xs.append(r.x1)
    if not xs:
        return 0.0
    return max(0.0, max(xs) - max_x) + 0.06 * max_x   # + slack for label text


def _measure_right(fig, ax, artists):
    """Right-most data-x reached by the given text artists, at current scale."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    mx = None
    for t in artists:
        bb = t.get_window_extent(renderer=r)
        for corner in ((bb.x0, bb.y0), (bb.x1, bb.y1), (bb.x1, bb.y0), (bb.x0, bb.y1)):
            dx = inv.transform(corner)[0]
            mx = dx if mx is None else max(mx, dx)
    return mx


def _set_equal_limits(ax, scene, pad=0.08):
    """Symmetric, content-fitting limits for equal-aspect layouts
    (circular / unrooted): no y flip, equal data scaling."""
    xmin, ymin, xmax, ymax = scene.bounds()
    dx = (xmax - xmin) or 1.0
    dy = (ymax - ymin) or 1.0
    ax.set_xlim(xmin - pad * dx, xmax + pad * dx)
    ax.set_ylim(ymin - pad * dy, ymax + pad * dy)


def _set_limits(ax, scene, max_x, right_pad):
    xmin, ymin, xmax, ymax = scene.bounds()
    dy = (ymax - ymin) or 1.0
    padx = 0.03 * (max_x or 1.0)
    ax.set_xlim(xmin - padx, max_x + right_pad)
    lo, hi = ymin - 0.04 * dy - 0.4, ymax + 0.04 * dy + 0.4
    if ax.yaxis_inverted():
        ax.set_ylim(hi, lo)
    else:
        ax.set_ylim(lo, hi)


# --------------------------------------------------------------------------
# plotly
# --------------------------------------------------------------------------
_DASH_PLOTLY = {None: "solid", "dash": "dash", "dot": "dot"}


def render_plotly(ctx: RenderContext, title: Optional[str] = None,
                  height: int = 700, width: Optional[int] = None):
    import plotly.graph_objects as go

    scene = ctx.scene
    polar = getattr(ctx.layout, "equal_aspect", ctx.layout.is_polar)
    fig = go.Figure()

    # plotly can't measure text, so estimate how far the aligned tracks must
    # shift right to clear the tip labels (proportional to the longest name).
    kind = getattr(ctx.layout, "kind", "rect")
    dx_align = 0.0
    if kind == "rect":
        maxlen = max((len(l.text) for l in scene.labels if l.role == "tiplab"),
                     default=0)
        dx_align = ctx.layout.max_x * min(0.8, 0.03 * maxlen)

    def shx(x, aligned):
        return x + dx_align if aligned else x

    # rasters (MSA tracks) -- go.Heatmap (categorical) so the axes are NOT
    # forced to a square pixel aspect the way go.Image would be.
    for r in scene.rasters:
        h, w = r.codes.shape[0], r.codes.shape[1]
        x0 = shx(r.x0, r.align)
        x1 = shx(r.x1, r.align)
        k = len(r.palette)
        cscale = []
        for i, hexc in enumerate(r.palette):
            cscale.append([i / k, hexc])
            cscale.append([(i + 1) / k, hexc])
        fig.add_trace(go.Heatmap(
            z=r.codes, x0=x0 + (x1 - x0) / (2 * w), dx=(x1 - x0) / w,
            y0=r.y0 + (r.y1 - r.y0) / (2 * h), dy=(r.y1 - r.y0) / h,
            zmin=-0.5, zmax=k - 0.5, colorscale=cscale, showscale=False,
            hoverinfo="skip",
        ))

    # polygons
    for poly in sorted(scene.polygons, key=lambda p: p.zorder):
        xs = [shx(pt[0], poly.align) for pt in poly.points] + [shx(poly.points[0][0], poly.align)]
        ys = [pt[1] for pt in poly.points] + [poly.points[0][1]]
        fig.add_trace(go.Scatter(
            x=xs, y=ys, fill="toself",
            fillcolor=_rgba(poly.facecolor, poly.alpha) if poly.facecolor and poly.facecolor != "none" else "rgba(0,0,0,0)",
            line=dict(color=poly.edgecolor or "rgba(0,0,0,0)", width=poly.width),
            mode="lines", hoverinfo="text" if poly.label else "skip",
            text=poly.label, showlegend=False,
        ))

    # paths -- group all branch segments into one trace with None breaks
    bx, by = [], []
    colored = {}
    for p in scene.paths:
        key = (p.color, p.width, p.dash)
        colored.setdefault(key, [[], []])
        for x, y in p.points:
            colored[key][0].append(x)
            colored[key][1].append(y)
        colored[key][0].append(None)
        colored[key][1].append(None)
    for (color, w, dash), (xs, ys) in colored.items():
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode="lines", showlegend=False, hoverinfo="skip",
            line=dict(color=color, width=w, dash=_DASH_PLOTLY.get(dash, "solid")),
        ))

    # markers
    if scene.markers:
        fig.add_trace(go.Scatter(
            x=[m.x for m in scene.markers],
            y=[m.y for m in scene.markers],
            mode="markers", showlegend=False,
            marker=dict(size=[m.size for m in scene.markers],
                        color=[m.color for m in scene.markers],
                        line=dict(width=0.5,
                                  color=[m.edgecolor or m.color for m in scene.markers])),
            text=[m.label or "" for m in scene.markers],
            hoverinfo="text",
        ))

    # labels via annotations (supports rotation)
    for lb in scene.labels:
        fig.add_annotation(
            x=shx(lb.x, lb.align), y=lb.y, text=lb.text, showarrow=False,
            font=dict(size=lb.size, color=lb.color),
            xanchor={"left": "left", "right": "right", "center": "center"}[lb.ha],
            yanchor={"top": "top", "bottom": "bottom", "center": "middle"}[lb.va],
            textangle=-lb.rotation,
        )

    # legends -> dummy traces so plotly draws a clickable legend
    _mk2plotly = {"o": "circle", "s": "square", "^": "triangle-up",
                  "D": "diamond", "v": "triangle-down", "P": "cross",
                  "X": "x", "*": "star", "<": "triangle-left", ">": "triangle-right"}
    for lt, entries in scene.legends:
        for e in entries:
            label, color = e[0], e[1]
            sym = _mk2plotly.get(e[2], "circle") if len(e) > 2 else "circle"
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode="markers", name=str(label),
                legendgroup=lt, legendgrouptitle_text=lt,
                marker=dict(size=8, color=color, symbol=sym), showlegend=True,
            ))

    # continuous colorbars (stacked on the right)
    for ci, (title, vmin, vmax, stops) in enumerate(scene.colorbars):
        cs = [[i / (len(stops) - 1), c] for i, c in enumerate(stops)]
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="markers", showlegend=False,
            marker=dict(colorscale=cs, cmin=vmin, cmax=vmax, color=[vmin],
                        colorbar=dict(title=title, len=0.4, thickness=12,
                                      x=1.02, y=0.4 - ci * 0.45, yanchor="top")),
        ))

    yaxis = dict(visible=False)
    if polar:
        yaxis["scaleanchor"] = "x"   # equal aspect
    else:
        yaxis["autorange"] = "reversed"
    fig.update_layout(
        title=title, height=height, width=width,
        xaxis=dict(visible=False), yaxis=yaxis,
        plot_bgcolor="white", margin=dict(l=10, r=10, t=40, b=10),
        legend=dict(itemsizing="constant"),
    )
    return fig


def _rgba(hex_color: str, alpha: float) -> str:
    h = hex_color.lstrip("#")
    if len(h) != 6:
        return hex_color
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"
