"""Backend-agnostic drawing primitives.

The whole point of this module is to let the *same* layout drive both a
static (matplotlib) and an interactive (plotly) renderer.  Layout and
geoms emit these dumb primitives in final display coordinates; a backend
only has to know how to draw a polyline, a marker, a text label and a
filled polygon.  Nothing here imports matplotlib or plotly.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

XY = Tuple[float, float]


@dataclass
class Path:
    """A polyline (used for branches, arcs, connectors)."""
    points: Sequence[XY]
    color: str = "black"
    width: float = 1.0
    dash: Optional[str] = None          # None | "dash" | "dot"
    zorder: float = 1.0
    opacity: float = 1.0
    align: bool = False                 # shift right to clear the tip-label band


@dataclass
class Marker:
    """A single point glyph (tip/node points)."""
    x: float
    y: float
    size: float = 6.0
    color: str = "black"
    edgecolor: Optional[str] = None
    marker: str = "o"
    zorder: float = 3.0
    opacity: float = 1.0
    label: Optional[str] = None         # hover text (interactive backend)


@dataclass
class Label:
    """A text label."""
    x: float
    y: float
    text: str
    size: float = 10.0
    color: str = "black"
    ha: str = "left"                    # left | center | right
    va: str = "center"                  # top | center | bottom
    rotation: float = 0.0               # degrees
    zorder: float = 4.0
    italic: bool = False
    align: bool = False                 # shift right to clear the tip-label band
    role: str = ""                      # "tiplab" labels drive alignment measurement


@dataclass
class Polygon:
    """A filled (and/or stroked) polygon: clade highlights, heatmap cells, wedges."""
    points: Sequence[XY]
    facecolor: Optional[str] = "none"
    edgecolor: Optional[str] = None
    alpha: float = 1.0
    width: float = 0.0
    zorder: float = 0.5
    label: Optional[str] = None         # hover text
    align: bool = False                 # shift right to clear the tip-label band
    rounded: bool = False               # draw as a soft rounded rectangle


@dataclass
class Raster:
    """A categorical image block (used for MSA tracks).

    Stored as an integer ``codes`` grid (H x W) plus a ``palette`` (hex per
    code), so matplotlib can imshow an RGB array and plotly can use a Heatmap
    (which, unlike go.Image, does not force a square pixel aspect ratio).
    Rendered as one raster, so a full alignment stays fast and the file small.
    """
    codes: object                        # numpy (H, W) int, values index palette
    palette: list                        # list of hex colours, indexed by code
    x0: float
    x1: float
    y0: float
    y1: float
    zorder: float = 2.0
    align: bool = False


@dataclass
class Scene:
    """Accumulates primitives plus a few rendering hints."""
    paths: List[Path] = field(default_factory=list)
    markers: List[Marker] = field(default_factory=list)
    labels: List[Label] = field(default_factory=list)
    polygons: List[Polygon] = field(default_factory=list)
    rasters: List[Raster] = field(default_factory=list)
    # legend entries: (title, [(label, color), ...])
    legends: List[Tuple[str, List[Tuple[str, str]]]] = field(default_factory=list)
    # continuous colorbars: (title, vmin, vmax, [hex stops low->high])
    colorbars: List[Tuple[str, float, float, List[str]]] = field(default_factory=list)

    def add(self, primitive) -> None:
        if isinstance(primitive, Path):
            self.paths.append(primitive)
        elif isinstance(primitive, Marker):
            self.markers.append(primitive)
        elif isinstance(primitive, Label):
            self.labels.append(primitive)
        elif isinstance(primitive, Polygon):
            self.polygons.append(primitive)
        elif isinstance(primitive, Raster):
            self.rasters.append(primitive)
        else:
            raise TypeError(f"unknown primitive {type(primitive)!r}")

    def add_legend(self, title: str, entries: List[Tuple[str, str]]) -> None:
        if entries:
            self.legends.append((title, entries))

    def add_colorbar(self, title: str, vmin: float, vmax: float,
                     stops: List[str]) -> None:
        if stops:
            self.colorbars.append((title, vmin, vmax, stops))

    def bounds(self) -> Tuple[float, float, float, float]:
        """(xmin, ymin, xmax, ymax) over every primitive."""
        xs: List[float] = []
        ys: List[float] = []
        for p in self.paths:
            for x, y in p.points:
                xs.append(x); ys.append(y)
        for poly in self.polygons:
            for x, y in poly.points:
                xs.append(x); ys.append(y)
        for m in self.markers:
            xs.append(m.x); ys.append(m.y)
        for lb in self.labels:
            xs.append(lb.x); ys.append(lb.y)
        for r in self.rasters:
            xs.extend([r.x0, r.x1]); ys.extend([r.y0, r.y1])
        if not xs:
            return (0.0, 0.0, 1.0, 1.0)
        return (min(xs), min(ys), max(xs), max(ys))
