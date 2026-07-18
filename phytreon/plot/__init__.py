"""Plotting: the :class:`TreeFigure` builder, elements and backends."""
from .figure import (
    TreeFigure,
    RenderContext,
    ColorScale,
    build_color_scale,
)
from .tangle import TangleFigure

__all__ = [
    "TreeFigure",
    "TangleFigure",
    "RenderContext",
    "ColorScale",
    "build_color_scale",
]
