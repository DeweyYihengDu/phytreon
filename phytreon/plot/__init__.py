"""Plotting: the :class:`TreeFigure` builder, elements and backends."""
from .figure import (
    TreeFigure,
    RenderContext,
    ColorScale,
    build_color_scale,
)
from .tangle import TangleFigure
from .densi import DensiTreeFigure

__all__ = [
    "TreeFigure",
    "TangleFigure",
    "DensiTreeFigure",
    "RenderContext",
    "ColorScale",
    "build_color_scale",
]
