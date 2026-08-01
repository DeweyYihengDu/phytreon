"""Plotting: the :class:`TreeFigure` builder, elements and backends."""
from .figure import (
    TreeFigure,
    RenderContext,
    ColorScale,
    build_color_scale,
)
from .tangle import TangleFigure
from .densi import DensiTreeFigure
from .network import SequenceNetwork
from .panels import PanelFigure, panels
from .splitnet import (SplitNetwork, neighbor_net,
                       neighbornet_ordering)

__all__ = [
    "TreeFigure",
    "TangleFigure",
    "DensiTreeFigure",
    "SequenceNetwork",
    "PanelFigure",
    "panels",
    "SplitNetwork", "neighbor_net", "neighbornet_ordering",
    "RenderContext",
    "ColorScale",
    "build_color_scale",
]
