"""Layout engines: topology -> display coordinates."""
from .base import Layout
from .rectangular import RectangularLayout, SlantedLayout, DendrogramLayout
from .circular import CircularLayout, CircularSlantedLayout, InwardCircularLayout
from .unrooted import EqualAngleLayout, DaylightLayout

#: registry used by ``TreeFigure(layout=...)``
LAYOUTS = {
    "rectangular": RectangularLayout,
    "rect": RectangularLayout,
    "phylogram": RectangularLayout,
    "slanted": SlantedLayout,
    "cladogram": SlantedLayout,
    "dendrogram": DendrogramLayout,
    "circular": CircularLayout,
    "fan": CircularLayout,
    "radial": CircularLayout,
    "circular_slanted": CircularSlantedLayout,
    "slanted_circular": CircularSlantedLayout,
    "fan_slanted": CircularSlantedLayout,
    "inward_circular": InwardCircularLayout,
    "inward": InwardCircularLayout,
    "unrooted": DaylightLayout,
    "daylight": DaylightLayout,
    "equal_angle": EqualAngleLayout,
}


def get_layout(name: str, **kwargs) -> Layout:
    key = name.lower()
    if key not in LAYOUTS:
        raise ValueError(
            f"unknown layout {name!r}; available: {sorted(set(LAYOUTS))}"
        )
    return LAYOUTS[key](**kwargs)


__all__ = ["Layout", "RectangularLayout", "SlantedLayout", "DendrogramLayout",
           "CircularLayout", "CircularSlantedLayout", "InwardCircularLayout",
           "EqualAngleLayout", "DaylightLayout", "LAYOUTS", "get_layout"]
