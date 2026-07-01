"""Colour palettes.

Two palette families, both returning plain hex strings so the backends stay
colour-agnostic:

* **categorical** colours are evenly spaced *hues* in the perceptual HCL
  (polar-CIELUV) space at fixed chroma/luminance (c=100, l=65). For n=3 this
  yields a balanced salmon / green / blue ``#F8766D #00BA38 #619CFF``.
* **continuous** scales default to a dark-blue -> light-blue gradient
  (``#132B43`` -> ``#56B1F7``), interpolated in linear light.
"""
from __future__ import annotations

import math
from typing import List, Tuple

# default continuous gradient (low -> high)
_BLUE_LOW = "#132B43"
_BLUE_HIGH = "#56B1F7"

# a few tasteful named qualitative palettes (ColorBrewer)
NAMED_PALETTES = {
    "set2": ["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854",
             "#ffd92f", "#e5c494", "#b3b3b3"],
    "dark2": ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
              "#e6ab02", "#a6761d", "#666666"],
    "tab10": ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
              "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"],
}


# --------------------------------------------------------------------------
# sRGB <-> linear helpers
# --------------------------------------------------------------------------
def _srgb_to_linear(c: float) -> float:
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def _linear_to_srgb(c: float) -> float:
    c = max(0.0, min(1.0, c))
    return c * 12.92 if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055


def _hex_to_rgb01(h: str) -> Tuple[float, float, float]:
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))  # type: ignore


def _rgb01_to_hex(rgb: Tuple[float, float, float]) -> str:
    return "#" + "".join(f"{int(round(max(0.0, min(1.0, c)) * 255)):02x}" for c in rgb)


def lerp_color(c1: str, c2: str, t: float) -> str:
    """Interpolate two hex colours in linear light (gamma-correct)."""
    a = [_srgb_to_linear(c) for c in _hex_to_rgb01(c1)]
    b = [_srgb_to_linear(c) for c in _hex_to_rgb01(c2)]
    mix = [a[i] + (b[i] - a[i]) * t for i in range(3)]
    return _rgb01_to_hex(tuple(_linear_to_srgb(c) for c in mix))


# --------------------------------------------------------------------------
# HCL (polar CIELUV) -> hex, the hue-wheel palette
# --------------------------------------------------------------------------
# D65 white point
_XN, _YN, _ZN = 95.047, 100.0, 108.883
_UN = 4 * _XN / (_XN + 15 * _YN + 3 * _ZN)
_VN = 9 * _YN / (_XN + 15 * _YN + 3 * _ZN)


def _hcl_to_hex(h_deg: float, c: float = 100.0, l: float = 65.0) -> str:
    if l <= 0:
        return "#000000"
    hr = math.radians(h_deg)
    u = c * math.cos(hr)
    v = c * math.sin(hr)
    # CIELUV -> XYZ
    if l > 8:
        Y = _YN * ((l + 16) / 116) ** 3
    else:
        Y = _YN * l / 903.3
    up = u / (13 * l) + _UN
    vp = v / (13 * l) + _VN
    X = Y * 9 * up / (4 * vp)
    Z = Y * (12 - 3 * up - 20 * vp) / (4 * vp)
    # XYZ (0..100) -> linear sRGB
    x, y, z = X / 100.0, Y / 100.0, Z / 100.0
    r = 3.2404542 * x - 1.5371385 * y - 0.4985314 * z
    g = -0.9692660 * x + 1.8760108 * y + 0.0415560 * z
    b = 0.0556434 * x - 0.2040259 * y + 1.0572252 * z
    return _rgb01_to_hex((_linear_to_srgb(r), _linear_to_srgb(g), _linear_to_srgb(b)))


def hue_palette(n: int, c: float = 100.0, l: float = 65.0,
                h_start: float = 15.0, h_end: float = 375.0) -> List[str]:
    """n evenly spaced HCL hues around the colour wheel."""
    if n <= 0:
        return []
    if n == 1:
        return [_hcl_to_hex(h_start, c, l)]
    step = (h_end - h_start) / n
    return [_hcl_to_hex((h_start + step * i) % 360, c, l) for i in range(n)]


def categorical_palette(n: int, name: str = "hue") -> List[str]:
    if name == "hue":
        return hue_palette(n)
    base = NAMED_PALETTES.get(name)
    if base is None:
        raise ValueError(f"unknown palette {name!r}; choose hue/"
                         + "/".join(NAMED_PALETTES))
    return [base[i % len(base)] for i in range(n)]


def continuous_mapper(vmin: float, vmax: float, cmap=None):
    """Return f(value)->hex for a continuous scale.

    ``cmap`` may be None (default blue gradient), a (low, high) hex pair, or
    the name of a matplotlib colormap.
    """
    span = (vmax - vmin) or 1.0

    if cmap is None:
        lo, hi = _BLUE_LOW, _BLUE_HIGH
    elif isinstance(cmap, (tuple, list)) and len(cmap) == 2:
        lo, hi = cmap
    else:
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
        colormap = cm.get_cmap(cmap)

        def mp(v, _lo=vmin, _span=span, _cm=colormap):
            if v is None:
                return "#cccccc"
            return mcolors.to_hex(_cm((float(v) - _lo) / _span))
        return mp

    def mp(v, _lo=vmin, _span=span, _a=lo, _b=hi):
        if v is None:
            return "#cccccc"
        return lerp_color(_a, _b, max(0.0, min(1.0, (float(v) - _lo) / _span)))
    return mp
