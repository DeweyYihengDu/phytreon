"""Several figures laid out as one multi-panel figure.

A grid of small trees answers a question a single big tree cannot: *do these
gene families tell the same story?* Eleven cramped panels sharing one colour
key let a reader compare eleven histories at a glance, which is why that
layout keeps appearing in comparative-genomics papers::

    pt.panels([pt.TreeFigure(t).titled(name) for name, t in gene_trees],
              ncols=4, share_legend=True).save("panel.pdf")

Any figure this package can draw works as a panel -- trees, tanglegrams,
DensiTree clouds, sequence networks -- because each one renders into an axes
it is handed.
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple


class PanelFigure:
    """A grid of figures rendered as one image.

    ``share_legend`` collects the panels' colour keys, draws each distinct one
    once beside the grid, and suppresses the per-panel copies. Panels of the
    same data type normally carry the *same* key, and repeating it in every
    cell wastes the space the panels need.
    """

    def __init__(self, figures: Sequence[object], *, ncols: int = 3,
                 share_legend: bool = True, titles: Optional[Sequence[str]] = None,
                 panel_size: Tuple[float, float] = (3.2, 2.8),
                 label_panels: bool = False):
        self.figures = list(figures)
        if not self.figures:
            raise ValueError("panels() needs at least one figure")
        self.ncols = max(1, int(ncols))
        self.share_legend = share_legend
        self.titles = list(titles) if titles is not None else None
        if self.titles is not None and len(self.titles) != len(self.figures):
            raise ValueError(
                f"got {len(self.titles)} titles for {len(self.figures)} panels")
        self.panel_size = panel_size
        #: label panels a, b, c ... in the corner, as a journal would
        self.label_panels = label_panels
        self.title: Optional[str] = None

    def titled(self, title: str) -> "PanelFigure":
        self.title = title
        return self

    @property
    def nrows(self) -> int:
        return -(-len(self.figures) // self.ncols)     # ceil

    # -- rendering -------------------------------------------------------
    def draw(self, backend: str = "mpl", figsize=None):
        if backend not in ("mpl", "matplotlib", "static"):
            raise ValueError(
                "panels render through matplotlib only; plotly has no "
                "equivalent multi-panel target here")
        import matplotlib.pyplot as plt
        from .backends import render_mpl

        pw, ph = self.panel_size
        legend_w = 1.8 if self.share_legend else 0.0
        if figsize is None:
            figsize = (self.ncols * pw + legend_w, self.nrows * ph)
        fig, axes = plt.subplots(self.nrows, self.ncols, figsize=figsize,
                                 squeeze=False)

        collected: List[Tuple[str, list]] = []
        swatches = {}
        for i, panel in enumerate(self.figures):
            ax = axes[i // self.ncols][i % self.ncols]
            ctx = panel._build()
            if self.share_legend:
                # hold the keys back and draw them once beside the grid
                collected += [e for e in ctx.scene.legends if e not in collected]
                swatches.update(ctx.scene.legend_swatch)
                ctx.scene.legends = []
                ctx.scene.colorbars = []
            title = (self.titles[i] if self.titles is not None
                     else getattr(panel, "title", None))
            render_mpl(ctx, title=title, ax=ax)
            if self.label_panels:
                ax.text(-0.02, 1.04, chr(ord("a") + i), transform=ax.transAxes,
                        fontsize=12, fontweight="bold", va="bottom", ha="right")

        for j in range(len(self.figures), self.nrows * self.ncols):
            axes[j // self.ncols][j % self.ncols].set_axis_off()

        extra = []
        if self.share_legend and collected:
            extra = self._draw_shared_legend(fig, axes[0][-1], collected,
                                             swatches)
        if self.title:
            fig.suptitle(self.title, fontsize=14)
        fig.tight_layout()
        fig._phytreon_extra_artists = extra or None
        return fig

    @staticmethod
    def _draw_shared_legend(fig, ax, entries, swatches):
        from matplotlib.legend import Legend
        from matplotlib.lines import Line2D
        from matplotlib.patches import Rectangle
        out = []
        y = 1.0
        for title, items in entries:
            handles, labels = [], []
            patch = swatches.get(title) == "patch"
            for e in items:
                lab, col = e[0], e[1]
                if patch:
                    handles.append(Rectangle((0, 0), 1, 1, facecolor=col,
                                             edgecolor="none"))
                else:
                    mk = e[2] if len(e) > 2 else "o"
                    handles.append(Line2D([0], [0], marker=mk, linestyle="None",
                                          markerfacecolor=col,
                                          markeredgecolor=col))
                labels.append(str(lab))
            leg = Legend(ax, handles, labels, title=title, loc="upper left",
                         bbox_to_anchor=(1.04, y), frameon=False,
                         fontsize=8, title_fontsize=9)
            leg._legend_box.align = "left"
            ax.add_artist(leg)
            out.append(leg)
            y -= 0.09 * (len(items) + 2)
        return out

    def save(self, path: str, dpi: int = 300, **kwargs) -> str:
        """Save the grid; SVG keeps editable text, as elsewhere in phytreon."""
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        ext = path.lower().rsplit(".", 1)[-1]
        fig = self.draw(**kwargs)
        extra = getattr(fig, "_phytreon_extra_artists", None)
        rc = {"svg.fonttype": "none"} if ext == "svg" else {}
        with mpl.rc_context(rc):
            fig.savefig(path, bbox_inches="tight", dpi=dpi,
                        bbox_extra_artists=extra)
        plt.close(fig)          # save() never hands the figure back -- see
        return path              # the same fix and note in plot/figure.py

    def show(self, **kwargs):
        import matplotlib.pyplot as plt
        fig = self.draw(**kwargs)
        plt.show()
        return fig


def panels(figures: Sequence[object], **kwargs) -> PanelFigure:
    """Lay several figures out as one multi-panel figure."""
    return PanelFigure(figures, **kwargs)
