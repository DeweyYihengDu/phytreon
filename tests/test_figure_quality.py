"""Do the figures come out publishable?

Not "does it render without raising" -- these check what a reviewer sees
first: names on top of each other, names on top of the data, text too small
to survive the journal's downscaling, a circle drawn as an oval, and strokes
too thin to print.

Two measurement details do most of the work. Overlap is taken off real glyph
ink rather than the layout box, and off *oriented* boxes -- a circular tree's
labels are rotated, and matplotlib's window extent for rotated text is the
axis-aligned envelope, nearly twice the ink in each direction for a
45-degree label, so comparing envelopes reports collisions between labels
that are nowhere near each other. And text sitting on a filled shape is not
by itself a fault, since a clade highlight is drawn behind its labels on
purpose; what matters is the luminance contrast between the two.
"""
import math

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath

import phytreon as pt


def _ink(text, size, weight, style):
    prop = FontProperties(size=size, weight=weight, style=style)
    try:
        ext = TextPath((0, 0), text, size=size, prop=prop).get_extents()
    except Exception:                                  # pragma: no cover
        return 0.6 * size * len(text), 0.7 * size
    return ext.width, ext.height


def _corners(t, fig):
    """The four corners of a text's ink box, in display pixels."""
    scale = fig.dpi / 72.0
    w, h = _ink(t.get_text(), t.get_fontsize(), t.get_fontweight(),
                t.get_fontstyle())
    w, h = w * scale, h * scale
    if w <= 0 or h <= 0:
        return None
    x, y = t.get_transform().transform(t.get_position())
    dx = {"left": 0.0, "center": -w / 2, "right": -w}.get(t.get_ha(), -w / 2)
    dy = {"bottom": 0.0, "baseline": 0.0, "top": -h}.get(t.get_va(), -h / 2)
    ang = math.radians(t.get_rotation())
    ca, sa = math.cos(ang), math.sin(ang)
    return [(x + ox * ca - oy * sa, y + ox * sa + oy * ca)
            for ox, oy in ((dx, dy), (dx + w, dy), (dx + w, dy + h), (dx, dy + h))]


def _bite(a, b):
    """How deeply two oriented rectangles interpenetrate, in pixels."""
    least = float("inf")
    for poly in (a, b):
        for i in range(len(poly)):
            (x0, y0), (x1, y1) = poly[i], poly[(i + 1) % len(poly)]
            nx, ny = -(y1 - y0), x1 - x0
            norm = math.hypot(nx, ny)
            if not norm:
                continue
            nx, ny = nx / norm, ny / norm
            pa = [p[0] * nx + p[1] * ny for p in a]
            pb = [p[0] * nx + p[1] * ny for p in b]
            gap = min(max(pa), max(pb)) - max(min(pa), min(pb))
            if gap <= 0:
                return 0.0
            least = min(least, gap)
    return least


def _all_text(fig):
    out = []
    for ax in fig.axes:
        out.extend(ax.texts)
        if ax.get_title():
            out.append(ax.title)
        legend = ax.get_legend()
        if legend is not None:
            out.extend(legend.get_texts())
            if legend.get_title() is not None:
                out.append(legend.get_title())
    out.extend(fig.texts)
    for art in getattr(fig, "_phytreon_extra_artists", []) or []:
        out.extend(getattr(art, "get_texts", lambda: [])())
        title = getattr(art, "get_title", lambda: None)()
        if title is not None:
            out.append(title)
    return [t for t in out if t.get_visible() and t.get_text().strip()]


def collisions(fig, bite=1.0):
    """Pairs of labels whose glyphs overlap by more than a hairline."""
    fig.canvas.draw()
    boxes = [(t, _corners(t, fig)) for t in _all_text(fig)]
    boxes = [(t, c) for t, c in boxes if c]
    hits = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            if _bite(boxes[i][1], boxes[j][1]) > bite:
                hits.append((boxes[i][0].get_text(), boxes[j][0].get_text()))
    return hits


@pytest.fixture(scope="module")
def big():
    import os
    import pandas as pd
    here = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "examples", "data")
    aln = os.path.join(here, "big16S_aligned.fasta")
    csv = os.path.join(here, "big16S_metadata.csv")
    if not (os.path.exists(aln) and os.path.exists(csv)):
        pytest.skip("example data not fetched")
    tree = pt.build_tree(aln, aligner="none", trim_kw=dict(max_gap=0.5),
                         method="nj", dist_model="k2p", root="midpoint")
    tree.join_data(pd.read_csv(csv), on="name")
    return tree


def _clean(fig):
    hits = collisions(fig)
    plt.close(fig)
    assert hits == [], "labels collide: %s" % hits[:5]


def test_a_rectangular_tree_of_a_hundred_taxa_prints_every_name(big):
    _clean(pt.TreeFigure(big).tip_points(color="phylum").tip_labels()
           .support_labels().scale_bar().titled("106 taxa").draw())


def test_a_circular_tree_of_a_hundred_taxa_prints_every_name(big):
    # the default canvas has to grow with the tip count for this to hold: the
    # names sit around a circumference, and a flat 8x8 ran them together past
    # about thirty tips
    _clean(pt.TreeFigure(big, layout="circular").tip_points(color="phylum")
           .tip_labels().titled("circular").draw())


def test_a_circular_tree_with_a_ring_still_prints_every_name(big):
    import pandas as pd
    import os
    here = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "examples", "data")
    meta = pd.read_csv(os.path.join(here, "big16S_metadata.csv"))
    _clean(pt.TreeFigure(big, layout="circular").tip_points(color="phylum")
           .tip_labels().ring(meta, columns=["domain"]).titled("ring").draw())


def test_thinning_the_labels_shrinks_the_canvas_instead_of_wasting_it():
    # asking for fewer names means fewer to seat, so the figure should not be
    # sized for names that are never drawn
    import os
    import pandas as pd
    here = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "examples", "data")
    aln = os.path.join(here, "big16S_aligned.fasta")
    if not os.path.exists(aln):
        pytest.skip("example data not fetched")
    tree = pt.build_tree(aln, aligner="none", trim_kw=dict(max_gap=0.5),
                         method="nj", dist_model="k2p", root="midpoint")
    tree.join_data(pd.read_csv(os.path.join(here, "big16S_metadata.csv")), on="name")
    full = pt.TreeFigure(tree, layout="circular").tip_labels().draw()
    thin = pt.TreeFigure(tree, layout="circular").tip_labels(max_labels=20).draw()
    assert full.get_size_inches()[0] > thin.get_size_inches()[0]
    plt.close(full)
    plt.close(thin)


def test_a_tanglegram_of_a_hundred_taxa_prints_every_name(big):
    _clean(pt.TangleFigure(big, big, titles=("a", "b")).untangle()
           .titled("106 taxa").draw())


def test_an_unrooted_tree_says_so_rather_than_piling_names_up(big):
    # it spaces tips by branch geometry, so two taxa a gene cannot separate
    # land on one spot and so do their names -- at any figure size
    with pytest.warns(UserWarning, match="max_labels"):
        fig = (pt.TreeFigure(big, layout="unrooted").tip_points(color="phylum")
               .tip_labels().draw())
    plt.close(fig)
    fig = (pt.TreeFigure(big, layout="unrooted").tip_points(color="phylum")
           .tip_labels(max_labels=20).draw())
    plt.close(fig)


def text_on_fill(fig):
    """Labels sitting on top of filled shapes.

    A track's name printed over the track's own data is unreadable *and* hides
    what it names, and no amount of text-versus-text checking sees it.
    """
    from matplotlib.patches import Patch
    from matplotlib.path import Path as MPath
    fig.canvas.draw()
    bad = []
    for ax in fig.axes:
        paths = []
        for p in ax.patches:
            if isinstance(p, Patch) and p.get_fill():
                try:
                    paths.append(p.get_path().transformed(p.get_transform()))
                except Exception:                      # pragma: no cover
                    pass
        for t in ax.texts:
            if not t.get_visible() or not t.get_text().strip():
                continue
            pts = _corners(t, fig)
            if not pts:
                continue
            box = MPath(pts + [pts[0]], closed=True)
            if any(p.intersects_path(box, filled=True) for p in paths):
                bad.append(t.get_text())
    return bad


def test_a_ring_name_is_not_printed_over_the_ring(big):
    # the names go in the fan opening, along the spoke: across the spoke a name
    # takes up its own width in angle, which is more than the opening leaves
    import os
    import pandas as pd
    here = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "examples", "data")
    meta = pd.read_csv(os.path.join(here, "big16S_metadata.csv"))
    fig = (pt.TreeFigure(big, layout="circular").tip_points(color="phylum")
           .tip_labels().ring(meta, columns=["domain"]).titled("ring").draw())
    on_data = text_on_fill(fig)
    plt.close(fig)
    assert "domain" not in on_data


def _luminance(rgb):
    def lin(c):
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (lin(v) for v in rgb[:3])
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def worst_contrast(fig):
    """The least readable text-on-fill pairing in the figure, as a ratio."""
    from matplotlib.colors import to_rgba
    from matplotlib.path import Path as MPath
    fig.canvas.draw()
    worst = 21.0
    for ax in fig.axes:
        shapes = []
        for p in ax.patches:
            if not (p.get_fill() and p.get_visible()):
                continue
            fc = to_rgba(p.get_facecolor())
            if fc[3] < 0.2:
                continue
            try:
                shapes.append((p.get_path().transformed(p.get_transform()), fc))
            except Exception:                          # pragma: no cover
                pass
        for t in ax.texts:
            pts = _corners(t, fig) if t.get_visible() else None
            if not pts:
                continue
            box = MPath(pts + [pts[0]], closed=True)
            tc = to_rgba(t.get_color())
            for path, fc in shapes:
                if path.intersects_path(box, filled=True):
                    lt, lf = _luminance(tc), _luminance(fc)
                    worst = min(worst, (max(lt, lf) + 0.05) / (min(lt, lf) + 0.05))
    return worst


def test_a_name_inside_a_coloured_block_still_reads():
    # it used to be white whatever the block, which is 2.2:1 against the paler
    # half of the palette -- below the 3:1 floor for large text
    tree = pt.Tree.from_newick("((P1:.1,P2:.1):.1,P3:.2);")
    arch = {"P1": [("wHTH", 60), ("ParB", 180)],
            "P2": [("DUF262", 210), ("ParBDB", 90)],
            "P3": [("TRD", 70), ("PUA", 95)]}
    fig = pt.TreeFigure(tree).tip_labels().domains(arch, labels=True).draw()
    worst = worst_contrast(fig)
    plt.close(fig)
    assert worst >= 4.5, "least readable block label is %.1f:1" % worst


def test_a_round_layout_is_not_drawn_as_an_oval(big):
    # scaling x and y differently turns a circular tree into an ellipse and
    # every ring sector into a different shape
    for layout in ("circular", "circular_slanted", "unrooted"):
        fig = pt.TreeFigure(big, layout=layout).tip_labels(max_labels=20).draw()
        for ax in fig.axes:
            if ax.get_aspect() not in ("equal", 1.0, 1):
                continue
            p0 = ax.transData.transform((0.0, 0.0))
            px = ax.transData.transform((1.0, 0.0))
            py = ax.transData.transform((0.0, 1.0))
            sx, sy = abs(px[0] - p0[0]), abs(py[1] - p0[1])
            assert abs(sx - sy) / max(sx, sy) < 0.01, layout
        plt.close(fig)


def test_no_default_figure_sets_text_too_small_to_print(big):
    # 5 pt is about the floor that survives a journal reducing a figure to one
    # column; anything under it is decoration, not information
    figs = [pt.TreeFigure(big).tip_points(color="phylum").tip_labels()
            .support_labels().scale_bar().draw(),
            pt.TreeFigure(big, layout="circular").tip_labels().draw()]
    for fig in figs:
        for ax in fig.axes:
            for t in list(ax.texts) + ([ax.title] if ax.get_title() else []):
                assert t.get_fontsize() >= 5.0, t.get_text()
        plt.close(fig)


def test_a_dendrogram_grows_sideways_not_downwards(big):
    # its leaves run along x, so it is the width that has to follow the tip
    # count -- it was being given the tall, narrow shape a rectangular tree
    # wants, which crushed 106 names into eight inches
    fig = pt.TreeFigure(big, layout="dendrogram").tip_labels().draw()
    w, h = fig.get_size_inches()
    assert w > h
    _clean(fig)


def test_network_edges_stay_thick_enough_to_print():
    # edge width carries the similarity, but scaled without a floor the weak
    # edges reach widths a press cannot hold -- so they fade instead
    net = pt.SequenceNetwork.from_pairs(
        [("s%d" % i, "s%d" % (i + 1), 0.05 + 0.15 * i) for i in range(5)])
    fig = net.draw()
    for ax in fig.axes:
        for ln in ax.lines:
            if ln.get_linestyle() not in ("None", "none", ""):
                assert ln.get_linewidth() >= 0.3, ln.get_linewidth()
    plt.close(fig)


def test_a_split_network_prints_every_name():
    import os
    here = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "examples", "data")
    aln = os.path.join(here, "tol_16S_aligned.fasta")
    if not os.path.exists(aln):
        pytest.skip("example data not fetched")
    seqs = pt.read_fasta(aln)
    names, mat = pt.infer.distance_matrix_model(
        pt.Alignment([n for n, _ in seqs], [s for _, s in seqs]), "k2p")
    _clean(pt.neighbor_net(names, mat, label_size=7)
           .titled("Neighbor-Net").draw(figsize=(11, 8)))
