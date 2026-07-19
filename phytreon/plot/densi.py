"""DensiTree: a whole set of trees drawn on top of each other.

The standard picture for topological *uncertainty* -- a posterior sample from
BEAST/MrBayes, a bootstrap set, or a collection of gene trees -- where a
single summary tree would hide the disagreement::

    trees = [pt.Tree.read(p) for p in paths]
    pt.DensiTreeFigure(trees).save("cloud.pdf")

Every tree is drawn translucent, so shared structure accumulates into dark,
confident edges while conflicting placements stay faint and fan out. It is the
same idea as ``ggdensitree``.

The trees only line up if their tips are in a common order, so each tree is
first rotated to match a reference (via :func:`phytreon.treeops.untangle`) --
rotation reorders children without touching topology or branch lengths, so
aligning the cloud never changes what any tree says.
"""
from __future__ import annotations

from typing import List, Optional, Sequence

from ..core.tree import Tree
from ..layout import get_layout
from ..scene import Path, Scene
from ..treeops import untangle as _untangle
from .figure import RenderContext, TreeFigure, _Renderable


class DensiTreeFigure(_Renderable):
    """Overlay a set of trees sharing the same taxa.

    ``consensus`` picks the tree whose tip order the others are rotated to
    match; by default the first one. ``align=False`` skips the rotation pass
    and draws the trees exactly as given (faster, but usually a tangle).

    Tips must be shared across the set; a tree missing some of the reference's
    taxa is still drawn, it simply contributes no edges for them.
    """

    def __init__(self, trees: Sequence[Tree], *, layout: str = "rectangular",
                 color: str = "#37618e", alpha: Optional[float] = None,
                 width: float = 0.6, align: bool = True,
                 consensus: int = 0, tip_labels: bool = True,
                 **layout_kwargs):
        self.trees = list(trees)
        if not self.trees:
            raise ValueError("DensiTreeFigure needs at least one tree")
        self.layout_name = layout
        self.layout_kwargs = layout_kwargs
        self.color = color
        self.width = width
        self.align = align
        self.consensus = consensus
        self.tip_labels = tip_labels
        self.title: Optional[str] = None
        # enough opacity that a lone tree still reads, faint enough that a big
        # sample does not saturate into a solid block
        self.alpha = alpha if alpha is not None \
            else max(0.04, min(0.8, 6.0 / len(self.trees)))

    # -- composition -----------------------------------------------------
    def titled(self, title: str) -> "DensiTreeFigure":
        self.title = title
        return self

    @property
    def reference(self) -> Tree:
        """The tree whose tip order the cloud is aligned to."""
        return self.trees[self.consensus]

    # -- building --------------------------------------------------------
    def _build(self) -> RenderContext:
        ref = self.reference
        if self.align:
            for tree in self.trees:
                if tree is not ref:
                    _untangle(ref, tree, fix="left")

        # the reference carries the labels and sets the frame; drawing it
        # through a normal TreeFigure keeps every element behaving as usual
        base = TreeFigure(ref, layout=self.layout_name, skeleton=False,
                          **self.layout_kwargs)
        if self.tip_labels:
            # every tree in the cloud ends its tips at a slightly different
            # depth, so labels are aligned in a column rather than chasing any
            # one tree's tip positions
            base.tip_labels(align=True)
        ctx = base._build()

        rows = {name: i for i, name in enumerate(ref.leaf_names())}
        for tree in self.trees:
            self._add_tree(ctx.scene, tree, rows, ctx.layout.max_x)
        return ctx

    def _add_tree(self, scene: Scene, tree: Tree, rows: dict,
                  ref_max_x: float) -> None:
        """Lay a single tree out on the reference's rows and stroke its edges."""
        layout = get_layout(self.layout_name, **self.layout_kwargs)
        layout.apply(tree)
        # rescale this tree's depth onto the reference's, so trees of different
        # total length still overlay. On a rectangular layout depth runs along
        # x alone; on a polar one it is the radius, and since a point is
        # (r cos a, r sin a), scaling the radius means scaling *both*
        # coordinates -- scaling x only would smear the tree into an ellipse.
        scale = (ref_max_x / layout.max_x) if layout.max_x else 1.0
        polar = getattr(layout, "is_polar", False)

        # nudge onto the reference's rows. For same-taxon trees both layouts
        # number the rows 0..n-1, so this is normally 0; it only bites when a
        # tree is missing taxa.
        shift: List[float] = []
        for leaf in tree.leaves():
            if leaf.name in rows:
                shift.append(rows[leaf.name] - leaf.data["_row"])
        dy = 0.0 if polar else (sum(shift) / len(shift) if shift else 0.0)

        for node in tree.traverse():
            if node.is_root:
                continue
            for pts in (layout.branch_path(node), layout.child_connector(node)):
                if not pts:
                    continue
                if polar:
                    moved = [(x * scale, y * scale) for x, y in pts]
                else:
                    moved = [(x * scale, y + dy) for x, y in pts]
                scene.add(Path(moved, color=self.color, width=self.width,
                               opacity=self.alpha, zorder=1))

    def _default_figsize(self, ctx: RenderContext = None):
        n = self.reference.n_leaves
        return (8.0, max(3.0, min(0.3 * n, 30.0)))
