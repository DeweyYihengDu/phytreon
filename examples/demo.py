"""End-to-end demo: every element of phytreon, both backends.

Run from the repo root:  python examples/demo.py
Outputs land in examples/out/.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import phytreon as pt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)


def demo_rectangular():
    tr = pt.datasets.primates()
    meta = pt.datasets.primates_metadata()
    tr.join_data(meta.reset_index(), on="name")

    apes = tr.get_mrca(["Human", "Gibbon"])
    p = (
        pt.TreeFigure(tr)
        .highlight(node=apes, fill="#cfe8f3")
        .branches(size=1.4)
        .tip_points(color="habitat", size=8)
        .tip_labels()
        .support_labels()
        .clade_label("Apes", node=apes)
        .titled("Primates - rectangular")
    )
    p.save(os.path.join(OUT, "rect.pdf"))
    p.save(os.path.join(OUT, "rect.png"))
    p.save(os.path.join(OUT, "rect.html"))
    print("[ok] rectangular -> rect.{pdf,png,html}")


def demo_circular():
    tr = pt.datasets.primates()
    meta = pt.datasets.primates_metadata()
    tr.join_data(meta.reset_index(), on="name")
    p = (
        pt.TreeFigure(tr, layout="circular", extent=320)
        .tip_points(color="habitat", size=10)
        .tip_labels(offset=0.01)
        .titled("Primates - circular")
    )
    p.save(os.path.join(OUT, "circular.pdf"))
    p.save(os.path.join(OUT, "circular.html"))
    print("[ok] circular -> circular.{pdf,html}")


def demo_heatmap():
    tr = pt.datasets.primates()
    meta = pt.datasets.primates_metadata()
    # numeric matrix for the heatmap
    mat = meta[["body_mass_kg"]].copy()
    mat["log_mass"] = (mat["body_mass_kg"] ** 0.5)
    p = (
        pt.TreeFigure(tr)
        .tip_labels(offset=0.01)
        .heatmap(mat, width=0.5, offset=0.12)
        .titled("Primates + trait heatmap")
    )
    p.save(os.path.join(OUT, "heatmap.pdf"))
    print("[ok] heatmap -> heatmap.pdf")


def demo_inference():
    names = ["A", "B", "C", "D", "E"]
    # symmetric distance matrix
    d = [
        [0, 5, 9, 9, 8],
        [5, 0, 10, 10, 9],
        [9, 10, 0, 8, 7],
        [9, 10, 8, 0, 3],
        [8, 9, 7, 3, 0],
    ]
    tr = pt.neighbor_joining(names, d)
    (pt.TreeFigure(tr.ladderize()).tip_labels(offset=0.05).tip_points()) \
        .save(os.path.join(OUT, "nj.pdf"))
    print("[ok] neighbor-joining -> nj.pdf  (newick: %s)" % tr.write())


def demo_ancestral():
    tr = pt.datasets.primates()
    trait = {
        "Human": "urban", "Chimp": "forest", "Gorilla": "forest",
        "Orangutan": "forest", "Gibbon": "forest",
        "Macaque": "savanna", "Baboon": "savanna",
    }
    pt.ace_ml(tr, trait)
    model = tr.root.data["_ace_model"]
    print("[ok] ACE (Mk-ML) rate=%.3f states=%s" % (model["rate"], model["states"]))
    # colour tips by observed state, internal nodes by reconstructed state
    for tip in tr.leaves():
        tip.data["state"] = trait[tip.name]
    p = (
        pt.TreeFigure(tr)
        .tip_labels(offset=0.01)
        .node_points(color="ace_state", size=11)
        .tip_points(color="state", size=9)
        .titled("Ancestral habitat (Mk-ML marginal)")
    )
    p.save(os.path.join(OUT, "ancestral.pdf"))
    print("[ok] ancestral reconstruction -> ancestral.pdf")


if __name__ == "__main__":
    demo_rectangular()
    demo_circular()
    demo_heatmap()
    demo_inference()
    demo_ancestral()
    print("\nAll demos written to:", OUT)
