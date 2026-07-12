"""Download a real single-cell CRISPR lineage-tracing demo dataset.

One real mouse tumor clonal sample (226 cells, 10 CRISPR integration
barcodes) from the Cassiopeia package (YosefLab/Cassiopeia, MIT license),
plus Cassiopeia's own published tree reconstruction for the same sample
(used as an external validation reference -- see examples/lineage_demo.py).
The published tree ships in a custom ete3 annotated-Newick format
(``)N|N|N|...`` per-node feature blocks instead of ``:branchlength``); this
strips that down to a plain bifurcating Newick tree. Writes:

    examples/data/lineage_alleletable.txt      raw allele table (unmodified)
    examples/data/lineage_reference_tree.nwk   published tree, plain Newick

Run:  python examples/data/fetch_lineage_data.py
"""
import os
import re
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "https://raw.githubusercontent.com/YosefLab/Cassiopeia/master/notebooks/data"
SAMPLE = "3432_NT_T1"


def _strip_ete3_annotations(newick: str) -> str:
    """Cassiopeia's tree.write() appends a ``|``-joined feature block after
    every node's closing paren (including single-leaf wraps like
    ``(CELLBARCODE)4|5|3|...``) instead of a standard ``:branchlength``.
    Strip those blocks and collapse the resulting redundant single-leaf
    parens down to a plain leaf name, leaving a bare topology-only Newick
    string (no branch lengths -- the published tree's are Cassiopeia-internal
    edit-distance units, not directly comparable to ours anyway; this demo
    only needs the topology for Robinson-Foulds comparison)."""
    stripped = re.sub(r"\)[0-9\-]+(?:\|[0-9\-]+)*", ")", newick.strip())
    prev = None
    while prev != stripped:
        prev = stripped
        stripped = re.sub(r"\(([A-Za-z0-9_.\-]+)\)", r"\1", stripped)
    return stripped


def main():
    allele_path = os.path.join(HERE, "lineage_alleletable.txt")
    tree_path = os.path.join(HERE, "lineage_reference_tree.nwk")

    print(f"downloading {SAMPLE} allele table ...")
    urllib.request.urlretrieve(f"{BASE}/{SAMPLE}_alleletable.txt", allele_path)
    n_rows = sum(1 for _ in open(allele_path, encoding="utf-8")) - 1
    print(f"  wrote {allele_path} ({n_rows} rows)")

    print(f"downloading {SAMPLE} published reference tree ...")
    with urllib.request.urlopen(f"{BASE}/{SAMPLE}_tree.processed.tree") as r:
        raw_newick = r.read().decode("utf-8")
    clean_newick = _strip_ete3_annotations(raw_newick)
    with open(tree_path, "w", encoding="utf-8") as f:
        f.write(clean_newick + "\n")

    from phytreon.core.tree import Tree
    t = Tree.from_newick(clean_newick)
    print(f"  wrote {tree_path} ({t.n_leaves} leaves, topology only)")


if __name__ == "__main__":
    main()
