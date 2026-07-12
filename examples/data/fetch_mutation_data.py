"""Download a real single-cell somatic-mutation demo dataset.

Hou, Y. et al. 2012, "Single-cell exome sequencing and monoclonal evolution
of a JAK2-negative myeloproliferative neoplasm", Cell 148:873-885 -- a real
clonal-evolution cancer dataset (18 genes, 58 cells), redistributed as
``dataHou18.csv``/``dataHou18.geneNames`` by the SCITE package
(cbg-ethz/SCITE, GPL3 license), which documents the encoding in its
README: ``0`` wild-type, ``1`` heterozygous mutation, ``2`` homozygous
mutation, ``3`` missing. Genes are rows / cells are columns in the raw
file, with no cell names.

phytreon's Camin-Sokal model treats every derived state as mutually
exclusive (it doesn't model a heterozygous -> homozygous progression at the
same site), so this deliberately collapses ``1``/``2`` into one "mutated"
label rather than silently mismodeling zygosity as two unrelated mutation
events -- a documented simplification, not a hidden one. Writes:

    examples/data/mutation_genotypes.csv   cells x genes, wild_type="WT"

Run:  python examples/data/fetch_mutation_data.py
"""
import os
import urllib.request

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "https://raw.githubusercontent.com/cbg-ethz/SCITE/master"

_RECODE = {"0": "WT", "1": "mutated", "2": "mutated", "3": None}


def main():
    genes_path = os.path.join(HERE, "_dataHou18.geneNames")
    matrix_path = os.path.join(HERE, "_dataHou18.csv")
    out_path = os.path.join(HERE, "mutation_genotypes.csv")

    print("downloading Hou et al. 2012 mutation matrix (via SCITE) ...")
    urllib.request.urlretrieve(f"{BASE}/dataHou18.csv", matrix_path)
    urllib.request.urlretrieve(f"{BASE}/dataHou18.geneNames", genes_path)

    genes = [line.strip() for line in open(genes_path, encoding="utf-8") if line.strip()]
    rows = [line.split() for line in open(matrix_path, encoding="utf-8") if line.strip()]
    assert len(rows) == len(genes), f"{len(rows)} rows but {len(genes)} gene names"
    n_cells = len(rows[0])

    # genes-as-rows/cells-as-columns -> cells-as-rows/genes-as-columns
    cells = [f"cell{i + 1}" for i in range(n_cells)]
    data = {gene: [_RECODE[rows[g][c]] for c in range(n_cells)]
           for g, gene in enumerate(genes)}
    df = pd.DataFrame(data, index=cells)
    df.to_csv(out_path)
    print(f"  wrote {out_path} ({df.shape[0]} cells x {df.shape[1]} genes)")

    os.remove(genes_path)
    os.remove(matrix_path)


if __name__ == "__main__":
    main()
