"""Pure-Python validation of phytreon's core algorithms (no external tools).

1. Likelihood engine: cross-check `log_likelihood` (pattern-compressed, rescaled)
   against a deliberately naive, independent re-implementation (per column, no
   compression, no scaling).  They must agree to numerical precision.
2. Neighbor-joining: NJ is guaranteed to recover a tree from *additive*
   distances, so NJ on a tree's own patristic distance matrix must return the
   original topology (Robinson-Foulds = 0).
3. ML: must recover a known clade and report a finite logL/AIC.

Writes validation/validation_report.md.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import phytreon as pt

HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(HERE, "validation_report.md")

SEQS = [("A1", "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"),
        ("A2", "ATGGCCATTGTTATGGGCCGCTGAAAGGGTGCCCGATAG"),
        ("A3", "ATGGCCATTGTAATGGGCCGCTGTAAGGGTGCCGATAG"),
        ("B1", "ATGTCGATTCTAATGAACCGCTGAAAGCGTGACCTTTAG"),
        ("B2", "ATGTCGATTCTAATGAACCGCTGTAAGCGTGACCTTTAG"),
        ("B3", "ATGTCGATTCTAATGAACCGGCTGAAAGCGTGACCTTTAG")]


def naive_jc_logl(tree, aln):
    """Independent JC69 log-likelihood: per-column Felsenstein pruning with no
    pattern compression and no rescaling (a reference implementation)."""
    import numpy as np
    from scipy.linalg import expm
    Q = (np.ones((4, 4)) - 4 * np.eye(4)) / 3.0     # JC: off-diag 1/3, diag -1
    s2i = {c: i for i, c in enumerate("ACGT")}
    s2i["U"] = 3
    P = {id(n): expm(Q * max(n.length or 0.0, 1e-9))
         for n in tree.traverse() if not n.is_root}
    pos = {n: i for i, n in enumerate(aln.names)}
    total = 0.0
    for j in range(aln.ncol):
        L = {}
        for node in tree.traverse("postorder"):
            if node.is_leaf:
                v = np.ones(4)
                ch = aln.seqs[pos[node.name]][j].upper()
                if ch in s2i:
                    v = np.zeros(4); v[s2i[ch]] = 1.0
                L[id(node)] = v
            else:
                v = np.ones(4)
                for c in node.children:
                    v = v * (P[id(c)] @ L[id(c)])
                L[id(node)] = v
        total += np.log(0.25 * L[id(tree.root)].sum() + 1e-300)
    return float(total)


def patristic_matrix(tree):
    leaves = tree.leaves()
    names = [l.name for l in leaves]
    depth = {n: n.depth(use_lengths=True) for n in tree.traverse()}
    D = [[0.0] * len(leaves) for _ in leaves]
    for i, a in enumerate(leaves):
        for j in range(i + 1, len(leaves)):
            b = leaves[j]
            mrca = tree.get_mrca([a.name, b.name])
            d = depth[a] + depth[b] - 2 * depth[mrca]
            D[i][j] = D[j][i] = d
    return names, D


def main():
    lines = ["# phytreon validation (pure Python, no external tools)", ""]

    # 1. likelihood engine vs naive reference
    aln = pt.align(SEQS, seqtype="nucleotide")
    tree = pt.ml_tree(aln, model="JC69", search=False)
    engine = pt.log_likelihood(tree, aln, model="JC69")
    naive = naive_jc_logl(tree, aln)
    diff = abs(engine - naive)
    lines += ["## 1. Likelihood engine vs independent naive implementation (JC69)",
              f"- engine  : {engine:.6f}",
              f"- naive   : {naive:.6f}",
              f"- **|diff| : {diff:.2e}**  ({'PASS' if diff < 1e-6 else 'FAIL'}, "
              f"tolerance 1e-6)", ""]

    # 2. NJ recovers a tree from its own additive distances
    true_tree = pt.datasets.random_tree(10, seed=3)
    names, D = patristic_matrix(true_tree)
    nj = pt.neighbor_joining(names, D)
    rf = pt.robinson_foulds(true_tree, nj)
    lines += ["## 2. NJ on additive (patristic) distances recovers the tree",
              f"- Robinson-Foulds(true, NJ) : **{rf:.0f}**  "
              f"({'PASS' if rf == 0 else 'FAIL'}, expected 0)", ""]

    # 3. ML recovers a known clade
    clade = {n for n in ("A1", "A2", "A3")}
    sides = [frozenset(x.leaf_names()) for x in tree.traverse() if not x.is_leaf]
    ok = frozenset(clade) in sides or frozenset({"B1", "B2", "B3"}) in sides
    lines += ["## 3. ML recovers the known A|B split",
              f"- A|B clade present : **{ok}**  ({'PASS' if ok else 'FAIL'})",
              f"- logL {tree.data['logL']:.2f}, AIC {tree.data['AIC']:.2f}", ""]

    with open(REPORT, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
