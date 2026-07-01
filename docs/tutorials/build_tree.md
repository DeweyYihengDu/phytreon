# Tutorial: from sequences to a tree

The whole pipeline — **align → trim → infer → bootstrap** — is one call, every
stage configurable.

```python
import phytreon as pt

# 1. Quick neighbour-joining tree (JC69-corrected distances, midpoint root)
tree = pt.build_tree(
    "examples/data/tol_16S.fasta",   # raw, unaligned sequences
    aligner="builtin",               # pure-Python progressive MSA
    trim_kw=dict(max_gap=0.5),       # the "cut" step
    method="nj",
    root="midpoint",
    bootstrap=200,                   # bipartition support
)
print(tree.write())                  # Newick
```

## Maximum likelihood

```python
ml = pt.build_tree("seqs.fasta", method="ml", ml_model="HKY85", ml_gamma=4)
print(ml.data["logL"], ml.data["AIC"], ml.data["gamma_shape"])

# pick the best model by AIC:
for row in pt.model_finder("seqs.fasta"):
    print(row)
```

## Parsimony

```python
mp = pt.build_tree("seqs.fasta", method="parsimony")
print(mp.data["parsimony_score"], mp.data["ci"], mp.data["ri"])
```

## Scaling up

The built-in aligner and ML are pure Python and target small/medium data.
For hundreds of taxa, plug in external engines:

```python
tree = pt.build_tree("big.fasta", aligner="mafft",
                     method="ml", ml_engine="iqtree")
```

## Comparing to a reference tree

```python
rf = pt.robinson_foulds(tree_a, tree_b)        # symmetric-difference distance
```
