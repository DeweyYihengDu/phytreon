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

### Protein sequences

Amino acid alignments use the same `build_tree(..., method="ml")` call, with
one of the empirical protein models (`"JTT"` / `"WAG"` / `"LG"`) passed as
`ml_model`. There is no automatic alphabet detection: `ml_model`'s default
stays `"HKY85"`, so a protein alignment with no explicit `ml_model` (or a
nucleotide model passed to protein data) raises a `ValueError` instead of
silently producing a meaningless result.

```python
ml = pt.build_tree("proteins.fasta", method="ml", ml_model="LG", bootstrap=100)
print(ml.data["logL"], ml.data["AIC"])

# model_finder ranks JTT/WAG/LG (instead of the nucleotide set) once it
# detects the alignment is protein:
for row in pt.model_finder("proteins.fasta"):
    print(row)
```

Distance-based methods (`method="nj"`/`"upgma"`) work on protein data too;
`dist_model`'s default (`"jc69"`) falls back to raw p-distance on protein
data unless you opt in to the protein-specific correction, `"poisson"`:

```python
nj = pt.build_tree("proteins.fasta", method="nj", dist_model="poisson")
```

## Parsimony

```python
mp = pt.build_tree("seqs.fasta", method="parsimony")
print(mp.data["parsimony_score"], mp.data["ci"], mp.data["ri"])
```

## From a distance or character matrix

A precomputed distance matrix skips alignment entirely:

```python
import pandas as pd
df = pd.read_csv("distances.csv", index_col=0)
tree = pt.neighbor_joining(list(df.index), df.values.tolist())
```

A discrete character/trait matrix (e.g. a 0/1 gene presence/absence table)
goes through `read_character_matrix` and into parsimony:

```python
aln = pt.read_character_matrix("genes.csv", taxa_col="name")
tree = pt.parsimony_tree(aln, search=True)
```

Any small set of hashable states per column works; missing values (`NaN`,
or an explicit `missing=` sentinel) are encoded as ambiguous so they never
force a false character change.

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
