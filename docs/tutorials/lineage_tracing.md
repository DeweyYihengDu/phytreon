# Tutorial: single-cell CRISPR lineage tracing

Cas9-based lineage tracing (GESTALT, Cassiopeia-style "allele table"
experiments, etc.) engineers cells to accumulate small indel "scars" at a
handful of genomic target sites as they divide. Cells sharing a scar
inherited it from a common ancestor -- the scar pattern is itself a
phylogenetic character matrix, just one with an **irreversible** evolutionary
model: a site, once cut, can never revert to unedited, and can't
spontaneously become a *different* edit either.

phytreon reconstructs the cell-division tree from raw allele-table data using
this irreversible ("Camin-Sokal") parsimony model -- a purely additive
capability alongside the existing reversible Fitch parsimony used for
discrete character/trait matrices.

## From a raw allele table to a tree

```python
import phytreon as pt

aln = pt.read_allele_table("alleletable.txt")   # cellBC/intBC/r1/r2/r3 columns
tree = pt.lineage_tree(aln, search=True)        # NNI hill-climbing

print(tree.data["camin_sokal_score"],       # reconstruction's total cost
      tree.data["min_possible_score"],      # topology-independent lower bound
      tree.data["excess_origins"])          # homoplasy the tree couldn't avoid
```

`read_allele_table` handles the two real correctness pitfalls in CRISPR
allele tables directly: allele dropout (a cell missing a row for some
integration site -- encoded as missing, not silently omitted) and sites that
happen to be edited in nearly every profiled cell (the true ancestral/uncut
state is still guaranteed to recode as character `"0"`, even though it's
never actually observed at that site).

## Beyond CRISPR: any irreversible mutation signal

A somatic mutation follows the exact same logic as a CRISPR scar: a mutated
gene doesn't spontaneously revert to wild-type, so cells sharing a mutation
are evidence of common descent. `read_mutation_matrix` covers this general
case -- a single gene, or a handful of genes -- without the CRISPR-specific
allele-table format:

```python
import pandas as pd

genotypes = pd.DataFrame({"TP53": ["R175H", "R175H", "R175H", "WT", "WT", "WT"]},
                         index=["A1", "A2", "A3", "B1", "B2", "B3"])
aln = pt.read_mutation_matrix(genotypes)        # wild_type="WT" by default
tree = pt.lineage_tree(aln, search=True)
```

Same guarantees as `read_allele_table`: the wild-type state always lands on
code `"0"`, even for a gene mutated in every profiled cell.

## Which mutation arose on which branch

A tree and a total cost tell you cells are related, but not *how the
lineage unfolded*. `reconstruct_ancestral_mutations` traces that back under
the same Camin-Sokal model, writing `node.data["mutations_acquired"]` --
the site/gene names that transitioned from ancestral to derived on that
node's incoming branch -- for every node in the tree:

```python
pt.reconstruct_ancestral_mutations(tree, aln, site_names=list(genotypes.columns))
for node in tree.traverse():
    if node.data["mutations_acquired"]:
        print(node.name or "(clade)", "->", node.data["mutations_acquired"])
```

`site_names` labels the alignment's columns in order (defaults to
`"site0"`, `"site1"`, ... if omitted) -- pass the same gene/site list you
used to build `aln`. Ties among equally-optimal reconstructions are broken
toward matching the parent's state, so an already-derived site doesn't get
reported as spuriously re-arising partway down a clade that simply
inherited it.

## Through the one-call pipeline

`build_tree(..., method="parsimony")` also accepts lineage data -- pass
`parsimony_model="camin_sokal"` for the irreversible model instead of the
default reversible Fitch parsimony:

```python
tree = pt.build_tree(aln, aligner="none", method="parsimony",
                     parsimony_model="camin_sokal", bootstrap=100)
```

## Validating against an independent reconstruction

If you have a reference tree for the same cells from another tool (or a
published reconstruction), compare topologies with the existing
`robinson_foulds()`:

```python
ref = pt.Tree.from_newick(open("reference.nwk").read())
rf = pt.robinson_foulds(tree, ref, normalized=True)   # 0 = identical
```

`robinson_foulds()` requires both trees to cover the same leaf set; use
`pt.prune_to_taxa(tree, shared_cells)` first if one tree has extra cells the
other doesn't.

## Scaling

The irreversible cost model is intrinsically pricier per NNI move than
reversible Fitch parsimony (a small cost-matrix computation per candidate
state rather than an O(1) bitmask op). On a real 226-cell, 30-site dataset
(`examples/lineage_demo.py`), scoring the NJ starting tree alone took ~15s;
the full NNI search took ~2.5 minutes. For much larger datasets, start with
`search=False` (branch/model fit on the NJ tree only) before committing to a
full topology search.
