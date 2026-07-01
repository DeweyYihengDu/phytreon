# Tutorial: comparative methods

## Ancestral state reconstruction

```python
import phytreon as pt
tr = pt.datasets.primates()
trait = {"Human": "urban", "Chimp": "forest", ...}

pt.ace_ml(tr, trait, model="ARD")          # ER | SYM | ARD Mk models
print(tr.root.data["_ace_model"]["rates"])

# visualise posterior probabilities as pies at internal nodes:
(pt.TreeFigure(tr).node_pies().tip_labels()).save("ace.pdf")
```

Continuous traits use Brownian-motion (independent contrasts):

```python
vals = pt.ace_continuous(tr, {"Human": 62.0, ...})
```

## Stochastic character mapping

```python
pt.stochastic_map(tr, trait, n=200, model="ARD")    # simulate histories
(pt.TreeFigure(tr).painted_branches()               # branches painted by state
    .tip_labels()).save("painted.pdf")
```

Each branch is split into segments proportional to the time spent in each
state. Node posteriors are stored in `node.data["ace_probs"]`.

## Time-scaled trees

```python
(pt.TreeFigure(dated_tree)                  # branch lengths = time
    .time_axis(geo=True, unit="Mya")        # geological period bands
    .tip_labels()).save("timetree.pdf")
```

## Reshaping

```python
pt.rotate(tr, node); pt.flip(tr, a, b)      # move branches
pt.collapse_low_support(tr, 70)             # weak edges -> polytomies
clusters = pt.cut_tree(tr, k=4)             # {tip: cluster_id}
pt.group_clade(tr, {node: "lineage A"})     # colour branches by lineage
```
