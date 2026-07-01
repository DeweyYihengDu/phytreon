# phytreon

**Phylogenetic trees and publication figures in Python** — a layered library
with a fluent figure-builder API and both static (matplotlib) and interactive
(plotly) backends.

```python
import phytreon as pt

tr = pt.datasets.primates()
(pt.TreeFigure(tr).tip_labels().support_labels()).save("tree.pdf")
```

## Install

```bash
pip install -e .                 # core
pip install -e .[interactive]    # + plotly backend
pip install -e .[dev]            # + pytest, plotly
```

## What's inside

- **core** — `Tree`/`Node` model, Newick/Nexus/PhyloXML I/O, metadata join
- **layout** — rectangular, slanted, dendrogram, circular, fan, radial,
  inward-circular, unrooted (equal-angle / equal-daylight)
- **infer** — NJ/UPGMA (model-corrected distances), native ML (JC69/K80/HKY85/
  GTR +Γ, NNI, AIC/`model_finder`), parsimony, bootstrap, built-in MSA, trimming
- **comparative** — ancestral states (parsimony / Mk-ML ER·SYM·ARD / Brownian),
  stochastic mapping (`stochastic_map`)
- **plot** — the `TreeFigure` builder: chain `.tip_labels()`, `.tip_points()`,
  `.heatmap()`, `.ring()` … onto a tree, then `.save()`

See the [tutorials](tutorials/build_tree.md) and the
[API reference](api.md). phytreon's core is validated in pure Python
(`validation/`): the likelihood engine matches an independent naive
implementation to machine precision, and NJ recovers a tree exactly from its
own additive distances.
