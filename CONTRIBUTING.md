# Contributing to phytreon

Thanks for your interest. phytreon is a pure-Python phylogenetics +
visualization library; contributions that keep it dependency-light and well
tested are very welcome.

## Setup

```bash
pip install -e .[dev]      # installs pytest
pytest -q                  # run the test suite
```

## Architecture (where things live)

- `phytreon/core/`        — `Tree`/`Node` model and I/O
- `phytreon/layout/`      — topology → display coordinates (subclass `Layout`)
- `phytreon/infer/`       — alignment, distances, ML, parsimony, bootstrap
- `phytreon/comparative/` — ancestral states, stochastic mapping
- `phytreon/plot/`        — TreeFigure builder, elements, matplotlib/plotly backends
- `phytreon/scene.py`     — backend-agnostic drawing primitives

The key invariant: **layout computes coordinates and emits scene primitives;
backends only render them.** A new element subclasses `plot.figure._Element`
and appends primitives; a new layout subclasses `layout.base.Layout`.

## Guidelines

- Add a test in `tests/` for any new feature or bug fix.
- Match the surrounding style; keep heavy deps (scipy/biopython/plotly)
  import-local where reasonable.
- Run `pytest -q` before opening a PR; CI must stay green.
- For phylogenetic methods, prefer adding a correctness check (e.g. recovering
  a known topology, or agreement with an independent reference as in
  `validation/`).
