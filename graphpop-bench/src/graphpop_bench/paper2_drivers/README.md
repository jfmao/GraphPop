# `paper2_drivers/` — Per-figure benchmark drivers for Paper 2

Each module orchestrates one figure's data pipeline + figure-gen
for **Paper 2** (kinship + ARG). The driver layer pulls together:

- competitor wrapper outputs (`graphpop_bench.competitors`),
- GraphPop's actual procedure output via Maven-driven Java test
  dumps (e.g. `-Dgraphpop.bench.dump.dir=<path>` triggers the
  hook in `BranchGrmProcedureTest`),
- frozen reference fixtures (`graphpop-procedures/src/test/
  resources/egrm_expected_*.json`),

joins them by sample pair, and emits a tidy panel CSV. A
sibling figure-gen module reads the CSV and writes a vector PDF
to `paper/paper2_kinship_arg/figures/` (a local-only path —
`paper/` is git-ignored per the existing private paper workflow).

## Shipped

### Fig 1d / Fig 1e — branch-GRM rel-err validation

- Driver:  `fig1de_panels.py` →
  `python -m graphpop_bench.paper2_drivers.fig1de_panels --help`
- Figures: `fig1de_figures.py` →
  `python -m graphpop_bench.paper2_drivers.fig1de_figures --help`
- Tests:   `graphpop-bench/tests/test_paper2_fig1de.py` (15
  unit tests on the join + rel-err logic; no Java / mvn / egrm
  dependency)

Per-paper plan + reproduction notes are in
`paper/paper2_kinship_arg/benchmarks/{PLAN_fig1de.md, README.md}`
(local-only).

#### Reproducing

```bash
# 1. Run the driver — invokes mvn (Java test dumps GraphPop's
#    branch_grm output as H1-schema TSV) + the H4 egrm wrapper,
#    then joins with the cached egrm_expected_*.json references.
python -m graphpop_bench.paper2_drivers.fig1de_panels \
    --java-dump-dir /tmp/fig1de/graphpop \
    --egrm-out-dir /tmp/fig1de/egrm_h4 \
    --output-dir /tmp/fig1de/panels

# 2. Render the figures from the panel CSVs.
python -m graphpop_bench.paper2_drivers.fig1de_figures \
    --fig1d-csv /tmp/fig1de/panels/fig1d_panel_data.csv \
    --fig1e-csv /tmp/fig1de/panels/fig1e_panel_data.csv \
    --fig1d-pdf /tmp/fig1de/fig1d.pdf \
    --fig1e-pdf /tmp/fig1de/fig1e.pdf
```

The driver gates on per-entry rel-err < 1e-6 vs the egrm
reference; non-zero exit on any panel exceeding the gate.

#### Last-run numbers (2026-05-11)

|        | Fig 1d (unconditional) | Fig 1e (pathway_half) |
|--------|------------------------|-----------------------|
| n pairs | 210                   | 210                   |
| max rel-err | 3.72e-10          | 3.70e-10              |
| max abs-diff | 4.33e-10         | 4.28e-10              |
| status | ✓ < 1e-6 gate          | ✓ < 1e-6 gate         |

Both panels clear the gate by ~4 orders of magnitude — the
GraphPop branch_grm procedure and `egrm.varGRM_C` implement
the same math to floating-point precision.

## Pending

- Fig 2 (posterior branch GRM)
- Fig 3 (annotation-conditional GRMs)
- Fig 4 (biobank-scale validation)
- Fig 5 (graph-native query plane)

Per the roadmap, drivers ship one figure at a time with an
explicit pause-and-review boundary between each.

## Optional install extras

The figure-gen modules require `matplotlib`. Install with the
`figures` extra:

```bash
pip install -e ".[figures]"
```
