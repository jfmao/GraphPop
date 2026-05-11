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

## Phase 1 Step I — closed (5 of 5 figures shipped)

| Figure | Driver module(s)                              | Tests | Status                |
|--------|-----------------------------------------------|-------|-----------------------|
| 1d/1e  | `fig1de_panels`, `fig1de_figures`             | 15    | shipped (Phase 1)     |
| 2      | `fig2_panels`, `fig2_figures`                 | 22    | shipped (Phase 1)     |
| 3a/3b  | `fig3ab_panels`, `fig3_figures`               | (with 3c) | shipped (Phase 1) |
| 3c     | `fig3c_panels`, `fig3_figures`                | 11    | shipped (Phase 1)     |
| 4a/4b  | `fig4_panels`, `fig4_figures`                 | 12    | shipped (P1 + P2.1)   |
| 5a/5d  | `fig5_panels`, `fig5_figures`                 | 9     | shipped (Phase 1)     |

Total tests: **161 green + 2 skipped** (skips gated by absent
ldsc.py / king binaries).

## Phase 2

| Sub-iteration | Status | Blocking dependency |
|---|---|---|
| P2.1 — Fig 4 N-extension to 1000 haploid | shipped | (none — msprime only) |
| P2.2 — Fig 5c LOC comparison | unblocked | (Python pipeline on simulated cohort) |
| P2.3 — Fig 5b cryptic-pair recall | gated | Neo4j ingest + 1000G chr22 |
| P2.4 — Fig 3e pathway-conditional on 1000G | gated | ingest + ARG inference tool |
| P2.5 — Fig 2d/2e SINGER posterior on 1000G | gated | SINGER install + 1000G EUR/YRI |
| P2.6 — Fig 3f HGDP cross-pop | gated | HGDP data + ingest |

### P2.1 — Fig 4 N-extension (2026-05-11)

Extended the Fig 4 scaling sweep to **n_diploid = 500
(1000 haploid)** — Phase 1 stopped at 250 (500 haploid).
~22 min wall-clock for the new cells (Phase 1 cells re-ran
deterministically from the same RNG seeds).

#### 5-point sweep, last-run wall-clock (s)

| n_hap | tskit  | egrm    | PLINK |
|-------|--------|---------|-------|
| 50    | 2.67   |  5.48   | 0.06  |
| 100   | 3.86   | 12.23   | 0.15  |
| 200   | 6.92   | 27.36   | 0.18  |
| 500   | 24.82  | 78.93   | 0.22  |
| 1000  | 88.99  | 182.47  | 0.49  |

Log-log slopes (50 → 1000): tskit 1.17, egrm 1.17, PLINK 0.61.
PLINK still wins on raw wall-clock at this scale, but its
RSS slope is steepening — PLINK RSS 22 MB → 83 MB across N
50 → 1000 (slope 0.44), tskit 110 → 188 (slope 0.18), egrm
180 → 247 (slope 0.10). Extrapolating, PLINK's memory
overtakes tskit at **~16k haploid** — the predicted
biobank-scale crossover. Confirming this empirically is a
follow-up sub-iteration (n_diploid ≥ 5000 would need a few
CPU-hours of additional sim + egrm-on-large-N).

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

## Optional install extras

The figure-gen modules require `matplotlib`. Install with the
`figures` extra:

```bash
pip install -e ".[figures]"
```
