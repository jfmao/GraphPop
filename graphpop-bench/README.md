# graphpop-bench

Cross-paper benchmarking + comparison harness for GraphPop.

Each GraphPop paper (Paper 2 kinship+ARG, Paper 3 recombination,
Paper 4 GNN, Paper 5 selection) ships with a head-to-head
comparison against the strongest existing baselines. This
package provides the shared infrastructure:

- **Profiling harness** — wraps any command under `/usr/bin/time
  -v` (Linux) or `resource.getrusage` fallback; captures wall-
  clock, RSS peak, user/system CPU; dumps a JSON receipt
  alongside the user's output for reproducibility.
- **Competitor wrappers** — one Python module per external tool
  (PLINK 2.0, KING-robust, tskit, egrm, s-LDSC, ADMIXTURE,
  Threads, SINGER, ARG-RHE). Each wrapper runs the competitor
  on the same input data + emits a GraphPop-schema-compatible
  output so paired comparison is mechanical.
- **CLI** — `graphpop-bench profile`, `graphpop-bench run`.

## Scope (v1)

Initial release covers **Phase 1 (simulation-only)** scope for
Paper 2: profiling + the Phase 1 wrapper subset (PLINK 2.0,
KING-robust, tskit branch_grm, egrm reference, s-LDSC).
**Phase 2 wrappers** (ARG-RHE, Threads-on-real-data, ADMIXTURE,
SINGER-on-1000G) are deferred until Phase 1 benchmarks are
complete and reviewable, per
`paper/paper2_kinship_arg/benchmark_plan.md` § 1.1.

## Install

```sh
pip install -e graphpop-bench                # core + profiling
pip install -e "graphpop-bench[compare]"     # also Python competitors
pip install -e "graphpop-bench[dev]"         # also pytest
```

External binaries (PLINK, KING, ADMIXTURE, etc.) install
separately; the wrappers shell out and skip cleanly when the
binary isn't on `PATH`.

## CLI sketch

```sh
# Profile any command
graphpop-bench profile ./output_dir -- sleep 0.1

# (Phase 1 wrappers will add subcommands; not in this release)
# graphpop-bench run plink_grm --input ...
# graphpop-bench run-figure fig1 --phase 1
```

## Companion artefacts

- `paper/paper2_kinship_arg/benchmark_plan.md` — per-figure
  benchmark recipes (Phase 1 + Phase 2).
- `docs/publication_benchmarking_plan.md` — Papers 2–5 strategic
  roadmap (cross-paper benchmarking framework).
- `IMPLEMENTATION_PLAN.md` (this folder) — the plan-before-act
  spec for this package.
