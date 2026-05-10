# graphpop-sim

Simulation integration for GraphPop (M12, Phase 6).

Three subcommands under one CLI:

- `graphpop-sim msprime` — orchestrate an msprime simulation from a
  YAML config; dumps `.trees` and (optionally) ingests it as a new
  `:ARGRun` via `graphpop-import.ARGIngester`.
- `graphpop-sim slim` — orchestrate a SLiM simulation (forward-time,
  selection-aware) from a curated template + parameter overrides.
- `graphpop-sim abc` — Rejection-ABC (Beaumont, Zhang & Balding
  2002): simulate N draws under a prior, compute summary statistics
  (π, θ_W, Tajima's D, F_ST), accept top-ε by Euclidean distance to
  observed, return posterior sample.

## Install

```sh
pip install graphpop-sim                      # msprime + ABC core
pip install "graphpop-sim[slim]"             # also SLiM (pyslim)
pip install "graphpop-sim[ingest]"           # also Neo4j ingest
```

## YAML config (msprime)

```yaml
n_diploid: 100
sequence_length: 1000000
recombination_rate: 1e-8
mutation_rate: 1e-8
demography:
  - {time: 0, Ne: 10000}
  - {time: 1000, Ne: 5000}
  - {time: 10000, Ne: 20000}
seed: 42
```

## Defaults

- Rejection-ABC fraction `epsilon = 0.05`
- Summary stats: π, θ_W, Tajima's D, mean F_ST across population pairs
- All operations deterministic given a seed
