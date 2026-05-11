# PLAN — KING-robust wrapper (Step H2)

**Status**: in progress, 2026-05-11.
**Driving rule**: plan-before-act.
**Companion artefacts**: `paper/paper2_kinship_arg/benchmark_plan.md`
(competitor table § 3; Fig 4a/b panels); `IMPLEMENTATION_PLAN.md`.
**Phase**: 1 (simulation-only); also Phase 2 Fig 5c.

## Goal

Wrap the KING-robust (Manichaikul et al. 2010, DOI
10.1093/bioinformatics/btq559) kinship inference so its output
is directly comparable with GraphPop's `relate.classify` /
`branch_grm` outputs. The wrapper:

1. Shells out to `king` (no version suffix in the binary name).
2. Auto-converts VCF → BED via PLINK if needed (KING accepts
   BED only).
3. Captures wall-clock + RAM via `profile_command`.
4. Parses KING's `.kin` (within-family) and `.kin0`
   (between-family) tab-separated output files.
5. Emits a normalised TSV (`king.tsv`) with columns
   `sample_a, sample_b, kinship` — one row per off-diagonal
   pair (KING does not emit self-pairs).
6. Writes a `receipt.json` alongside.

## Important note on the statistic

**KING reports phi (kinship coefficient), not a GRM entry.**

- PLINK GRM: standardised genotype covariance; diagonals
  ≈ 1 + inbreeding; off-diagonals 0 for unrelated, ~0.5 for
  parent-child.
- KING phi: kinship coefficient (Manichaikul 2010); phi(self) =
  0.5 by definition, phi(parent, child) = 0.25, phi(full-sib) =
  0.25, phi(unrelated) ≈ 0.

These are not directly comparable on the same axis. The
normalised TSV preserves whatever the tool computes; figure-
generation code is responsible for any conversion. For Paper 2
benchmarks, KING is the **classical-kinship reference** while
PLINK GRM is the **GRM reference**.

## Architecture

```
src/graphpop_bench/competitors/king.py
├── KingResult                # dataclass
├── KingRunner                # is_available(), run(...)
├── parse_kin_files(...)      # standalone (.kin + .kin0 parser)
└── _build_normalised_tsv()   # parsed rows → schema TSV
```

CLI: `graphpop-bench run king --input <vcf-or-bfile> --output <dir>`.

## KING output formats

### `.kin0` (between-family pairs)

```
FID1  IID1  FID2  IID2  N_SNP  HetHet  IBS0  HetConc  HomIBS0  Kinship
F1    A     F2    B     1000   0.123   0.001 ...      ...      0.0001
```

Tab-separated text. The `Kinship` column is phi. Column count
varies slightly across KING versions (2.x adds more columns over
time); the parser keys by header name, not position.

### `.kin` (within-family pairs)

Same general shape, single FID column. Within-family pairs in
simulated cohorts are typically empty (all samples are in their
own family) — so this file may be empty.

## Output schema (normalised)

`king.tsv`:

```
sample_a    sample_b    kinship
S0001       S0002       0.2503
S0001       S0003       0.0001
S0002       S0003       0.0002
...
```

- One row per **off-diagonal** sample pair (no self-pairs).
- `kinship` is the KING phi from the `Kinship` column.
- IID-only (FID dropped after parsing).

## Tests strategy

Parser tests run without KING; integration test skips when
`king` (or PLINK for VCF inputs) is missing.

## Out of scope (deferred)

- KING's `.seg` IBD-segment output — superseded by M11 IBD logic.
- `.con` condensed-identity-descriptors parsing.
- Sparse-output modes (`--related --degree 3`) — supported via
  `mode="related"` + `related_degree=N` but tested only in v1
  smoke; full sparse-vs-dense semantics deferred.
- KING-homo (non-robust mode) — not used; we ship `--kinship`
  which is the standard.

## Execution order

1. Write `king.py` (parser first; runner second).
2. Write `tests/test_king.py` (parser unit + integration).
3. Extend `cli.py` with the `run king` subcommand.
4. Re-install + run tests.
5. Commit + push.
