# PLAN — H5: s-LDSC partitioned-heritability wrapper

**Date:** 2026-05-11
**Owner:** GraphPop (Paper 2 Phase 1)
**Discipline:** plan-before-act (per
`paper/paper2_kinship_arg/literature_survey.md` + roadmap)

## Goal

Wrap **stratified LD Score Regression** (Finucane et al. 2015,
*Nature Genetics*; tool: `ldsc` at github.com/bulik/ldsc) as
the **annotation-stratified heritability competitor** for
Paper 2's Fig 3 panels:

- **Fig 3a** — pathway-restricted h² (GraphPop branch_grm with
  `restrict_to_pathway` predicate vs s-LDSC pathway annotation).
- **Fig 3b** — LoF-class h² (GraphPop with `mutation_filter` vs
  s-LDSC functional category annotation).

Paper 2's G2 (annotation-conditional GRM) is novel **vis-à-vis
the graph-native composition story**, not vis-à-vis s-LDSC's
ability to partition heritability — s-LDSC has done that since
2015. The benchmark question is **can GraphPop's composable
predicate recover the same per-category point estimates s-LDSC
does, with comparable standard errors, on a unified data model
that also supports kinship + ARG queries.**

## Inputs / outputs

`ldsc.py --h2` takes:
- `--h2 <sumstats>` — gzipped GWAS summary statistics
  (LDSC `.sumstats.gz` schema: `SNP, A1, A2, N, Z, …`).
- `--ref-ld-chr <prefix>` — per-chromosome LD-score files
  for stratified annotations.
- `--w-ld-chr <prefix>` — regression-weight LD-scores.
- `--overlap-annot` — partition with overlapping categories.
- `--frqfile-chr <prefix>` — allele frequencies.
- `--out <out_prefix>` — output prefix.

Produces:
- `<prefix>.results` — per-category partitioned-h² table.
- `<prefix>.log` — total h², intercept, λ_gc, sample-size,
  warnings.

**Output directory** of the wrapper contains:

| File | Contents |
|------|----------|
| `ldsc.results` | raw s-LDSC per-category table (passthrough) |
| `ldsc.log` | raw s-LDSC log (passthrough) |
| `partition_h2.tsv` | normalised per-category: `category, prop_snps, prop_h2, prop_h2_se, enrichment, enrichment_se, enrichment_p` |
| `total_h2.tsv` | normalised total: `total_h2, total_h2_se, intercept, intercept_se, lambda_gc, n_snps` |
| `receipt.json` | standard graphpop-bench receipt |

This wrapper does **not** emit a pairwise-kinship TSV (the
schema of H1/H3/H4 doesn't apply — s-LDSC produces partitioned
variance components, not pairwise kinship). The normalised
schema here is per-category-h².

## API surface

```python
class SLdscRunner:
    def __init__(self, ldsc_binary: str | None = None): ...

    @staticmethod
    def is_available() -> bool: ...

    @staticmethod
    def detect_version(binary: str) -> str | None: ...

    def run(
        self,
        sumstats: str | Path,
        output_dir: str | Path,
        *,
        ref_ld_chr: str | Path,
        w_ld_chr: str | Path,
        frqfile_chr: str | Path | None = None,
        overlap_annot: bool = True,
        extra_args: Sequence[str] | None = None,
        seed: int | None = None,
        graphpop_commit: str | None = None,
    ) -> SLdscResult: ...


@dataclass
class SLdscResult:
    categories: List[CategoryRow]   # parsed .results rows
    total_h2: float
    total_h2_se: float
    intercept: float
    intercept_se: float
    lambda_gc: float | None
    n_snps: int | None
    output_dir: Path
    partition_tsv: Path
    total_tsv: Path
    profiling: ProfilingResult


@dataclass
class CategoryRow:
    category: str
    prop_snps: float
    prop_h2: float
    prop_h2_se: float
    enrichment: float
    enrichment_se: float
    enrichment_p: float
```

## Binary detection (mirrors H1 / H2)

LDSC ships as a Python script `ldsc.py`. We detect via:

1. `shutil.which("ldsc.py")` (preferred — works with the
   conda-forge or bioconda installs).
2. Failing that, `shutil.which("ldsc")`.

The wrapper invokes `<python> <ldsc.py-path> ...` if the
detected binary's first line is a shebang we don't trust; in
practice modern installs make `ldsc.py` directly executable, so
we go via PATH like H1 / H2.

If not found, `is_available() → False`; the integration test
skips. Parser tests run from hand-crafted fixtures and need no
binary.

## Parser

`.results` is a header-row TSV from LDSC. The exact column
order has been stable across LDSC versions since 2018:

```
Category Prop._SNPs Prop._h2 Prop._h2_std_error Enrichment \
    Enrichment_std_error Enrichment_p Coefficient \
    Coefficient_std_error Coefficient_z-score
```

The wrapper:

- Parses the header row to map column names → indices (defensive
  against minor LDSC-version reorderings).
- Yields `CategoryRow` per data line; skips blank / malformed
  rows silently.
- Returns the list keyed by the `Category` column.

`.log` parsing — small grep-style extraction for:

- `Total Observed scale h2: <h2> (<se>)`
- `Intercept: <intercept> (<se>)`
- `Lambda GC: <lambda>` (when present; some runs omit)
- `<N> SNPs remain` (regression-window SNP count)

These are stable single-line patterns; if a pattern is missing
we record `None` rather than fail — the wrapper preserves the
raw `.log` for diagnostics.

## Subprocess invocation

Standard `profile_command(cmd, output_dir)` path, exactly like
H1 / H2. The `cmd` list:

```python
[ldsc_binary, "--h2", str(sumstats),
 "--ref-ld-chr", str(ref_ld_chr),
 "--w-ld-chr", str(w_ld_chr),
 "--out", str(out_prefix),
 "--overlap-annot"]
# + ["--frqfile-chr", str(frqfile_chr)] when given
# + extra_args
```

Non-zero exit → `RuntimeError` with stderr tail (same pattern
as H1 / H2 / H3 / H4). The `.log` is also preserved so failing
runs are diagnosable from disk artefacts.

## Test plan

- `test_parse_results_header_only` — empty body returns `[]`.
- `test_parse_results_row_canonical_format` — three-category
  fixture with realistic LDSC numbers; verify each
  `CategoryRow` field.
- `test_parse_results_handles_extra_columns` — newer LDSC
  versions may add `Coefficient_z-score` etc.; parser tolerates.
- `test_parse_results_missing_required_column_raises` — drop
  `Enrichment` from header → `ValueError`.
- `test_parse_log_extracts_total_h2_intercept_lambda` — log
  with all four metrics present.
- `test_parse_log_handles_missing_lambda` — log without the GC
  line returns `lambda_gc=None`.
- `test_build_partition_tsv` — schema row count + header.
- `test_build_total_tsv` — single-row TSV with column names.
- `test_is_available_returns_bool`.
- `test_detect_version_returns_none_for_nonexistent_binary`.
- `test_run_raises_clean_error_when_ldsc_missing`.
- `test_run_rejects_missing_required_args` (no `ref_ld_chr` →
  `TypeError` from the dataclass-like call).
- **Integration (skip-if-missing)**: only if `ldsc.py` is on
  PATH; the s-LDSC reference 1000G-EUR LD-scores are a large
  download we will NOT pull into CI. The test loops a tiny
  hand-crafted fixture; on every machine that lacks LDSC it
  skips. This is the same shape as H2's "skip if KING refuses
  this tiny fixture" guard.

Expected suite: 73 → ~85 tests passing.

## Out of scope (deferred)

- Generating sumstats from simulated phenotypes (that's the
  job of `graphpop-sim` + Fig 3 driver scripts in Step I, not
  this wrapper).
- LDSC reference panel construction (1000G EUR baseline LD
  scores). The user supplies these as paths to the wrapper.
- `ldsc.py --rg` (genetic correlation) — not needed for Paper 2.
- Direct comparison against MTAG, S-PrediXcan, etc. — outside
  Paper 2's claim space.

## CLI

```
graphpop-bench run s_ldsc \
    --sumstats traits.sumstats.gz \
    --ref-ld-chr ref/baseline. \
    --w-ld-chr weights/weights.hm3_noMHC. \
    --frqfile-chr freq/1000G.EUR.QC. \
    --output ./out
```

## Methodology / discipline notes

- Plan written **before** any code edit (discipline rule).
- Phase 1 plan: parser + skip-if-missing integration test only.
  Actual head-to-head Fig 3 panel execution uses driver scripts
  in Step I; the wrapper's job is to *make those drivers
  trivial*.
- The wrapper does not emit the H1-schema pairwise TSV because
  s-LDSC produces variance components, not pairwise kinship —
  schema divergence is correct here.
- Mirrors the H1/H2 binary-detection pattern (not the
  H3/H4 Python-subprocess pattern) because `ldsc.py` is a
  shippable executable, not a Python library we import.
