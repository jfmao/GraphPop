# PLAN — PLINK 2.0 GRM wrapper (Step H1)

**Status**: in progress, 2026-05-11.
**Driving rule**: plan-before-act.
**Companion artefacts**: `paper/paper2_kinship_arg/benchmark_plan.md`
(competitor table § 3; Fig 4a/b panels); `IMPLEMENTATION_PLAN.md`
(package architecture).
**Phase**: 1 (simulation-only); also used in Phase 2 Fig 5c.

## Goal

Wrap the PLINK 2.0 GRM (`--make-grm-bin`) computation so its
output can be compared head-to-head with GraphPop's
`branch_grm` output on identical inputs. The wrapper:

1. Shells out to `plink2` (or `plink` v1.9 if v2 is missing).
2. Captures wall-clock + RAM via the `profile_command` harness.
3. Parses the binary GRM (`.grm.bin` lower-triangle float32) +
   sample ID list (`.grm.id`).
4. Emits a normalised TSV (`plink_grm.tsv`) with columns
   matching the GraphPop kinship-output schema:
   `sample_a, sample_b, kinship` — one row per unordered
   sample pair (including diagonal i=j).
5. Writes the profiling receipt (`receipt.json`) alongside.

## Architecture

```
src/graphpop_bench/competitors/plink_grm.py
├── PlinkGrmResult           # dataclass
├── PlinkGrmRunner           # is_available(), run(...)
├── parse_grm_bin(...)       # standalone parser (testable w/o PLINK)
└── _build_normalised_tsv()  # GRM array → schema TSV
```

CLI integration: extend `cli.py` with
`graphpop-bench run plink_grm --input <vcf-or-bfile> --output <dir>`.

## Input formats

PLINK accepts:

- `--vcf <vcf.gz>` for VCF input.
- `--bfile <prefix>` for binary PED (.bed/.bim/.fam).
- `--pfile <prefix>` for PLINK 2.0 .pgen format.

v1 wrapper: support `--vcf` and `--bfile`. Other formats deferred.

## PLINK GRM binary format

`.grm.bin`: `float32` lower-triangle including diagonal, packed
row-major:

```
GRM[0,0]
GRM[1,0] GRM[1,1]
GRM[2,0] GRM[2,1] GRM[2,2]
...
```

For `n` samples there are `n*(n+1)/2` float32 entries. Sample
order from `.grm.id` (two-column FID/IID TSV).

## Output schema (normalised)

`plink_grm.tsv`:

```
sample_a    sample_b    kinship
S0001       S0001       1.00
S0001       S0002       0.51
S0002       S0002       0.99
...
```

- One row per unordered pair, **including diagonal** (sample
  with itself).
- IID-only (FID is ignored after the .grm.id is read; we keep
  whichever GraphPop uses as `sampleId`).
- `kinship` is the float32 GRM entry from PLINK.

This schema matches the output of (future) GraphPop kinship
exports for comparison.

## Wrapper API

```python
@dataclass
class PlinkGrmResult:
    grm: np.ndarray              # shape (n, n) symmetric
    sample_ids: list[str]
    output_dir: Path
    normalised_tsv: Path
    profiling: ProfilingResult

class PlinkGrmRunner:
    def __init__(self, plink_binary: str | None = None):
        ...

    @staticmethod
    def is_available() -> bool:
        """plink2 (preferred) or plink on PATH."""

    @staticmethod
    def detect_version(binary: str) -> str | None:
        """Run `<binary> --version` and return the version string."""

    def run(
        self,
        input_path: Path,           # .vcf(.gz) or BED prefix
        output_dir: Path,
        *,
        input_kind: str = "auto",   # "vcf" | "bfile" | "auto"
        extra_args: list[str] | None = None,
        seed: int | None = None,
        graphpop_commit: str | None = None,
    ) -> PlinkGrmResult:
        ...
```

## Tests (`tests/test_plink_grm.py`)

1. **Parser unit tests** (runnable WITHOUT PLINK):
   - Build a 5×5 symmetric matrix in numpy; pack to
     `.grm.bin` format; write `.grm.id`; verify the parser
     reconstructs the full matrix.
   - Test missing `.grm.id` raises `FileNotFoundError`.
   - Test malformed `.grm.bin` (wrong byte count) raises a
     clear error.

2. **TSV-emission unit tests**:
   - Given a hand-crafted 3×3 GRM, verify the emitted TSV has
     6 rows (3 diagonal + 3 off-diagonal), correct columns,
     correct values.

3. **`is_available()` smoke**:
   - Returns a bool. Does not throw if PLINK is missing.

4. **Integration test** (`pytest.mark.skipif` when PLINK
   absent):
   - Build a tiny hand-crafted VCF (5 samples, 20 SNPs).
   - Run wrapper; verify TSV emitted, GRM is symmetric,
     diagonals near 1.0.

## Skip strategy

If `plink2` (preferred) or `plink` (fallback for v1.9) is not
on `PATH`, the integration test is skipped via
`pytest.mark.skipif(not PlinkGrmRunner.is_available(), reason=...)`.
The parser + TSV emission tests run regardless — they don't
need the binary.

## Success criterion (Step H1)

- `pytest graphpop-bench/tests/test_plink_grm.py` green.
- Integration test passes when PLINK is on PATH; skips cleanly
  when not.
- `graphpop-bench run plink_grm --help` works.
- `parse_grm_bin()` round-trips a hand-crafted GRM bit-for-bit.

## Out of scope (deferred)

- `.king` / `--make-king` integration — separate H2 wrapper
  (KING-robust).
- Sparse GRM (`--make-grm-list`) — Phase 2.
- Multi-chromosome merging — wrapper accepts a single input
  for v1.
- VCF preprocessing (filtering, MAF, missingness) — caller's
  responsibility; the wrapper passes `extra_args` through to
  PLINK.

## Execution order

1. Write `plink_grm.py` (parser first; runner second).
2. Write `tests/test_plink_grm.py`.
3. Extend `cli.py` with the `run plink_grm` subcommand.
4. Re-install + run tests.
5. Commit + push.
