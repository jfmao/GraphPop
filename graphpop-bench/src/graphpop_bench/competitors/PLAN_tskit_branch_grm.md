# PLAN — tskit branch_grm wrapper (Step H3)

**Status**: in progress, 2026-05-11.
**Driving rule**: plan-before-act.
**Companion artefacts**: `paper/paper2_kinship_arg/benchmark_plan.md`
(Fig 1d/e, Fig 4a/b/d, sensitivity ED); `IMPLEMENTATION_PLAN.md`.
**Phase**: 1.

## Goal

Wrap tskit's `TreeSequence.genetic_relatedness_matrix(...,
mode='branch')` — the branch GRM implementation shipped with
tskit by Tang & Chiang 2025 *Genetics* — so it's a first-class
competitor in the benchmark harness on the same axes as PLINK
2.0 GRM (H1) and KING-robust (H2). The wrapper:

1. Accepts a tskit `.trees` file (no VCF/BED conversion needed).
2. Spawns a Python subprocess that imports tskit and computes
   the GRM. Subprocess isolation gives clean RSS measurement
   (no pytest/click baseline contamination), matching the
   H1/H2 instrumentation pattern.
3. Persists the resulting (n, n) matrix + sample list as
   `grm.npy` + `grm.ids` in the output dir.
4. Reads them back; emits the normalised TSV with schema
   `sample_a, sample_b, kinship` (one row per unordered pair
   *including diagonal* — matches H1 schema).
5. Writes a `receipt.json` alongside.

## Architecture

```
src/graphpop_bench/competitors/tskit_branch_grm.py
├── TskitBranchGrmResult        # dataclass
├── TskitBranchGrmRunner        # is_available(), run(...)
├── _inner_compute_and_dump()   # subprocess entry; loads .trees,
│                                 # writes grm.npy + grm.ids
├── _load_grm_files()           # round-trip the inner output
└── _build_normalised_tsv()     # GRM array → schema TSV
```

CLI: `graphpop-bench run tskit_branch_grm --input <ts.trees>
--output <dir>`.

The module is *both* importable as a Python wrapper AND
executable via `python -m
graphpop_bench.competitors.tskit_branch_grm <input> <output>
[<mode>]` for the subprocess path.

## Statistic equivalence

tskit `mode='branch'` GRM is the same statistic that
GraphPop's `kinship.branch_grm` computes. M4.1 already
validates GraphPop's implementation to rel-err < 10⁻⁶ vs
`egrm.varGRM`. The tskit wrapper is the reference for
wall-clock + RSS comparison on Paper 2 Fig 4a/b.

## Output schema (normalised)

`tskit_branch_grm.tsv`:

```
sample_a    sample_b    kinship
0           0           1.234
0           1           0.456
1           1           1.567
...
```

- One row per unordered pair (i ≤ j) **including diagonal**
  — matches H1 schema for direct join with PLINK output.
- `sample_a` / `sample_b` are tskit sample node-ids (integers,
  printed as strings).
- `kinship` is the float64 GRM entry from tskit.

## Wrapper API

```python
@dataclass
class TskitBranchGrmResult:
    grm: np.ndarray              # (n, n)
    sample_ids: list[str]
    output_dir: Path
    normalised_tsv: Path
    profiling: ProfilingResult

class TskitBranchGrmRunner:
    def __init__(self, python_binary: str | None = None):
        ...

    @staticmethod
    def is_available() -> bool:
        """tskit is importable in the current environment."""

    @staticmethod
    def detect_version() -> str | None:
        """tskit.__version__ if importable."""

    def run(
        self,
        input_path: Path,            # .trees
        output_dir: Path,
        *,
        mode: str = "branch",        # "branch" | "site"
        seed: int | None = None,
        graphpop_commit: str | None = None,
    ) -> TskitBranchGrmResult:
        ...
```

## Subprocess entry signature

```
python -m graphpop_bench.competitors.tskit_branch_grm \
    <input.trees> <output_dir> [<mode>]
```

Default mode: `branch`. Writes:

- `<output_dir>/grm.npy` — float64 (n, n) matrix.
- `<output_dir>/grm.ids` — newline-separated sample node ids.

## Tests (`tests/test_tskit_branch_grm.py`)

1. **`is_available()` / `detect_version()`** — bool + string
   when tskit is importable; no throw when missing.
2. **TSV-emission unit test** (no tskit): hand-crafted 3×3
   GRM → 6-row TSV with correct header + values.
3. **Inner-compute test** (skip if tskit missing): invoke
   `_inner_compute_and_dump` on the 20-sample fixture at
   `graphpop-procedures/src/test/resources/egrm_fixture_20samples.trees`;
   verify `grm.npy` shape (20, 20), symmetric, diagonals > 0.
4. **End-to-end runner test** (skip if tskit missing):
   `TskitBranchGrmRunner.run(...)` on the same fixture;
   verify receipt JSON, normalised TSV with 210 rows
   (20 + 20\*19/2), symmetric GRM in the dataclass.
5. **Module-as-script test** (skip if tskit missing): run
   `python -m graphpop_bench.competitors.tskit_branch_grm`
   via subprocess; verify exit 0, outputs present.

## Skip strategy

- `pytest.mark.skipif(not TskitBranchGrmRunner.is_available(),
  reason=...)` on every tskit-dependent test. The local
  `graphmana` env has tskit (already used by graphpop-sim and
  the M4.1 fixtures); CI without it skips cleanly.

## Success criterion (Step H3)

- `pytest graphpop-bench/tests/test_tskit_branch_grm.py` green
  on `graphmana`.
- `graphpop-bench run tskit_branch_grm --help` works.
- Full pipeline on the 20-sample fixture produces a symmetric
  20×20 GRM, normalised TSV with 210 rows, receipt JSON.

## Out of scope (deferred)

- `mode='site'` benchmarking — supported via parameter but not
  exercised in v1 tests.
- `genetic_relatedness_vector` (Algorithm V matvec) — separate
  H3b later if Fig 4a needs matvec-only benchmarks.
- Multi-window outputs — v1 uses `windows=None` (whole-genome).
- Sample subsetting — v1 always uses all `ts.samples()`.
- VCF/BED inputs — tskit reads `.trees` natively; users
  convert via `tsinfer + tsdate` separately.

## Execution order

1. Write `tskit_branch_grm.py` (subprocess inner + wrapper).
2. Write `tests/test_tskit_branch_grm.py`.
3. Update `competitors/__init__.py` exports.
4. Extend `cli.py` with the `run tskit_branch_grm` subcommand.
5. Re-install + run tests.
6. Commit + push.
