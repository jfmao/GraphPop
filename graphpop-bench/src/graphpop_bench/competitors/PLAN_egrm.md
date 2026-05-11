# PLAN — H4: egrm reference wrapper

**Date:** 2026-05-11
**Owner:** GraphPop (Paper 2 Phase 1)
**Discipline:** plan-before-act (per
`paper/paper2_kinship_arg/literature_survey.md` + roadmap)

## Goal

Wrap the **egrm** package (Fan, Mancuso & Chiang 2022,
*PLOS Genetics*; `pip install egrm`, currently 0.1) as a
first-class competitor in `graphpop-bench`. The package exposes
`varGRM_C(trees, ...) -> (egrm, vargrm, total_mu)` where `egrm`
is the N×N double-centred eGRM matrix that GraphPop's M4.1
branch-GRM procedure was originally validated against (relative
error < 10⁻⁶ per the test suite at
`graphpop-procedures/src/test/java/.../KinshipProcedureTest.java`).

For Paper 2 this wrapper produces the **reference ground-truth
matrix** that Fig 1d/1e (branch-GRM rel-err) compares GraphPop's
output against, and the head-to-head competitor for the Fig 4a/b
sim-scaling timing panels (alongside H1 PLINK, H2 KING,
H3 tskit).

## Inputs / outputs

**Input**: a tskit `.trees` file path.
**Output directory** contains:

| File | Contents |
|------|----------|
| `egrm.npy` | float64 (N, N) eGRM matrix |
| `vargrm.npy` | float64 (N, N) varGRM matrix (when `var=True`) |
| `egrm.ids` | newline-separated sample node ids (matches H3) |
| `egrm.tsv` | normalised TSV: `sample_a, sample_b, kinship` (upper triangle + diagonal — matches H1 / H3 schema) |
| `receipt.json` | standard graphpop-bench receipt |

## API surface

```python
class EgrmRunner:
    def __init__(self, python_binary: str | None = None): ...

    @staticmethod
    def is_available() -> bool: ...

    @staticmethod
    def detect_version() -> str | None: ...

    def run(
        self,
        input_path: str | Path,
        output_dir: str | Path,
        *,
        compute_var: bool = True,
        rlim: float = 0.0,
        alim: float | None = None,   # None → math.inf
        left: float = 0.0,
        right: float | None = None,  # None → math.inf
        seed: int | None = None,
        graphpop_commit: str | None = None,
    ) -> EgrmResult: ...


@dataclass
class EgrmResult:
    egrm: np.ndarray            # (N, N)
    vargrm: np.ndarray | None   # (N, N) when compute_var=True
    total_mu: float             # egrm's mu normalisation scalar
    sample_ids: list[str]
    output_dir: Path
    normalised_tsv: Path
    profiling: ProfilingResult
```

Module is **both importable AND runnable**:

```
python -m graphpop_bench.competitors.egrm \
    <input.trees> <output_dir> [<var_flag>] [<rlim>] [<alim>] \
    [<left>] [<right>]
```

`var_flag` ∈ {`var`, `novar`}; the four numeric args accept
`inf` (alim, right defaults).

## Subprocess pattern (mirrors H3)

The inner subprocess loads the `.trees` file, calls
`egrm.varGRM_C(trees, ...)`, writes `egrm.npy` (+ optional
`vargrm.npy`) + `egrm.ids` + a small `meta.json` containing
`total_mu`. The outer wrapper:

1. Builds the cmd `python -m graphpop_bench.competitors.egrm ...`.
2. Calls `profile_command` for RSS + wall-clock + CPU.
3. Reads back the artefacts via `_load_egrm_files`.
4. Emits the H1-schema TSV (one row per (i ≤ j) pair, including
   diagonal — matches PLINK GRM + tskit branch GRM rows).
5. Writes the `receipt.json` with `tool="egrm"`,
   `tool_version=egrm.__version__`, `n_samples`, `total_mu`,
   `compute_var`.

RSS isolation matters: `egrm.varGRM_C` allocates two N×N C
matrices internally, and a tqdm bar — running it in-process
contaminates pytest / click memory measurements just like the
tskit branch_grm path.

## Normalised TSV schema

`sample_a, sample_b, kinship`, where `kinship` is the **eGRM
matrix entry G_ij**. Upper triangle + diagonal. Identical to
H1 (PLINK GRM) and H3 (tskit branch GRM) schemas — Fig 1d/1e
join code can mechanically pair any two competitor TSVs.

The `vargrm.npy` is kept on disk for downstream uncertainty
panels (Fig 2 G1 posterior MSE comparison) but does NOT enter
the TSV — the TSV is point-estimate kinship only.

## Test plan

- `test_load_egrm_files_roundtrip` — synthetic .npy + .ids
  round-trip including the varGRM file.
- `test_load_egrm_files_missing_*` — fail-fast on missing files.
- `test_build_normalised_tsv_*` — H1-schema diagonal + dimensions.
- `test_is_available_returns_bool` — env-portable.
- `test_detect_version_matches_available` — `'0.1'` when present.
- `test_run_rejects_bad_args` — argument validation.
- `test_run_on_20_sample_fixture` (skip-if-missing) — uses
  the same fixture as H3 (`egrm_fixture_20samples.trees`).
  Asserts:
  - egrm shape (20, 20), symmetric, finite.
  - Row sums ≈ 0 (double-centred property).
  - varGRM shape (20, 20), symmetric, finite, all-positive
    diagonal.
  - TSV row count = 210 (= 20·21/2).
  - Receipt records `tool == "egrm"`, `compute_var == True`,
    `total_mu > 0`.
- `test_run_compute_var_false` — varGRM file absent, but egrm
  + TSV still produced.
- `test_run_subprocess_failure_surfaces_runtime_error` — bogus
  input path causes subprocess to exit non-zero.

Expected suite: 54 → ~63 tests passing.

## Out of scope (deferred)

- `egrm.mTMRCA_C` (mean-TMRCA matrix) — not part of Paper 2's
  benchmark plan; skip the wrapper for now.
- Custom `g` allele-frequency-scaling function — leave at the
  package default (standard GRM scaling).
- Custom `Gmap` (recombination map) — leave at the package
  default (uniform map); Paper 3 may revisit.
- Configurable per-population stratification — not part of the
  egrm API; that's Paper 2's `branch_grm_by_ancestry`
  territory.

## CLI

`graphpop-bench run egrm --input <.trees> --output <dir>
[--no-var] [--rlim 0] [--alim inf] [--left 0] [--right inf]
[--seed N] [--graphpop-commit SHA]`.

## Mathematical note (for the paper's appendix)

The egrm "variance" matrix `vargrm` returned by `varGRM_C` is
the Poisson-empirical-Bayes shrinkage variance from
Fan et al. 2022 § Methods. For Paper 2's G1 panel (posterior
MSE), we compare this against GraphPop's `branch_grm_posterior`
Welford-aggregated variance across SINGER posterior samples —
the two are conceptually different (Fan's is per-tree-prior
Poisson uncertainty; ours is per-ARG-posterior epistemic
uncertainty). The wrapper preserves both so the figure code
makes that distinction explicit.

## Methodology / discipline notes

- Plan written **before** any code edit (discipline rule).
- Pattern mirrors H3 exactly to keep the cross-wrapper code
  surface small.
- `egrm.varGRM_C` is C-accelerated; use it (not the pure
  Python `varGRM`) by default.
- The subprocess writes a `meta.json` (one line, ~50 bytes) so
  the outer can record `total_mu` in the receipt without
  re-running egrm.
