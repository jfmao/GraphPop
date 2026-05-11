# graphpop-bench — Implementation Plan

**Status**: scaffold initial draft, 2026-05-11.
**Driving rule**: plan-before-act (per
`feedback_paper_discipline.md`).
**Companion artefacts**: `paper/paper2_kinship_arg/benchmark_plan.md`
(per-figure benchmark recipes); the personal plan at
`~/.claude/plans/for-next-stage-of-peppy-noodle.md` Step G.

## Goal

Build the cross-paper benchmarking infrastructure that Papers 2–5
will share. v1 ships only the **Phase 1 (simulation-only)** scope
required by Paper 2; Phase 2 wrappers (ARG-RHE, Threads-on-real-
data, ADMIXTURE, SINGER-on-1000G) deferred until Phase 1
benchmarks are complete and reviewable.

## Package layout

```
graphpop-bench/
├── pyproject.toml          # deps: click, numpy, pyyaml; [compare] = competitor extras
├── README.md
├── IMPLEMENTATION_PLAN.md  # this file
├── src/graphpop_bench/
│   ├── __init__.py
│   ├── profiling.py        # /usr/bin/time -v subprocess wrapper
│   ├── competitors/        # one module per external tool
│   │   ├── __init__.py
│   │   └── (wrappers added in Step H1+)
│   └── cli.py              # graphpop-bench CLI entry
└── tests/
    └── test_profiling.py   # smoke tests with simple commands
```

Same package-layout pattern as `graphpop-sim`, `graphpop-gnn`,
`graphpop-pedigree` (existing sister packages on `develop`).

## What Step G produces (this scaffold)

1. `pyproject.toml` with light deps (click, numpy, pyyaml).
2. `README.md` — usage sketch + Phase 1/2 scope note.
3. `src/graphpop_bench/__init__.py` — package marker + version.
4. `src/graphpop_bench/profiling.py` — the core utility:
   - `ProfilingResult` dataclass: wall-clock, RSS peak,
     user/system CPU, exit code, stdout/stderr.
   - `profile_command(cmd, output_dir, seed=None)` — runs a
     command under `/usr/bin/time -v`, parses output, writes
     a receipt JSON next to the user-requested output dir.
   - `ReceiptWriter` — dumps GraphPop SHA + competitor version
     + hardware fingerprint alongside every run (per
     `benchmark_plan.md` § 8).
5. `src/graphpop_bench/competitors/__init__.py` — stub for the
   subsequent `Step H` wrappers.
6. `src/graphpop_bench/cli.py` — `graphpop-bench profile <cmd>`
   entry (single command initially; subcommands added per
   wrapper).
7. `tests/test_profiling.py` — smoke tests that profile a
   trivial command (`sleep 0.1`, `python -c "pass"`) and assert
   the dataclass fields populate. Skip the `/usr/bin/time -v`
   path on macOS (BSD time has different flags); use a
   pure-Python `resource.getrusage(RUSAGE_CHILDREN)` fallback.

## What Step G does NOT produce (deferred to H1+)

- PLINK 2.0 / KING / tskit / egrm / s-LDSC wrappers.
- Per-figure benchmark scripts.
- Figure-generation scripts.

## Core design decisions

### Profiling backend

Two paths, chosen at runtime:

- **Linux** (the user's deploy target per CLAUDE.md): use
  `/usr/bin/time -v -o <tmp_file> <cmd>` and parse the verbose
  output. Captures everything we need (Max RSS, wall clock,
  user/system CPU).
- **Other** (or `/usr/bin/time` unavailable): fall back to
  Python's `subprocess.run` + `resource.getrusage(RUSAGE_CHILDREN)`
  delta. Less detailed but portable; tests use this path.

### Receipt format

Per benchmark run, alongside the user's output TSV, dump a
`receipt.json`:

```json
{
  "graphpop_commit": "abc1234",
  "tool": "graphpop",
  "tool_version": "0.2.0.dev0",
  "seed": 42,
  "wall_clock_s": 12.3,
  "rss_peak_mb": 1024.5,
  "user_cpu_s": 11.8,
  "system_cpu_s": 0.5,
  "exit_code": 0,
  "host_fingerprint": "linux-x86_64-64gb",
  "timestamp_iso": "2026-05-11T..."
}
```

### CLI surface

Initial CLI: `graphpop-bench profile <output_dir> -- <cmd> <args...>`.
The `--` separator demarcates the profiled command from
`graphpop-bench` options. Future Step H1+ adds:

- `graphpop-bench run <competitor> <input>` — dispatches to a
  registered competitor wrapper.
- `graphpop-bench run-figure <fig_id> --phase 1` — runs every
  benchmark for a given Paper 2 figure.

These are documented in this plan but **not built in Step G**.

### Dependency footprint

Hard deps (always installed):
- `click >= 8.0`
- `numpy >= 1.22`
- `pyyaml >= 6.0`

Optional `[compare]` extras (per Step H1+ as wrappers land):
- tskit, egrm (Python competitors — runtime imports).

External binaries (not packaged):
- PLINK 2.0, KING, ADMIXTURE, Threads, SINGER, s-LDSC, ARG-RHE.
  Wrappers shell out; tests skip when the binary is missing.

## Success criterion (Step G)

- Package installs via `pip install -e graphpop-bench`.
- `graphpop-bench --help` works.
- `tests/test_profiling.py` passes locally.
- `profile_command(["sleep", "0.1"])` returns a `ProfilingResult`
  with `wall_clock_s ≥ 0.1` and `exit_code == 0`.

## Out of scope (Step G)

- Any competitor wrapper (Step H1+).
- Per-figure benchmark scripts (Step I).
- Figure-generation scripts (Step J).
- Real-data ETL (Phase 2).

## Risks (Step G)

- `/usr/bin/time -v` parsing is fragile (Linux GNU-time format
  varies slightly across distros). Mitigation: regex-based
  parser with key-value fallback; tests on the user's machine.
- Receipt schema may need fields we don't anticipate (GPU
  utilisation, disk I/O). Mitigation: add fields incrementally
  as wrappers reveal needs.

## Execution order

1. Write `pyproject.toml`, `README.md`, package skeleton (G2).
2. Write `profiling.py` + smoke tests (G3).
3. Install + verify + commit (G4).

Each is a small, well-defined edit; no per-step plan-before-act
needed beyond this file.
