"""Paper 2 Fig 3c driver — composed-predicate rel-err.

Mirrors the Fig 1d/1e pipeline for the predicate-composition
panel (G2 novelty claim). For each of the four conditional
fixtures the M4.1 unit tests already validate against, this
driver:

1. Invokes the Java test with `-Dgraphpop.bench.dump.dir` to
   dump GraphPop's actual procedure output as H1-schema TSV.
2. Loads the frozen `conditional_egrm` reference JSON (these
   ARE Paper 2's predicate-algebra ground truth).
3. Joins by sample pair, computes per-entry rel-err, emits
   tidy panel CSVs for the figure renderer.

Headline panel: `branch_grm_composed_pathway_time.tsv` vs
`egrm_expected_composed_pathway_time.json` — the
`restrict_to_pathway × time_window` predicate composition.

Reuses `fig1de_panels` helpers verbatim; the only new logic is
multi-predicate orchestration + a cross-predicate summary CSV.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

from .fig1de_panels import (
    JoinedRow,
    join_pairs,
    load_egrm_json_matrix,
    load_pair_tsv,
    matrix_to_pair_map,
    summarise,
    write_panel_csv,
)


# Mapping: predicate label → (Java TSV filename, JSON reference filename).
PREDICATES: Dict[str, tuple[str, str]] = {
    "composed_pathway_time": (
        "branch_grm_composed_pathway_time.tsv",
        "egrm_expected_composed_pathway_time.json",
    ),
    "pathway_half": (
        "branch_grm_pathway_half.tsv",
        "egrm_expected_pathway_half.json",
    ),
    "time_window": (
        "branch_grm_time_window.tsv",
        "egrm_expected_time_window.json",
    ),
    "consequence_missense": (
        "branch_grm_consequence_missense.tsv",
        "egrm_expected_consequence_missense.json",
    ),
}


# ---------------------------------------------------------------------------
# Integration steps
# ---------------------------------------------------------------------------

def run_java_dump(
    repo_root: Path, java_dump_dir: Path,
    *, mvn_binary: str = "mvn",
) -> None:
    """Invoke the Maven test that dumps GraphPop's branch_grm TSVs
    (unconditional + 3 conditional fixtures)."""
    java_dump_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        mvn_binary, "test",
        "-Dtest=BranchGrmProcedureTest",
        f"-Dgraphpop.bench.dump.dir={java_dump_dir}",
        "-q",
    ]
    cwd = repo_root / "graphpop-procedures"
    print(f"[fig3c] invoking (cwd={cwd}): {' '.join(cmd)}",
          file=sys.stderr)
    r = subprocess.run(cmd, cwd=cwd)
    if r.returncode != 0:
        raise RuntimeError(
            f"mvn test exited with code {r.returncode}")


# ---------------------------------------------------------------------------
# Per-predicate join
# ---------------------------------------------------------------------------

@dataclass
class PredicateResult:
    label: str
    rows: List[JoinedRow]
    summary: dict
    panel_csv: Path

    def as_summary_row(self) -> dict:
        return {"predicate": self.label, **self.summary}


def join_one_predicate(
    label: str,
    java_dump_dir: Path,
    json_reference_dir: Path,
    output_dir: Path,
) -> PredicateResult:
    """Load Java TSV + JSON reference for one predicate; emit
    `{label}_panel_data.csv` and return summary stats."""
    tsv_name, json_name = PREDICATES[label]
    tsv_path = java_dump_dir / tsv_name
    json_path = json_reference_dir / json_name
    sample_ids, matrix = load_egrm_json_matrix(json_path)
    ref_pairs = matrix_to_pair_map(sample_ids, matrix)
    gp_pairs = load_pair_tsv(tsv_path)
    rows = join_pairs(ref_pairs, gp_pairs)
    panel_csv = output_dir / f"fig3c_{label}_panel_data.csv"
    write_panel_csv(rows, panel_csv)
    summary = summarise(rows)
    return PredicateResult(
        label=label, rows=rows, summary=summary,
        panel_csv=panel_csv)


def write_summary_csv(
    results: List[PredicateResult], path: Path,
) -> None:
    """Cross-predicate one-row-per-predicate summary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "predicate", "n_pairs",
            "max_rel_err", "max_abs_diff",
            "mean_rel_err", "mean_abs_diff",
        ])
        for r in results:
            s = r.summary
            w.writerow([
                r.label,
                s["n_pairs"],
                f"{s['max_rel_err']:.6g}",
                f"{s['max_abs_diff']:.6g}",
                f"{s['mean_rel_err']:.6g}",
                f"{s['mean_abs_diff']:.6g}",
            ])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_REPO_ROOT = Path("/mnt/data/GraphPop")
DEFAULT_REF_DIR = (
    DEFAULT_REPO_ROOT / "graphpop-procedures/src/test/resources"
)
DEFAULT_OUTPUT = (
    DEFAULT_REPO_ROOT
    / "paper/paper2_kinship_arg/benchmarks/fig3_out"
)


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--repo-root", type=Path, default=DEFAULT_REPO_ROOT)
    p.add_argument("--json-reference-dir", type=Path,
                   default=DEFAULT_REF_DIR)
    p.add_argument("--java-dump-dir", type=Path, required=True,
                   help="Where the Java test writes branch_grm TSVs")
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--predicates", nargs="+",
                   default=list(PREDICATES.keys()),
                   choices=list(PREDICATES.keys()))
    p.add_argument("--skip-mvn", action="store_true",
                   help="Skip mvn invocation; reuse existing dumps")
    p.add_argument("--mvn-binary", default="mvn")
    p.add_argument("--max-rel-err-gate", type=float, default=1e-6,
                   help="Exit non-zero if any predicate exceeds gate")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_mvn:
        run_java_dump(args.repo_root, args.java_dump_dir,
                      mvn_binary=args.mvn_binary)

    results: List[PredicateResult] = []
    worst = 0.0
    for label in args.predicates:
        r = join_one_predicate(
            label, args.java_dump_dir, args.json_reference_dir,
            args.output_dir)
        results.append(r)
        worst = max(worst, r.summary["max_rel_err"])
        print(f"[fig3c:{label}] {r.summary} → {r.panel_csv}",
              file=sys.stderr)

    summary_csv = args.output_dir / "fig3c_summary.csv"
    write_summary_csv(results, summary_csv)
    print(f"[fig3c] summary → {summary_csv}", file=sys.stderr)

    if worst >= args.max_rel_err_gate:
        print(
            f"[fig3c] FAIL: max rel-err {worst:.3e} ≥ gate "
            f"{args.max_rel_err_gate:.3e}", file=sys.stderr)
        return 2
    print(
        f"[fig3c] OK: max rel-err {worst:.3e} < gate "
        f"{args.max_rel_err_gate:.3e}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
