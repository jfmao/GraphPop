"""Paper 2 Fig 1d/1e driver — branch-GRM rel-err panel data.

Pipeline:
1. Invoke `mvn -pl graphpop-procedures test -Dtest=BranchGrmProcedureTest
   -Dgraphpop.bench.dump.dir=<java-dump-dir>` to produce GraphPop's
   actual branch_grm procedure output as H1-schema TSVs (the Java
   test only writes the TSV when the system property is set, after
   the assertion passes).
2. Invoke the H4 egrm wrapper (`graphpop-bench run egrm`) on the
   same `.trees` fixture to produce a fresh egrm.tsv reference,
   with profiling receipt.
3. Load the frozen JSON references:
   - egrm_expected_20samples.json     → Fig 1d ground truth
   - egrm_expected_pathway_half.json  → Fig 1e ground truth
4. Join {GraphPop output} × {JSON reference} by canonical sample
   pair `(min(a,b), max(a,b))`; compute per-entry rel-err.
5. Emit fig1d_panel_data.csv + fig1e_panel_data.csv.

The pure-logic helpers (load_*, join_pairs, compute_relerr) are
unit-tested in `tests/test_fig1de_driver.py`; the mvn + egrm
wrapper invocations are integration-only.
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

PairKey = Tuple[int, int]


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_egrm_json_matrix(path: Path) -> tuple[List[int], List[List[float]]]:
    """Read an `egrm_expected_*.json` fixture into (sample_ids, matrix).

    Schema variants produced by `build_egrm_fixture.py`:
      - Unconditional fixture: keys include `sample_ids` (list of ints).
      - Conditional fixtures (pathway_half, time_window, etc.): no
        `sample_ids` key, only `n_samples` + `matrix`.  We synthesise
        IDs as `range(n_samples)` in that case — this matches the
        sample indexing the Java test + H4 wrapper both use.
    """
    if not path.exists():
        raise FileNotFoundError(f"egrm JSON reference not found: {path}")
    blob = json.loads(path.read_text())
    matrix = [[float(x) for x in row] for row in blob["matrix"]]
    n_from_matrix = len(matrix)
    if "sample_ids" in blob:
        sample_ids = [int(s) for s in blob["sample_ids"]]
    elif "n_samples" in blob:
        sample_ids = list(range(int(blob["n_samples"])))
    else:
        sample_ids = list(range(n_from_matrix))
    n = len(sample_ids)
    if n_from_matrix != n or any(len(row) != n for row in matrix):
        raise ValueError(
            f"matrix shape mismatch in {path}: {n_from_matrix} rows "
            f"vs {n} sample_ids")
    return sample_ids, matrix


def load_pair_tsv(path: Path) -> Dict[PairKey, float]:
    """Read an H1-schema TSV (`sample_a, sample_b, kinship`) into a
    pair-keyed map. Pairs are stored under the canonical key
    `(min(a,b), max(a,b))` so look-ups are order-invariant.
    """
    if not path.exists():
        raise FileNotFoundError(f"TSV not found: {path}")
    out: Dict[PairKey, float] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "kinship" not in reader.fieldnames:
            raise ValueError(
                f"TSV {path} missing 'kinship' column; "
                f"header: {reader.fieldnames!r}")
        for row in reader:
            a = int(row["sample_a"])
            b = int(row["sample_b"])
            key = (a, b) if a <= b else (b, a)
            out[key] = float(row["kinship"])
    return out


# ---------------------------------------------------------------------------
# Pure-logic helpers — join + rel-err
# ---------------------------------------------------------------------------

def matrix_to_pair_map(
    sample_ids: List[int], matrix: List[List[float]],
) -> Dict[PairKey, float]:
    """Flatten an (n, n) matrix to a canonical-pair-keyed map (upper
    triangle including diagonal).
    """
    out: Dict[PairKey, float] = {}
    n = len(sample_ids)
    for i in range(n):
        for j in range(i, n):
            key = (sample_ids[i], sample_ids[j])
            if key[0] > key[1]:
                key = (key[1], key[0])
            out[key] = matrix[i][j]
    return out


@dataclass
class JoinedRow:
    """One per-pair joined record (rel-err + absolute diff)."""
    sample_a: int
    sample_b: int
    egrm_ref: float
    graphpop: float
    abs_diff: float
    rel_err: float
    egrm_h4: float | None = None  # optional fresh-wrapper reading

    def as_csv_row(self) -> list[str]:
        return [
            str(self.sample_a),
            str(self.sample_b),
            f"{self.egrm_ref:.10g}",
            f"{self.graphpop:.10g}",
            f"{self.abs_diff:.6g}",
            f"{self.rel_err:.6g}",
            ("" if self.egrm_h4 is None else f"{self.egrm_h4:.10g}"),
        ]


def join_pairs(
    egrm_ref: Dict[PairKey, float],
    graphpop: Dict[PairKey, float],
    *,
    egrm_h4: Dict[PairKey, float] | None = None,
) -> List[JoinedRow]:
    """Inner-join two pair maps; emit `JoinedRow`s with rel_err.

    Raises if a key is in `egrm_ref` but missing from `graphpop`
    (this would indicate either a broken Java dump or a corrupted
    fixture — a silent join would hide the bug).
    """
    rows: List[JoinedRow] = []
    for key in sorted(egrm_ref.keys()):
        if key not in graphpop:
            raise KeyError(
                f"pair {key!r} present in egrm reference but missing "
                f"from graphpop output; cannot join")
        x = egrm_ref[key]
        y = graphpop[key]
        abs_diff = abs(x - y)
        denom = max(abs(x), 1e-12)
        rel_err = abs_diff / denom
        h4 = None if egrm_h4 is None else egrm_h4.get(key)
        rows.append(JoinedRow(
            sample_a=key[0], sample_b=key[1],
            egrm_ref=x, graphpop=y,
            abs_diff=abs_diff, rel_err=rel_err,
            egrm_h4=h4,
        ))
    return rows


def write_panel_csv(rows: Iterable[JoinedRow], path: Path) -> None:
    """Emit the tidy panel CSV for downstream figure-gen scripts."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "sample_a", "sample_b",
            "egrm_ref", "graphpop_branch_grm",
            "abs_diff", "rel_err", "egrm_h4",
        ])
        for row in rows:
            w.writerow(row.as_csv_row())


def summarise(rows: List[JoinedRow]) -> dict:
    """Return max/mean rel-err + abs-diff for stderr printing."""
    if not rows:
        return {"n_pairs": 0,
                "max_rel_err": 0.0, "max_abs_diff": 0.0,
                "mean_rel_err": 0.0, "mean_abs_diff": 0.0}
    rel = [r.rel_err for r in rows]
    abs_ = [r.abs_diff for r in rows]
    return {
        "n_pairs": len(rows),
        "max_rel_err": max(rel),
        "max_abs_diff": max(abs_),
        "mean_rel_err": sum(rel) / len(rel),
        "mean_abs_diff": sum(abs_) / len(abs_),
    }


# ---------------------------------------------------------------------------
# Integration steps — subprocess invocations
# ---------------------------------------------------------------------------

def run_java_dump(
    repo_root: Path, java_dump_dir: Path,
    *, mvn_binary: str = "mvn",
) -> None:
    """Invoke the Maven test that dumps GraphPop's branch_grm TSVs.

    Runs from `<repo_root>/graphpop-procedures` (no parent reactor
    pom exists at repo root, so `-pl` doesn't apply).
    """
    java_dump_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        mvn_binary,
        "test",
        "-Dtest=BranchGrmProcedureTest",
        f"-Dgraphpop.bench.dump.dir={java_dump_dir}",
        "-q",
    ]
    cwd = repo_root / "graphpop-procedures"
    print(f"[fig1de] invoking (cwd={cwd}): {' '.join(cmd)}", file=sys.stderr)
    r = subprocess.run(cmd, cwd=cwd)
    if r.returncode != 0:
        raise RuntimeError(
            f"mvn test exited with code {r.returncode}; see Maven "
            "output above for details")


def run_egrm_wrapper(
    fixture_trees: Path, egrm_out_dir: Path,
    *, graphpop_bench_binary: str = "graphpop-bench",
) -> None:
    """Invoke the H4 egrm wrapper on the fixture .trees file."""
    egrm_out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        graphpop_bench_binary, "run", "egrm",
        "--input", str(fixture_trees),
        "--output", str(egrm_out_dir),
    ]
    print(f"[fig1de] invoking: {' '.join(cmd)}", file=sys.stderr)
    r = subprocess.run(cmd)
    if r.returncode != 0:
        raise RuntimeError(
            f"graphpop-bench run egrm exited with code {r.returncode}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

DEFAULT_REPO_ROOT = Path("/mnt/data/GraphPop")
DEFAULT_FIXTURE = (
    DEFAULT_REPO_ROOT
    / "graphpop-procedures/src/test/resources/egrm_fixture_20samples.trees"
)
DEFAULT_REF_JSON_UNCONDITIONAL = (
    DEFAULT_REPO_ROOT
    / "graphpop-procedures/src/test/resources/egrm_expected_20samples.json"
)
DEFAULT_REF_JSON_PATHWAY = (
    DEFAULT_REPO_ROOT
    / "graphpop-procedures/src/test/resources/egrm_expected_pathway_half.json"
)


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--repo-root", type=Path, default=DEFAULT_REPO_ROOT)
    p.add_argument("--fixture-trees", type=Path, default=DEFAULT_FIXTURE,
                   help="Path to egrm_fixture_20samples.trees")
    p.add_argument("--ref-json-1d", type=Path,
                   default=DEFAULT_REF_JSON_UNCONDITIONAL)
    p.add_argument("--ref-json-1e", type=Path,
                   default=DEFAULT_REF_JSON_PATHWAY)
    p.add_argument("--java-dump-dir", type=Path, required=True,
                   help="Where the Java test writes branch_grm TSVs")
    p.add_argument("--egrm-out-dir", type=Path, required=True,
                   help="Where the H4 egrm wrapper writes outputs")
    p.add_argument("--output-dir", type=Path, required=True,
                   help="Where to write fig1d/fig1e panel CSVs")
    p.add_argument("--skip-mvn", action="store_true",
                   help="Skip the mvn invocation; reuse existing dumps")
    p.add_argument("--skip-egrm", action="store_true",
                   help="Skip the egrm wrapper invocation; "
                        "JSON reference is still used")
    p.add_argument("--mvn-binary", default="mvn")
    p.add_argument("--graphpop-bench-binary", default="graphpop-bench")
    p.add_argument("--max-rel-err-gate", type=float, default=1e-6,
                   help="Exit non-zero if any panel exceeds this rel-err")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)

    if not args.skip_mvn:
        run_java_dump(args.repo_root, args.java_dump_dir,
                      mvn_binary=args.mvn_binary)
    if not args.skip_egrm:
        run_egrm_wrapper(args.fixture_trees, args.egrm_out_dir,
                         graphpop_bench_binary=args.graphpop_bench_binary)

    # Fig 1d
    ref_ids_1d, ref_mat_1d = load_egrm_json_matrix(args.ref_json_1d)
    ref_pairs_1d = matrix_to_pair_map(ref_ids_1d, ref_mat_1d)
    gp_pairs_1d = load_pair_tsv(
        args.java_dump_dir / "branch_grm_unconditional.tsv")
    egrm_h4_pairs: Dict[PairKey, float] | None = None
    h4_tsv = args.egrm_out_dir / "egrm.tsv"
    if h4_tsv.exists():
        # The H4 wrapper writes pairs keyed by *node-ID string*; the
        # JSON reference uses 0-based integer sample IDs. Both happen
        # to coincide for `ts.samples()` on this fixture (haplotype
        # nodes 0..19), so the canonical-pair-keyed map indexes the
        # same way.
        egrm_h4_pairs = load_pair_tsv(h4_tsv)
    rows_1d = join_pairs(ref_pairs_1d, gp_pairs_1d, egrm_h4=egrm_h4_pairs)
    panel_1d = args.output_dir / "fig1d_panel_data.csv"
    write_panel_csv(rows_1d, panel_1d)
    summary_1d = summarise(rows_1d)
    print(f"[fig1d] {summary_1d} → {panel_1d}", file=sys.stderr)

    # Fig 1e
    ref_ids_1e, ref_mat_1e = load_egrm_json_matrix(args.ref_json_1e)
    ref_pairs_1e = matrix_to_pair_map(ref_ids_1e, ref_mat_1e)
    gp_pairs_1e = load_pair_tsv(
        args.java_dump_dir / "branch_grm_pathway_half.tsv")
    rows_1e = join_pairs(ref_pairs_1e, gp_pairs_1e)
    panel_1e = args.output_dir / "fig1e_panel_data.csv"
    write_panel_csv(rows_1e, panel_1e)
    summary_1e = summarise(rows_1e)
    print(f"[fig1e] {summary_1e} → {panel_1e}", file=sys.stderr)

    # Validation gate.
    worst = max(summary_1d["max_rel_err"], summary_1e["max_rel_err"])
    if worst >= args.max_rel_err_gate:
        print(
            f"[fig1de] FAIL: max rel-err {worst:.3e} ≥ gate "
            f"{args.max_rel_err_gate:.3e}", file=sys.stderr)
        return 2
    print(
        f"[fig1de] OK: max rel-err {worst:.3e} < gate "
        f"{args.max_rel_err_gate:.3e}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
