#!/usr/bin/env python3
"""
Benchmark all GraphPop procedures across three region sizes.

Metrics captured per test:
  - Wall-clock time (cold run + warm runs, median of warm)
  - Peak memory delta (GC before → raw snapshot after → delta)
  - Throughput (rows / second)

Output:
  - benchmarks/results/<label>.jsonl   (JSON Lines, one record per test)
  - Summary table to stdout

Usage:
    python scripts/benchmark-all.py [--label baseline] [--runs 3]
"""

import argparse
import json
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "benchmarks" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CYPHER_SHELL = "cypher-shell"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"

# Region definitions (chr22)
REGIONS = {
    "small":  {"chr": "chr22", "start": 16000000, "end": 16100000, "desc": "~1K variants"},
    "medium": {"chr": "chr22", "start": 16000000, "end": 16500000, "desc": "~5K variants"},
    "large":  {"chr": "chr22", "start": 16000000, "end": 17000000, "desc": "~10K variants"},
}

# Procedure definitions: name, cypher template, needs_write (if it writes data)
PROCEDURES = [
    {
        "name": "diversity",
        "cypher": "CALL graphpop.diversity('{chr}', {start}, {end}, '{pop}') YIELD pi RETURN count(*);",
        "pop": "EUR",
        "write": False,
    },
    {
        "name": "divergence",
        "cypher": "CALL graphpop.divergence('{chr}', {start}, {end}, '{pop}', '{pop2}') YIELD dxy RETURN count(*);",
        "pop": "EUR", "pop2": "AFR",
        "write": False,
    },
    {
        "name": "sfs",
        "cypher": "CALL graphpop.sfs('{chr}', {start}, {end}, '{pop}') YIELD n_variants RETURN count(*);",
        "pop": "EUR",
        "write": False,
    },
    {
        "name": "joint_sfs",
        "cypher": "CALL graphpop.joint_sfs('{chr}', {start}, {end}, '{pop}', '{pop2}') YIELD n_variants RETURN count(*);",
        "pop": "EUR", "pop2": "AFR",
        "write": False,
    },
    {
        "name": "genome_scan",
        "cypher": (
            "CALL graphpop.genome_scan('{chr}', '{pop}', 100000, 50000, "
            "{{pop2: '{pop2}', run_id: 'bench'}}) "
            "YIELD window_id RETURN count(*);"
        ),
        "pop": "EUR", "pop2": "AFR",
        "write": True,
    },
    {
        "name": "ld",
        "cypher": "CALL graphpop.ld('{chr}', {start}, {end}, '{pop}', 100000, 0.2) YIELD variant1 RETURN count(*);",
        "pop": "EUR",
        "write": True,
    },
    {
        "name": "ihs",
        "cypher": (
            "CALL graphpop.ihs('{chr}', '{pop}', "
            "{{start: {start}, end: {end}, min_af: 0.05}}) "
            "YIELD variantId RETURN count(*);"
        ),
        "pop": "EUR",
        "write": True,
    },
    {
        "name": "xpehh",
        "cypher": (
            "CALL graphpop.xpehh('{chr}', '{pop}', '{pop2}', "
            "{{start: {start}, end: {end}}}) "
            "YIELD variantId RETURN count(*);"
        ),
        "pop": "EUR", "pop2": "AFR",
        "write": True,
    },
]


def run_cypher(query: str, timeout: int = 600) -> tuple[str, float]:
    """Execute a Cypher query via cypher-shell and return (output, elapsed_seconds)."""
    cmd = [
        CYPHER_SHELL,
        "-u", NEO4J_USER,
        "-p", NEO4J_PASS,
        "--format", "plain",
        query,
    ]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    elapsed = time.perf_counter() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"cypher-shell failed: {proc.stderr.strip()}")
    return proc.stdout.strip(), elapsed


def parse_heap_output(out: str) -> dict:
    """Parse the plain-format heap stats output into a dict."""
    lines = [l for l in out.split("\n") if l.strip()]
    if len(lines) >= 2:
        vals = lines[-1].split(",")
        if len(vals) >= 4:
            return {
                "total_mb": int(vals[0].strip()),
                "used_mb": int(vals[1].strip()),
                "free_mb": int(vals[2].strip()),
                "max_mb": int(vals[3].strip()),
            }
    return {}


def force_gc() -> dict:
    """Force GC and return clean heap stats (baseline for memory delta)."""
    out, _ = run_cypher(
        "CALL graphpop.gc() YIELD total_mb, used_mb, free_mb, max_mb "
        "RETURN total_mb, used_mb, free_mb, max_mb;")
    return parse_heap_output(out)


def get_heap_stats() -> dict:
    """Get raw heap stats (no GC — captures transient allocations)."""
    out, _ = run_cypher(
        "CALL graphpop.heap_stats() YIELD total_mb, used_mb, free_mb, max_mb "
        "RETURN total_mb, used_mb, free_mb, max_mb;")
    return parse_heap_output(out)


def cleanup_bench_artifacts():
    """Remove any benchmark-created nodes/edges to keep state clean."""
    try:
        run_cypher("MATCH (w:GenomicWindow) WHERE w.run_id = 'bench' DETACH DELETE w;", timeout=60)
    except Exception:
        pass
    try:
        run_cypher("MATCH ()-[r:LD]->() DELETE r;", timeout=120)
    except Exception:
        pass


def format_cypher(proc: dict, region: dict) -> str:
    """Fill in the Cypher template with region and population parameters."""
    return proc["cypher"].format(
        chr=region["chr"],
        start=region["start"],
        end=region["end"],
        pop=proc.get("pop", "EUR"),
        pop2=proc.get("pop2", "AFR"),
    )


def run_benchmark(label: str, n_runs: int, skip_write_regions: list[str] | None = None):
    """Run the full benchmark suite.

    Run 0 is the cold run (first execution after cleanup, page cache cold).
    Runs 1..n are warm runs.  Median is computed from warm runs only.
    Peak memory delta = max(heap_after.used_mb - heap_before.used_mb) across runs.
    """
    outfile = RESULTS_DIR / f"{label}.jsonl"
    results = []

    total_tests = 0
    for proc in PROCEDURES:
        for rname in REGIONS:
            if skip_write_regions and proc["write"] and rname in skip_write_regions:
                continue
            total_tests += 1

    # n_runs is warm runs; we always do 1 cold + n_runs warm = n_runs+1 total
    total_runs = n_runs + 1

    print(f"\n{'='*72}")
    print(f" GraphPop Benchmark Suite — {label}")
    print(f" {total_tests} tests × {total_runs} runs (1 cold + {n_runs} warm)")
    print(f"{'='*72}\n")

    test_num = 0
    for proc in PROCEDURES:
        for rname, region in REGIONS.items():
            if skip_write_regions and proc["write"] and rname in skip_write_regions:
                continue

            test_num += 1
            query = format_cypher(proc, region)
            cold_time = None
            warm_timings = []
            row_counts = []
            mem_deltas = []  # heap used_mb delta per run

            print(f"[{test_num}/{total_tests}] {proc['name']:12s} {rname:8s} ", end="", flush=True)

            for run in range(total_runs):
                # Clean up write artifacts before each run
                if proc["write"]:
                    cleanup_bench_artifacts()

                # Force GC to get clean baseline
                heap_before = force_gc()
                try:
                    out, elapsed = run_cypher(query, timeout=600)

                    # Raw snapshot immediately after (captures transient allocations)
                    heap_after = get_heap_stats()

                    # Memory delta
                    if heap_before and heap_after:
                        delta = heap_after.get("used_mb", 0) - heap_before.get("used_mb", 0)
                        mem_deltas.append(max(delta, 0))

                    # Parse row count
                    lines = [l for l in out.split("\n") if l.strip()]
                    if lines:
                        try:
                            row_counts.append(int(lines[-1].strip().strip('"')))
                        except ValueError:
                            row_counts.append(0)

                    if run == 0:
                        cold_time = elapsed
                        print("c", end="", flush=True)
                    else:
                        warm_timings.append(elapsed)
                        print(".", end="", flush=True)

                except Exception as e:
                    print(f"E({e})", end="", flush=True)
                    if run == 0:
                        cold_time = float("nan")
                    else:
                        warm_timings.append(float("nan"))

            # Compute median of warm runs
            valid_warm = [t for t in warm_timings if t == t]
            median_warm = statistics.median(valid_warm) if valid_warm else float("nan")

            # Peak memory delta across all runs
            peak_mem_delta = max(mem_deltas) if mem_deltas else 0

            # Throughput: rows / warm-median seconds
            rows = row_counts[0] if row_counts else 0
            throughput = rows / median_warm if median_warm == median_warm and median_warm > 0 else 0.0

            record = {
                "label": label,
                "procedure": proc["name"],
                "region": rname,
                "region_desc": region["desc"],
                "chr": region["chr"],
                "start": region["start"],
                "end": region["end"],
                "population": proc.get("pop", "EUR"),
                "n_warm_runs": n_runs,
                "cold_s": cold_time,
                "warm_timings_s": warm_timings,
                "median_warm_s": median_warm,
                "row_count": rows,
                "throughput_rows_per_s": round(throughput, 1),
                "peak_mem_delta_mb": peak_mem_delta,
                "mem_deltas_mb": mem_deltas,
            }
            results.append(record)

            cold_str = f"{cold_time:.2f}" if cold_time and cold_time == cold_time else "?"
            print(f"  warm={median_warm:7.2f}s  cold={cold_str}s  "
                  f"mem=+{peak_mem_delta}MB  "
                  f"({rows} rows, {throughput:.0f}/s)")

            # Clean up after write procedures
            if proc["write"]:
                cleanup_bench_artifacts()

    # Write JSONL
    with open(outfile, "w") as f:
        for r in results:
            f.write(json.dumps(r) + "\n")
    print(f"\nResults written to {outfile}")

    # Print summary table
    print_summary(results)

    return results


def print_summary(results: list[dict]):
    """Print a formatted summary table."""
    print(f"\n{'Procedure':<14s} {'Region':<10s} {'Cold (s)':<10s} {'Warm (s)':<10s} "
          f"{'Rows':<8s} {'Rows/s':<10s} {'PeakΔMB':<10s}")
    print("-" * 75)
    for r in results:
        cold = r["cold_s"]
        warm = r["median_warm_s"]
        cold_str = f"{cold:.3f}" if cold and cold == cold else "FAIL"
        warm_str = f"{warm:.3f}" if warm == warm else "FAIL"
        tp_str = f"{r['throughput_rows_per_s']:.0f}" if r['throughput_rows_per_s'] > 0 else "—"
        print(f"{r['procedure']:<14s} {r['region']:<10s} {cold_str:<10s} {warm_str:<10s} "
              f"{r['row_count']:<8d} {tp_str:<10s} +{r['peak_mem_delta_mb']:<9d}")


def main():
    parser = argparse.ArgumentParser(description="GraphPop benchmark suite")
    parser.add_argument("--label", default="baseline", help="Label for this benchmark run")
    parser.add_argument("--runs", type=int, default=3, help="Number of warm runs per test")
    parser.add_argument("--skip-large-write", action="store_true",
                        help="Skip large region for write procedures (ld, ihs, xpehh, genome_scan)")
    args = parser.parse_args()

    skip = ["large"] if args.skip_large_write else None

    # Verify connectivity
    try:
        run_cypher("RETURN 1;", timeout=10)
    except Exception as e:
        print(f"ERROR: Cannot connect to Neo4j: {e}", file=sys.stderr)
        sys.exit(1)

    run_benchmark(args.label, args.runs, skip_write_regions=skip)


if __name__ == "__main__":
    main()
