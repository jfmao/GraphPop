#!/usr/bin/env python3
"""
Compare two GraphPop benchmark runs and produce a formatted table with:
  - Wall-clock speedup (warm median)
  - Cold vs warm ratio
  - Peak memory delta comparison
  - Throughput comparison

Usage:
    python scripts/benchmark-compare.py benchmarks/results/baseline.jsonl benchmarks/results/optimized.jsonl
"""

import argparse
import json
import sys
from pathlib import Path


def load_results(path: str) -> dict[tuple[str, str], dict]:
    """Load JSONL results into a dict keyed by (procedure, region)."""
    results = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)
            key = (r["procedure"], r["region"])
            results[key] = r
    return results


def safe_val(record, key, default=float("nan")):
    """Get a numeric value from a record, returning default if missing/NaN."""
    if record is None:
        return default
    v = record.get(key, default)
    if v is None:
        return default
    return v


def main():
    parser = argparse.ArgumentParser(description="Compare two benchmark runs")
    parser.add_argument("before", help="Path to baseline JSONL")
    parser.add_argument("after", help="Path to optimized JSONL")
    args = parser.parse_args()

    before = load_results(args.before)
    after = load_results(args.after)

    # Collect all keys in order
    all_keys = []
    for k in before:
        if k not in all_keys:
            all_keys.append(k)
    for k in after:
        if k not in all_keys:
            all_keys.append(k)

    region_order = {"small": 0, "medium": 1, "large": 2}
    all_keys.sort(key=lambda k: (k[0], region_order.get(k[1], 99)))

    # ---- Wall-clock comparison ----
    print(f"\n{'='*85}")
    print(f" WALL-CLOCK TIME (warm median)")
    print(f"{'='*85}")
    print(f"{'Procedure':<14s} {'Region':<10s} {'Before (s)':<12s} {'After (s)':<12s} "
          f"{'Speedup':<10s} {'Δ Time':<10s}")
    print("-" * 85)

    current_proc = None
    total_before = 0.0
    total_after = 0.0
    n_compared = 0

    for key in all_keys:
        proc, region = key
        if proc != current_proc:
            if current_proc is not None:
                print()
            current_proc = proc

        b = before.get(key)
        a = after.get(key)

        # Support both old format (median_s) and new format (median_warm_s)
        b_time = safe_val(b, "median_warm_s", safe_val(b, "median_s"))
        a_time = safe_val(a, "median_warm_s", safe_val(a, "median_s"))

        b_str = f"{b_time:.3f}" if b_time == b_time else "—"
        a_str = f"{a_time:.3f}" if a_time == a_time else "—"

        if b_time == b_time and a_time == a_time and a_time > 0:
            speedup = b_time / a_time
            delta = b_time - a_time
            speedup_str = f"{speedup:.2f}x"
            delta_str = f"-{delta:.3f}s" if delta > 0 else f"+{abs(delta):.3f}s"
            total_before += b_time
            total_after += a_time
            n_compared += 1
        else:
            speedup_str = "—"
            delta_str = "—"

        print(f"{proc:<14s} {region:<10s} {b_str:<12s} {a_str:<12s} "
              f"{speedup_str:<10s} {delta_str:<10s}")

    if n_compared > 0 and total_after > 0:
        overall_speedup = total_before / total_after
        print(f"\n{'-'*85}")
        print(f"{'TOTAL':<14s} {'':10s} {total_before:<12.3f} {total_after:<12.3f} "
              f"{overall_speedup:.2f}x      -{total_before - total_after:.3f}s")

    # ---- Peak memory comparison ----
    has_mem = any(
        before.get(k, {}).get("peak_mem_delta_mb") is not None or
        after.get(k, {}).get("peak_mem_delta_mb") is not None
        for k in all_keys
    )
    if has_mem:
        print(f"\n{'='*70}")
        print(f" PEAK MEMORY DELTA (MB above GC baseline)")
        print(f"{'='*70}")
        print(f"{'Procedure':<14s} {'Region':<10s} {'Before':<12s} {'After':<12s} {'Δ Memory':<12s}")
        print("-" * 60)

        current_proc = None
        for key in all_keys:
            proc, region = key
            if proc != current_proc:
                if current_proc is not None:
                    print()
                current_proc = proc

            b = before.get(key)
            a = after.get(key)

            b_mem = safe_val(b, "peak_mem_delta_mb", None)
            a_mem = safe_val(a, "peak_mem_delta_mb", None)

            b_str = f"+{int(b_mem)}MB" if b_mem is not None else "—"
            a_str = f"+{int(a_mem)}MB" if a_mem is not None else "—"

            if b_mem is not None and a_mem is not None:
                delta = a_mem - b_mem
                delta_str = f"{int(delta):+d}MB"
            else:
                delta_str = "—"

            print(f"{proc:<14s} {region:<10s} {b_str:<12s} {a_str:<12s} {delta_str:<12s}")

    # ---- Throughput comparison ----
    has_tp = any(
        before.get(k, {}).get("throughput_rows_per_s", 0) > 0 or
        after.get(k, {}).get("throughput_rows_per_s", 0) > 0
        for k in all_keys
    )
    if has_tp:
        print(f"\n{'='*70}")
        print(f" THROUGHPUT (rows/sec)")
        print(f"{'='*70}")
        print(f"{'Procedure':<14s} {'Region':<10s} {'Before':<12s} {'After':<12s} {'Speedup':<10s}")
        print("-" * 60)

        current_proc = None
        for key in all_keys:
            proc, region = key
            if proc != current_proc:
                if current_proc is not None:
                    print()
                current_proc = proc

            b = before.get(key)
            a = after.get(key)

            b_tp = safe_val(b, "throughput_rows_per_s", 0)
            a_tp = safe_val(a, "throughput_rows_per_s", 0)

            b_str = f"{b_tp:.0f}/s" if b_tp > 0 else "—"
            a_str = f"{a_tp:.0f}/s" if a_tp > 0 else "—"

            if b_tp > 0 and a_tp > 0:
                speedup_str = f"{a_tp / b_tp:.2f}x"
            else:
                speedup_str = "—"

            print(f"{proc:<14s} {region:<10s} {b_str:<12s} {a_str:<12s} {speedup_str:<10s}")

    # ---- Cold vs Warm ----
    has_cold = any(
        after.get(k, {}).get("cold_s") is not None
        for k in all_keys
    )
    if has_cold:
        print(f"\n{'='*60}")
        print(f" COLD vs WARM (cache effect, 'after' run only)")
        print(f"{'='*60}")
        print(f"{'Procedure':<14s} {'Region':<10s} {'Cold (s)':<12s} {'Warm (s)':<12s} {'Ratio':<10s}")
        print("-" * 60)

        current_proc = None
        for key in all_keys:
            proc, region = key
            if proc != current_proc:
                if current_proc is not None:
                    print()
                current_proc = proc

            a = after.get(key)
            if not a:
                continue

            cold = safe_val(a, "cold_s")
            warm = safe_val(a, "median_warm_s", safe_val(a, "median_s"))

            cold_str = f"{cold:.3f}" if cold == cold else "—"
            warm_str = f"{warm:.3f}" if warm == warm else "—"

            if cold == cold and warm == warm and warm > 0:
                ratio_str = f"{cold / warm:.2f}x"
            else:
                ratio_str = "—"

            print(f"{proc:<14s} {region:<10s} {cold_str:<12s} {warm_str:<12s} {ratio_str:<10s}")


if __name__ == "__main__":
    main()
