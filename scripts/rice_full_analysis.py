#!/usr/bin/env python3
"""Rice 3K Genome — Full Per-Population Analysis.

Runs comprehensive per-population statistics across all chromosomes (Phase 1)
and targeted XP-EHH for key population pairs (Phase 2).

Usage:
    python scripts/rice_full_analysis.py --phase 1
    python scripts/rice_full_analysis.py --phase 2
    python scripts/rice_full_analysis.py --phase all

Features:
    - Resume support: skips already-completed (pop, chr, proc) combos
    - Progress reporting with ETA
    - Error handling: logs failures and continues
    - JSON output: rice_full_results.json
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from neo4j import GraphDatabase

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"
NEO4J_DB = "neo4j"

OUTPUT_FILE = Path("rice_full_results.json")
LOG_FILE = Path("rice_full_analysis.log")

# Configure logging to file (unbuffered) + stderr
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%H:%M:%S",
    handlers=[
        logging.FileHandler(LOG_FILE, mode="a"),
        logging.StreamHandler(sys.stderr),
    ],
)
log = logging.getLogger("rice")

# 12 populations (skip 'na' with 1 sample)
POPULATIONS = [
    "XI-adm", "XI-3", "GJ-trp", "GJ-tmp", "XI-2", "XI-1A", "XI-1B",
    "cA-Aus", "GJ-sbtrp", "admix", "GJ-adm", "cB-Bas",
]

CHROMOSOMES = [
    "Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6",
    "Chr7", "Chr8", "Chr9", "Chr10", "Chr11", "Chr12",
]

# XP-EHH pairs ordered by ascending memory (smallest first)
XPEHH_PAIRS = [
    ("cA-Aus", "cB-Bas"),       # 277 samples, ~1.4 GB
    ("XI-1A", "cA-Aus"),        # 410 samples, ~2.1 GB
    ("GJ-tmp", "cA-Aus"),       # 488 samples, ~2.4 GB
    ("GJ-tmp", "XI-1A"),        # 496 samples, ~2.5 GB
    ("GJ-trp", "XI-1A"),        # 581 samples, ~2.9 GB
    ("GJ-tmp", "GJ-trp"),       # 659 samples, ~3.3 GB
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _set_output_file(path: Path) -> None:
    global OUTPUT_FILE
    OUTPUT_FILE = path.resolve()


def load_results() -> dict:
    """Load existing results file for resume support."""
    if OUTPUT_FILE.exists():
        with open(OUTPUT_FILE) as f:
            return json.load(f)
    return {"phase1": {}, "phase2": {}}


def save_results(results: dict) -> None:
    """Atomically save results to JSON."""
    tmp = OUTPUT_FILE.with_suffix(".json.tmp")
    with open(tmp, "w") as f:
        json.dump(results, f, indent=2)
    tmp.replace(OUTPUT_FILE)


def fmt_time(seconds: float) -> str:
    """Format seconds as HH:MM:SS."""
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    s = int(seconds % 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


def get_chr_lengths(session) -> dict[str, int]:
    """Get all chromosome lengths from DB."""
    recs = session.run(
        "MATCH (c:Chromosome) RETURN c.chromosomeId AS chr, c.length AS len"
    )
    return {r["chr"]: r["len"] for r in recs}


_RETRYABLE_ERRORS = (
    "MemoryPoolOutOfMemoryError",
    "DatabaseUnavailable",
    "connection refused",
    "SessionExpired",
    "Connection refused",
    "ConnectionResetError",
)
_MAX_RETRIES = 3
_BASE_DELAY = 5  # seconds


def run_cypher(driver, cypher: str) -> tuple[list[dict], float, str]:
    """Execute Cypher, return (records_as_dicts, wall_sec, error_str).

    Retries up to 3 times with exponential backoff for transient errors
    (memory exhaustion, connection issues, session expiry).
    Creates a fresh session for each attempt to handle SessionExpired errors.
    """
    t0 = time.time()
    last_error = ""
    for attempt in range(_MAX_RETRIES + 1):
        try:
            with driver.session(database=NEO4J_DB) as session:
                records = []
                for rec in session.run(cypher):
                    records.append({k: _serialize(rec[k]) for k in rec.keys()})
                return records, time.time() - t0, ""
        except Exception as e:
            last_error = str(e)
            if attempt < _MAX_RETRIES and any(
                err in last_error for err in _RETRYABLE_ERRORS
            ):
                delay = _BASE_DELAY * (2 ** attempt)
                log.warning(
                    "  Retryable error (attempt %d/%d, retry in %ds): %s",
                    attempt + 1, _MAX_RETRIES, delay, last_error[:100],
                )
                time.sleep(delay)
            else:
                return [], time.time() - t0, last_error
    return [], time.time() - t0, last_error


def _serialize(v):
    """Make a value JSON-serializable."""
    if isinstance(v, (int, float, str, bool, type(None))):
        return v
    if isinstance(v, (list, tuple)):
        return [_serialize(x) for x in v]
    return str(v)


# ---------------------------------------------------------------------------
# Phase 1: Per-Population Statistics
# ---------------------------------------------------------------------------

def phase1_procedures(chrom: str, pop: str, chr_len: int) -> list[tuple[str, str]]:
    """Return list of (proc_name, cypher) for Phase 1."""
    return [
        ("diversity", f"CALL graphpop.diversity('{chrom}', 1, {chr_len}, '{pop}')"),
        ("sfs", f"CALL graphpop.sfs('{chrom}', 1, {chr_len}, '{pop}', false)"),
        ("pop_summary", f"CALL graphpop.pop_summary('{chrom}', '{pop}')"),
        ("roh", f"CALL graphpop.roh('{chrom}', '{pop}')"),
        ("ihs", f"CALL graphpop.ihs('{chrom}', '{pop}')"),
        ("nsl", f"CALL graphpop.nsl('{chrom}', '{pop}')"),
        ("garud_h", f"CALL graphpop.garud_h('{chrom}', '{pop}', 100000, 50000)"),
    ]


def is_phase1_done(results: dict, pop: str, chrom: str, proc: str) -> bool:
    """Check if a Phase 1 result already exists (no error)."""
    return (
        pop in results.get("phase1", {})
        and chrom in results["phase1"][pop]
        and proc in results["phase1"][pop][chrom]
        and not results["phase1"][pop][chrom][proc].get("error")
    )


def _summarize_phase1(proc_name: str, records: list[dict]) -> dict:
    """Build summary entry for a Phase 1 procedure result."""
    if proc_name in ("diversity", "pop_summary") and records:
        return {"result": records[0]}
    if proc_name == "sfs" and records:
        return {"n_bins": len(records), "result": records}
    if proc_name == "roh" and records:
        total_len = sum(r.get("total_length", 0) for r in records)
        n_samples_with_roh = len(set(r.get("sampleId", "") for r in records))
        return {
            "n_segments": len(records),
            "total_roh_bp": total_len,
            "n_samples_with_roh": n_samples_with_roh,
        }
    if proc_name in ("ihs", "nsl") and records:
        return {"n_variants": len(records)}
    if proc_name == "garud_h" and records:
        return {"n_windows": len(records), "result": records}
    return {"n_records": len(records)}


def _status_str(entry: dict, wall_sec: float) -> str:
    """One-line status string for a completed procedure."""
    s = f"OK ({wall_sec:.1f}s"
    if "n_variants" in entry:
        s += f", {entry['n_variants']} variants"
    elif "n_windows" in entry:
        s += f", {entry['n_windows']} windows"
    elif "n_segments" in entry:
        s += f", {entry['n_segments']} ROH segments"
    return s + ")"


def _run_phase1_chr(
    driver, pop: str, chrom: str, chr_len: int, done_procs: set[str],
) -> list[tuple[str, dict]]:
    """Run all Phase 1 procedures for one pop×chr. Thread-safe."""
    out: list[tuple[str, dict]] = []
    procs = phase1_procedures(chrom, pop, chr_len)

    for proc_name, cypher in procs:
        if proc_name in done_procs:
            continue

        records, wall_sec, error = run_cypher(driver, cypher)

        entry: dict = {"wall_sec": round(wall_sec, 2)}
        if error:
            entry["error"] = error
            log.info(f"  {pop}/{chrom}/{proc_name}: ERROR {error[:80]}")
        else:
            entry.update(_summarize_phase1(proc_name, records))
            log.info(
                f"  {pop}/{chrom}/{proc_name}: "
                f"{_status_str(entry, wall_sec)}"
            )
        out.append((proc_name, entry))
    return out


def run_phase1(driver, results: dict, chr_lens: dict[str, int],
               workers: int = 2) -> None:
    """Run Phase 1: 7 procedures × 12 pops × 12 chromosomes.

    Parallelizes chromosomes within each population. Each thread gets its
    own Neo4j session and runs the 7 procedures sequentially for one
    chromosome. Populations run sequentially.
    """
    total = len(POPULATIONS) * len(CHROMOSOMES) * 7
    done = 0
    skipped = 0
    t_start = time.time()
    lock = threading.Lock()

    for pop in POPULATIONS:
        results.setdefault("phase1", {}).setdefault(pop, {})

        # Identify pending chromosomes and already-done procs
        pending: list[tuple[str, int, set]] = []  # (chrom, chr_len, done_procs)
        for chrom in CHROMOSOMES:
            results["phase1"][pop].setdefault(chrom, {})
            chr_len = chr_lens.get(chrom, 50_000_000)
            done_procs = set()
            for proc_name, _ in phase1_procedures(chrom, pop, chr_len):
                if is_phase1_done(results, pop, chrom, proc_name):
                    done_procs.add(proc_name)
                    done += 1
                    skipped += 1
            n_pending = 7 - len(done_procs)
            if n_pending > 0:
                pending.append((chrom, chr_len, done_procs))

        if not pending:
            log.info(f"  {pop}: all 12 chr done, skipping")
            continue

        n_pending_total = sum(7 - len(dp) for _, _, dp in pending)
        elapsed = time.time() - t_start
        log.info(
            f"  {pop}: {len(pending)} chr to run ({n_pending_total} procs, "
            f"{workers} workers, elapsed {fmt_time(elapsed)})"
        )

        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {}
            for chrom, chr_len, done_procs in pending:
                f = pool.submit(
                    _run_phase1_chr, driver, pop, chrom, chr_len, done_procs,
                )
                futures[f] = chrom

            for future in as_completed(futures):
                chrom = futures[future]
                try:
                    proc_results = future.result()
                except Exception as e:
                    log.info(f"  {pop}/{chrom}: THREAD ERROR {e}")
                    continue

                with lock:
                    for proc_name, entry in proc_results:
                        results["phase1"][pop][chrom][proc_name] = entry
                        done += 1
                    save_results(results)

                elapsed = time.time() - t_start
                completed_new = done - skipped
                if completed_new > 0:
                    eta = elapsed / completed_new * (total - done)
                    eta_str = fmt_time(eta)
                else:
                    eta_str = "..."
                log.info(
                    f"  [{done}/{total}] {pop}/{chrom} complete "
                    f"(elapsed {fmt_time(elapsed)}, ETA {eta_str})"
                )

    elapsed_total = time.time() - t_start
    log.info(
        f"Phase 1 complete: {done} calls ({skipped} skipped) "
        f"in {fmt_time(elapsed_total)}"
    )


# ---------------------------------------------------------------------------
# Phase 2: XP-EHH for Key Population Pairs
# ---------------------------------------------------------------------------


def pair_key(pop1: str, pop2: str) -> str:
    return f"{pop1}_vs_{pop2}"


def is_phase2_done(results: dict, pkey: str, chrom: str) -> bool:
    return (
        pkey in results.get("phase2", {})
        and chrom in results["phase2"][pkey]
        and not results["phase2"][pkey][chrom].get("error")
    )


def _run_xpehh_one(driver, pop1: str, pop2: str, chrom: str) -> dict:
    """Run a single XP-EHH call in its own session. Thread-safe."""
    cypher = f"CALL graphpop.xpehh('{chrom}', '{pop1}', '{pop2}')"
    records, wall_sec, error = run_cypher(driver, cypher)

    entry: dict = {"wall_sec": round(wall_sec, 2)}
    if error:
        entry["error"] = error
    else:
        entry["n_variants"] = len(records)
        top_hits = [
            r for r in records
            if abs(r.get("xpehh", 0) or 0) > 3
        ]
        entry["n_top_hits"] = len(top_hits)
        entry["top_hits"] = top_hits
    return entry


def run_phase2(driver, results: dict, workers: int = 3) -> None:
    """Run Phase 2: XP-EHH for 6 key pairs × 12 chromosomes.

    Parallelizes chromosomes within each pair (each thread gets its own
    Neo4j session). Pairs run sequentially to bound peak JVM heap usage.
    """
    total = len(XPEHH_PAIRS) * len(CHROMOSOMES)
    done = 0
    skipped = 0
    t_start = time.time()
    lock = threading.Lock()

    for pop1, pop2 in XPEHH_PAIRS:
        pkey = pair_key(pop1, pop2)
        results.setdefault("phase2", {}).setdefault(pkey, {})

        # Collect pending chromosomes for this pair
        pending_chroms = []
        for chrom in CHROMOSOMES:
            if is_phase2_done(results, pkey, chrom):
                done += 1
                skipped += 1
            else:
                pending_chroms.append(chrom)

        if not pending_chroms:
            log.info(f"  {pkey}: all 12 chr already done, skipping")
            continue

        log.info(
            f"  {pkey}: {len(pending_chroms)} chr to run "
            f"({workers} workers)"
        )

        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {}
            for chrom in pending_chroms:
                f = pool.submit(_run_xpehh_one, driver, pop1, pop2, chrom)
                futures[f] = chrom

            for future in as_completed(futures):
                chrom = futures[future]
                try:
                    entry = future.result()
                except Exception as e:
                    entry = {"wall_sec": 0, "error": str(e)}

                done += 1
                elapsed = time.time() - t_start
                completed_new = done - skipped
                if completed_new > 0:
                    eta = elapsed / completed_new * (total - done)
                    eta_str = fmt_time(eta)
                else:
                    eta_str = "..."

                if entry.get("error"):
                    log.info(
                        f"[{done}/{total}] {pkey}/{chrom}: "
                        f"ERROR {entry['error'][:80]} "
                        f"(elapsed {fmt_time(elapsed)}, ETA {eta_str})"
                    )
                else:
                    log.info(
                        f"[{done}/{total}] {pkey}/{chrom}: "
                        f"OK ({entry['wall_sec']:.1f}s, "
                        f"{entry['n_variants']} variants, "
                        f"{entry['n_top_hits']} |XP-EHH|>3) "
                        f"(elapsed {fmt_time(elapsed)}, ETA {eta_str})"
                    )

                with lock:
                    results["phase2"][pkey][chrom] = entry
                    save_results(results)

    elapsed_total = time.time() - t_start
    log.info(
        f"Phase 2 complete: {done} calls ({skipped} skipped) "
        f"in {fmt_time(elapsed_total)}"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rice 3K Genome — Full Per-Population Analysis"
    )
    parser.add_argument(
        "--phase", choices=["1", "2", "all"], default="all",
        help="Which phase to run (default: all)",
    )
    parser.add_argument(
        "--output", default=None,
        help=f"Output JSON file (default: {OUTPUT_FILE})",
    )
    parser.add_argument(
        "--workers", type=int, default=3,
        help="Concurrent chromosomes for Phase 2 XP-EHH (default: 3)",
    )
    args = parser.parse_args()

    if args.output is not None:
        _set_output_file(Path(args.output))

    results = load_results()

    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))

    with driver.session(database=NEO4J_DB) as session:
        chr_lens = get_chr_lengths(session)
        log.info(f"Loaded {len(chr_lens)} chromosome lengths")

    if args.phase in ("1", "all"):
        log.info(f"=== Phase 1: Per-Population Statistics ({args.workers} workers) ===")
        run_phase1(driver, results, chr_lens, workers=args.workers)

    if args.phase in ("2", "all"):
        log.info(f"=== Phase 2: XP-EHH Key Pairs ({args.workers} workers) ===")
        run_phase2(driver, results, workers=args.workers)

    driver.close()
    save_results(results)
    log.info(f"Results saved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
