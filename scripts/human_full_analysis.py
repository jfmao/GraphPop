#!/usr/bin/env python3
"""1000 Genomes Project — Full Per-Population Analysis.

Runs comprehensive population genetics statistics across all 22 autosomes
for 26 sub-populations (Phase 1), targeted XP-EHH for key pairs (Phase 2),
ancestral allele analyses (Phase 3), and population-specific analyses (Phase 4).

Usage:
    python scripts/human_full_analysis.py --phase 1
    python scripts/human_full_analysis.py --phase 2
    python scripts/human_full_analysis.py --phase 3
    python scripts/human_full_analysis.py --phase 4
    python scripts/human_full_analysis.py --phase all
    python scripts/human_full_analysis.py --phase 1 --workers 4

Features:
    - Resume support: skips already-completed (pop, chr, proc) combos
    - Progress reporting with ETA
    - Error handling: logs failures and continues
    - JSON output: human_full_results.json
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path

from neo4j import GraphDatabase
import os

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = os.environ.get("GRAPHPOP_USER", "neo4j")
NEO4J_PASS = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")
NEO4J_DB = "neo4j"

OUTPUT_FILE = Path("human_full_results.json")
LOG_FILE = Path("human_full_analysis.log")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%H:%M:%S",
    handlers=[
        logging.FileHandler(LOG_FILE, mode="a"),
        logging.StreamHandler(sys.stderr),
    ],
)
log = logging.getLogger("human")

# 26 sub-populations (1000G)
POPULATIONS = [
    # AFR (7)
    "YRI", "LWK", "GWD", "MSL", "ESN", "ACB", "ASW",
    # EUR (5)
    "CEU", "TSI", "FIN", "GBR", "IBS",
    # EAS (5)
    "CHB", "JPT", "CHS", "CDX", "KHV",
    # SAS (5)
    "GIH", "PJL", "BEB", "STU", "ITU",
    # AMR (4)
    "MXL", "PUR", "CLM", "PEL",
]

CHROMOSOMES = [f"chr{i}" for i in range(1, 23)]

# XP-EHH pairs (18 key pairs: cross-continental + within-continental)
XPEHH_PAIRS = [
    # Cross-continental (5)
    ("YRI", "CEU"),
    ("YRI", "CHB"),
    ("CEU", "CHB"),
    ("YRI", "GIH"),
    ("CEU", "JPT"),
    # Within AFR (3)
    ("YRI", "LWK"),
    ("YRI", "GWD"),
    ("ESN", "GWD"),
    # Within EUR (3)
    ("CEU", "TSI"),
    ("CEU", "FIN"),
    ("GBR", "IBS"),
    # Within EAS (3)
    ("CHB", "JPT"),
    ("CHB", "CHS"),
    ("CHB", "CDX"),
    # Within SAS (2)
    ("GIH", "BEB"),
    ("PJL", "STU"),
    # Within AMR (2)
    ("MXL", "PUR"),
    ("CLM", "PEL"),
]

# PBS configurations: (focal, sister, outgroup)
PBS_CONFIGS = [
    ("CEU", "GBR", "YRI"),      # Northern European-specific
    ("YRI", "LWK", "CEU"),      # West African-specific
    ("CHB", "JPT", "YRI"),      # Han Chinese-specific
    ("GIH", "BEB", "YRI"),      # Gujarati-specific
    ("PEL", "MXL", "CHB"),      # Peruvian-specific
    ("FIN", "CEU", "CHB"),      # Finnish-specific
]

# Joint SFS key pairs (same as XPEHH_PAIRS)
JSFS_PAIRS = XPEHH_PAIRS

# Genome scan: key populations
GSCAN_KEY_POPS = ["YRI", "CEU", "CHB", "GIH", "PEL", "FIN"]

WINDOW_SIZE = 100_000
WINDOW_STEP = 50_000

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def load_results() -> dict:
    if OUTPUT_FILE.exists():
        with open(OUTPUT_FILE) as f:
            return json.load(f)
    return {}


def save_results(results: dict) -> None:
    tmp = OUTPUT_FILE.with_suffix(".json.tmp")
    with open(tmp, "w") as f:
        json.dump(results, f, indent=2)
    tmp.replace(OUTPUT_FILE)


def fmt_time(seconds: float) -> str:
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    s = int(seconds % 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


def get_chr_lengths(session) -> dict[str, int]:
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
    if isinstance(v, (int, float, str, bool, type(None))):
        return v
    if isinstance(v, (list, tuple)):
        return [_serialize(x) for x in v]
    if isinstance(v, dict):
        return {k: _serialize(val) for k, val in v.items()}
    return str(v)


# ---------------------------------------------------------------------------
# Phase 1: Per-Population Statistics
# ---------------------------------------------------------------------------


IHS_REGION_SIZE = 50_000_000  # 50 Mb per region to avoid transaction memory OOM


def phase1_procedures(chrom: str, pop: str, chr_len: int) -> list[tuple[str, str]]:
    """Return list of (proc_name, cypher) for Phase 1.

    Note: 'ihs' is listed here for resume tracking but handled specially
    by _run_phase1_chr via _run_ihs_one() for per-region splitting.
    """
    return [
        ("diversity", f"CALL graphpop.diversity('{chrom}', 1, {chr_len}, '{pop}')"),
        ("sfs", f"CALL graphpop.sfs('{chrom}', 1, {chr_len}, '{pop}', false)"),
        ("pop_summary", f"CALL graphpop.pop_summary('{chrom}', '{pop}')"),
        ("roh", f"CALL graphpop.roh('{chrom}', '{pop}')"),
        ("ihs", None),  # handled by _run_ihs_one
        ("nsl", f"CALL graphpop.nsl('{chrom}', '{pop}')"),
        ("garud_h", f"CALL graphpop.garud_h('{chrom}', '{pop}', {WINDOW_SIZE}, {WINDOW_STEP})"),
    ]


def is_done(results: dict, section: str, *keys: str) -> bool:
    d = results.get(section, {})
    for k in keys:
        if not isinstance(d, dict):
            return False
        d = d.get(k, {})
    return bool(d) and not (isinstance(d, dict) and d.get("error"))


def _summarize_phase1(proc_name: str, records: list[dict]) -> dict:
    if proc_name in ("diversity", "pop_summary") and records:
        return {"result": records[0]}
    if proc_name == "sfs" and records:
        return {"n_bins": len(records), "result": records}
    if proc_name == "roh" and records:
        total_len = sum(r.get("total_length", 0) for r in records)
        n_with_roh = sum(1 for r in records if r.get("n_roh", 0) > 0)
        froh_values = [r.get("froh", 0) for r in records if r.get("n_roh", 0) > 0]
        return {
            "n_samples": len(records),
            "n_samples_with_roh": n_with_roh,
            "total_roh_bp": total_len,
            "mean_froh": sum(froh_values) / len(froh_values) if froh_values else 0.0,
        }
    if proc_name in ("ihs", "nsl") and records:
        return {"n_variants": len(records)}
    if proc_name == "garud_h" and records:
        return {"n_windows": len(records), "result": records}
    return {"n_records": len(records)}


def _run_ihs_one(driver, pop: str, chrom: str,
                  chr_len: int = 250_000_000) -> dict:
    """Run iHS for one population on one chromosome.

    Splits into ~50 Mb regions (each a separate Neo4j transaction) to avoid
    transaction memory exhaustion on large chromosomes.  Per-region calls use
    unstd_only mode; AF-bin standardization is done in Python afterwards.
    """
    t0 = time.time()
    prop_unstd = f"ihs_unstd_{pop}"
    prop_std = f"ihs_{pop}"
    AF_BINS = 20

    # Build region list
    regions = []
    pos = 1
    while pos <= chr_len:
        end = min(pos + IHS_REGION_SIZE - 1, chr_len)
        regions.append((pos, end))
        pos = end + 1

    # Phase 1: per-region unstandardized iHS (each in its own transaction)
    total_variants = 0
    errors = []
    for rstart, rend in regions:
        cypher = (
            f"CALL graphpop.ihs('{chrom}', '{pop}', "
            f"{{start: {rstart}, end: {rend}, unstd_only: true}})"
        )
        try:
            with driver.session(database=NEO4J_DB) as session:
                for rec in session.run(cypher):
                    if rec["variantId"] == "SUMMARY":
                        total_variants += rec["pos"]
            log.info("    %s %s region %d-%d: %d variants so far",
                     chrom, pop, rstart, rend, total_variants)
        except Exception as e:
            errors.append(f"region {rstart}-{rend}: {str(e)[:120]}")
            log.warning("    %s %s region %d-%d: ERROR %s",
                        chrom, pop, rstart, rend, str(e)[:80])

    entry: dict = {"wall_sec": round(time.time() - t0, 2)}
    if errors:
        entry["error"] = "; ".join(errors)
        entry["n_variants"] = total_variants
        return entry
    if total_variants == 0:
        entry.update(n_variants=0)
        return entry

    # Phase 2: AF-bin standardization in Python
    try:
        # Resolve the population's index in per-variant arrays (pop_ids order)
        with driver.session(database=NEO4J_DB) as session:
            pop_ids_rec = session.run(
                "MATCH (v:Variant) WHERE v.chr = $chr AND v.pop_ids IS NOT NULL "
                "RETURN v.pop_ids AS pids LIMIT 1",
                chr=chrom,
            ).single()
        pop_idx = 0
        if pop_ids_rec and pop_ids_rec["pids"] is not None:
            pids = list(pop_ids_rec["pids"])
            if pop in pids:
                pop_idx = pids.index(pop)

        # Read all unstandardized values + per-population AF
        with driver.session(database=NEO4J_DB) as session:
            rows = [
                (rec["vid"], rec["pos"], rec["unstd"], rec["af"])
                for rec in session.run(
                    f"MATCH (v:Variant) WHERE v.chr = '{chrom}' "
                    f"AND v.{prop_unstd} IS NOT NULL "
                    f"RETURN v.variantId AS vid, v.pos AS pos, "
                    f"v.{prop_unstd} AS unstd, "
                    f"v.af[{pop_idx}] AS af"
                )
            ]

        n = len(rows)
        if n == 0:
            entry.update(n_variants=0)
            return entry

        # AF-bin standardization (same as Java IHSProcedure)
        bins: dict[int, list[int]] = {}
        for i, (vid, pos, unstd, af) in enumerate(rows):
            b = min(int(af * AF_BINS), AF_BINS - 1) if af is not None else 0
            bins.setdefault(b, []).append(i)

        ihs_std = [0.0] * n
        for bin_indices in bins.values():
            vals = [rows[i][2] for i in bin_indices]
            mean = sum(vals) / len(vals)
            if len(vals) > 1:
                var_sum = sum((v - mean) ** 2 for v in vals)
                std = (var_sum / (len(vals) - 1)) ** 0.5
            else:
                std = 1.0
            for i in bin_indices:
                ihs_std[i] = (rows[i][2] - mean) / std if std > 0 else 0.0

        # Write standardized values back in batches
        BATCH = 10_000
        for i in range(0, n, BATCH):
            batch = [
                {"vid": rows[j][0], "ihs": ihs_std[j]}
                for j in range(i, min(i + BATCH, n))
            ]
            with driver.session(database=NEO4J_DB) as session:
                session.run(
                    "UNWIND $batch AS row "
                    "MATCH (v:Variant {variantId: row.vid}) "
                    f"SET v.{prop_std} = row.ihs",
                    batch=batch,
                )

        entry["n_variants"] = n
        log.info("    %s %s: standardized %d variants (%d AF bins)",
                 chrom, pop, n, len(bins))
    except Exception as e:
        entry["error"] = f"standardization: {e}"
        entry["n_variants"] = total_variants

    entry["wall_sec"] = round(time.time() - t0, 2)
    return entry


def _run_phase1_chr(
    driver, pop: str, chrom: str, chr_len: int, done_procs: set[str],
) -> list[tuple[str, dict]]:
    out: list[tuple[str, dict]] = []
    procs = phase1_procedures(chrom, pop, chr_len)
    for proc_name, cypher in procs:
        if proc_name in done_procs:
            continue

        # iHS: use per-region splitting to avoid transaction memory OOM
        if proc_name == "ihs":
            entry = _run_ihs_one(driver, pop, chrom, chr_len)
            if entry.get("error"):
                log.info("  %s/%s/ihs: ERROR %s", pop, chrom,
                         entry["error"][:80])
            else:
                log.info("  %s/%s/ihs: OK (%.1fs, %d variants)",
                         pop, chrom, entry["wall_sec"],
                         entry.get("n_variants", 0))
            out.append((proc_name, entry))
            continue

        records, wall_sec, error = run_cypher(driver, cypher)
        entry: dict = {"wall_sec": round(wall_sec, 2)}
        if error:
            entry["error"] = error
            log.info("  %s/%s/%s: ERROR %s", pop, chrom, proc_name, error[:80])
        else:
            entry.update(_summarize_phase1(proc_name, records))
            log.info("  %s/%s/%s: OK (%.1fs)", pop, chrom, proc_name, wall_sec)
        out.append((proc_name, entry))
    return out


def run_phase1(driver, results: dict, chr_lens: dict[str, int],
               workers: int = 2) -> None:
    total = len(POPULATIONS) * len(CHROMOSOMES) * 7
    done = 0
    skipped = 0
    t_start = time.time()
    lock = threading.Lock()

    for pop in POPULATIONS:
        results.setdefault("phase1", {}).setdefault(pop, {})
        pending: list[tuple[str, int, set]] = []
        for chrom in CHROMOSOMES:
            results["phase1"][pop].setdefault(chrom, {})
            chr_len = chr_lens.get(chrom, 250_000_000)
            done_procs = set()
            for proc_name, _ in phase1_procedures(chrom, pop, chr_len):
                if is_done(results, "phase1", pop, chrom, proc_name):
                    done_procs.add(proc_name)
                    done += 1
                    skipped += 1
            if len(done_procs) < 7:
                pending.append((chrom, chr_len, done_procs))

        if not pending:
            log.info("  %s: all %d chr done, skipping", pop, len(CHROMOSOMES))
            continue

        n_pending = sum(7 - len(dp) for _, _, dp in pending)
        elapsed = time.time() - t_start
        log.info(
            "  %s: %d chr to run (%d procs, %d workers, elapsed %s)",
            pop, len(pending), n_pending, workers, fmt_time(elapsed),
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
                    log.info("  %s/%s: THREAD ERROR %s", pop, chrom, e)
                    continue

                with lock:
                    for proc_name, entry in proc_results:
                        results["phase1"][pop][chrom][proc_name] = entry
                        done += 1
                    save_results(results)

                elapsed = time.time() - t_start
                completed_new = done - skipped
                if completed_new > 0:
                    eta_str = fmt_time(elapsed / completed_new * (total - done))
                else:
                    eta_str = "..."
                log.info(
                    "  [%d/%d] %s/%s complete (elapsed %s, ETA %s)",
                    done, total, pop, chrom, fmt_time(elapsed), eta_str,
                )

    log.info(
        "Phase 1 complete: %d calls (%d skipped) in %s",
        done, skipped, fmt_time(time.time() - t_start),
    )


# ---------------------------------------------------------------------------
# Phase 2: Pairwise Statistics (XP-EHH + Divergence)
# ---------------------------------------------------------------------------


def pair_key(pop1: str, pop2: str) -> str:
    return f"{pop1}_vs_{pop2}"


XPEHH_REGION_SIZE = 50_000_000  # 50 Mb per region to avoid transaction memory OOM


def _run_xpehh_one(driver, pop1: str, pop2: str, chrom: str,
                    chr_len: int = 250_000_000) -> dict:
    """Run XP-EHH for one pair on one chromosome.

    Splits into ~50 Mb regions (each a separate Neo4j transaction) to avoid
    transaction memory exhaustion on large chromosomes.  Per-region calls use
    unstd_only mode; global standardization is done in Python afterwards.
    """
    t0 = time.time()
    prop_unstd = f"xpehh_unstd_{pop1}_{pop2}"
    prop_std = f"xpehh_{pop1}_{pop2}"

    # Build region list
    regions = []
    pos = 1
    while pos <= chr_len:
        end = min(pos + XPEHH_REGION_SIZE - 1, chr_len)
        regions.append((pos, end))
        pos = end + 1

    # Phase 1: per-region unstandardized XP-EHH (each in its own transaction)
    total_variants = 0
    errors = []
    for rstart, rend in regions:
        cypher = (
            f"CALL graphpop.xpehh('{chrom}', '{pop1}', '{pop2}', "
            f"{{start: {rstart}, end: {rend}, unstd_only: true}})"
        )
        try:
            with driver.session(database=NEO4J_DB) as session:
                for rec in session.run(cypher):
                    if rec["variantId"] == "SUMMARY":
                        total_variants += rec["pos"]
            log.info("    %s %s/%s region %d-%d: %d variants so far",
                     chrom, pop1, pop2, rstart, rend, total_variants)
        except Exception as e:
            errors.append(f"region {rstart}-{rend}: {str(e)[:120]}")
            log.warning("    %s %s/%s region %d-%d: ERROR %s",
                        chrom, pop1, pop2, rstart, rend, str(e)[:80])

    entry: dict = {"wall_sec": round(time.time() - t0, 2)}
    if errors:
        entry["error"] = "; ".join(errors)
        entry["n_variants"] = total_variants
        return entry
    if total_variants == 0:
        entry.update(n_variants=0, n_top_hits=0, top_hits=[])
        return entry

    # Phase 2: global standardization in Python
    try:
        # Read all unstandardized values (lightweight read — no packed arrays)
        with driver.session(database=NEO4J_DB) as session:
            rows = [
                (rec["vid"], rec["pos"], rec["unstd"])
                for rec in session.run(
                    f"MATCH (v:Variant) WHERE v.chr = '{chrom}' "
                    f"AND v.{prop_unstd} IS NOT NULL "
                    f"RETURN v.variantId AS vid, v.pos AS pos, "
                    f"v.{prop_unstd} AS unstd"
                )
            ]

        n = len(rows)
        if n == 0:
            entry.update(n_variants=0, n_top_hits=0, top_hits=[])
            return entry

        # Compute global mean/std
        vals = [r[2] for r in rows]
        mean = sum(vals) / n
        var_sum = sum((v - mean) ** 2 for v in vals)
        std = (var_sum / (n - 1)) ** 0.5 if n > 1 else 1.0

        # Write standardized values back in batches
        top_hits = []
        BATCH = 10_000
        for i in range(0, n, BATCH):
            batch = []
            for vid, pos, unstd in rows[i:i + BATCH]:
                xpehh = (unstd - mean) / std if std > 0 else 0.0
                batch.append({"vid": vid, "xpehh": xpehh})
                if abs(xpehh) > 3:
                    top_hits.append({"variantId": vid, "pos": pos,
                                     "xpehh": round(xpehh, 4)})
            with driver.session(database=NEO4J_DB) as session:
                session.run(
                    "UNWIND $batch AS row "
                    "MATCH (v:Variant {variantId: row.vid}) "
                    f"SET v.{prop_std} = row.xpehh",
                    batch=batch,
                )

        entry["n_variants"] = n
        entry["n_top_hits"] = len(top_hits)
        entry["top_hits"] = sorted(
            top_hits, key=lambda x: abs(x["xpehh"]), reverse=True
        )[:100]
        log.info("    %s %s/%s: standardized %d variants (mean=%.4f, std=%.4f, %d |XP-EHH|>3)",
                 chrom, pop1, pop2, n, mean, std, len(top_hits))
    except Exception as e:
        entry["error"] = f"standardization: {e}"
        entry["n_variants"] = total_variants

    entry["wall_sec"] = round(time.time() - t0, 2)
    return entry


def run_phase2(driver, results: dict, chr_lens: dict[str, int],
               workers: int = 3) -> None:
    """Phase 2A: XP-EHH for 18 key pairs × 22 chromosomes."""
    total = len(XPEHH_PAIRS) * len(CHROMOSOMES)
    done = 0
    skipped = 0
    t_start = time.time()
    lock = threading.Lock()

    for pop1, pop2 in XPEHH_PAIRS:
        pkey = pair_key(pop1, pop2)
        results.setdefault("phase2_xpehh", {}).setdefault(pkey, {})

        pending_chroms = []
        for chrom in CHROMOSOMES:
            if is_done(results, "phase2_xpehh", pkey, chrom):
                done += 1
                skipped += 1
            else:
                pending_chroms.append(chrom)

        if not pending_chroms:
            log.info("  %s: all %d chr done, skipping", pkey, len(CHROMOSOMES))
            continue

        log.info("  %s: %d chr to run (%d workers)", pkey, len(pending_chroms), workers)

        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {}
            for chrom in pending_chroms:
                f = pool.submit(_run_xpehh_one, driver, pop1, pop2, chrom,
                                chr_lens.get(chrom, 250_000_000))
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
                eta_str = fmt_time(elapsed / completed_new * (total - done)) if completed_new > 0 else "..."

                if entry.get("error"):
                    log.info(
                        "[%d/%d] %s/%s: ERROR %s (elapsed %s, ETA %s)",
                        done, total, pkey, chrom, entry["error"][:80],
                        fmt_time(elapsed), eta_str,
                    )
                else:
                    log.info(
                        "[%d/%d] %s/%s: OK (%.1fs, %d variants, %d |XP-EHH|>3) "
                        "(elapsed %s, ETA %s)",
                        done, total, pkey, chrom, entry["wall_sec"],
                        entry["n_variants"], entry["n_top_hits"],
                        fmt_time(elapsed), eta_str,
                    )

                with lock:
                    results["phase2_xpehh"][pkey][chrom] = entry
                    save_results(results)

    log.info("Phase 2A (XP-EHH) complete: %d calls (%d skipped) in %s",
             done, skipped, fmt_time(time.time() - t_start))

    # Phase 2B: Pairwise divergence (all 325 pairs)
    run_divergence(driver, results, chr_lens, workers)


def run_divergence(driver, results: dict, chr_lens: dict[str, int],
                   workers: int = 6) -> None:
    """Phase 2B: Pairwise divergence for all 325 pairs × 22 chromosomes."""
    all_pairs = list(combinations(POPULATIONS, 2))
    total = len(all_pairs) * len(CHROMOSOMES)
    done = 0
    skipped = 0
    t_start = time.time()
    lock = threading.Lock()

    results.setdefault("phase2_divergence", {})

    def run_one(pop1, pop2, chrom):
        key = f"{pop1}|{pop2}|{chrom}"
        t0 = time.time()
        chr_len = chr_lens.get(chrom, 250_000_000)
        cypher = (
            f"CALL graphpop.divergence('{chrom}', 1, {chr_len}, "
            f"'{pop1}', '{pop2}')"
        )
        try:
            with driver.session(database=NEO4J_DB) as s:
                rec = s.run(cypher).single()
                result = {k: _serialize(rec[k]) for k in rec.keys()} if rec else {}
            return key, {"result": result, "wall_sec": time.time() - t0}
        except Exception as e:
            return key, {"error": str(e), "wall_sec": time.time() - t0}

    existing = results["phase2_divergence"]

    for i, (pop1, pop2) in enumerate(all_pairs):
        tasks = []
        for chrom in CHROMOSOMES:
            key = f"{pop1}|{pop2}|{chrom}"
            if key in existing and "error" not in existing[key]:
                done += 1
                skipped += 1
            else:
                tasks.append((pop1, pop2, chrom))

        if not tasks:
            continue

        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(run_one, *t): t for t in tasks}
            for fut in as_completed(futures):
                key, val = fut.result()
                with lock:
                    existing[key] = val
                    done += 1

        if (i + 1) % 10 == 0 or i == len(all_pairs) - 1:
            elapsed = time.time() - t_start
            log.info(
                "  [%d/%d pairs] %s vs %s: done (%d/%d total, elapsed %s)",
                i + 1, len(all_pairs), pop1, pop2, done, total, fmt_time(elapsed),
            )
            save_results(results)

    log.info("Phase 2B (divergence) complete: %d calls (%d skipped) in %s",
             done, skipped, fmt_time(time.time() - t_start))


# ---------------------------------------------------------------------------
# Phase 3: Ancestral Allele Analyses
# ---------------------------------------------------------------------------


def run_phase3(driver, results: dict, chr_lens: dict[str, int],
               workers: int = 4) -> None:
    """Phase 3: Fay & Wu's H (diversity re-run), polarized SFS, ROH HMM, H genome scan."""
    lock = threading.Lock()

    # 3A: Fay & Wu's H (diversity already includes it if ancestral alleles present)
    log.info("=== Phase 3A: Fay & Wu's H (diversity) ===")
    results.setdefault("phase3_fay_wu", {})
    existing = results["phase3_fay_wu"]
    n_total = len(POPULATIONS) * len(CHROMOSOMES)
    n_done = sum(1 for k in existing if "error" not in existing[k])

    def run_diversity(pop, chrom):
        key = f"{pop}|{chrom}"
        if key in existing and "error" not in existing[key]:
            return key, existing[key]
        t0 = time.time()
        chr_len = chr_lens.get(chrom, 250_000_000)
        cypher = f"CALL graphpop.diversity('{chrom}', 1, {chr_len}, '{pop}')"
        try:
            with driver.session(database=NEO4J_DB) as s:
                rec = s.run(cypher).single()
                result = {k: _serialize(rec[k]) for k in rec.keys()} if rec else {}
            return key, {"result": result, "wall_sec": time.time() - t0}
        except Exception as e:
            return key, {"error": str(e), "wall_sec": time.time() - t0}

    for pop in POPULATIONS:
        tasks = [(pop, c) for c in CHROMOSOMES
                 if f"{pop}|{c}" not in existing or "error" in existing.get(f"{pop}|{c}", {})]
        if not tasks:
            continue
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(run_diversity, *t): t for t in tasks}
            for fut in as_completed(futures):
                key, val = fut.result()
                with lock:
                    existing[key] = val
                    n_done += 1
        log.info("  %s: done (%d/%d)", pop, n_done, n_total)
        save_results(results)

    # 3B: Polarized (unfolded) SFS
    log.info("=== Phase 3B: Polarized SFS ===")
    results.setdefault("phase3_usfs", {})
    existing_usfs = results["phase3_usfs"]
    n_total = len(POPULATIONS) * len(CHROMOSOMES)
    n_done = sum(1 for k in existing_usfs if "error" not in existing_usfs[k])

    def run_sfs(pop, chrom):
        key = f"{pop}|{chrom}"
        if key in existing_usfs and "error" not in existing_usfs[key]:
            return key, existing_usfs[key]
        t0 = time.time()
        chr_len = chr_lens.get(chrom, 250_000_000)
        cypher = f"CALL graphpop.sfs('{chrom}', 1, {chr_len}, '{pop}', false)"
        try:
            with driver.session(database=NEO4J_DB) as s:
                rec = s.run(cypher).single()
                result = {k: _serialize(rec[k]) for k in rec.keys()} if rec else {}
            return key, {"result": result, "wall_sec": time.time() - t0}
        except Exception as e:
            return key, {"error": str(e), "wall_sec": time.time() - t0}

    for pop in POPULATIONS:
        tasks = [(pop, c) for c in CHROMOSOMES
                 if f"{pop}|{c}" not in existing_usfs or "error" in existing_usfs.get(f"{pop}|{c}", {})]
        if not tasks:
            continue
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(run_sfs, *t): t for t in tasks}
            for fut in as_completed(futures):
                key, val = fut.result()
                with lock:
                    existing_usfs[key] = val
                    n_done += 1
        log.info("  %s: done (%d/%d)", pop, n_done, n_total)
        save_results(results)

    # 3C: ROH (HMM method)
    log.info("=== Phase 3C: ROH HMM ===")
    results.setdefault("phase3_roh_hmm", {})
    existing_roh = results["phase3_roh_hmm"]
    n_total = len(POPULATIONS) * len(CHROMOSOMES)
    n_done = sum(1 for k in existing_roh if "error" not in existing_roh[k])

    def run_roh(pop, chrom):
        key = f"{pop}|{chrom}"
        if key in existing_roh and "error" not in existing_roh[key]:
            return key, existing_roh[key]
        t0 = time.time()
        cypher = f"CALL graphpop.roh('{chrom}', '{pop}')"
        try:
            with driver.session(database=NEO4J_DB) as s:
                records = s.run(cypher).data()
            total_len = sum(r.get("total_length", 0) for r in records)
            froh_values = [r.get("froh", 0) for r in records if r.get("n_roh", 0) > 0]
            result = {
                "n_samples": len(records),
                "n_samples_with_roh": sum(1 for r in records if r.get("n_roh", 0) > 0),
                "total_roh_bp": total_len,
                "mean_froh": sum(froh_values) / len(froh_values) if froh_values else 0.0,
                "max_froh": max(froh_values) if froh_values else 0.0,
                "per_sample": [
                    {"sampleId": r["sampleId"], "n_roh": r["n_roh"],
                     "total_length": r["total_length"], "froh": r["froh"],
                     "mean_length": r["mean_length"], "max_length": r["max_length"]}
                    for r in records if r.get("n_roh", 0) > 0
                ][:50],
            }
            return key, {**result, "wall_sec": time.time() - t0}
        except Exception as e:
            return key, {"error": str(e), "wall_sec": time.time() - t0}

    for pop in POPULATIONS:
        tasks = [(pop, c) for c in CHROMOSOMES
                 if f"{pop}|{c}" not in existing_roh or "error" in existing_roh.get(f"{pop}|{c}", {})]
        if not tasks:
            log.info("  %s: all cached", pop)
            continue
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(run_roh, *t): t for t in tasks}
            for fut in as_completed(futures):
                key, val = fut.result()
                with lock:
                    existing_roh[key] = val
                    n_done += 1
        log.info("  %s: done (%d/%d)", pop, n_done, n_total)
        save_results(results)

    # 3D: Genome scan for Fay & Wu's H (6 key pops)
    log.info("=== Phase 3D: Fay & Wu's H genome scan ===")
    results.setdefault("phase3_hwscan", {})
    existing_hw = results["phase3_hwscan"]
    n_total = len(GSCAN_KEY_POPS) * len(CHROMOSOMES)
    n_done = sum(1 for k in existing_hw if "error" not in existing_hw[k])

    # genome_scan is WRITE mode — run sequentially
    for pop in GSCAN_KEY_POPS:
        for chrom in CHROMOSOMES:
            key = f"{pop}|{chrom}"
            if key in existing_hw and "error" not in existing_hw[key]:
                n_done += 1
                continue
            t0 = time.time()
            cypher = (
                f"CALL graphpop.genome_scan('{chrom}', '{pop}', "
                f"{WINDOW_SIZE}, {WINDOW_STEP}, {{}})"
            )
            try:
                with driver.session(database=NEO4J_DB) as s:
                    records = s.run(cypher).data()
                    windows = [{k: _serialize(v) for k, v in r.items()} for r in records]
                existing_hw[key] = {"result": windows, "wall_sec": time.time() - t0}
            except Exception as e:
                existing_hw[key] = {"error": str(e), "wall_sec": time.time() - t0}
            n_done += 1
            log.info("  [%d/%d] %s/%s: %d windows (%.1fs)",
                     n_done, n_total, pop, chrom,
                     len(existing_hw[key].get("result", [])),
                     existing_hw[key].get("wall_sec", 0))
            if n_done % 22 == 0:
                save_results(results)

    save_results(results)
    log.info("Phase 3 complete")


# ---------------------------------------------------------------------------
# Phase 4: Population-Specific Analyses
# ---------------------------------------------------------------------------


def run_phase4(driver, results: dict, chr_lens: dict[str, int],
               workers: int = 4) -> None:
    """Phase 4: PBS genome scans + Joint SFS for key pairs."""
    lock = threading.Lock()

    # 4A: PBS genome scans
    log.info("=== Phase 4A: PBS genome scans ===")
    results.setdefault("phase4_pbs", {})
    existing_pbs = results["phase4_pbs"]
    n_total = len(PBS_CONFIGS) * len(CHROMOSOMES)
    n_done = sum(1 for k in existing_pbs if "error" not in existing_pbs[k])

    # genome_scan is WRITE mode — run sequentially
    for focal, sister, outgroup in PBS_CONFIGS:
        for chrom in CHROMOSOMES:
            key = f"{focal}|{sister}|{outgroup}|{chrom}"
            if key in existing_pbs and "error" not in existing_pbs[key]:
                n_done += 1
                continue
            t0 = time.time()
            cypher = (
                f"CALL graphpop.genome_scan('{chrom}', '{focal}', "
                f"{WINDOW_SIZE}, {WINDOW_STEP}, "
                f"{{pop2: '{sister}', pop3: '{outgroup}'}})"
            )
            try:
                with driver.session(database=NEO4J_DB) as s:
                    records = s.run(cypher).data()
                    windows = [{k: _serialize(v) for k, v in r.items()} for r in records]
                existing_pbs[key] = {"result": windows, "wall_sec": time.time() - t0}
            except Exception as e:
                existing_pbs[key] = {"error": str(e), "wall_sec": time.time() - t0}
            n_done += 1
            log.info("  [%d/%d] PBS %s (vs %s, out %s) %s: %d windows (%.1fs)",
                     n_done, n_total, focal, sister, outgroup, chrom,
                     len(existing_pbs[key].get("result", [])),
                     existing_pbs[key].get("wall_sec", 0))
            if n_done % 22 == 0:
                save_results(results)

    save_results(results)

    # 4B: Joint SFS for key pairs
    log.info("=== Phase 4B: Joint SFS ===")
    results.setdefault("phase4_jsfs", {})
    existing_jsfs = results["phase4_jsfs"]
    n_total = len(JSFS_PAIRS) * len(CHROMOSOMES)
    n_done = sum(1 for k in existing_jsfs if "error" not in existing_jsfs[k])

    def run_jsfs(pop1, pop2, chrom):
        key = f"{pop1}|{pop2}|{chrom}"
        if key in existing_jsfs and "error" not in existing_jsfs[key]:
            return key, existing_jsfs[key]
        t0 = time.time()
        chr_len = chr_lens.get(chrom, 250_000_000)
        cypher = (
            f"CALL graphpop.joint_sfs('{chrom}', 1, {chr_len}, "
            f"'{pop1}', '{pop2}')"
        )
        try:
            with driver.session(database=NEO4J_DB) as s:
                records = []
                for rec in s.run(cypher):
                    records.append({k: _serialize(rec[k]) for k in rec.keys()})
            return key, {"n_cells": len(records), "wall_sec": time.time() - t0, "result": records}
        except Exception as e:
            return key, {"error": str(e), "wall_sec": time.time() - t0}

    for pop1, pop2 in JSFS_PAIRS:
        tasks = [(pop1, pop2, c) for c in CHROMOSOMES
                 if f"{pop1}|{pop2}|{c}" not in existing_jsfs
                 or "error" in existing_jsfs.get(f"{pop1}|{pop2}|{c}", {})]
        if not tasks:
            continue
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(run_jsfs, *t): t for t in tasks}
            for fut in as_completed(futures):
                key, val = fut.result()
                with lock:
                    existing_jsfs[key] = val
                    n_done += 1
        log.info("  %s vs %s: done (%d/%d)", pop1, pop2, n_done, n_total)
        save_results(results)

    log.info("Phase 4 complete")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def _set_output_file(path: Path) -> None:
    global OUTPUT_FILE
    OUTPUT_FILE = path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="1000 Genomes Project — Full Per-Population Analysis"
    )
    parser.add_argument(
        "--phase", choices=["1", "2", "3", "4", "all"], default="all",
        help="Which phase to run (default: all)",
    )
    parser.add_argument(
        "--output", default=None,
        help="Output JSON file (default: human_full_results.json)",
    )
    parser.add_argument(
        "--workers", type=int, default=3,
        help="Concurrent chromosomes per population (default: 3)",
    )
    args = parser.parse_args()

    if args.output is not None:
        _set_output_file(Path(args.output))

    results = load_results()
    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))

    with driver.session(database=NEO4J_DB) as session:
        chr_lens = get_chr_lengths(session)
        log.info("Loaded %d chromosome lengths", len(chr_lens))

    if args.phase in ("1", "all"):
        log.info("=== Phase 1: Per-Population Statistics (%d workers) ===", args.workers)
        run_phase1(driver, results, chr_lens, workers=args.workers)

    if args.phase in ("2", "all"):
        log.info("=== Phase 2: Pairwise Statistics (%d workers) ===", args.workers)
        run_phase2(driver, results, chr_lens, workers=args.workers)

    if args.phase in ("3", "all"):
        log.info("=== Phase 3: Ancestral Allele Analyses (%d workers) ===", args.workers)
        run_phase3(driver, results, chr_lens, workers=args.workers)

    if args.phase in ("4", "all"):
        log.info("=== Phase 4: Population-Specific Analyses (%d workers) ===", args.workers)
        run_phase4(driver, results, chr_lens, workers=args.workers)

    driver.close()
    save_results(results)
    log.info("All results saved to %s", OUTPUT_FILE)


if __name__ == "__main__":
    main()
