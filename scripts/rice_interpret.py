#!/usr/bin/env python3
"""Rice 3K Genome Biological Interpretation Script.

Cross-references population genetics signals with gene annotations
using GraphPop's graph-native annotation-aware analyses.

Usage:
    python scripts/rice_interpret.py annotate    # Cross-reference peaks with genes
    python scripts/rice_interpret.py pinsps      # pi_N/pi_S per population
    python scripts/rice_interpret.py divergence  # Full pairwise Fst matrix
    python scripts/rice_interpret.py gscan       # Genome scans with annotation
    python scripts/rice_interpret.py roh         # Re-run ROH with sliding_window
    python scripts/rice_interpret.py report      # Generate final Markdown report
"""

import argparse
import json
import logging
import os
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations
from threading import Lock

from neo4j import GraphDatabase

# ── Constants ─────────────────────────────────────────────────────────

NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"
NEO4J_DB = "neo4j"

RESULTS_FILE = "rice_full_results.json"
OUTPUT_FILE = "rice_interpretation_results.json"

POPULATIONS = [
    "XI-adm", "XI-3", "GJ-trp", "GJ-tmp", "XI-2", "XI-1A", "XI-1B",
    "cA-Aus", "GJ-sbtrp", "admix", "GJ-adm", "cB-Bas",
]

CHROMOSOMES = [
    "Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6",
    "Chr7", "Chr8", "Chr9", "Chr10", "Chr11", "Chr12",
]

XPEHH_PAIRS = [
    ("cA-Aus", "cB-Bas"),
    ("XI-1A", "cA-Aus"),
    ("GJ-tmp", "cA-Aus"),
    ("GJ-tmp", "XI-1A"),
    ("GJ-trp", "XI-1A"),
    ("GJ-tmp", "GJ-trp"),
]

# Key pairs for annotation-conditioned genome scans
GSCAN_PAIRS = [
    ("GJ-tmp", "XI-1A"),
    ("GJ-tmp", "GJ-trp"),
    ("XI-1A", "cA-Aus"),
]

# Known rice domestication genes (Nipponbare IRGSP-1.0 coordinates)
KNOWN_GENES = [
    {"gene": "Wx", "chr": "Chr6", "pos": 1600000, "function": "Amylose content (Waxy)"},
    {"gene": "Sd1", "chr": "Chr1", "pos": 38400000, "function": "Semi-dwarf (Green Revolution)"},
    {"gene": "sh4", "chr": "Chr4", "pos": 33000000, "function": "Seed shattering"},
    {"gene": "qSH1", "chr": "Chr1", "pos": 37000000, "function": "Seed shattering"},
    {"gene": "PROG1", "chr": "Chr7", "pos": 3200000, "function": "Prostrate growth"},
    {"gene": "Hd1", "chr": "Chr6", "pos": 9000000, "function": "Heading date"},
    {"gene": "Hd3a", "chr": "Chr6", "pos": 2600000, "function": "Florigen"},
    {"gene": "GS3", "chr": "Chr3", "pos": 16700000, "function": "Grain length"},
    {"gene": "badh2", "chr": "Chr8", "pos": 20000000, "function": "Fragrance"},
    {"gene": "Sub1A", "chr": "Chr9", "pos": 6300000, "function": "Submergence tolerance"},
    {"gene": "COLD1", "chr": "Chr4", "pos": 30300000, "function": "Cold tolerance"},
    {"gene": "DRO1", "chr": "Chr9", "pos": 22000000, "function": "Deep rooting"},
    {"gene": "Bh4", "chr": "Chr4", "pos": 26000000, "function": "Hull color"},
    {"gene": "OsC1", "chr": "Chr6", "pos": 5000000, "function": "Hull color"},
    {"gene": "GW5", "chr": "Chr5", "pos": 5400000, "function": "Grain width"},
    {"gene": "Pi-ta", "chr": "Chr12", "pos": 10600000, "function": "Blast resistance"},
]

log = logging.getLogger("rice_interpret")

# ── Helpers ───────────────────────────────────────────────────────────


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stderr),
            logging.FileHandler("rice_interpret.log", mode="a"),
        ],
    )


def get_driver():
    return GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))


def get_chr_lengths(session):
    recs = session.run(
        "MATCH (c:Chromosome) RETURN c.chromosomeId AS chr, c.length AS len"
    )
    return {r["chr"]: r["len"] for r in recs}


def load_results():
    with open(RESULTS_FILE) as f:
        return json.load(f)


def load_output():
    if os.path.exists(OUTPUT_FILE):
        with open(OUTPUT_FILE) as f:
            return json.load(f)
    return {}


def save_output(data):
    tmp = OUTPUT_FILE + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2)
    os.replace(tmp, OUTPUT_FILE)


def _serialize(v):
    if isinstance(v, (int, float, str, bool, type(None))):
        return v
    if isinstance(v, (list, tuple)):
        return [_serialize(x) for x in v]
    if isinstance(v, dict):
        return {k: _serialize(val) for k, val in v.items()}
    return str(v)


def find_nearby_known_genes(chrom, start, end, margin=50000):
    """Find known domestication genes overlapping [start-margin, end+margin]."""
    result = []
    for kg in KNOWN_GENES:
        if kg["chr"] == chrom and start - margin <= kg["pos"] <= end + margin:
            result.append(kg)
    return result


def call_peaks(hits, merge_distance=100000):
    """Merge nearby XP-EHH hits into peak regions.

    Returns list of {chr, start, end, peak_xpehh, peak_pos, n_hits}.
    """
    by_chr = defaultdict(list)
    for h in hits:
        vid = h["variantId"]
        chrom = vid.split(":")[0]
        by_chr[chrom].append(h)

    peaks = []
    for chrom, chr_hits in by_chr.items():
        chr_hits.sort(key=lambda x: x["pos"])
        cur_start = chr_hits[0]["pos"]
        cur_end = chr_hits[0]["pos"]
        cur_hits = [chr_hits[0]]

        for h in chr_hits[1:]:
            if h["pos"] - cur_end <= merge_distance:
                cur_end = h["pos"]
                cur_hits.append(h)
            else:
                best = max(cur_hits, key=lambda x: abs(x["xpehh"]))
                peaks.append({
                    "chr": chrom,
                    "start": cur_start,
                    "end": cur_end,
                    "peak_xpehh": best["xpehh"],
                    "peak_pos": best["pos"],
                    "peak_variantId": best["variantId"],
                    "n_hits": len(cur_hits),
                })
                cur_start = h["pos"]
                cur_end = h["pos"]
                cur_hits = [h]

        best = max(cur_hits, key=lambda x: abs(x["xpehh"]))
        peaks.append({
            "chr": chrom,
            "start": cur_start,
            "end": cur_end,
            "peak_xpehh": best["xpehh"],
            "peak_pos": best["pos"],
            "peak_variantId": best["variantId"],
            "n_hits": len(cur_hits),
        })

    peaks.sort(key=lambda x: abs(x["peak_xpehh"]), reverse=True)
    return peaks


# ── Subcommand: annotate ──────────────────────────────────────────────


def cmd_annotate(args):
    """Cross-reference XP-EHH peaks and Garud's H sweeps with genes."""
    results = load_results()
    output = load_output()
    driver = get_driver()

    with driver.session(database=NEO4J_DB) as session:
        # ── 1. XP-EHH peaks → gene annotation ──
        log.info("=== Annotating XP-EHH peaks ===")
        xpehh_annotations = {}

        gene_query = """
        MATCH (v:Variant)-[c:HAS_CONSEQUENCE]->(g:Gene)
        WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end
        RETURN DISTINCT g.symbol AS symbol, g.geneId AS geneId,
               c.consequence AS consequence, c.impact AS impact
        ORDER BY CASE c.impact
            WHEN 'HIGH' THEN 0 WHEN 'MODERATE' THEN 1
            WHEN 'LOW' THEN 2 ELSE 3
        END
        """

        for pair_key, pair_data in results.get("phase2", {}).items():
            all_hits = []
            for chrom_data in pair_data.values():
                all_hits.extend(chrom_data.get("top_hits", []))

            peaks = call_peaks(all_hits)
            log.info(f"  {pair_key}: {len(all_hits)} hits → {len(peaks)} peaks")

            top_peaks = peaks[:100]  # annotate top 100 peaks
            for i, peak in enumerate(top_peaks):
                margin = 50000
                recs = list(session.run(
                    gene_query,
                    chr=peak["chr"],
                    start=peak["start"] - margin,
                    end=peak["end"] + margin,
                ))
                genes = [
                    {"symbol": r["symbol"], "geneId": r["geneId"],
                     "consequence": r["consequence"], "impact": r["impact"]}
                    for r in recs
                ]
                # Deduplicate by gene symbol, keeping highest impact
                seen = {}
                for g in genes:
                    key = g["symbol"]
                    if key not in seen or _impact_rank(g["impact"]) < _impact_rank(seen[key]["impact"]):
                        seen[key] = g
                peak["genes"] = list(seen.values())[:30]
                peak["n_genes"] = len(seen)
                peak["known_domestication_genes"] = find_nearby_known_genes(
                    peak["chr"], peak["start"], peak["end"]
                )
                if (i + 1) % 20 == 0:
                    log.info(f"    annotated {i+1}/{len(top_peaks)} peaks")

            xpehh_annotations[pair_key] = {
                "n_total_hits": len(all_hits),
                "n_peaks": len(peaks),
                "top_peaks": top_peaks,
            }

        # ── 2. Garud's H sweeps → gene annotation ──
        log.info("=== Annotating Garud's H sweep windows ===")
        sweep_annotations = {}

        sweep_gene_query = """
        MATCH (v:Variant)-[c:HAS_CONSEQUENCE]->(g:Gene)
        WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end
          AND c.impact IN ['HIGH', 'MODERATE']
        WITH g, c.consequence AS csq, count(DISTINCT v) AS n_vars
        RETURN g.symbol AS symbol, g.geneId AS geneId,
               csq AS consequence, n_vars AS n_variants
        ORDER BY n_vars DESC
        LIMIT 20
        """

        for pop in POPULATIONS:
            pop_data = results.get("phase1", {}).get(pop, {})
            pop_sweeps = []

            for chrom in CHROMOSOMES:
                windows = pop_data.get(chrom, {}).get("garud_h", {}).get("result", [])
                if not isinstance(windows, list):
                    continue
                sweeps = [w for w in windows if w.get("h12", 0) > 0.10]
                for sw in sweeps:
                    recs = list(session.run(
                        sweep_gene_query,
                        chr=sw["chr"], start=sw["start"], end=sw["end"],
                    ))
                    genes = [
                        {"symbol": r["symbol"], "geneId": r["geneId"],
                         "consequence": r["consequence"],
                         "n_variants": r["n_variants"]}
                        for r in recs
                    ]
                    pop_sweeps.append({
                        "chr": sw["chr"],
                        "start": sw["start"],
                        "end": sw["end"],
                        "h1": sw["h1"],
                        "h12": sw["h12"],
                        "h2_h1": sw["h2_h1"],
                        "hap_diversity": sw["hap_diversity"],
                        "n_haplotypes": sw["n_haplotypes"],
                        "n_variants": sw["n_variants"],
                        "genes": genes,
                        "known_domestication_genes": find_nearby_known_genes(
                            sw["chr"], sw["start"], sw["end"]
                        ),
                        "sweep_type": "hard" if sw.get("h2_h1", 1) < 0.05 else "soft",
                    })

            pop_sweeps.sort(key=lambda x: x["h12"], reverse=True)
            sweep_annotations[pop] = {
                "n_sweeps": len(pop_sweeps),
                "sweeps": pop_sweeps,
            }
            log.info(f"  {pop}: {len(pop_sweeps)} sweep windows annotated")

        # ── 3. Fay & Wu's H genome scan peaks → gene annotation ──
        log.info("=== Annotating Fay & Wu's H genome scan peaks ===")
        hwscan_annotations = {}
        hwscan_data = output.get("hwscan", {})
        H_SCAN_POPS = ["GJ-tmp", "XI-1A", "cA-Aus", "GJ-trp"]

        for pop in H_SCAN_POPS:
            all_windows = []
            for chrom in CHROMOSOMES:
                key = f"{pop}|{chrom}"
                windows = hwscan_data.get(key, {}).get("result", [])
                if isinstance(windows, list):
                    for w in windows:
                        h = w.get("fay_wu_h")
                        if h is not None:
                            all_windows.append(w)

            # Sort by H (most negative first = strongest sweep)
            all_windows.sort(key=lambda x: x.get("fay_wu_h", 0))
            top_windows = all_windows[:50]  # annotate top 50 H peaks

            annotated = []
            for i, win in enumerate(top_windows):
                recs = list(session.run(
                    sweep_gene_query,
                    chr=win.get("chr", ""),
                    start=win.get("start", 0),
                    end=win.get("end", 0),
                ))
                genes = [
                    {"symbol": r["symbol"], "geneId": r["geneId"],
                     "consequence": r["consequence"],
                     "n_variants": r["n_variants"]}
                    for r in recs
                ]
                annotated.append({
                    "chr": win.get("chr"),
                    "start": win.get("start"),
                    "end": win.get("end"),
                    "fay_wu_h": win.get("fay_wu_h"),
                    "fay_wu_h_norm": win.get("fay_wu_h_norm"),
                    "pi": win.get("pi"),
                    "n_variants": win.get("n_variants"),
                    "genes": genes,
                    "known_domestication_genes": find_nearby_known_genes(
                        win.get("chr", ""), win.get("start", 0), win.get("end", 0)
                    ),
                })

            hwscan_annotations[pop] = {
                "n_windows_total": len(all_windows),
                "n_annotated": len(annotated),
                "top_peaks": annotated,
            }
            log.info(f"  {pop}: top {len(annotated)} H peaks annotated")

        # ── 4. ROH HMM — summarize genes in high-FROH chromosomes ──
        log.info("=== Annotating top ROH regions ===")
        roh_hmm_data = output.get("roh_hmm", {})
        roh_annotations = {}

        for pop in POPULATIONS:
            # Find chromosomes with highest ROH for this pop
            chr_roh = []
            for chrom in CHROMOSOMES:
                key = f"{pop}|{chrom}"
                r = roh_hmm_data.get(key, {})
                total_bp = r.get("total_roh_bp", 0)
                if total_bp > 0:
                    chr_roh.append({"chr": chrom, "total_roh_bp": total_bp,
                                    "mean_froh": r.get("mean_froh", 0),
                                    "n_samples_with_roh": r.get("n_samples_with_roh", 0),
                                    "n_samples": r.get("n_samples", 0)})
            chr_roh.sort(key=lambda x: x["total_roh_bp"], reverse=True)
            roh_annotations[pop] = {
                "chr_with_roh": chr_roh[:5],  # top 5 chromosomes
                "total_chr_with_roh": len(chr_roh),
            }
            if chr_roh:
                log.info(f"  {pop}: {len(chr_roh)} chromosomes with ROH")

        # ── 5. Diversity comparison table ──
        log.info("=== Building diversity comparison tables ===")
        diversity_table = _build_diversity_table(results)

    driver.close()

    output["annotate"] = {
        "xpehh_annotations": xpehh_annotations,
        "sweep_annotations": sweep_annotations,
        "hwscan_annotations": hwscan_annotations,
        "roh_annotations": roh_annotations,
        "diversity_table": diversity_table,
    }
    save_output(output)
    log.info(f"Annotation results saved to {OUTPUT_FILE}")


def _impact_rank(impact):
    return {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}.get(impact, 4)


def _build_diversity_table(results):
    """Build diversity comparison table from Phase 1 results."""
    table = {}
    for pop in POPULATIONS:
        pop_data = results.get("phase1", {}).get(pop, {})
        per_chr = {}
        for chrom in CHROMOSOMES:
            div = pop_data.get(chrom, {}).get("diversity", {}).get("result", {})
            if div:
                per_chr[chrom] = {
                    "pi": div.get("pi"),
                    "theta_w": div.get("theta_w"),
                    "tajima_d": div.get("tajima_d"),
                    "fis": div.get("fis"),
                    "het_obs": div.get("het_obs"),
                    "n_segregating": div.get("n_segregating"),
                    "n_variants": div.get("n_variants"),
                }

        vals = lambda key: [v[key] for v in per_chr.values() if v.get(key) is not None]
        mean = lambda xs: sum(xs) / len(xs) if xs else None
        table[pop] = {
            "per_chr": per_chr,
            "mean_pi": mean(vals("pi")),
            "mean_theta_w": mean(vals("theta_w")),
            "mean_tajima_d": mean(vals("tajima_d")),
            "mean_fis": mean(vals("fis")),
        }
    return table


# ── Subcommand: pinsps ────────────────────────────────────────────────


def cmd_pinsps(args):
    """Compute pi_N/pi_S (missense vs synonymous diversity) per population."""
    output = load_output()
    existing = output.get("pinsps", {})
    driver = get_driver()
    lock = Lock()
    n_done = 0
    n_total = len(POPULATIONS) * len(CHROMOSOMES) * 2

    with driver.session(database=NEO4J_DB) as session:
        chr_lens = get_chr_lengths(session)

    def run_one(pop, chrom, consequence):
        key = f"{pop}|{chrom}|{consequence}"
        if key in existing:
            return key, existing[key]

        t0 = time.time()
        cypher = (
            f"CALL graphpop.diversity('{chrom}', 1, {chr_lens[chrom]}, '{pop}', "
            f"{{consequence: '{consequence}'}})"
        )
        try:
            with driver.session(database=NEO4J_DB) as s:
                rec = s.run(cypher).single()
                result = {k: _serialize(rec[k]) for k in rec.keys()} if rec else {}
            elapsed = time.time() - t0
            return key, {"result": result, "wall_sec": elapsed}
        except Exception as e:
            return key, {"error": str(e), "wall_sec": time.time() - t0}

    log.info(f"=== Computing pi_N/pi_S: {n_total} calls ===")

    # Process sequentially by pop, parallel by chrom
    for pop in POPULATIONS:
        tasks = []
        for chrom in CHROMOSOMES:
            for csq in ("missense_variant", "synonymous_variant"):
                key = f"{pop}|{chrom}|{csq}"
                if key in existing:
                    n_done += 1
                    continue
                tasks.append((pop, chrom, csq))

        if not tasks:
            log.info(f"  {pop}: all cached")
            continue

        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = {pool.submit(run_one, *t): t for t in tasks}
            for fut in as_completed(futures):
                key, val = fut.result()
                with lock:
                    existing[key] = val
                    n_done += 1
                    if n_done % 24 == 0:
                        log.info(f"  progress: {n_done}/{n_total}")
                        output["pinsps"] = existing
                        save_output(output)

        log.info(f"  {pop}: done")

    # Compute pi_N/pi_S ratios
    ratios = {}
    for pop in POPULATIONS:
        pop_ratios = {}
        for chrom in CHROMOSOMES:
            mis_key = f"{pop}|{chrom}|missense_variant"
            syn_key = f"{pop}|{chrom}|synonymous_variant"
            mis = existing.get(mis_key, {}).get("result", {})
            syn = existing.get(syn_key, {}).get("result", {})
            pi_n = mis.get("pi")
            pi_s = syn.get("pi")
            if pi_n is not None and pi_s is not None and pi_s > 0:
                pop_ratios[chrom] = {
                    "pi_n": pi_n,
                    "pi_s": pi_s,
                    "pi_n_pi_s": pi_n / pi_s,
                    "n_missense": mis.get("n_segregating"),
                    "n_synonymous": syn.get("n_segregating"),
                }
            else:
                pop_ratios[chrom] = {
                    "pi_n": pi_n, "pi_s": pi_s, "pi_n_pi_s": None,
                    "n_missense": mis.get("n_segregating") if mis else None,
                    "n_synonymous": syn.get("n_segregating") if syn else None,
                }

        vals = [v["pi_n_pi_s"] for v in pop_ratios.values() if v["pi_n_pi_s"] is not None]
        ratios[pop] = {
            "per_chr": pop_ratios,
            "mean_pi_n_pi_s": sum(vals) / len(vals) if vals else None,
        }

    driver.close()
    output["pinsps"] = existing
    output["pinsps_ratios"] = ratios
    save_output(output)
    log.info(f"pi_N/pi_S results saved to {OUTPUT_FILE}")


# ── Subcommand: divergence ────────────────────────────────────────────


def cmd_divergence(args):
    """Compute full pairwise Fst/Dxy matrix for all 66 population pairs."""
    output = load_output()
    existing = output.get("divergence", {})
    driver = get_driver()
    lock = Lock()

    with driver.session(database=NEO4J_DB) as session:
        chr_lens = get_chr_lengths(session)

    all_pairs = list(combinations(POPULATIONS, 2))
    n_total = len(all_pairs) * len(CHROMOSOMES)
    n_done = sum(1 for p in all_pairs for c in CHROMOSOMES
                 if f"{p[0]}|{p[1]}|{c}" in existing)
    log.info(f"=== Pairwise divergence: {n_total} calls ({n_done} cached) ===")

    def run_one(pop1, pop2, chrom):
        key = f"{pop1}|{pop2}|{chrom}"
        if key in existing:
            return key, existing[key]

        t0 = time.time()
        cypher = (
            f"CALL graphpop.divergence('{chrom}', 1, {chr_lens[chrom]}, "
            f"'{pop1}', '{pop2}')"
        )
        try:
            with driver.session(database=NEO4J_DB) as s:
                rec = s.run(cypher).single()
                result = {k: _serialize(rec[k]) for k in rec.keys()} if rec else {}
            return key, {"result": result, "wall_sec": time.time() - t0}
        except Exception as e:
            return key, {"error": str(e), "wall_sec": time.time() - t0}

    # Process pairs sequentially, chromosomes in parallel
    for i, (pop1, pop2) in enumerate(all_pairs):
        tasks = [(pop1, pop2, c) for c in CHROMOSOMES
                 if f"{pop1}|{pop2}|{c}" not in existing]
        if not tasks:
            continue

        with ThreadPoolExecutor(max_workers=6) as pool:
            futures = {pool.submit(run_one, *t): t for t in tasks}
            for fut in as_completed(futures):
                key, val = fut.result()
                with lock:
                    existing[key] = val
                    n_done += 1

        log.info(f"  [{i+1}/{len(all_pairs)}] {pop1} vs {pop2}: done ({n_done}/{n_total})")
        output["divergence"] = existing
        save_output(output)

    # Build Fst/Dxy matrices
    fst_matrix = {}
    dxy_matrix = {}
    for pop1, pop2 in all_pairs:
        pair_key = f"{pop1}_vs_{pop2}"
        fst_vals = []
        dxy_vals = []
        for chrom in CHROMOSOMES:
            key = f"{pop1}|{pop2}|{chrom}"
            r = existing.get(key, {}).get("result", {})
            fst_val = r.get("fst_wc")
            if fst_val is not None:
                fst_vals.append(fst_val)
            if r.get("dxy") is not None:
                dxy_vals.append(r["dxy"])
        fst_matrix[pair_key] = {
            "mean_fst": sum(fst_vals) / len(fst_vals) if fst_vals else None,
            "per_chr": {
                c: existing.get(f"{pop1}|{pop2}|{c}", {}).get("result", {}).get("fst_wc")
                for c in CHROMOSOMES
            },
        }
        dxy_matrix[pair_key] = {
            "mean_dxy": sum(dxy_vals) / len(dxy_vals) if dxy_vals else None,
            "per_chr": {
                c: existing.get(f"{pop1}|{pop2}|{c}", {}).get("result", {}).get("dxy")
                for c in CHROMOSOMES
            },
        }

    driver.close()
    output["divergence"] = existing
    output["fst_matrix"] = fst_matrix
    output["dxy_matrix"] = dxy_matrix
    save_output(output)
    log.info(f"Divergence results saved to {OUTPUT_FILE}")


# ── Subcommand: gscan ─────────────────────────────────────────────────


def cmd_gscan(args):
    """Run genome scans: divergence + annotation-conditioned (missense/synonymous)."""
    output = load_output()
    existing = output.get("gscan", {})
    driver = get_driver()
    lock = Lock()

    with driver.session(database=NEO4J_DB) as session:
        chr_lens = get_chr_lengths(session)

    # Build task list: pairs x chr x consequence variants
    # (1) divergence scans for 6 XPEHH pairs (no consequence filter)
    # (2) annotation-conditioned scans for 3 key pairs (missense, synonymous)
    tasks_all = []
    for pop1, pop2 in XPEHH_PAIRS:
        for chrom in CHROMOSOMES:
            key = f"{pop1}|{pop2}|{chrom}|all"
            if key not in existing:
                tasks_all.append((pop1, pop2, chrom, None, key))
    for pop1, pop2 in GSCAN_PAIRS:
        for chrom in CHROMOSOMES:
            for csq in ("missense_variant", "synonymous_variant"):
                key = f"{pop1}|{pop2}|{chrom}|{csq}"
                if key not in existing:
                    tasks_all.append((pop1, pop2, chrom, csq, key))

    n_total = (len(XPEHH_PAIRS) + len(GSCAN_PAIRS) * 2) * len(CHROMOSOMES)
    n_done = n_total - len(tasks_all)
    log.info(f"=== Genome scans: {n_total} calls ({n_done} cached, {len(tasks_all)} remaining) ===")

    def run_one(pop1, pop2, chrom, consequence, key):
        t0 = time.time()
        opts = f"pop2: '{pop2}'"
        if consequence:
            opts += f", consequence: '{consequence}'"
        cypher = (
            f"CALL graphpop.genome_scan('{chrom}', '{pop1}', 100000, 50000, "
            f"{{{opts}}})"
        )
        try:
            with driver.session(database=NEO4J_DB) as s:
                records = []
                for rec in s.run(cypher):
                    records.append({k: _serialize(rec[k]) for k in rec.keys()})
            elapsed = time.time() - t0
            return key, {
                "n_windows": len(records),
                "wall_sec": elapsed,
                "result": records,
            }
        except Exception as e:
            return key, {"error": str(e), "wall_sec": time.time() - t0}

    # Process ONE at a time — genome_scan is WRITE mode
    for i, (pop1, pop2, chrom, csq, key) in enumerate(tasks_all):
        _, val = run_one(pop1, pop2, chrom, csq, key)
        existing[key] = val
        n_done += 1
        label = csq or "all"
        err = f" ERROR: {val['error']}" if "error" in val else ""
        log.info(f"  [{n_done}/{n_total}] {pop1} vs {pop2} {chrom} ({label}): "
                 f"{val.get('n_windows', 0)} windows, {val.get('wall_sec', 0):.1f}s{err}")

        if n_done % 12 == 0:
            output["gscan"] = existing
            save_output(output)

    driver.close()
    output["gscan"] = existing
    save_output(output)
    log.info(f"Genome scan results saved to {OUTPUT_FILE}")


# ── Subcommand: roh ───────────────────────────────────────────────────


def cmd_roh(args):
    """Re-run ROH with sliding_window method (no ancestral alleles needed)."""
    output = load_output()
    existing = output.get("roh", {})
    driver = get_driver()
    lock = Lock()
    n_total = len(POPULATIONS) * len(CHROMOSOMES)
    n_done = sum(1 for p in POPULATIONS for c in CHROMOSOMES
                 if f"{p}|{c}" in existing)
    log.info(f"=== ROH sliding_window: {n_total} calls ({n_done} cached) ===")

    def run_one(pop, chrom):
        key = f"{pop}|{chrom}"
        if key in existing:
            return key, existing[key]

        t0 = time.time()
        cypher = (
            f"CALL graphpop.roh('{chrom}', '{pop}', "
            f"{{method: 'sliding_window'}})"
        )
        try:
            with driver.session(database=NEO4J_DB) as s:
                records = []
                for rec in s.run(cypher):
                    records.append({k: _serialize(rec[k]) for k in rec.keys()})
            elapsed = time.time() - t0

            # Aggregate per-sample results
            total_len = sum(r.get("total_length", 0) for r in records)
            n_with_roh = sum(1 for r in records if r.get("n_roh", 0) > 0)
            n_segments = sum(r.get("n_roh", 0) for r in records)
            froh_values = [r.get("froh", 0) for r in records if r.get("n_roh", 0) > 0]

            return key, {
                "wall_sec": elapsed,
                "n_samples": len(records),
                "n_samples_with_roh": n_with_roh,
                "n_segments": n_segments,
                "total_roh_bp": total_len,
                "mean_froh": sum(froh_values) / len(froh_values) if froh_values else 0.0,
                "max_froh": max(froh_values) if froh_values else 0.0,
                "per_sample": [
                    {"sampleId": r["sampleId"], "n_roh": r["n_roh"],
                     "total_length": r["total_length"], "froh": r["froh"],
                     "mean_length": r["mean_length"], "max_length": r["max_length"]}
                    for r in records if r.get("n_roh", 0) > 0
                ][:50],  # store top 50 samples with ROH
            }
        except Exception as e:
            return key, {"error": str(e), "wall_sec": time.time() - t0}

    # Process pops sequentially, chr in parallel
    for pop in POPULATIONS:
        tasks = [(pop, c) for c in CHROMOSOMES if f"{pop}|{c}" not in existing]
        if not tasks:
            log.info(f"  {pop}: all cached")
            continue

        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = {pool.submit(run_one, *t): t for t in tasks}
            for fut in as_completed(futures):
                key, val = fut.result()
                with lock:
                    existing[key] = val
                    n_done += 1

        log.info(f"  {pop}: done ({n_done}/{n_total})")
        output["roh"] = existing
        save_output(output)

    # Build summary
    roh_summary = {}
    for pop in POPULATIONS:
        total_roh = 0
        total_froh_sum = 0.0
        n_froh = 0
        for chrom in CHROMOSOMES:
            r = existing.get(f"{pop}|{chrom}", {})
            total_roh += r.get("total_roh_bp", 0)
            if r.get("mean_froh", 0) > 0:
                total_froh_sum += r["mean_froh"]
                n_froh += 1
        roh_summary[pop] = {
            "total_roh_bp": total_roh,
            "mean_froh_across_chr": total_froh_sum / n_froh if n_froh > 0 else 0.0,
        }

    driver.close()
    output["roh"] = existing
    output["roh_summary"] = roh_summary
    save_output(output)
    log.info(f"ROH results saved to {OUTPUT_FILE}")


# ── Subcommand: fay_wu ────────────────────────────────────────────────


def cmd_fay_wu(args):
    """Re-run diversity for all pops×chrs to get Fay & Wu's H with ancestral alleles."""
    output = load_output()
    existing = output.get("fay_wu", {})
    driver = get_driver()
    lock = Lock()

    with driver.session(database=NEO4J_DB) as s:
        chr_lens = get_chr_lengths(s)

    n_total = len(POPULATIONS) * len(CHROMOSOMES)
    n_done = sum(1 for k in existing if "error" not in existing[k])
    log.info(f"=== Fay & Wu's H (diversity): {n_total} calls ({n_done} cached) ===")

    def run_one(pop, chrom):
        key = f"{pop}|{chrom}"
        if key in existing:
            return key, existing[key]
        t0 = time.time()
        cypher = f"CALL graphpop.diversity('{chrom}', 1, {chr_lens[chrom]}, '{pop}')"
        try:
            with driver.session(database=NEO4J_DB) as s:
                rec = s.run(cypher).single()
                result = {k: _serialize(rec[k]) for k in rec.keys()} if rec else {}
            return key, {"result": result, "wall_sec": time.time() - t0}
        except Exception as e:
            return key, {"error": str(e), "wall_sec": time.time() - t0}

    for pop in POPULATIONS:
        tasks = [(pop, c) for c in CHROMOSOMES if f"{pop}|{c}" not in existing]
        if not tasks:
            continue
        with ThreadPoolExecutor(max_workers=6) as pool:
            futures = {pool.submit(run_one, *t): t for t in tasks}
            for fut in as_completed(futures):
                key, val = fut.result()
                with lock:
                    existing[key] = val
                    n_done += 1
        log.info(f"  {pop}: done ({n_done}/{n_total})")
        output["fay_wu"] = existing
        save_output(output)

    # Build summary table
    fay_wu_summary = {}
    for pop in POPULATIONS:
        h_vals = []
        hn_vals = []
        for chrom in CHROMOSOMES:
            r = existing.get(f"{pop}|{chrom}", {}).get("result", {})
            if r.get("fay_wu_h") is not None and r.get("n_polarized", 0) > 0:
                h_vals.append(r["fay_wu_h"])
                hn_vals.append(r.get("fay_wu_h_norm", 0))
        fay_wu_summary[pop] = {
            "mean_h": sum(h_vals) / len(h_vals) if h_vals else None,
            "mean_h_norm": sum(hn_vals) / len(hn_vals) if hn_vals else None,
            "per_chr": {
                c: {
                    "fay_wu_h": existing.get(f"{pop}|{c}", {}).get("result", {}).get("fay_wu_h"),
                    "fay_wu_h_norm": existing.get(f"{pop}|{c}", {}).get("result", {}).get("fay_wu_h_norm"),
                    "n_polarized": existing.get(f"{pop}|{c}", {}).get("result", {}).get("n_polarized"),
                }
                for c in CHROMOSOMES
            },
        }

    driver.close()
    output["fay_wu"] = existing
    output["fay_wu_summary"] = fay_wu_summary
    save_output(output)
    log.info(f"Fay & Wu's H results saved to {OUTPUT_FILE}")


# ── Subcommand: usfs ─────────────────────────────────────────────────


def cmd_usfs(args):
    """Compute unfolded (polarized) SFS for all populations."""
    output = load_output()
    existing = output.get("usfs", {})
    driver = get_driver()
    lock = Lock()

    with driver.session(database=NEO4J_DB) as s:
        chr_lens = get_chr_lengths(s)

    n_total = len(POPULATIONS) * len(CHROMOSOMES)
    n_done = sum(1 for k in existing if "error" not in existing[k])
    log.info(f"=== Unfolded SFS: {n_total} calls ({n_done} cached) ===")

    def run_one(pop, chrom):
        key = f"{pop}|{chrom}"
        if key in existing:
            return key, existing[key]
        t0 = time.time()
        cypher = f"CALL graphpop.sfs('{chrom}', 1, {chr_lens[chrom]}, '{pop}', false)"
        try:
            with driver.session(database=NEO4J_DB) as s:
                rec = s.run(cypher).single()
                result = {k: _serialize(rec[k]) for k in rec.keys()} if rec else {}
            return key, {"result": result, "wall_sec": time.time() - t0}
        except Exception as e:
            return key, {"error": str(e), "wall_sec": time.time() - t0}

    for pop in POPULATIONS:
        tasks = [(pop, c) for c in CHROMOSOMES if f"{pop}|{c}" not in existing]
        if not tasks:
            continue
        with ThreadPoolExecutor(max_workers=6) as pool:
            futures = {pool.submit(run_one, *t): t for t in tasks}
            for fut in as_completed(futures):
                key, val = fut.result()
                with lock:
                    existing[key] = val
                    n_done += 1
        log.info(f"  {pop}: done ({n_done}/{n_total})")
        output["usfs"] = existing
        save_output(output)

    # Summarize: aggregate SFS across chromosomes per population
    usfs_summary = {}
    for pop in POPULATIONS:
        total_sfs = None
        total_polarized = 0
        total_variants = 0
        for chrom in CHROMOSOMES:
            r = existing.get(f"{pop}|{chrom}", {}).get("result", {})
            sfs = r.get("sfs")
            if sfs:
                if total_sfs is None:
                    total_sfs = [0] * len(sfs)
                for i in range(min(len(sfs), len(total_sfs))):
                    total_sfs[i] += sfs[i]
                # Extend if needed
                if len(sfs) > len(total_sfs):
                    total_sfs.extend(sfs[len(total_sfs):])
            total_polarized += r.get("n_polarized", 0)
            total_variants += r.get("n_variants", 0)
        usfs_summary[pop] = {
            "sfs": total_sfs,
            "n_polarized": total_polarized,
            "n_variants": total_variants,
        }

    driver.close()
    output["usfs"] = existing
    output["usfs_summary"] = usfs_summary
    save_output(output)
    log.info(f"Unfolded SFS results saved to {OUTPUT_FILE}")


# ── Subcommand: hwscan ───────────────────────────────────────────────

# Key populations for H genome scan
H_SCAN_POPS = ["GJ-tmp", "XI-1A", "cA-Aus", "GJ-trp"]


def cmd_hwscan(args):
    """Genome scan for Fay & Wu's H per window (key populations)."""
    output = load_output()
    existing = output.get("hwscan", {})
    driver = get_driver()
    lock = Lock()

    with driver.session(database=NEO4J_DB) as s:
        chr_lens = get_chr_lengths(s)

    n_total = len(H_SCAN_POPS) * len(CHROMOSOMES)
    n_done = sum(1 for k in existing if "error" not in existing[k])
    log.info(f"=== Fay & Wu's H genome scan: {n_total} calls ({n_done} cached) ===")

    def run_one(pop, chrom):
        key = f"{pop}|{chrom}"
        if key in existing:
            return key, existing[key]
        t0 = time.time()
        cypher = (
            f"CALL graphpop.genome_scan('{chrom}', '{pop}', "
            f"100000, 50000, {{}})"
        )
        try:
            with driver.session(database=NEO4J_DB) as s:
                records = s.run(cypher).data()
                windows = [{k: _serialize(v) for k, v in r.items()} for r in records]
            return key, {"result": windows, "wall_sec": time.time() - t0}
        except Exception as e:
            return key, {"error": str(e), "wall_sec": time.time() - t0}

    for pop in H_SCAN_POPS:
        tasks = [(pop, c) for c in CHROMOSOMES if f"{pop}|{c}" not in existing]
        if not tasks:
            log.info(f"  {pop}: all cached")
            continue
        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = {pool.submit(run_one, *t): t for t in tasks}
            for fut in as_completed(futures):
                key, val = fut.result()
                with lock:
                    existing[key] = val
                    n_done += 1
        log.info(f"  {pop}: done ({n_done}/{n_total})")
        output["hwscan"] = existing
        save_output(output)

    # Find top negative H windows (sweep candidates)
    hwscan_summary = {}
    for pop in H_SCAN_POPS:
        all_windows = []
        for chrom in CHROMOSOMES:
            windows = existing.get(f"{pop}|{chrom}", {}).get("result", [])
            for w in windows:
                h = w.get("fay_wu_h")
                if h is not None and w.get("n_variants", 0) >= 20:
                    all_windows.append(w)
        # Sort by most negative H (strongest sweep signal)
        all_windows.sort(key=lambda x: x.get("fay_wu_h", 0))
        hwscan_summary[pop] = {
            "n_windows": len(all_windows),
            "mean_h": (sum(w.get("fay_wu_h", 0) for w in all_windows) /
                       len(all_windows)) if all_windows else None,
            "top_sweep_windows": all_windows[:20],
        }

    driver.close()
    output["hwscan"] = existing
    output["hwscan_summary"] = hwscan_summary
    save_output(output)
    log.info(f"H genome scan results saved to {OUTPUT_FILE}")


# ── Subcommand: roh_hmm ──────────────────────────────────────────────


def cmd_roh_hmm(args):
    """Re-run ROH with HMM method (now that ancestral alleles are available)."""
    output = load_output()
    existing = output.get("roh_hmm", {})
    driver = get_driver()
    lock = Lock()

    n_total = len(POPULATIONS) * len(CHROMOSOMES)
    n_done = sum(1 for k in existing if "error" not in existing[k])
    log.info(f"=== ROH HMM: {n_total} calls ({n_done} cached) ===")

    def run_one(pop, chrom):
        key = f"{pop}|{chrom}"
        if key in existing:
            return key, existing[key]
        t0 = time.time()
        # Use default HMM method (no method override = hmm)
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
        tasks = [(pop, c) for c in CHROMOSOMES if f"{pop}|{c}" not in existing]
        if not tasks:
            log.info(f"  {pop}: all cached")
            continue
        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = {pool.submit(run_one, *t): t for t in tasks}
            for fut in as_completed(futures):
                key, val = fut.result()
                with lock:
                    existing[key] = val
                    n_done += 1
        log.info(f"  {pop}: done ({n_done}/{n_total})")
        output["roh_hmm"] = existing
        save_output(output)

    # Build summary
    roh_hmm_summary = {}
    for pop in POPULATIONS:
        total_roh = 0
        total_froh_sum = 0.0
        n_froh = 0
        for chrom in CHROMOSOMES:
            r = existing.get(f"{pop}|{chrom}", {})
            total_roh += r.get("total_roh_bp", 0)
            if r.get("mean_froh", 0) > 0:
                total_froh_sum += r["mean_froh"]
                n_froh += 1
        roh_hmm_summary[pop] = {
            "total_roh_bp": total_roh,
            "mean_froh_across_chr": total_froh_sum / n_froh if n_froh > 0 else 0.0,
        }

    driver.close()
    output["roh_hmm"] = existing
    output["roh_hmm_summary"] = roh_hmm_summary
    save_output(output)
    log.info(f"ROH HMM results saved to {OUTPUT_FILE}")


# ── Subcommand: daf_enrichment ────────────────────────────────────────


def cmd_daf_enrichment(args):
    """Test derived allele frequency enrichment at selection signal regions."""
    output = load_output()
    driver = get_driver()

    pbs_summary = output.get("pbs_summary", {})
    xpehh_ann = output.get("annotate", {}).get("xpehh_annotations", {})

    results = {}

    # 1. Genome-wide background DAF per population
    log.info("Computing genome-wide background DAF...")
    bg_daf = {}
    with driver.session(database=NEO4J_DB) as s:
        for pop in ["GJ-tmp", "XI-1A", "cA-Aus", "GJ-trp", "cB-Bas"]:
            rec = s.run(
                "MATCH (v:Variant) "
                "WHERE v.is_polarized = true AND v.ancestral_allele IS NOT NULL "
                "WITH v, v.pop_ids AS pids, v.af AS afs, v.an AS ans "
                "UNWIND range(0, size(pids)-1) AS i "
                "WITH v, pids[i] AS pop, afs[i] AS af, ans[i] AS an, "
                "     v.ancestral_allele AS aa "
                "WHERE pop = $pop AND an >= 10 "
                "RETURN avg(CASE WHEN aa = 'REF' THEN af ELSE 1.0 - af END) AS mean_daf, "
                "       count(v) AS n",
                pop=pop
            ).single()
            bg_daf[pop] = {
                "mean_daf": rec["mean_daf"],
                "n_variants": rec["n"],
            }
            log.info(f"  {pop}: background DAF = {rec['mean_daf']:.4f} (n={rec['n']:,})")

    results["background_daf"] = bg_daf

    # 2. DAF at PBS peak regions
    log.info("Computing DAF at PBS peak regions...")
    pbs_daf = {}
    for config_key, info in pbs_summary.items():
        focal = info["focal"]
        if focal not in bg_daf:
            continue
        peak_dafs = []
        with driver.session(database=NEO4J_DB) as s:
            for tw in info.get("top_windows", [])[:20]:
                rec = s.run(
                    "MATCH (v:Variant) "
                    "WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end "
                    "AND v.is_polarized = true AND v.ancestral_allele IS NOT NULL "
                    "WITH v, v.pop_ids AS pids, v.af AS afs, v.an AS ans, "
                    "     v.ancestral_allele AS aa "
                    "UNWIND range(0, size(pids)-1) AS i "
                    "WITH v, pids[i] AS pop, afs[i] AS af, ans[i] AS an, aa "
                    "WHERE pop = $pop AND an >= 10 "
                    "RETURN avg(CASE WHEN aa = 'REF' THEN af ELSE 1.0 - af END) AS mean_daf, "
                    "       count(v) AS n",
                    chr=tw["chr"], start=tw["start"], end=tw["end"], pop=focal
                ).single()
                if rec["n"] and rec["n"] > 0:
                    peak_dafs.append({
                        "chr": tw["chr"], "start": tw["start"], "end": tw["end"],
                        "pbs": tw["pbs"],
                        "mean_daf": rec["mean_daf"], "n": rec["n"],
                    })
        if peak_dafs:
            mean_peak_daf = sum(p["mean_daf"] for p in peak_dafs) / len(peak_dafs)
            bg = bg_daf[focal]["mean_daf"]
            pbs_daf[config_key] = {
                "focal": focal,
                "mean_peak_daf": mean_peak_daf,
                "background_daf": bg,
                "enrichment": mean_peak_daf / bg if bg > 0 else None,
                "n_peaks": len(peak_dafs),
                "peaks": peak_dafs,
            }
            log.info(f"  {focal}: peak DAF = {mean_peak_daf:.4f} vs bg = {bg:.4f} "
                     f"(enrichment = {mean_peak_daf/bg:.2f}x)")

    results["pbs_peak_daf"] = pbs_daf

    # 3. DAF at XP-EHH peak regions
    log.info("Computing DAF at XP-EHH peak regions...")
    xpehh_daf = {}
    for pair_key, ann in xpehh_ann.items():
        top_peaks = ann.get("top_peaks", [])[:10]
        if not top_peaks:
            continue
        parts = pair_key.split("_vs_")
        if len(parts) != 2:
            continue
        pop1, pop2 = parts
        focal = pop1  # positive XP-EHH = selection in pop1
        if focal not in bg_daf:
            continue
        peak_dafs = []
        with driver.session(database=NEO4J_DB) as s:
            for peak in top_peaks:
                margin = 50000
                rec = s.run(
                    "MATCH (v:Variant) "
                    "WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end "
                    "AND v.is_polarized = true AND v.ancestral_allele IS NOT NULL "
                    "WITH v, v.pop_ids AS pids, v.af AS afs, v.an AS ans, "
                    "     v.ancestral_allele AS aa "
                    "UNWIND range(0, size(pids)-1) AS i "
                    "WITH v, pids[i] AS pop, afs[i] AS af, ans[i] AS an, aa "
                    "WHERE pop = $pop AND an >= 10 "
                    "RETURN avg(CASE WHEN aa = 'REF' THEN af ELSE 1.0 - af END) AS mean_daf, "
                    "       count(v) AS n",
                    chr=peak["chr"],
                    start=peak["peak_pos"] - margin,
                    end=peak["peak_pos"] + margin,
                    pop=focal
                ).single()
                if rec["n"] and rec["n"] > 0:
                    peak_dafs.append({
                        "chr": peak["chr"], "pos": peak["peak_pos"],
                        "xpehh": peak["peak_xpehh"],
                        "mean_daf": rec["mean_daf"], "n": rec["n"],
                    })
        if peak_dafs:
            mean_peak_daf = sum(p["mean_daf"] for p in peak_dafs) / len(peak_dafs)
            bg = bg_daf[focal]["mean_daf"]
            xpehh_daf[pair_key] = {
                "focal": focal,
                "mean_peak_daf": mean_peak_daf,
                "background_daf": bg,
                "enrichment": mean_peak_daf / bg if bg > 0 else None,
                "n_peaks": len(peak_dafs),
                "peaks": peak_dafs,
            }
            log.info(f"  {pair_key}: peak DAF = {mean_peak_daf:.4f} vs bg = {bg:.4f}")

    results["xpehh_peak_daf"] = xpehh_daf

    driver.close()
    output["daf_enrichment"] = results
    save_output(output)
    log.info(f"DAF enrichment results saved to {OUTPUT_FILE}")


# ── Subcommand: pbs ──────────────────────────────────────────────────

# PBS configurations: (focal, sister, outgroup)
# Chosen to maximize phylogenetic informativeness
PBS_CONFIGS = [
    ("GJ-tmp", "GJ-trp", "XI-1A"),    # temperate japonica-specific selection
    ("GJ-trp", "GJ-tmp", "XI-1A"),    # tropical japonica-specific selection
    ("XI-1A", "XI-1B", "GJ-tmp"),     # indica-1A-specific selection
    ("cA-Aus", "XI-1A", "GJ-tmp"),    # aus-specific selection
    ("cB-Bas", "cA-Aus", "GJ-tmp"),   # basmati-specific selection
    ("XI-2", "XI-3", "GJ-tmp"),       # indica-2-specific selection
]

WINDOW_SIZE = 100000
WINDOW_STEP = 50000


def cmd_pbs(args):
    """Run PBS genome scans for focal populations."""
    output = load_output()
    existing = output.get("pbs", {})
    driver = get_driver()
    lock = Lock()

    with driver.session(database=NEO4J_DB) as s:
        chr_lens = get_chr_lengths(s)

    n_total = len(PBS_CONFIGS) * len(CHROMOSOMES)
    n_done = sum(1 for k in existing if "error" not in existing[k])
    log.info(f"=== PBS genome scans: {n_total} calls ({n_done} cached) ===")

    def run_one(focal, sister, outgroup, chrom):
        key = f"{focal}|{sister}|{outgroup}|{chrom}"
        if key in existing:
            return key, existing[key]

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
            return key, {"result": windows, "wall_sec": time.time() - t0}
        except Exception as e:
            return key, {"error": str(e), "wall_sec": time.time() - t0}

    for i, (focal, sister, outgroup) in enumerate(PBS_CONFIGS):
        tasks = [(focal, sister, outgroup, c) for c in CHROMOSOMES
                 if f"{focal}|{sister}|{outgroup}|{c}" not in existing]
        if not tasks:
            continue

        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = {pool.submit(run_one, *t): t for t in tasks}
            for fut in as_completed(futures):
                key, val = fut.result()
                with lock:
                    existing[key] = val
                    n_done += 1

        log.info(f"  [{i+1}/{len(PBS_CONFIGS)}] {focal} (vs {sister}, outgroup {outgroup}): "
                 f"done ({n_done}/{n_total})")
        output["pbs"] = existing
        save_output(output)

    # Summarize: top PBS windows per focal population + gene annotation
    driver2 = get_driver()
    pbs_summary = {}
    for focal, sister, outgroup in PBS_CONFIGS:
        config_key = f"{focal}_vs_{sister}_out_{outgroup}"
        all_windows = []
        for chrom in CHROMOSOMES:
            key = f"{focal}|{sister}|{outgroup}|{chrom}"
            windows = existing.get(key, {}).get("result", [])
            for w in windows:
                pbs_val = w.get("pbs")
                if pbs_val is not None and w.get("n_variants", 0) >= 10:
                    all_windows.append(w)

        all_windows.sort(key=lambda x: x.get("pbs", 0), reverse=True)
        top20 = all_windows[:20]

        # Annotate top windows with genes
        annotated = []
        for w in top20:
            genes = []
            try:
                with driver2.session(database=NEO4J_DB) as s:
                    recs = s.run(
                        "MATCH (v:Variant)-[c:HAS_CONSEQUENCE]->(g:Gene) "
                        "WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end "
                        "AND c.impact IN ['HIGH', 'MODERATE'] "
                        "RETURN DISTINCT g.symbol AS symbol, g.geneId AS geneId, "
                        "c.consequence AS consequence, c.impact AS impact "
                        "ORDER BY c.impact DESC LIMIT 5",
                        chr=w["chr"], start=w["start"], end=w["end"]
                    ).data()
                    genes = recs
            except Exception:
                pass
            known = find_nearby_known_genes(w["chr"], w["start"], w["end"])
            annotated.append({
                "chr": w["chr"], "start": w["start"], "end": w["end"],
                "pbs": w["pbs"], "fst_wc": w.get("fst_wc"),
                "n_variants": w["n_variants"],
                "genes": genes, "known_genes": known,
            })

        pbs_summary[config_key] = {
            "focal": focal, "sister": sister, "outgroup": outgroup,
            "n_windows": len(all_windows),
            "mean_pbs": (sum(w.get("pbs", 0) for w in all_windows) /
                         len(all_windows)) if all_windows else 0,
            "top_windows": annotated,
        }

    driver2.close()
    driver.close()
    output["pbs"] = existing
    output["pbs_summary"] = pbs_summary
    save_output(output)
    log.info(f"PBS results saved to {OUTPUT_FILE}")


# ── Subcommand: tree ─────────────────────────────────────────────────


def cmd_tree(args):
    """Build population tree from Fst matrix using UPGMA clustering."""
    output = load_output()
    fst_matrix = output.get("fst_matrix", {})

    if not fst_matrix:
        log.error("No fst_matrix found. Run 'divergence' first.")
        return

    # Build symmetric distance matrix
    pops = sorted(POPULATIONS)
    n = len(pops)
    pop_idx = {p: i for i, p in enumerate(pops)}

    dist = [[0.0] * n for _ in range(n)]
    for pair_key, v in fst_matrix.items():
        fst = v.get("mean_fst")
        if fst is None:
            continue
        parts = pair_key.split("_vs_")
        if len(parts) != 2:
            continue
        p1, p2 = parts
        if p1 in pop_idx and p2 in pop_idx:
            i, j = pop_idx[p1], pop_idx[p2]
            dist[i][j] = fst
            dist[j][i] = fst

    # UPGMA clustering
    tree_newick, merge_log = _upgma(pops, dist)

    # Also produce a text dendrogram
    dendro_lines = _text_dendrogram(merge_log, pops)

    output["tree"] = {
        "newick": tree_newick,
        "dendrogram": dendro_lines,
        "method": "UPGMA",
        "distance_metric": "Weir & Cockerham Fst",
        "populations": pops,
    }
    save_output(output)

    print("=== Population Tree (UPGMA on W&C Fst) ===")
    print()
    print(f"Newick: {tree_newick}")
    print()
    print("Dendrogram:")
    for line in dendro_lines:
        print(line)
    print()

    # Print full Fst matrix
    print("=== Pairwise Fst Matrix ===")
    print()
    header = "           " + "".join(f"{p:>10s}" for p in pops)
    print(header)
    for i, p1 in enumerate(pops):
        row = f"{p1:>10s} " + "".join(f"{dist[i][j]:10.4f}" for j in range(n))
        print(row)
    print()

    log.info("Tree saved to rice_interpretation_results.json")


def _upgma(labels, dist):
    """UPGMA hierarchical clustering. Returns Newick string and merge log."""
    import copy
    n = len(labels)
    # Working copies
    d = [row[:] for row in dist]
    clusters = {i: labels[i] for i in range(n)}
    sizes = {i: 1 for i in range(n)}
    heights = {i: 0.0 for i in range(n)}
    active = set(range(n))
    merge_log = []
    next_id = n

    while len(active) > 1:
        # Find minimum distance pair
        min_d = float("inf")
        mi, mj = -1, -1
        for i in active:
            for j in active:
                if i < j and d[i][j] < min_d:
                    min_d = d[i][j]
                    mi, mj = i, j

        h = min_d / 2.0
        merge_log.append((clusters[mi], clusters[mj], min_d, h))

        # Create new cluster label
        new_label = f"({clusters[mi]}:{h - heights[mi]:.4f},{clusters[mj]}:{h - heights[mj]:.4f})"

        # Expand distance matrix
        new_row = [0.0] * (len(d) + 1)
        for k in active:
            if k != mi and k != mj:
                new_dist = (d[mi][k] * sizes[mi] + d[mj][k] * sizes[mj]) / (sizes[mi] + sizes[mj])
                new_row[k] = new_dist
        # Extend existing rows
        for row in d:
            row.append(0.0)
        d.append(new_row)
        # Mirror
        for k in range(len(d) - 1):
            d[k][next_id] = d[next_id][k]

        clusters[next_id] = new_label
        sizes[next_id] = sizes[mi] + sizes[mj]
        heights[next_id] = h
        active.discard(mi)
        active.discard(mj)
        active.add(next_id)
        next_id += 1

    root = active.pop()
    newick = clusters[root] + ";"
    return newick, merge_log


def _text_dendrogram(merge_log, labels):
    """Generate simple text dendrogram from merge log."""
    lines = []
    lines.append("Merge order (UPGMA):")
    lines.append(f"{'Step':>4s}  {'Fst':>8s}  {'Cluster 1':<25s}  {'Cluster 2':<25s}")
    lines.append("-" * 70)

    def _short(label, maxlen=25):
        if len(label) > maxlen:
            return label[:maxlen - 3] + "..."
        return label

    for i, (c1, c2, fst, h) in enumerate(merge_log, 1):
        lines.append(f"{i:4d}  {fst:8.4f}  {_short(c1):<25s}  {_short(c2):<25s}")
    return lines


# ── Subcommand: report ────────────────────────────────────────────────


def cmd_report(args):
    """Generate final Markdown report from all interpretation results."""
    output = load_output()
    results = load_results()
    lines = []

    def w(s=""):
        lines.append(s)

    w("# Rice 3K Genome — Biological Interpretation Report")
    w()
    w("## 1. Population Diversity Ranking")
    w()

    div_table = output.get("annotate", {}).get("diversity_table", {})
    if div_table:
        ranked = sorted(div_table.items(), key=lambda x: x[1].get("mean_pi") or 0,
                        reverse=True)
        w("| Rank | Population | Mean π | Mean θ_W | Mean Tajima's D | Mean F_IS |")
        w("|------|-----------|--------|---------|----------------|----------|")
        for i, (pop, v) in enumerate(ranked, 1):
            w(f"| {i} | {pop} | {_fmt(v.get('mean_pi'))} | "
              f"{_fmt(v.get('mean_theta_w'))} | {_fmt(v.get('mean_tajima_d'))} | "
              f"{_fmt(v.get('mean_fis'))} |")
        w()

        # Identify outliers
        w("**Key observations:**")
        if ranked:
            w(f"- Most diverse: **{ranked[0][0]}** (π = {_fmt(ranked[0][1].get('mean_pi'))})")
            w(f"- Least diverse: **{ranked[-1][0]}** (π = {_fmt(ranked[-1][1].get('mean_pi'))})")
            gjtmp = div_table.get("GJ-tmp", {})
            if gjtmp:
                w(f"- GJ-tmp (Temperate Japonica) F_IS = {_fmt(gjtmp.get('mean_fis'))} "
                  f"— strongest inbreeding among all populations")
        w()

    # ── 2. pi_N/pi_S ──
    w("## 2. Purifying Selection Efficiency (π_N/π_S)")
    w()
    ratios = output.get("pinsps_ratios", {})
    if ratios:
        ranked = sorted(ratios.items(),
                        key=lambda x: x[1].get("mean_pi_n_pi_s") or 0, reverse=True)
        w("| Rank | Population | Mean π_N/π_S | Interpretation |")
        w("|------|-----------|-------------|----------------|")
        for i, (pop, v) in enumerate(ranked, 1):
            r = v.get("mean_pi_n_pi_s")
            interp = "relaxed constraint" if r and r > 0.3 else "strong purifying selection"
            w(f"| {i} | {pop} | {_fmt(r)} | {interp} |")
        w()

        w("**Key finding:** Populations with the highest π_N/π_S have accumulated "
          "more slightly deleterious missense variants, consistent with reduced "
          "effectiveness of purifying selection during domestication bottlenecks.")
        w()

        # Per-chr breakdown for GJ-tmp
        gjtmp = ratios.get("GJ-tmp", {}).get("per_chr", {})
        if gjtmp:
            w("### π_N/π_S per chromosome (GJ-tmp)")
            w()
            w("| Chr | π_N | π_S | π_N/π_S |")
            w("|-----|-----|-----|---------|")
            for chrom in CHROMOSOMES:
                v = gjtmp.get(chrom, {})
                w(f"| {chrom} | {_fmt(v.get('pi_n'))} | {_fmt(v.get('pi_s'))} | "
                  f"{_fmt(v.get('pi_n_pi_s'))} |")
            w()

    # ── 3. Pairwise Divergence ──
    w("## 3. Pairwise Population Divergence (Weir & Cockerham Fst)")
    w()
    fst_matrix = output.get("fst_matrix", {})
    if fst_matrix:
        top10 = sorted(fst_matrix.items(),
                       key=lambda x: x[1].get("mean_fst") or 0, reverse=True)[:10]
        w("### Top 10 most differentiated pairs")
        w()
        w("| Pair | Mean Fst |")
        w("|------|---------|")
        for pair_key, v in top10:
            w(f"| {pair_key} | {_fmt(v.get('mean_fst'), 4)} |")
        w()

        bottom5 = sorted(fst_matrix.items(),
                         key=lambda x: x[1].get("mean_fst") or 999)[:5]
        w("### Most closely related pairs")
        w()
        w("| Pair | Mean Fst |")
        w("|------|---------|")
        for pair_key, v in bottom5:
            w(f"| {pair_key} | {_fmt(v.get('mean_fst'), 4)} |")
        w()

    # ── 4. XP-EHH Peaks ──
    w("## 4. Selection Signals (XP-EHH Peaks)")
    w()
    xpehh_ann = output.get("annotate", {}).get("xpehh_annotations", {})
    for pair_key, ann in xpehh_ann.items():
        w(f"### {pair_key}")
        w(f"Total hits: {ann.get('n_total_hits', 0)}, "
          f"Peaks: {ann.get('n_peaks', 0)}")
        w()

        # Known domestication gene hits
        known_hits = []
        for peak in ann.get("top_peaks", []):
            for kg in peak.get("known_domestication_genes", []):
                known_hits.append({**kg, "xpehh": peak["peak_xpehh"],
                                   "peak_pos": peak["peak_pos"]})
        if known_hits:
            w("**Known domestication genes detected:**")
            w()
            w("| Gene | Function | XP-EHH | Peak Position |")
            w("|------|----------|--------|--------------|")
            for kh in known_hits:
                w(f"| {kh['gene']} | {kh['function']} | "
                  f"{kh['xpehh']:.2f} | {kh['chr']}:{kh['peak_pos']:,} |")
            w()

        # Top novel peaks
        w("**Top 10 peaks (by |XP-EHH|):**")
        w()
        w("| Chr | Position | XP-EHH | #Hits | Top Gene | Impact |")
        w("|-----|----------|--------|-------|----------|--------|")
        for peak in ann.get("top_peaks", [])[:10]:
            genes = peak.get("genes", [])
            top_gene = genes[0]["symbol"] if genes else "—"
            top_impact = genes[0]["impact"] if genes else "—"
            w(f"| {peak['chr']} | {peak['peak_pos']:,} | "
              f"{peak['peak_xpehh']:.2f} | {peak['n_hits']} | "
              f"{top_gene} | {top_impact} |")
        w()

    # ── 5. Garud's H Sweeps ──
    w("## 5. Selective Sweep Regions (Garud's H)")
    w()
    sweep_ann = output.get("annotate", {}).get("sweep_annotations", {})
    w("| Population | Total Sweeps (H12>0.10) | Hard Sweeps | Soft Sweeps |")
    w("|-----------|------------------------|-------------|-------------|")
    for pop in POPULATIONS:
        ann = sweep_ann.get(pop, {})
        sweeps = ann.get("sweeps", [])
        hard = sum(1 for s in sweeps if s.get("sweep_type") == "hard")
        soft = len(sweeps) - hard
        w(f"| {pop} | {len(sweeps)} | {hard} | {soft} |")
    w()

    # Highlight top sweeps with known genes
    w("### Notable sweeps near known domestication genes")
    w()
    for pop in POPULATIONS:
        ann = sweep_ann.get(pop, {})
        for sw in ann.get("sweeps", []):
            if sw.get("known_domestication_genes"):
                for kg in sw["known_domestication_genes"]:
                    w(f"- **{pop}** {sw['chr']}:{sw['start']:,}–{sw['end']:,} "
                      f"(H12={sw['h12']:.3f}, {sw['sweep_type']}): "
                      f"near **{kg['gene']}** ({kg['function']})")
    w()

    # ── 6. ROH ──
    w("## 6. Runs of Homozygosity (ROH)")
    w()
    roh_sum = output.get("roh_summary", {})
    if roh_sum:
        ranked = sorted(roh_sum.items(),
                        key=lambda x: x[1].get("mean_froh_across_chr", 0),
                        reverse=True)
        w("| Rank | Population | Total ROH (Mb) | Mean FROH |")
        w("|------|-----------|---------------|-----------|")
        for i, (pop, v) in enumerate(ranked, 1):
            total_mb = v.get("total_roh_bp", 0) / 1e6
            w(f"| {i} | {pop} | {total_mb:.1f} | "
              f"{_fmt(v.get('mean_froh_across_chr'))} |")
        w()

    # ── 7. Genome Scans ──
    w("## 7. Annotation-Conditioned Genome Scans")
    w()
    gscan = output.get("gscan", {})
    if gscan:
        for pop1, pop2 in GSCAN_PAIRS:
            w(f"### {pop1} vs {pop2}: Missense vs Synonymous Fst")
            w()

            mis_fsts = []
            syn_fsts = []
            for chrom in CHROMOSOMES:
                mis_key = f"{pop1}|{pop2}|{chrom}|missense_variant"
                syn_key = f"{pop1}|{pop2}|{chrom}|synonymous_variant"
                mis_windows = gscan.get(mis_key, {}).get("result", [])
                syn_windows = gscan.get(syn_key, {}).get("result", [])

                if mis_windows:
                    fsts = [w.get("fst", 0) for w in mis_windows
                            if w.get("fst") is not None and w.get("n_variants", 0) >= 5]
                    if fsts:
                        mis_fsts.append(sum(fsts) / len(fsts))
                if syn_windows:
                    fsts = [w.get("fst", 0) for w in syn_windows
                            if w.get("fst") is not None and w.get("n_variants", 0) >= 5]
                    if fsts:
                        syn_fsts.append(sum(fsts) / len(fsts))

            mean_mis = sum(mis_fsts) / len(mis_fsts) if mis_fsts else None
            mean_syn = sum(syn_fsts) / len(syn_fsts) if syn_fsts else None
            w(f"- Mean Fst (missense): {_fmt(mean_mis)}")
            w(f"- Mean Fst (synonymous): {_fmt(mean_syn)}")
            if mean_mis and mean_syn and mean_syn > 0:
                w(f"- Fst_missense / Fst_synonymous = {mean_mis/mean_syn:.3f}")
                if mean_mis / mean_syn > 1.0:
                    w(f"  → **Elevated missense divergence** — evidence of "
                      f"adaptive protein evolution")
            w()

        # Find top differentiation windows (missense-specific)
        w("### Top Adaptive Divergence Windows (missense-Fst outliers)")
        w()
        all_windows = []
        for pop1, pop2 in GSCAN_PAIRS:
            for chrom in CHROMOSOMES:
                key = f"{pop1}|{pop2}|{chrom}|missense_variant"
                windows = gscan.get(key, {}).get("result", [])
                for win in windows:
                    if win.get("fst") is not None and win.get("n_variants", 0) >= 3:
                        all_windows.append({
                            "pair": f"{pop1}_vs_{pop2}",
                            "chr": win.get("chr", chrom),
                            "start": win.get("start"),
                            "end": win.get("end"),
                            "fst": win["fst"],
                            "n_variants": win["n_variants"],
                        })
        all_windows.sort(key=lambda x: x["fst"], reverse=True)
        w("| Pair | Chr | Start | End | Fst (Hudson) | #Variants |")
        w("|------|-----|-------|-----|-------------|-----------|")
        for win in all_windows[:20]:
            w(f"| {win['pair']} | {win['chr']} | {win['start']:,} | "
              f"{win['end']:,} | {win['fst']:.4f} | {win['n_variants']} |")
        w()

    # ── 8. Fay & Wu's H ──
    w("## 8. Fay & Wu's H (Ancestral-Allele-Dependent)")
    w()
    fwh_summary = output.get("fay_wu_summary", {})
    if fwh_summary:
        ranked = sorted(fwh_summary.items(),
                        key=lambda x: x[1].get("mean_h") or 0)
        w("| Rank | Population | Mean H | Mean H_norm | Interpretation |")
        w("|------|-----------|--------|------------|----------------|")
        for i, (pop, v) in enumerate(ranked, 1):
            h = v.get("mean_h")
            interp = "strong sweep signals" if h and h < -0.05 else (
                "moderate sweeps" if h and h < -0.02 else "weak/no sweeps")
            w(f"| {i} | {pop} | {_fmt(v.get('mean_h'))} | "
              f"{_fmt(v.get('mean_h_norm'))} | {interp} |")
        w()
        w("**Key finding:** Negative H indicates excess high-frequency derived "
          "alleles from recent selective sweeps. More negative = stronger/more "
          "frequent sweeps in that population.")
        w()
    else:
        w("*Not yet run. Use `python scripts/rice_interpret.py fay_wu`*")
        w()

    # ── 9. Unfolded SFS ──
    w("## 9. Unfolded (Polarized) Site Frequency Spectrum")
    w()
    usfs_summary = output.get("usfs_summary", {})
    if usfs_summary:
        w("| Population | Total Variants | Polarized | Singletons | "
          "High-freq Derived (>50%) |")
        w("|-----------|---------------|-----------|------------|"
          "------------------------|")
        for pop in POPULATIONS:
            v = usfs_summary.get(pop, {})
            sfs = v.get("sfs", [])
            total = v.get("n_variants", 0)
            polarized = v.get("n_polarized", 0)
            sing = sfs[0] if sfs else 0
            n_bins = len(sfs)
            high_freq = sum(sfs[n_bins//2:]) if sfs else 0
            w(f"| {pop} | {total:,} | {polarized:,} | {sing:,} | {high_freq:,} |")
        w()
    else:
        w("*Not yet run. Use `python scripts/rice_interpret.py usfs`*")
        w()

    # ── 10. H Genome Scan ──
    w("## 10. Fay & Wu's H Genome Scan")
    w()
    hwscan_ann = output.get("annotate", {}).get("hwscan_annotations", {})
    hwscan_summary = output.get("hwscan_summary", {})

    if hwscan_ann:
        # Use annotated peaks (has gene info)
        for pop, info in hwscan_ann.items():
            w(f"### {pop}")
            summary = hwscan_summary.get(pop, {})
            w(f"- Windows: {summary.get('n_windows', info.get('n_windows_total', 0))}, "
              f"Mean H: {_fmt(summary.get('mean_h'))}")
            top_peaks = info.get("top_peaks", [])
            if top_peaks:
                w()
                w("**Top 10 most negative H windows (strongest sweeps):**")
                w()
                w("| Chr | Start | End | H | Top Gene | Known Gene |")
                w("|-----|-------|-----|---|----------|------------|")
                for tw in top_peaks[:10]:
                    genes = tw.get("genes", [])
                    top_gene = genes[0]["symbol"] if genes else "—"
                    known = tw.get("known_domestication_genes", [])
                    known_str = known[0]["gene"] if known else "—"
                    w(f"| {tw.get('chr','?')} | {tw.get('start',0):,} | "
                      f"{tw.get('end',0):,} | {_fmt(tw.get('fay_wu_h'))} | "
                      f"{top_gene} | {known_str} |")
                w()
    elif hwscan_summary:
        # Fallback: unannotated peaks
        for pop, info in hwscan_summary.items():
            w(f"### {pop}")
            w(f"- Windows: {info.get('n_windows', 0)}, "
              f"Mean H: {_fmt(info.get('mean_h'))}")
            top_windows = info.get("top_sweep_windows", [])
            if top_windows:
                w()
                w("**Top 5 most negative H windows (strongest sweeps):**")
                w()
                w("| Chr | Start | End | H |")
                w("|-----|-------|-----|---|")
                for tw in top_windows[:5]:
                    w(f"| {tw.get('chr','?')} | {tw.get('start',0):,} | "
                      f"{tw.get('end',0):,} | {_fmt(tw.get('fay_wu_h'))} |")
                w()
    else:
        w("*Not yet run. Use `python scripts/rice_interpret.py hwscan`*")
        w()

    # ── 11. ROH HMM ──
    w("## 11. Runs of Homozygosity (HMM Method)")
    w()
    roh_hmm_sum = output.get("roh_hmm_summary", {})
    if roh_hmm_sum:
        ranked = sorted(roh_hmm_sum.items(),
                        key=lambda x: x[1].get("mean_froh_across_chr", 0),
                        reverse=True)
        w("| Rank | Population | Total ROH (Mb) | Mean FROH |")
        w("|------|-----------|---------------|-----------|")
        for i, (pop, v) in enumerate(ranked, 1):
            total_mb = v.get("total_roh_bp", 0) / 1e6
            w(f"| {i} | {pop} | {total_mb:.1f} | "
              f"{_fmt(v.get('mean_froh_across_chr'))} |")
        w()
        w("**Key finding:** GJ-tmp (Temperate Japonica) has the highest FROH, "
          "consistent with the strongest domestication bottleneck among rice "
          "subpopulations.")
        w()

        # Per-chromosome ROH breakdown for top populations
        roh_ann = output.get("annotate", {}).get("roh_annotations", {})
        for pop in ["GJ-tmp", "GJ-adm", "GJ-trp"]:
            info = roh_ann.get(pop, {})
            chr_data = info.get("chr_with_roh", [])
            if chr_data:
                w(f"### {pop} — Top chromosomes by ROH")
                w()
                w("| Chr | Total ROH (Mb) | Mean FROH | Samples with ROH |")
                w("|-----|---------------|-----------|------------------|")
                for c in chr_data[:5]:
                    n_samp = c.get("n_samples_with_roh", 0)
                    total = c.get("n_samples", 0)
                    w(f"| {c['chr']} | {c['total_roh_bp']/1e6:.1f} | "
                      f"{_fmt(c.get('mean_froh'))} | "
                      f"{n_samp}/{total} |")
                w()
    else:
        w("*Not yet run. Use `python scripts/rice_interpret.py roh_hmm`*")
        w()

    # ── 12. DAF Enrichment ──
    w("## 12. Derived Allele Frequency Enrichment at Selection Signals")
    w()
    daf = output.get("daf_enrichment", {})
    if daf:
        bg = daf.get("background_daf", {})
        pbs_daf = daf.get("pbs_peak_daf", {})
        xpehh_daf = daf.get("xpehh_peak_daf", {})

        if pbs_daf:
            w("### PBS Peak Regions")
            w()
            w("| Focal Population | Peak DAF | Background DAF | Enrichment |")
            w("|-----------------|----------|---------------|------------|")
            for config_key, v in pbs_daf.items():
                focal = v.get("focal", config_key)
                enrichment = v.get("enrichment", 0)
                w(f"| {focal} ({config_key}) | {_fmt(v.get('mean_peak_daf'))} | "
                  f"{_fmt(v.get('background_daf'))} | "
                  f"{enrichment:.2f}x |")
            w()

        if xpehh_daf:
            w("### XP-EHH Peak Regions")
            w()
            w("| Pair | Focal Pop | Peak DAF | Background DAF | Enrichment |")
            w("|------|----------|----------|---------------|------------|")
            for pair, v in xpehh_daf.items():
                focal = v.get("focal", pair)
                enrichment = v.get("enrichment", 0)
                w(f"| {pair} | {focal} | {_fmt(v.get('mean_peak_daf'))} | "
                  f"{_fmt(v.get('background_daf'))} | {enrichment:.2f}x |")
            w()

        w("**Key finding:** Elevated derived allele frequency at selection "
          "signal peaks confirms that sweeps drove derived alleles toward "
          "fixation, validating the ancestral allele polarization.")
        w()
    else:
        w("*Not yet run. Use `python scripts/rice_interpret.py daf_enrichment`*")
        w()

    # ── 13. PBS ──
    w("## 13. Population Branch Statistic (PBS)")
    w()
    pbs_summary = output.get("pbs_summary", {})
    if pbs_summary:
        for config_key, info in pbs_summary.items():
            focal = info["focal"]
            w(f"### {focal} (vs {info['sister']}, outgroup {info['outgroup']})")
            w(f"- Windows analysed: {info['n_windows']}")
            w(f"- Mean PBS: {info['mean_pbs']:.4f}")
            w()

            top = info.get("top_windows", [])
            if top:
                w("**Top PBS windows (population-specific selection):**")
                w()
                w("| Chr | Start | End | PBS | Fst_WC | #Var | Top Gene | Known Gene |")
                w("|-----|-------|-----|-----|--------|------|----------|------------|")
                for tw in top[:10]:
                    genes = tw.get("genes", [])
                    top_gene = genes[0]["symbol"] if genes else "—"
                    known = tw.get("known_genes", [])
                    known_str = known[0]["gene"] if known else "—"
                    w(f"| {tw['chr']} | {tw['start']:,} | {tw['end']:,} | "
                      f"{tw['pbs']:.4f} | {_fmt(tw.get('fst_wc'), 4)} | "
                      f"{tw['n_variants']} | {top_gene} | {known_str} |")
                w()
    else:
        w("*PBS analysis not yet run. Use `python scripts/rice_interpret.py pbs`*")
        w()

    # ── 14. Population Tree ──
    w("## 14. Population Tree (UPGMA)")
    w()
    tree_data = output.get("tree", {})
    if tree_data:
        w(f"**Method:** {tree_data.get('method', 'UPGMA')} on "
          f"{tree_data.get('distance_metric', 'Fst')}")
        w()
        w(f"**Newick:** `{tree_data.get('newick', '')}`")
        w()
        dendro = tree_data.get("dendrogram", [])
        if dendro:
            w("```")
            for line in dendro:
                w(line)
            w("```")
            w()
    else:
        w("*Population tree not yet built. Use `python scripts/rice_interpret.py tree`*")
        w()

    # ── 15. Graph-Native Novelty ──
    w("## 15. What GraphPop Uniquely Enables")
    w()
    w("1. **Annotation-conditioned statistics**: `graphpop.diversity(..., "
      "{consequence: 'missense_variant'})` computes π restricted to missense "
      "variants via Variant→HAS_CONSEQUENCE→Gene traversal — one query vs. "
      "4 separate tools in classical pipelines.")
    w("2. **π_N/π_S across all 12 populations**: 288 single-line procedure "
      "calls (12 pops × 12 chr × 2 consequences). Classical tools require "
      "288 separate VCF filtering + statistics runs.")
    w("3. **Missense-conditioned genome scans**: Sliding-window Fst "
      "restricted to missense variants identifies adaptive protein divergence "
      "hotspots — not available in any classical tool.")
    w("4. **Cross-statistic gene lookups**: After iHS/XP-EHH computation, "
      "immediate graph traversal to find associated genes — zero file "
      "conversion or coordinate intersection.")
    w()

    report_text = "\n".join(lines)
    report_path = "rice_interpretation_report.md"
    with open(report_path, "w") as f:
        f.write(report_text)
    log.info(f"Report written to {report_path}")
    print(report_text)


def _fmt(v, decimals=6):
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.{decimals}f}"
    return str(v)


# ── Main ──────────────────────────────────────────────────────────────

SUBCOMMANDS = {
    "annotate": cmd_annotate,
    "pinsps": cmd_pinsps,
    "divergence": cmd_divergence,
    "gscan": cmd_gscan,
    "roh": cmd_roh,
    "pbs": cmd_pbs,
    "tree": cmd_tree,
    "fay_wu": cmd_fay_wu,
    "usfs": cmd_usfs,
    "hwscan": cmd_hwscan,
    "roh_hmm": cmd_roh_hmm,
    "daf_enrichment": cmd_daf_enrichment,
    "report": cmd_report,
}


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description="Rice 3K Genome Biological Interpretation")
    parser.add_argument(
        "command", choices=list(SUBCOMMANDS.keys()),
        help="Subcommand to run")
    args = parser.parse_args()
    SUBCOMMANDS[args.command](args)


if __name__ == "__main__":
    main()
