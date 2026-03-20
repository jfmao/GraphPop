#!/usr/bin/env python3
"""1000 Genomes Biological Interpretation Script.

Cross-references population genetics signals with gene annotations
using GraphPop's graph-native annotation-aware analyses.

Usage:
    python scripts/human_interpret.py annotate       # Cross-reference peaks with genes
    python scripts/human_interpret.py pinsps         # pi_N/pi_S per population
    python scripts/human_interpret.py divergence     # Full pairwise Fst matrix
    python scripts/human_interpret.py gscan          # Genome scans with annotation
    python scripts/human_interpret.py pbs            # PBS genome scans
    python scripts/human_interpret.py tree           # UPGMA population tree
    python scripts/human_interpret.py fay_wu         # Fay & Wu's H
    python scripts/human_interpret.py usfs           # Unfolded SFS
    python scripts/human_interpret.py hwscan         # H genome scan
    python scripts/human_interpret.py roh_hmm        # ROH HMM
    python scripts/human_interpret.py daf_enrichment # DAF at peaks
    python scripts/human_interpret.py report         # Generate Markdown report
    python scripts/human_interpret.py hap_stats      # Joint nSL+iHS+XP-EHH from hap cache
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

import numpy as np
import genome_scan_numpy

# ── Constants ─────────────────────────────────────────────────────────

NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"
NEO4J_DB = "neo4j"

RESULTS_FILE = "human_full_results.json"
OUTPUT_FILE = "human_interpretation_results.json"

POPULATIONS = [
    "YRI", "LWK", "GWD", "MSL", "ESN", "ACB", "ASW",
    "CEU", "TSI", "FIN", "GBR", "IBS",
    "CHB", "JPT", "CHS", "CDX", "KHV",
    "GIH", "PJL", "BEB", "STU", "ITU",
    "MXL", "PUR", "CLM", "PEL",
]

SUPERPOPS = {
    "AFR": ["YRI", "LWK", "GWD", "MSL", "ESN", "ACB", "ASW"],
    "EUR": ["CEU", "TSI", "FIN", "GBR", "IBS"],
    "EAS": ["CHB", "JPT", "CHS", "CDX", "KHV"],
    "SAS": ["GIH", "PJL", "BEB", "STU", "ITU"],
    "AMR": ["MXL", "PUR", "CLM", "PEL"],
}

CHROMOSOMES = [f"chr{i}" for i in range(1, 23)]

XPEHH_PAIRS = [
    ("YRI", "CEU"), ("YRI", "CHB"), ("CEU", "CHB"),
    ("YRI", "GIH"), ("CEU", "JPT"),
    ("YRI", "LWK"), ("YRI", "GWD"), ("ESN", "GWD"),
    ("CEU", "TSI"), ("CEU", "FIN"), ("GBR", "IBS"),
    ("CHB", "JPT"), ("CHB", "CHS"), ("CHB", "CDX"),
    ("GIH", "BEB"), ("PJL", "STU"),
    ("MXL", "PUR"), ("CLM", "PEL"),
]

GSCAN_PAIRS = [
    ("YRI", "CEU"), ("CEU", "CHB"), ("YRI", "CHB"),
]

PBS_CONFIGS = [
    ("CEU", "GBR", "YRI"),
    ("YRI", "LWK", "CEU"),
    ("CHB", "JPT", "YRI"),
    ("GIH", "BEB", "YRI"),
    ("PEL", "MXL", "CHB"),
    ("FIN", "CEU", "CHB"),
]

H_SCAN_POPS = ["YRI", "CEU", "CHB", "GIH", "PEL", "FIN"]

# Known human selection genes
KNOWN_GENES = [
    {"gene": "LCT", "chr": "chr2", "pos": 136608646, "function": "Lactase persistence"},
    {"gene": "SLC24A5", "chr": "chr15", "pos": 48426484, "function": "Skin pigmentation (light)"},
    {"gene": "SLC45A2", "chr": "chr5", "pos": 33951693, "function": "Skin pigmentation (light)"},
    {"gene": "EDAR", "chr": "chr2", "pos": 109513601, "function": "Hair/teeth morphology (EAS)"},
    {"gene": "EPAS1", "chr": "chr2", "pos": 46588376, "function": "High-altitude adaptation (Tibetan)"},
    {"gene": "DARC", "chr": "chr1", "pos": 159204875, "function": "Duffy blood group (malaria resist.)"},
    {"gene": "HBB", "chr": "chr11", "pos": 5227002, "function": "Sickle cell (malaria resist.)"},
    {"gene": "G6PD", "chr": "chrX", "pos": 154531391, "function": "G6PD deficiency (malaria resist.)"},
    {"gene": "HERC2", "chr": "chr15", "pos": 28356859, "function": "Eye color (blue eyes, EUR)"},
    {"gene": "OCA2", "chr": "chr15", "pos": 28000020, "function": "Skin/eye pigmentation"},
    {"gene": "MC1R", "chr": "chr16", "pos": 89919709, "function": "Red hair/fair skin"},
    {"gene": "KITLG", "chr": "chr12", "pos": 88886568, "function": "Skin pigmentation"},
    {"gene": "TYRP1", "chr": "chr9", "pos": 12699268, "function": "Skin pigmentation"},
    {"gene": "ADH1B", "chr": "chr4", "pos": 99318162, "function": "Alcohol metabolism (EAS)"},
    {"gene": "ALDH2", "chr": "chr12", "pos": 111803962, "function": "Alcohol metabolism (EAS)"},
    {"gene": "FADS1", "chr": "chr11", "pos": 61567099, "function": "Fatty acid desaturation"},
    {"gene": "FADS2", "chr": "chr11", "pos": 61595006, "function": "Fatty acid desaturation"},
    {"gene": "ABCC11", "chr": "chr16", "pos": 48224287, "function": "Earwax type (EAS)"},
    {"gene": "TLR1", "chr": "chr4", "pos": 38799307, "function": "Innate immunity"},
    {"gene": "TLR6", "chr": "chr4", "pos": 38830268, "function": "Innate immunity"},
    {"gene": "HLA-A", "chr": "chr6", "pos": 29942470, "function": "Adaptive immunity (MHC)"},
    {"gene": "HLA-B", "chr": "chr6", "pos": 31353872, "function": "Adaptive immunity (MHC)"},
    {"gene": "HLA-DRB1", "chr": "chr6", "pos": 32578769, "function": "Adaptive immunity (MHC)"},
    {"gene": "LARGE1", "chr": "chr22", "pos": 33532074, "function": "Lassa fever resistance (AFR)"},
]

WINDOW_SIZE = 100_000
WINDOW_STEP = 50_000

log = logging.getLogger("human_interpret")

# ── Helpers ───────────────────────────────────────────────────────────


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stderr),
            logging.FileHandler("human_interpret.log", mode="a"),
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


def find_nearby_known_genes(chrom, start, end, margin=100000):
    result = []
    for kg in KNOWN_GENES:
        if kg["chr"] == chrom and start - margin <= kg["pos"] <= end + margin:
            result.append(kg)
    return result


def call_peaks(hits, merge_distance=200000):
    """Merge nearby hits into peak regions."""
    by_chr = defaultdict(list)
    for h in hits:
        vid = h.get("variantId", "")
        chrom = vid.split(":")[0] if ":" in vid else ""
        by_chr[chrom].append(h)

    peaks = []
    for chrom, chr_hits in by_chr.items():
        chr_hits.sort(key=lambda x: x.get("pos", 0))
        cur_start = chr_hits[0].get("pos", 0)
        cur_end = cur_start
        cur_hits = [chr_hits[0]]

        for h in chr_hits[1:]:
            if h.get("pos", 0) - cur_end <= merge_distance:
                cur_end = h["pos"]
                cur_hits.append(h)
            else:
                best = max(cur_hits, key=lambda x: abs(x.get("xpehh", 0)))
                peaks.append({
                    "chr": chrom,
                    "start": cur_start,
                    "end": cur_end,
                    "peak_xpehh": best.get("xpehh", 0),
                    "peak_pos": best.get("pos", 0),
                    "peak_variantId": best.get("variantId", ""),
                    "n_hits": len(cur_hits),
                })
                cur_start = h.get("pos", 0)
                cur_end = cur_start
                cur_hits = [h]

        best = max(cur_hits, key=lambda x: abs(x.get("xpehh", 0)))
        peaks.append({
            "chr": chrom,
            "start": cur_start,
            "end": cur_end,
            "peak_xpehh": best.get("xpehh", 0),
            "peak_pos": best.get("pos", 0),
            "peak_variantId": best.get("variantId", ""),
            "n_hits": len(cur_hits),
        })

    peaks.sort(key=lambda x: abs(x["peak_xpehh"]), reverse=True)
    return peaks


def _impact_rank(impact):
    return {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}.get(impact, 4)


# ── Subcommand: annotate ──────────────────────────────────────────────


def cmd_annotate(args):
    """Cross-reference XP-EHH peaks and Garud's H sweeps with genes."""
    results = load_results()
    output = load_output()
    driver = get_driver()

    with driver.session(database=NEO4J_DB) as session:
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

        # 1. XP-EHH peaks → gene annotation
        log.info("=== Annotating XP-EHH peaks ===")
        xpehh_annotations = {}

        for pair_key_str, pair_data in results.get("phase2_xpehh", {}).items():
            all_hits = []
            for chrom_data in pair_data.values():
                if isinstance(chrom_data, dict):
                    all_hits.extend(chrom_data.get("top_hits", []))

            if not all_hits:
                continue

            peaks = call_peaks(all_hits)
            log.info("  %s: %d hits -> %d peaks", pair_key_str, len(all_hits), len(peaks))

            top_peaks = peaks[:100]
            for i, peak in enumerate(top_peaks):
                margin = 100000
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
                seen = {}
                for g in genes:
                    key = g["symbol"]
                    if key not in seen or _impact_rank(g["impact"]) < _impact_rank(seen[key]["impact"]):
                        seen[key] = g
                peak["genes"] = list(seen.values())[:30]
                peak["n_genes"] = len(seen)
                peak["known_selection_genes"] = find_nearby_known_genes(
                    peak["chr"], peak["start"], peak["end"]
                )
                if (i + 1) % 20 == 0:
                    log.info("    annotated %d/%d peaks", i + 1, len(top_peaks))

            xpehh_annotations[pair_key_str] = {
                "n_total_hits": len(all_hits),
                "n_peaks": len(peaks),
                "top_peaks": top_peaks,
            }

        # 2. Garud's H sweeps
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
                sweeps = [w for w in windows if w.get("h12", 0) > 0.05]
                for sw in sweeps:
                    recs = list(session.run(
                        sweep_gene_query,
                        chr=sw.get("chr", chrom), start=sw.get("start", 0),
                        end=sw.get("end", 0),
                    ))
                    genes = [
                        {"symbol": r["symbol"], "geneId": r["geneId"],
                         "consequence": r["consequence"],
                         "n_variants": r["n_variants"]}
                        for r in recs
                    ]
                    pop_sweeps.append({
                        "chr": sw.get("chr", chrom),
                        "start": sw.get("start"),
                        "end": sw.get("end"),
                        "h1": sw.get("h1"),
                        "h12": sw.get("h12"),
                        "h2_h1": sw.get("h2_h1"),
                        "hap_diversity": sw.get("hap_diversity"),
                        "n_haplotypes": sw.get("n_haplotypes"),
                        "n_variants": sw.get("n_variants"),
                        "genes": genes,
                        "known_selection_genes": find_nearby_known_genes(
                            sw.get("chr", chrom), sw.get("start", 0), sw.get("end", 0)
                        ),
                        "sweep_type": "hard" if sw.get("h2_h1", 1) < 0.05 else "soft",
                    })

            pop_sweeps.sort(key=lambda x: x.get("h12", 0), reverse=True)
            sweep_annotations[pop] = {
                "n_sweeps": len(pop_sweeps),
                "sweeps": pop_sweeps,
            }
            log.info("  %s: %d sweep windows annotated", pop, len(pop_sweeps))

        # 3. Diversity comparison table
        log.info("=== Building diversity comparison tables ===")
        diversity_table = _build_diversity_table(results)

    driver.close()

    output["annotate"] = {
        "xpehh_annotations": xpehh_annotations,
        "sweep_annotations": sweep_annotations,
        "diversity_table": diversity_table,
    }
    save_output(output)
    log.info("Annotation results saved to %s", OUTPUT_FILE)


def _build_diversity_table(results):
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
    """Compute pi_N/pi_S per population — numpy bypass (same .npz cache as gscan)."""
    output   = load_output()
    existing = output.get("pinsps", {})

    n_total = len(POPULATIONS) * len(CHROMOSOMES) * 2
    n_done  = sum(1 for k in existing if "error" not in existing.get(k, {}))
    log.info("=== Computing pi_N/pi_S: %d calls (%d cached) ===", n_total, n_done)

    with open(os.path.join(GSCAN_CACHE_DIR, "pop_meta.json")) as f:
        pop_ids = json.load(f)["pop_ids"]
    pop_idx_map = {p: pop_ids.index(p) for p in POPULATIONS}

    for csq in ("missense_variant", "synonymous_variant"):
        chroms_needed = sorted({
            c for pop in POPULATIONS for c in CHROMOSOMES
            if f"{pop}|{c}|{csq}" not in existing
            or "error" in existing.get(f"{pop}|{c}|{csq}", {})
        })
        if not chroms_needed:
            log.info("  %s: all cached", csq)
            continue

        _ensure_variant_cache(chroms_needed, GSCAN_CACHE_DIR)
        for chrom in chroms_needed:
            t0 = time.time()
            new = _compute_diversity_chrom(chrom, pop_idx_map, csq, GSCAN_CACHE_DIR)
            existing.update(new)
            n_done += len(POPULATIONS)
            log.info("  %s %s: %d pops (%.1fs)", csq[:3], chrom,
                     len(POPULATIONS), time.time() - t0)
        output["pinsps"] = existing
        save_output(output)

    # Compute ratios
    ratios = {}
    for pop in POPULATIONS:
        pop_ratios = {}
        for chrom in CHROMOSOMES:
            mis = existing.get(f"{pop}|{chrom}|missense_variant",  {}).get("result", {})
            syn = existing.get(f"{pop}|{chrom}|synonymous_variant", {}).get("result", {})
            pi_n = mis.get("pi")
            pi_s = syn.get("pi")
            pop_ratios[chrom] = {
                "pi_n": pi_n, "pi_s": pi_s,
                "pi_n_pi_s": pi_n / pi_s if pi_n is not None and pi_s and pi_s > 0 else None,
                "n_missense":   mis.get("n_segregating"),
                "n_synonymous": syn.get("n_segregating"),
            }
        vals = [v["pi_n_pi_s"] for v in pop_ratios.values() if v["pi_n_pi_s"] is not None]
        ratios[pop] = {
            "per_chr":        pop_ratios,
            "mean_pi_n_pi_s": sum(vals) / len(vals) if vals else None,
        }

    output["pinsps"] = existing
    output["pinsps_ratios"] = ratios
    save_output(output)
    log.info("pi_N/pi_S results saved to %s", OUTPUT_FILE)


# ── Subcommand: divergence ────────────────────────────────────────────


def _divergence_one_chrom(chrom, all_pairs, pop_idx, cache_dir):
    """Compute divergence stats for all population pairs on one chromosome.

    Loads the .npz cache once, then computes fst_hudson, fst_wc, dxy, da,
    n_variants for every pair using the same formulas as DivergenceProcedure.java.
    Returns dict {key: {"result": {...}, "wall_sec": float}}.
    """
    t0 = time.time()
    with np.load(os.path.join(cache_dir, f"{chrom}.npz")) as data:
        # Load once — slice per-pair below (float64 for precision)
        ac_all  = data["ac"].astype(np.float64)    # (M, 26)
        an_all  = data["an"].astype(np.float64)
        af_all  = data["af"].astype(np.float64)
        het_all = data["het_count"].astype(np.float64)

    results = {}
    for pop1, pop2 in all_pairs:
        i1, i2 = pop_idx[pop1], pop_idx[pop2]
        ac1, an1, af1, het1 = ac_all[:, i1], an_all[:, i1], af_all[:, i1], het_all[:, i1]
        ac2, an2, af2, het2 = ac_all[:, i2], an_all[:, i2], af_all[:, i2], het_all[:, i2]

        valid = (an1 >= 2) & (an2 >= 2)
        n_var = int(valid.sum())
        key = f"{pop1}|{pop2}|{chrom}"

        if n_var == 0:
            results[key] = {"result": {
                "fst_hudson": 0.0, "fst_wc": 0.0,
                "dxy": 0.0, "da": 0.0, "n_variants": 0,
            }, "wall_sec": 0.0}
            continue

        with np.errstate(invalid="ignore", divide="ignore"):
            # Hudson Fst (ratio-of-sums estimator)
            hw1 = 2.0 * af1 * (1.0 - af1) * an1 / (an1 - 1.0)
            hw2 = 2.0 * af2 * (1.0 - af2) * an2 / (an2 - 1.0)
            hb  = (af1 - af2)**2 - hw1 / (2.0 * an1) - hw2 / (2.0 * an2)
            fst_num = float(np.where(valid, hb, 0.0).sum())
            fst_den = float(np.where(valid, hb + (hw1 + hw2) / 2.0, 0.0).sum())

            # Dxy and pi-within
            dxy_sum = float(np.where(valid,
                af1 * (1.0 - af2) + af2 * (1.0 - af1), 0.0).sum())
            pi1_sum = float(np.where(valid & (an1 > 1),
                2.0 * af1 * (1.0 - af1) * an1 / (an1 - 1.0), 0.0).sum())
            pi2_sum = float(np.where(valid & (an2 > 1),
                2.0 * af2 * (1.0 - af2) * an2 / (an2 - 1.0), 0.0).sum())

        # W&C Fst (reuse _wc_components from genome_scan_numpy)
        a, b, c = genome_scan_numpy._wc_components(
            ac1, an1, het1, ac2, an2, het2, valid)
        wc_denom = float(a.sum() + b.sum() + c.sum())

        dxy = dxy_sum / n_var
        results[key] = {"result": {
            "fst_hudson": fst_num / fst_den if fst_den > 0 else 0.0,
            "fst_wc":     float(a.sum()) / wc_denom if wc_denom > 0 else 0.0,
            "dxy":        dxy,
            "da":         dxy - (pi1_sum + pi2_sum) / (2.0 * n_var),
            "n_variants": n_var,
        }, "wall_sec": time.time() - t0}

    return results


def cmd_divergence(args):
    """Compute full pairwise Fst/Dxy matrix — numpy bypass (same .npz cache as gscan)."""
    output   = load_output()
    existing = output.get("divergence", {})

    all_pairs = list(combinations(POPULATIONS, 2))
    n_total   = len(all_pairs) * len(CHROMOSOMES)
    n_cached  = sum(1 for p1, p2 in all_pairs for c in CHROMOSOMES
                    if f"{p1}|{p2}|{c}" in existing
                    and "error" not in existing.get(f"{p1}|{p2}|{c}", {}))
    log.info("=== Pairwise divergence: %d pairs × %d chroms (%d cached) ===",
             len(all_pairs), len(CHROMOSOMES), n_cached)

    # Find chromosomes that have at least one missing pair
    chroms_needed = sorted({
        c for p1, p2 in all_pairs for c in CHROMOSOMES
        if f"{p1}|{p2}|{c}" not in existing
        or "error" in existing.get(f"{p1}|{p2}|{c}", {})
    })

    if chroms_needed:
        _ensure_variant_cache(chroms_needed, GSCAN_CACHE_DIR)

        with open(os.path.join(GSCAN_CACHE_DIR, "pop_meta.json")) as f:
            pop_ids = json.load(f)["pop_ids"]
        pop_idx = {p: i for i, p in enumerate(pop_ids)}

        for chrom in chroms_needed:
            t0 = time.time()
            new = _divergence_one_chrom(chrom, all_pairs, pop_idx, GSCAN_CACHE_DIR)
            existing.update(new)
            log.info("  %s: %d pairs (%.1fs)", chrom, len(all_pairs), time.time() - t0)
            output["divergence"] = existing
            save_output(output)
    else:
        log.info("  all cached")

    # Build Fst / Dxy matrices
    fst_matrix = {}
    dxy_matrix = {}
    for pop1, pop2 in all_pairs:
        pk = f"{pop1}_vs_{pop2}"
        fst_vals, dxy_vals = [], []
        for chrom in CHROMOSOMES:
            r = existing.get(f"{pop1}|{pop2}|{chrom}", {}).get("result", {})
            if r.get("fst_wc") is not None:
                fst_vals.append(r["fst_wc"])
            if r.get("dxy") is not None:
                dxy_vals.append(r["dxy"])
        fst_matrix[pk] = {"mean_fst": sum(fst_vals) / len(fst_vals) if fst_vals else None}
        dxy_matrix[pk] = {"mean_dxy": sum(dxy_vals) / len(dxy_vals) if dxy_vals else None}

    output["divergence"] = existing
    output["fst_matrix"] = fst_matrix
    output["dxy_matrix"] = dxy_matrix
    save_output(output)
    log.info("Divergence results saved to %s", OUTPUT_FILE)


# ── Subcommand: gscan ─────────────────────────────────────────────────

GSCAN_CACHE_DIR = "data/variants_cache"
GSCAN_WORKERS   = 8   # parallel chromosomes; each uses ~1.2 GB RAM peak


def _ensure_variant_cache(chromosomes, cache_dir=GSCAN_CACHE_DIR):
    """Build any missing .npz cache files before the parallel scan starts."""
    missing = [
        c for c in chromosomes
        if not os.path.exists(os.path.join(cache_dir, f"{c}.npz"))
    ]
    meta_missing = not os.path.exists(os.path.join(cache_dir, "pop_meta.json"))

    if not missing and not meta_missing:
        return

    log.info("Variant cache incomplete — exporting %d chromosome(s) + meta...",
             len(missing) + int(meta_missing))
    import subprocess, sys as _sys
    script = os.path.join(os.path.dirname(__file__), "export_variant_cache.py")
    cmd = [_sys.executable, script,
           "--cache-dir", cache_dir,
           "--workers", "4"]
    if missing:
        cmd += ["--chr"] + missing
    ret = subprocess.run(cmd, check=False)
    if ret.returncode != 0:
        raise RuntimeError(
            "export_variant_cache.py failed — check export_variant_cache.log"
        )


def cmd_gscan(args):
    """Run genome scans: XPEHH pairs + annotation-conditioned, via numpy cache."""
    output   = load_output()
    existing = output.get("gscan", {})
    lock     = Lock()

    # Build task list (resume: skip already-completed keys)
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
    n_done  = n_total - len(tasks_all)
    log.info("=== Genome scans (numpy): %d calls (%d cached, %d remaining) ===",
             n_total, n_done, len(tasks_all))

    if not tasks_all:
        log.info("All gscan calls cached — nothing to do.")
        return

    # Ensure .npz cache exists for every chromosome we'll scan
    chroms_needed = sorted({t[2] for t in tasks_all})
    _ensure_variant_cache(chroms_needed, GSCAN_CACHE_DIR)

    def run_one(pop1, pop2, chrom, consequence, key):
        t0 = time.time()
        try:
            windows = genome_scan_numpy.scan(
                chrom, pop1, pop2,
                consequence=consequence,
                window_size=WINDOW_SIZE,
                step_size=WINDOW_STEP,
                cache_dir=GSCAN_CACHE_DIR,
            )
            return key, {
                "n_windows": len(windows),
                "wall_sec":  time.time() - t0,
                "result":    windows,
            }
        except Exception as e:
            return key, {"error": str(e), "wall_sec": time.time() - t0}

    # Parallel execution — no DB WRITE locks, safe to run all chromosomes at once
    with ThreadPoolExecutor(max_workers=GSCAN_WORKERS) as pool:
        futures = {pool.submit(run_one, *t): t for t in tasks_all}
        for fut in as_completed(futures):
            pop1, pop2, chrom, consequence, key = futures[fut]
            val = fut.result()[1]
            label = consequence or "all"
            with lock:
                existing[key] = val
                n_done += 1
                log.info("  [%d/%d] %s vs %s %s (%s): %d windows, %.1fs",
                         n_done, n_total, pop1, pop2, chrom, label,
                         val.get("n_windows", 0), val.get("wall_sec", 0))
                if n_done % 22 == 0:
                    output["gscan"] = existing
                    save_output(output)

    output["gscan"] = existing
    save_output(output)
    log.info("Genome scan results saved to %s", OUTPUT_FILE)


# ── Subcommand: pbs ──────────────────────────────────────────────────


def cmd_pbs(args):
    """Run PBS genome scans for focal populations — numpy bypass (same .npz cache as gscan)."""
    output   = load_output()
    existing = output.get("pbs", {})
    lock     = Lock()

    tasks_all = [
        (focal, sister, outgroup, chrom,
         f"{focal}|{sister}|{outgroup}|{chrom}")
        for focal, sister, outgroup in PBS_CONFIGS
        for chrom in CHROMOSOMES
        if (f"{focal}|{sister}|{outgroup}|{chrom}" not in existing
            or "error" in existing.get(f"{focal}|{sister}|{outgroup}|{chrom}", {}))
    ]

    n_total = len(PBS_CONFIGS) * len(CHROMOSOMES)
    n_done  = n_total - len(tasks_all)
    log.info("=== PBS genome scans: %d calls (%d cached) ===", n_total, n_done)

    if tasks_all:
        chroms_needed = sorted({t[3] for t in tasks_all})
        _ensure_variant_cache(chroms_needed, GSCAN_CACHE_DIR)

        def run_one(focal, sister, outgroup, chrom, key):
            t0 = time.time()
            windows = genome_scan_numpy.scan(
                chrom, focal, pop2=sister, pop3=outgroup,
                window_size=WINDOW_SIZE, step_size=WINDOW_STEP,
                cache_dir=GSCAN_CACHE_DIR,
            )
            return {"result": windows, "wall_sec": time.time() - t0}

        with ThreadPoolExecutor(max_workers=GSCAN_WORKERS) as pool:
            futures = {pool.submit(run_one, *t): t for t in tasks_all}
            for fut in as_completed(futures):
                focal, sister, outgroup, chrom, key = futures[fut]
                try:
                    val = fut.result()
                except Exception as e:
                    log.error("  %s|%s|%s %s FAILED: %s", focal, sister, outgroup, chrom, e)
                    val = {"error": str(e)}
                with lock:
                    existing[key] = val
                    n_done += 1
                    log.info("  [%d/%d] PBS %s|%s|%s %s: %d windows (%.1fs)",
                             n_done, n_total, focal, sister, outgroup, chrom,
                             len(val.get("result", [])), val.get("wall_sec", 0))
                    if n_done % 22 == 0:
                        output["pbs"] = existing
                        save_output(output)

    # Summarize top PBS windows per config — gene annotation via Neo4j (READ only)
    driver = get_driver()
    pbs_summary = {}
    for focal, sister, outgroup in PBS_CONFIGS:
        config_key = f"{focal}_vs_{sister}_out_{outgroup}"
        all_windows = []
        for chrom in CHROMOSOMES:
            key = f"{focal}|{sister}|{outgroup}|{chrom}"
            windows = existing.get(key, {}).get("result", [])
            for w in windows:
                if w.get("pbs") is not None and w.get("n_variants", 0) >= 10:
                    all_windows.append(w)

        all_windows.sort(key=lambda x: x.get("pbs", 0), reverse=True)
        top20 = all_windows[:20]

        annotated = []
        for w in top20:
            genes = []
            try:
                with driver.session(database=NEO4J_DB) as s:
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

    driver.close()
    output["pbs"] = existing
    output["pbs_summary"] = pbs_summary
    save_output(output)
    log.info("PBS results saved to %s", OUTPUT_FILE)


# ── Subcommand: tree ─────────────────────────────────────────────────


def cmd_tree(args):
    """Build population tree from Fst matrix using UPGMA."""
    output = load_output()
    fst_matrix = output.get("fst_matrix", {})
    if not fst_matrix:
        log.error("No fst_matrix found. Run 'divergence' first.")
        return

    pops = sorted(POPULATIONS)
    n = len(pops)
    pop_idx = {p: i for i, p in enumerate(pops)}
    dist = [[0.0] * n for _ in range(n)]

    for pair_key_str, v in fst_matrix.items():
        fst = v.get("mean_fst")
        if fst is None:
            continue
        parts = pair_key_str.split("_vs_")
        if len(parts) != 2:
            continue
        p1, p2 = parts
        if p1 in pop_idx and p2 in pop_idx:
            i, j = pop_idx[p1], pop_idx[p2]
            dist[i][j] = fst
            dist[j][i] = fst

    tree_newick, merge_log = _upgma(pops, dist)
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
    print("=== Pairwise Fst Matrix ===")
    print()
    header = "     " + "".join(f"{p:>7s}" for p in pops)
    print(header)
    for i, p1 in enumerate(pops):
        row = f"{p1:>4s} " + "".join(f"{dist[i][j]:7.4f}" for j in range(n))
        print(row)


def _upgma(labels, dist):
    n = len(labels)
    d = [row[:] for row in dist]
    clusters = {i: labels[i] for i in range(n)}
    sizes = {i: 1 for i in range(n)}
    heights = {i: 0.0 for i in range(n)}
    active = set(range(n))
    merge_log = []
    next_id = n

    while len(active) > 1:
        min_d = float("inf")
        mi, mj = -1, -1
        for i in active:
            for j in active:
                if i < j and d[i][j] < min_d:
                    min_d = d[i][j]
                    mi, mj = i, j

        h = min_d / 2.0
        merge_log.append((clusters[mi], clusters[mj], min_d, h))
        new_label = f"({clusters[mi]}:{h - heights[mi]:.4f},{clusters[mj]}:{h - heights[mj]:.4f})"

        new_row = [0.0] * (len(d) + 1)
        for k in active:
            if k != mi and k != mj:
                new_dist = (d[mi][k] * sizes[mi] + d[mj][k] * sizes[mj]) / (sizes[mi] + sizes[mj])
                new_row[k] = new_dist
        for row in d:
            row.append(0.0)
        d.append(new_row)
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
    lines = []
    lines.append("Merge order (UPGMA):")
    lines.append(f"{'Step':>4s}  {'Fst':>8s}  {'Cluster 1':<30s}  {'Cluster 2':<30s}")
    lines.append("-" * 80)
    for i, (c1, c2, fst, h) in enumerate(merge_log, 1):
        s1 = c1[:27] + "..." if len(c1) > 30 else c1
        s2 = c2[:27] + "..." if len(c2) > 30 else c2
        lines.append(f"{i:4d}  {fst:8.4f}  {s1:<30s}  {s2:<30s}")
    return lines


# ── Shared diversity helper (fay_wu / pinsps numpy bypass) ───────────


def _compute_diversity_chrom(chrom, pop_idx_map, consequence, cache_dir):
    """Compute whole-chromosome diversity stats for multiple populations.

    Loads the .npz cache once per chromosome, then iterates over populations.
    Returns {key: {"result": {...}, "wall_sec": float}} where
    key = f"{pop}|{chrom}" (consequence=None) or f"{pop}|{chrom}|{consequence}".

    Computes the same fields as graphpop.diversity:
    pi, theta_w, tajima_d, fay_wu_h, fay_wu_h_norm, het_exp, het_obs, fis,
    n_variants, n_segregating, n_polarized.
    """
    t0 = time.time()
    with np.load(os.path.join(cache_dir, f"{chrom}.npz")) as data:
        anc = data["ancestral_allele"]               # (M,) int8
        if consequence == "missense_variant":
            csq = data["has_missense"].copy()
        elif consequence == "synonymous_variant":
            csq = data["has_synonymous"].copy()
        else:
            csq = np.ones(len(anc), dtype=bool)
        # Load full (M,26) arrays; slice per-pop below
        ac_all  = data["ac"]                         # int32
        an_all  = data["an"]                         # int32
        af_all  = data["af"]                         # float32
        het_all = data["het_count"]                  # int32

    results = {}
    for pop, idx in pop_idx_map.items():
        key = f"{pop}|{chrom}" if consequence is None else f"{pop}|{chrom}|{consequence}"

        valid = csq & (an_all[:, idx] >= 2)
        valid_idx = np.where(valid)[0]
        n_var = len(valid_idx)

        if n_var == 0:
            results[key] = {"result": {
                "pi": 0.0, "theta_w": 0.0, "tajima_d": 0.0,
                "fay_wu_h": 0.0, "fay_wu_h_norm": 0.0,
                "het_exp": 0.0, "het_obs": 0.0, "fis": 0.0,
                "n_variants": 0, "n_segregating": 0, "n_polarized": 0,
            }, "wall_sec": 0.0}
            continue

        ac  = ac_all[:, idx].astype(np.float64)
        an  = an_all[:, idx].astype(np.float64)
        af  = af_all[:, idx].astype(np.float64)
        het = het_all[:, idx].astype(np.float64)

        # n = haploid AN of first valid variant (matches Java behaviour)
        n    = int(an_all[valid_idx[0], idx])
        a_n  = genome_scan_numpy._harmonic(n - 1)  if n > 1 else 0.0
        a_n2 = genome_scan_numpy._harmonic2(n - 1) if n > 1 else 0.0

        with np.errstate(invalid="ignore", divide="ignore"):
            pi_arr = np.where(valid & (an > 1),
                              2.0 * af * (1.0 - af) * an / (an - 1.0), 0.0)

        pi_sum = float(pi_arr.sum())
        n_seg  = int((valid & (ac > 0) & (ac < an)).sum())

        # Tajima's D
        if n_seg > 0 and a_n > 0 and n >= 4:
            theta_w_raw = n_seg / a_n
            d   = pi_sum - theta_w_raw
            b1  = (n + 1) / (3.0 * (n - 1))
            b2  = 2.0 * (n**2 + n + 3) / (9.0 * n * (n - 1))
            c1  = b1 - 1.0 / a_n
            c2  = b2 - (n + 2) / (a_n * n) + a_n2 / (a_n**2)
            e1  = c1 / a_n
            e2  = c2 / (a_n**2 + a_n2)
            var_d = e1 * n_seg + e2 * n_seg * (n_seg - 1)
            tajima_d     = d / var_d**0.5 if var_d > 0 else 0.0
            theta_w_site = theta_w_raw / n_var
        else:
            tajima_d = theta_w_site = 0.0

        # Heterozygosity / Fis
        he_sum = float(np.where(valid, 2.0 * af * (1.0 - af), 0.0).sum())
        with np.errstate(invalid="ignore", divide="ignore"):
            ho_sum = float(np.where(valid & (an > 0), het / (an / 2.0), 0.0).sum())
        het_exp = he_sum / n_var
        het_obs = ho_sum / n_var
        fis = 1.0 - het_obs / het_exp if het_exp > 0 else 0.0

        # Fay & Wu's H
        derived = np.where(anc == 1, ac,
                  np.where(anc == 2, an - ac, -1.0))
        pol_valid = valid & (derived >= 0) & (derived > 0) & (derived < an)
        n_pol = int((valid & (anc > 0)).sum())  # ALL polarized sites (matches Java)

        with np.errstate(invalid="ignore", divide="ignore"):
            theta_h_arr  = np.where(pol_valid, 2.0 * derived**2 / (an * (an - 1.0)), 0.0)
            pi_polar_arr = np.where(pol_valid, pi_arr, 0.0)

        pi_pol_sum = float(pi_polar_arr.sum())
        th_sum     = float(theta_h_arr.sum())
        H          = pi_pol_sum - th_sum
        fay_wu_h   = H / n_pol if n_pol > 0 else 0.0

        # Normalized H (Zeng et al. 2006)
        if n_pol > 0 and n >= 4 and a_n > 0:
            theta = n_pol / a_n
            b_n   = genome_scan_numpy._harmonic2(n - 1)
            nD    = float(n)
            term1 = theta * (nD - 2.0) / (6.0 * (nD - 1.0))
            term2 = theta**2 * (
                (18.0 * nD**2 * (3.0 * nD + 2.0) * b_n
                 - (88.0 * nD**3 + 9.0 * nD**2 - 13.0 * nD + 6.0))
                / (9.0 * nD * (nD - 1.0)**2)
            )
            var_h = term1 + term2
            fay_wu_h_norm = H / var_h**0.5 if var_h > 0 else 0.0
        else:
            fay_wu_h_norm = 0.0

        results[key] = {"result": {
            "pi":            pi_sum / n_var,
            "theta_w":       theta_w_site,
            "tajima_d":      tajima_d,
            "fay_wu_h":      fay_wu_h,
            "fay_wu_h_norm": fay_wu_h_norm,
            "het_exp":       het_exp,
            "het_obs":       het_obs,
            "fis":           fis,
            "n_variants":    n_var,
            "n_segregating": n_seg,
            "n_polarized":   n_pol,
        }, "wall_sec": time.time() - t0}

    return results


# ── Subcommand: fay_wu ────────────────────────────────────────────────


def cmd_fay_wu(args):
    """Compute Fay & Wu's H — numpy bypass (same .npz cache as gscan)."""
    output   = load_output()
    existing = output.get("fay_wu", {})

    chroms_needed = sorted({
        c for pop in POPULATIONS for c in CHROMOSOMES
        if f"{pop}|{c}" not in existing
        or "error" in existing.get(f"{pop}|{c}", {})
    })

    n_total = len(POPULATIONS) * len(CHROMOSOMES)
    n_done  = n_total - len(POPULATIONS) * len(chroms_needed)
    log.info("=== Fay & Wu's H: %d calls (%d cached) ===", n_total, n_done)

    if chroms_needed:
        _ensure_variant_cache(chroms_needed, GSCAN_CACHE_DIR)
        with open(os.path.join(GSCAN_CACHE_DIR, "pop_meta.json")) as f:
            pop_ids = json.load(f)["pop_ids"]
        pop_idx_map = {p: pop_ids.index(p) for p in POPULATIONS}

        for chrom in chroms_needed:
            t0 = time.time()
            new = _compute_diversity_chrom(chrom, pop_idx_map, None, GSCAN_CACHE_DIR)
            existing.update(new)
            log.info("  %s: %d pops (%.1fs)", chrom, len(POPULATIONS), time.time() - t0)
            output["fay_wu"] = existing
            save_output(output)
    else:
        log.info("  all cached")

    fay_wu_summary = {}
    for pop in POPULATIONS:
        h_vals, hn_vals = [], []
        for chrom in CHROMOSOMES:
            r = existing.get(f"{pop}|{chrom}", {}).get("result", {})
            if r.get("fay_wu_h") is not None and r.get("n_polarized", 0) > 0:
                h_vals.append(r["fay_wu_h"])
                hn_vals.append(r.get("fay_wu_h_norm", 0))
        fay_wu_summary[pop] = {
            "mean_h":      sum(h_vals)  / len(h_vals)  if h_vals  else None,
            "mean_h_norm": sum(hn_vals) / len(hn_vals) if hn_vals else None,
        }

    output["fay_wu"] = existing
    output["fay_wu_summary"] = fay_wu_summary
    save_output(output)
    log.info("Fay & Wu's H results saved to %s", OUTPUT_FILE)


# ── Subcommand: usfs ─────────────────────────────────────────────────


def _compute_usfs_chrom(chrom, pop_idx_map, cache_dir):
    """Compute unfolded SFS for multiple populations from one .npz file.

    Matches SFSProcedure.java (unfolded=true): histogram of derived allele
    counts indexed 0..max_AN. Unpolarized variants use raw AC as count.
    Returns {f"{pop}|{chrom}": {"result": {...}, "wall_sec": float}}.
    """
    t0 = time.time()
    with np.load(os.path.join(cache_dir, f"{chrom}.npz")) as data:
        anc    = data["ancestral_allele"]    # (M,) int8
        ac_all = data["ac"]                  # (M,26) int32
        an_all = data["an"]                  # (M,26) int32

    results = {}
    for pop, idx in pop_idx_map.items():
        key = f"{pop}|{chrom}"
        ac_pop = ac_all[:, idx].astype(np.int64)
        an_pop = an_all[:, idx].astype(np.int64)

        valid  = an_pop >= 2
        n_var  = int(valid.sum())
        if n_var == 0:
            results[key] = {"result": {
                "sfs": [], "n_variants": 0, "max_ac": 0, "n_polarized": 0,
            }, "wall_sec": 0.0}
            continue

        max_an  = int(an_pop[valid].max())
        # Derived allele count: REF-ancestral → derived=ac; ALT-ancestral → derived=an−ac
        derived = np.where(anc == 1, ac_pop,
                  np.where(anc == 2, an_pop - ac_pop,
                           ac_pop))            # unpolarized: use raw ac
        n_pol   = int(((anc == 1) | (anc == 2))[valid].sum())

        counts  = derived[valid].clip(0, max_an).astype(np.int64)
        sfs     = np.bincount(counts, minlength=max_an + 1)

        results[key] = {"result": {
            "sfs":         sfs.tolist(),
            "n_variants":  n_var,
            "max_ac":      max_an,
            "n_polarized": n_pol,
        }, "wall_sec": time.time() - t0}

    return results


def cmd_usfs(args):
    """Compute unfolded (polarized) SFS — numpy bypass (same .npz cache as gscan)."""
    output   = load_output()
    existing = output.get("usfs", {})

    chroms_needed = sorted({
        c for pop in POPULATIONS for c in CHROMOSOMES
        if f"{pop}|{c}" not in existing
        or "error" in existing.get(f"{pop}|{c}", {})
    })

    n_total = len(POPULATIONS) * len(CHROMOSOMES)
    n_done  = n_total - len(POPULATIONS) * len(chroms_needed)
    log.info("=== Unfolded SFS: %d calls (%d cached) ===", n_total, n_done)

    if chroms_needed:
        _ensure_variant_cache(chroms_needed, GSCAN_CACHE_DIR)
        with open(os.path.join(GSCAN_CACHE_DIR, "pop_meta.json")) as f:
            pop_ids = json.load(f)["pop_ids"]
        pop_idx_map = {p: pop_ids.index(p) for p in POPULATIONS}

        for chrom in chroms_needed:
            t0 = time.time()
            new = _compute_usfs_chrom(chrom, pop_idx_map, GSCAN_CACHE_DIR)
            existing.update(new)
            log.info("  %s: %d pops (%.1fs)", chrom, len(POPULATIONS), time.time() - t0)
            output["usfs"] = existing
            save_output(output)
    else:
        log.info("  all cached")

    usfs_summary = {}
    for pop in POPULATIONS:
        total_sfs, total_pol, total_var = None, 0, 0
        for chrom in CHROMOSOMES:
            r = existing.get(f"{pop}|{chrom}", {}).get("result", {})
            sfs = r.get("sfs")
            if sfs:
                if total_sfs is None:
                    total_sfs = [0] * len(sfs)
                for i in range(min(len(sfs), len(total_sfs))):
                    total_sfs[i] += sfs[i]
                if len(sfs) > len(total_sfs):
                    total_sfs.extend(sfs[len(total_sfs):])
            total_pol += r.get("n_polarized", 0)
            total_var += r.get("n_variants", 0)
        usfs_summary[pop] = {"sfs": total_sfs, "n_polarized": total_pol, "n_variants": total_var}

    output["usfs"] = existing
    output["usfs_summary"] = usfs_summary
    save_output(output)
    log.info("Unfolded SFS results saved to %s", OUTPUT_FILE)


# ── Subcommand: hwscan ───────────────────────────────────────────────


def cmd_hwscan(args):
    """Fay & Wu's H genome scan — numpy bypass (same .npz cache as gscan)."""
    output   = load_output()
    existing = output.get("hwscan", {})
    lock     = Lock()

    # Build task list (skip cached)
    tasks_all = [
        (pop, chrom, f"{pop}|{chrom}")
        for pop in H_SCAN_POPS
        for chrom in CHROMOSOMES
        if f"{pop}|{chrom}" not in existing or "error" in existing.get(f"{pop}|{chrom}", {})
    ]

    n_total = len(H_SCAN_POPS) * len(CHROMOSOMES)
    n_done  = n_total - len(tasks_all)
    log.info("=== Fay & Wu's H genome scan: %d calls (%d cached) ===", n_total, n_done)

    if tasks_all:
        chroms_needed = sorted({t[1] for t in tasks_all})
        _ensure_variant_cache(chroms_needed, GSCAN_CACHE_DIR)

        def run_one(pop, chrom, key):
            t0 = time.time()
            windows = genome_scan_numpy.scan(
                chrom, pop,
                window_size=WINDOW_SIZE, step_size=WINDOW_STEP,
                cache_dir=GSCAN_CACHE_DIR,
            )
            return {"result": windows, "wall_sec": time.time() - t0}

        with ThreadPoolExecutor(max_workers=GSCAN_WORKERS) as pool:
            futures = {pool.submit(run_one, *t): t for t in tasks_all}
            for fut in as_completed(futures):
                pop, chrom, key = futures[fut]
                try:
                    val = fut.result()
                except Exception as e:
                    log.error("  %s/%s FAILED: %s", pop, chrom, e)
                    val = {"error": str(e)}
                with lock:
                    existing[key] = val
                    n_done += 1
                    log.info("  [%d/%d] %s/%s: %d windows (%.1fs)",
                             n_done, n_total, pop, chrom,
                             len(val.get("result", [])), val.get("wall_sec", 0))
                    if n_done % 22 == 0:
                        output["hwscan"] = existing
                        save_output(output)

    hwscan_summary = {}
    for pop in H_SCAN_POPS:
        all_windows = []
        for chrom in CHROMOSOMES:
            windows = existing.get(f"{pop}|{chrom}", {}).get("result", [])
            for w in windows:
                h = w.get("fay_wu_h")
                if h is not None and w.get("n_variants", 0) >= 20:
                    all_windows.append(w)
        all_windows.sort(key=lambda x: x.get("fay_wu_h", 0))
        hwscan_summary[pop] = {
            "n_windows": len(all_windows),
            "mean_h": (sum(w.get("fay_wu_h", 0) for w in all_windows) /
                       len(all_windows)) if all_windows else None,
            "top_sweep_windows": all_windows[:20],
        }

    output["hwscan"] = existing
    output["hwscan_summary"] = hwscan_summary
    save_output(output)
    log.info("H genome scan results saved to %s", OUTPUT_FILE)


# ── Subcommand: roh_hmm ──────────────────────────────────────────────


def cmd_roh_hmm(args):
    """ROH with HMM method."""
    output = load_output()
    existing = output.get("roh_hmm", {})
    driver = get_driver()
    lock = Lock()

    n_total = len(POPULATIONS) * len(CHROMOSOMES)
    n_done = sum(1 for k in existing if "error" not in existing[k])
    log.info("=== ROH HMM: %d calls (%d cached) ===", n_total, n_done)

    def run_one(pop, chrom):
        key = f"{pop}|{chrom}"
        if key in existing and "error" not in existing[key]:
            return key, existing[key]
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

    # graphpop.roh is Mode.READ — submit all tasks to one global pool
    all_tasks = [
        (pop, c)
        for pop in POPULATIONS
        for c in CHROMOSOMES
        if f"{pop}|{c}" not in existing or "error" in existing.get(f"{pop}|{c}", {})
    ]

    if not all_tasks:
        log.info("  all cached")
    else:
        with ThreadPoolExecutor(max_workers=8) as pool:
            futures = {pool.submit(run_one, *t): t for t in all_tasks}
            for fut in as_completed(futures):
                pop, chrom = futures[fut]
                key = f"{pop}|{chrom}"
                _, val = fut.result()
                with lock:
                    existing[key] = val
                    n_done += 1
                    log.info("  [%d/%d] %s/%s: %d samples with ROH (%.1fs)",
                             n_done, n_total, pop, chrom,
                             val.get("n_samples_with_roh", 0), val.get("wall_sec", 0))
                    if n_done % 22 == 0:
                        output["roh_hmm"] = existing
                        save_output(output)

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
    log.info("ROH HMM results saved to %s", OUTPUT_FILE)


# ── Subcommand: daf_enrichment ────────────────────────────────────────


def _bg_daf_from_cache(key_pops, cache_dir=GSCAN_CACHE_DIR):
    """Compute genome-wide background DAF per population from .npz variant cache.

    Replaces the full-table Neo4j scan (~15 min/pop) with a numpy pass over
    the pre-exported .npz files (~2-3 s total for all populations).

    ANC encoding (from export_variant_cache.py): 0=unknown, 1=REF, 2=ALT.
    DAF = af   when ancestral==REF (derived allele is ALT)
    DAF = 1-af when ancestral==ALT (derived allele is REF)
    """
    with open(os.path.join(cache_dir, "pop_meta.json")) as f:
        meta = json.load(f)
    pop_ids = meta["pop_ids"]
    pop_idx = {p: i for i, p in enumerate(pop_ids)}

    sum_daf = {p: 0.0 for p in key_pops}
    cnt     = {p: 0   for p in key_pops}

    for chrom in CHROMOSOMES:
        with np.load(os.path.join(cache_dir, f"{chrom}.npz")) as data:
            anc = data["ancestral_allele"]   # (M,) int8: 0/1/2
            pol = anc > 0                    # polarized mask
            af_all = data["af"]              # (M, 26) float32
            an_all = data["an"]              # (M, 26) int32
            for pop in key_pops:
                idx  = pop_idx[pop]
                mask = pol & (an_all[:, idx] >= 10)
                af_p = af_all[mask, idx]
                anc_p = anc[mask]
                daf  = np.where(anc_p == 1, af_p, 1.0 - af_p)
                sum_daf[pop] += float(daf.sum())
                cnt[pop]     += int(mask.sum())

    result = {}
    for pop in key_pops:
        n = cnt[pop]
        result[pop] = {
            "mean_daf":   float(sum_daf[pop] / n) if n > 0 else None,
            "n_variants": n,
        }
        log.info("  %s: background DAF = %.4f (n=%d)",
                 pop, result[pop]["mean_daf"] or 0, n)
    return result


def cmd_daf_enrichment(args):
    """Test derived allele frequency enrichment at selection signal regions."""
    output = load_output()
    driver = get_driver()

    pbs_summary = output.get("pbs_summary", {})
    xpehh_ann = output.get("annotate", {}).get("xpehh_annotations", {})
    results = {}

    KEY_POPS = ["YRI", "CEU", "CHB", "GIH", "PEL", "FIN"]

    log.info("Computing genome-wide background DAF (numpy cache bypass)...")
    bg_daf = _bg_daf_from_cache(KEY_POPS, GSCAN_CACHE_DIR)
    results["background_daf"] = bg_daf

    # PBS peak DAF
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
                "enrichment": mean_peak_daf / bg if bg and bg > 0 else None,
                "n_peaks": len(peak_dafs),
            }
            log.info("  %s: peak DAF = %.4f vs bg = %.4f",
                     focal, mean_peak_daf, bg or 0)

    results["pbs_peak_daf"] = pbs_daf

    # XP-EHH peak DAF
    log.info("Computing DAF at XP-EHH peak regions...")
    xpehh_daf = {}
    for pair_key_str, ann in xpehh_ann.items():
        top_peaks = ann.get("top_peaks", [])[:10]
        if not top_peaks:
            continue
        parts = pair_key_str.split("_vs_")
        if len(parts) != 2:
            continue
        focal = parts[0]
        if focal not in bg_daf:
            continue
        peak_dafs = []
        with driver.session(database=NEO4J_DB) as s:
            for peak in top_peaks:
                margin = 100000
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
            xpehh_daf[pair_key_str] = {
                "focal": focal,
                "mean_peak_daf": mean_peak_daf,
                "background_daf": bg,
                "enrichment": mean_peak_daf / bg if bg and bg > 0 else None,
                "n_peaks": len(peak_dafs),
            }

    results["xpehh_peak_daf"] = xpehh_daf

    driver.close()
    output["daf_enrichment"] = results
    save_output(output)
    log.info("DAF enrichment results saved to %s", OUTPUT_FILE)


# ── Subcommand: report ────────────────────────────────────────────────


def cmd_report(args):
    """Generate final Markdown report."""
    output = load_output()
    results = load_results()
    lines = []

    def w(s=""):
        lines.append(s)

    w("# 1000 Genomes Project — Biological Interpretation Report")
    w()
    w("## 1. Population Diversity Ranking")
    w()

    div_table = output.get("annotate", {}).get("diversity_table", {})
    if div_table:
        ranked = sorted(div_table.items(),
                        key=lambda x: x[1].get("mean_pi") or 0, reverse=True)
        w("| Rank | Population | Continent | Mean pi | Mean theta_W | Mean Tajima's D | Mean F_IS |")
        w("|------|-----------|-----------|---------|-------------|----------------|----------|")
        for i, (pop, v) in enumerate(ranked, 1):
            continent = _pop_continent(pop)
            w(f"| {i} | {pop} | {continent} | {_fmt(v.get('mean_pi'))} | "
              f"{_fmt(v.get('mean_theta_w'))} | {_fmt(v.get('mean_tajima_d'))} | "
              f"{_fmt(v.get('mean_fis'))} |")
        w()
        if ranked:
            w("**Key observations:**")
            w(f"- Most diverse: **{ranked[0][0]}** (pi = {_fmt(ranked[0][1].get('mean_pi'))})")
            w(f"- Least diverse: **{ranked[-1][0]}** (pi = {_fmt(ranked[-1][1].get('mean_pi'))})")
            w("- African populations expected to show highest diversity (Out-of-Africa)")
            w()

    # 2. pi_N/pi_S
    w("## 2. Purifying Selection Efficiency (pi_N/pi_S)")
    w()
    ratios = output.get("pinsps_ratios", {})
    if ratios:
        ranked = sorted(ratios.items(),
                        key=lambda x: x[1].get("mean_pi_n_pi_s") or 0, reverse=True)
        w("| Rank | Population | Continent | Mean pi_N/pi_S |")
        w("|------|-----------|-----------|---------------|")
        for i, (pop, v) in enumerate(ranked, 1):
            w(f"| {i} | {pop} | {_pop_continent(pop)} | {_fmt(v.get('mean_pi_n_pi_s'))} |")
        w()
        w("**Expectation:** Non-African populations should show elevated pi_N/pi_S "
          "due to reduced efficiency of purifying selection during Out-of-Africa bottleneck.")
        w()

    # 3. Pairwise Divergence
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
        for pk, v in top10:
            w(f"| {pk} | {_fmt(v.get('mean_fst'), 4)} |")
        w()

        bottom5 = sorted(fst_matrix.items(),
                         key=lambda x: x[1].get("mean_fst") or 999)[:5]
        w("### Most closely related pairs")
        w()
        w("| Pair | Mean Fst |")
        w("|------|---------|")
        for pk, v in bottom5:
            w(f"| {pk} | {_fmt(v.get('mean_fst'), 4)} |")
        w()

    # 4. XP-EHH Peaks
    w("## 4. Selection Signals (XP-EHH Peaks)")
    w()
    xpehh_ann = output.get("annotate", {}).get("xpehh_annotations", {})
    for pk, ann in xpehh_ann.items():
        w(f"### {pk}")
        w(f"Total hits: {ann.get('n_total_hits', 0)}, Peaks: {ann.get('n_peaks', 0)}")
        w()

        known_hits = []
        for peak in ann.get("top_peaks", []):
            for kg in peak.get("known_selection_genes", []):
                known_hits.append({**kg, "xpehh": peak["peak_xpehh"],
                                   "peak_pos": peak["peak_pos"]})
        if known_hits:
            w("**Known selection genes detected:**")
            w()
            w("| Gene | Function | XP-EHH | Peak Position |")
            w("|------|----------|--------|--------------|")
            for kh in known_hits:
                w(f"| {kh['gene']} | {kh['function']} | "
                  f"{kh['xpehh']:.2f} | {kh['chr']}:{kh['peak_pos']:,} |")
            w()

        w("**Top 10 peaks (by |XP-EHH|):**")
        w()
        w("| Chr | Position | XP-EHH | #Hits | Top Gene | Impact |")
        w("|-----|----------|--------|-------|----------|--------|")
        for peak in ann.get("top_peaks", [])[:10]:
            genes = peak.get("genes", [])
            top_gene = genes[0]["symbol"] if genes else "-"
            top_impact = genes[0]["impact"] if genes else "-"
            w(f"| {peak['chr']} | {peak['peak_pos']:,} | "
              f"{peak['peak_xpehh']:.2f} | {peak['n_hits']} | "
              f"{top_gene} | {top_impact} |")
        w()

    # 5. Garud's H
    w("## 5. Selective Sweep Regions (Garud's H)")
    w()
    sweep_ann = output.get("annotate", {}).get("sweep_annotations", {})
    if sweep_ann:
        w("| Population | Continent | Sweeps (H12>0.05) | Hard | Soft |")
        w("|-----------|-----------|-------------------|------|------|")
        for pop in POPULATIONS:
            ann = sweep_ann.get(pop, {})
            sweeps = ann.get("sweeps", [])
            hard = sum(1 for s in sweeps if s.get("sweep_type") == "hard")
            w(f"| {pop} | {_pop_continent(pop)} | {len(sweeps)} | {hard} | {len(sweeps) - hard} |")
        w()

        w("### Notable sweeps near known selection genes")
        w()
        for pop in POPULATIONS:
            ann = sweep_ann.get(pop, {})
            for sw in ann.get("sweeps", []):
                if sw.get("known_selection_genes"):
                    for kg in sw["known_selection_genes"]:
                        w(f"- **{pop}** {sw['chr']}:{sw['start']:,}-{sw['end']:,} "
                          f"(H12={sw['h12']:.3f}, {sw['sweep_type']}): "
                          f"near **{kg['gene']}** ({kg['function']})")
        w()

    # 6. ROH
    w("## 6. Runs of Homozygosity (ROH)")
    w()
    roh_sum = output.get("roh_hmm_summary", {})
    if roh_sum:
        ranked = sorted(roh_sum.items(),
                        key=lambda x: x[1].get("mean_froh_across_chr", 0), reverse=True)
        w("| Rank | Population | Continent | Total ROH (Mb) | Mean FROH |")
        w("|------|-----------|-----------|---------------|-----------|")
        for i, (pop, v) in enumerate(ranked, 1):
            total_mb = v.get("total_roh_bp", 0) / 1e6
            w(f"| {i} | {pop} | {_pop_continent(pop)} | {total_mb:.1f} | "
              f"{_fmt(v.get('mean_froh_across_chr'))} |")
        w()
        w("**Expectation:** Populations with recent bottlenecks (PEL, FIN) or "
          "consanguinity practices should show elevated FROH.")
        w()

    # 7. Fay & Wu's H
    w("## 7. Fay & Wu's H")
    w()
    fwh_summary = output.get("fay_wu_summary", {})
    if fwh_summary:
        ranked = sorted(fwh_summary.items(), key=lambda x: x[1].get("mean_h") or 0)
        w("| Rank | Population | Continent | Mean H | Mean H_norm |")
        w("|------|-----------|-----------|--------|------------|")
        for i, (pop, v) in enumerate(ranked, 1):
            w(f"| {i} | {pop} | {_pop_continent(pop)} | {_fmt(v.get('mean_h'))} | "
              f"{_fmt(v.get('mean_h_norm'))} |")
        w()

    # 8. PBS
    w("## 8. Population Branch Statistic (PBS)")
    w()
    pbs_summary = output.get("pbs_summary", {})
    if pbs_summary:
        for config_key, info in pbs_summary.items():
            w(f"### {info['focal']} (vs {info['sister']}, outgroup {info['outgroup']})")
            w(f"- Windows: {info['n_windows']}, Mean PBS: {info['mean_pbs']:.4f}")
            w()
            top = info.get("top_windows", [])
            if top:
                w("**Top 10 PBS windows:**")
                w()
                w("| Chr | Start | End | PBS | Top Gene | Known Gene |")
                w("|-----|-------|-----|-----|----------|------------|")
                for tw in top[:10]:
                    genes = tw.get("genes", [])
                    top_gene = genes[0]["symbol"] if genes else "-"
                    known = tw.get("known_genes", [])
                    known_str = known[0]["gene"] if known else "-"
                    w(f"| {tw['chr']} | {tw['start']:,} | {tw['end']:,} | "
                      f"{tw['pbs']:.4f} | {top_gene} | {known_str} |")
                w()

    # 9. Population Tree
    w("## 9. Population Tree (UPGMA)")
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
        w("**Expected topology:** (((EUR, SAS), EAS), (AMR)), AFR)")
        w()

    # 10. DAF Enrichment
    w("## 10. Derived Allele Frequency Enrichment at Selection Signals")
    w()
    daf = output.get("daf_enrichment", {})
    if daf:
        pbs_daf = daf.get("pbs_peak_daf", {})
        if pbs_daf:
            w("### PBS Peak Regions")
            w()
            w("| Focal | Peak DAF | Background DAF | Enrichment |")
            w("|-------|----------|---------------|------------|")
            for ck, v in pbs_daf.items():
                w(f"| {v.get('focal', ck)} | {_fmt(v.get('mean_peak_daf'))} | "
                  f"{_fmt(v.get('background_daf'))} | "
                  f"{v.get('enrichment', 0):.2f}x |")
            w()

    # 11. Verification Checkpoints
    w("## 11. Verification Checkpoints")
    w()
    w("| Checkpoint | Expected | Status |")
    w("|-----------|----------|--------|")
    w("| Population tree topology | (((EUR,SAS),EAS),AMR),AFR | Check tree output |")
    w("| LCT region (chr2:136Mb) | Strong EUR selection | Check XP-EHH YRI_vs_CEU |")
    w("| SLC24A5 (chr15:48Mb) | EUR-specific sweep | Check Garud's H CEU |")
    w("| EDAR (chr2:109Mb) | EAS-specific sweep | Check XP-EHH YRI_vs_CHB |")
    w("| AFR highest diversity | pi_AFR > pi_others | Check diversity table |")
    w("| AMR highest FROH | Bottleneck + founder | Check ROH table |")
    w()

    # 12. GraphPop Unique Value
    w("## 12. What GraphPop Uniquely Enables")
    w()
    w("1. **Annotation-conditioned statistics**: `graphpop.diversity(..., "
      "{consequence: 'missense_variant'})` computes pi restricted to missense "
      "variants via Variant->HAS_CONSEQUENCE->Gene traversal.")
    w("2. **pi_N/pi_S across 26 populations**: 1,144 single-line procedure "
      "calls (26 pops x 22 chr x 2 consequences).")
    w("3. **Missense-conditioned genome scans**: Sliding-window Fst restricted "
      "to missense variants identifies adaptive protein divergence hotspots.")
    w("4. **Cross-statistic gene lookups**: After iHS/XP-EHH computation, "
      "immediate graph traversal to find associated genes.")
    w("5. **Full 1000G at scale**: 88M variants, 3,202 samples, 26 populations, "
      "12 procedures, all in a single graph database.")
    w()

    report_text = "\n".join(lines)
    report_path = "human_interpretation_report.md"
    with open(report_path, "w") as f:
        f.write(report_text)
    log.info("Report written to %s", report_path)
    print(report_text)


def _fmt(v, decimals=6):
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:.{decimals}f}"
    return str(v)


def _pop_continent(pop):
    for continent, pops in SUPERPOPS.items():
        if pop in pops:
            return continent
    return "?"


# ── Subcommand: hap_stats ────────────────────────────────────────────

HAP_CACHE_DIR = "data/hap_cache"


def _load_hap(pop, chrom, hap_cache_dir=HAP_CACHE_DIR):
    """Load haplotype columns for one population from the bit-packed cache."""
    import json as _json
    meta_path = os.path.join(hap_cache_dir, "sample_meta.json")
    npz_path  = os.path.join(hap_cache_dir, f"{chrom}.npz")
    with open(meta_path) as f:
        meta = _json.load(f)
    pop_ids   = meta["pop_ids"]
    n_samples = meta["n_samples"]
    sample_idx = np.array([i for i, p in enumerate(pop_ids) if p == pop], dtype=np.int32)
    d          = np.load(npz_path, mmap_mode="r")
    pos, hap_packed = np.array(d["pos"]), np.array(d["hap"])
    hap_full = np.unpackbits(hap_packed, axis=1, bitorder="big")[:, :2 * n_samples]
    hap_idx = np.empty(2 * len(sample_idx), dtype=np.int32)
    hap_idx[0::2] = 2 * sample_idx
    hap_idx[1::2] = 2 * sample_idx + 1
    return pos, hap_full[:, hap_idx].astype(np.int8)


def cmd_hap_stats(args):
    """
    Compute nSL + iHS + XP-EHH from haplotype cache jointly.

    Loads pop1 + pop2 haplotype arrays ONCE per chromosome, then
    computes all three statistics sharing the pre-loaded arrays.

    Output stored under "hap_stats" key in human_interpretation_results.json:
      key = "{pop1}|{pop2}|{chrom}"
      value = {
        "nsl_{pop1}":          {n_variants, nsl_mean, nsl_std, time_s},
        "ihs_{pop1}":          {n_variants, ihs_mean, ihs_std, time_s},
        "xpehh_{pop1}_{pop2}": {n_variants, xpehh_mean, xpehh_std, time_s},
        "load_time_s": float,
        "total_time_s": float,
      }
    """
    import allel

    output   = load_output()
    existing = output.setdefault("hap_stats", {})
    lock     = Lock()

    available_chroms = [
        c for c in CHROMOSOMES
        if os.path.exists(os.path.join(HAP_CACHE_DIR, f"{c}.npz"))
    ]
    if not available_chroms:
        log.error("No hap_cache found in %s — run export_hap_cache.py first", HAP_CACHE_DIR)
        return

    tasks = []
    for pop1, pop2 in XPEHH_PAIRS:
        for chrom in available_chroms:
            key = f"{pop1}|{pop2}|{chrom}"
            if key not in existing:
                tasks.append((pop1, pop2, chrom, key))
    log.info("hap_stats: %d tasks pending (%d already done)",
             len(tasks), len(existing))

    if not tasks:
        log.info("hap_stats: all tasks cached — nothing to do.")
        return

    def run_one(pop1, pop2, chrom, key):
        t0 = time.time()
        try:
            # Load both populations once
            t_load0 = time.perf_counter()
            pos1, hap1 = _load_hap(pop1, chrom)
            pos2, hap2 = _load_hap(pop2, chrom)
            load_s = time.perf_counter() - t_load0

            # Shared MAF filter
            h1   = allel.HaplotypeArray(hap1)
            h2   = allel.HaplotypeArray(hap2)
            ac1  = h1.count_alleles()
            af1  = ac1.to_frequencies()[:, 1]
            flt  = ac1.is_segregating() & (af1 >= 0.05) & (af1 <= 0.95)
            h1f  = h1[flt]; h2f = h2[flt]
            pos_flt = pos1[flt]; ac1f = ac1[flt]

            # nSL (pop1)
            t1 = time.perf_counter()
            nsl_raw = allel.nsl(h1f)
            try:
                nsl_std, _ = allel.standardize_by_allele_count(
                    nsl_raw, ac1f[:, 1], diagnostics=False)
            except Exception:
                nsl_std = (nsl_raw - np.nanmean(nsl_raw)) / max(np.nanstd(nsl_raw), 1e-9)
            nsl_time  = time.perf_counter() - t1
            nsl_valid = ~np.isnan(nsl_std)

            # iHS (pop1)
            t2 = time.perf_counter()
            ihs_raw = allel.ihs(h1f, pos_flt, include_edges=True)
            try:
                ihs_std, _ = allel.standardize_by_allele_count(
                    ihs_raw, ac1f[:, 1], diagnostics=False)
            except Exception:
                ihs_std = (ihs_raw - np.nanmean(ihs_raw)) / max(np.nanstd(ihs_raw), 1e-9)
            ihs_time  = time.perf_counter() - t2
            ihs_valid = ~np.isnan(ihs_std)

            # XP-EHH (pop1 vs pop2)
            t3 = time.perf_counter()
            xp_raw  = allel.xpehh(h1f, h2f, pos_flt, include_edges=True)
            xp_time = time.perf_counter() - t3
            xp_valid = ~np.isnan(xp_raw)

            result = {
                f"nsl_{pop1}": {
                    "n_variants": int(nsl_valid.sum()),
                    "nsl_mean":   float(np.nanmean(nsl_std)),
                    "nsl_std":    float(np.nanstd(nsl_std)),
                    "time_s":     nsl_time,
                },
                f"ihs_{pop1}": {
                    "n_variants": int(ihs_valid.sum()),
                    "ihs_mean":   float(np.nanmean(ihs_std)),
                    "ihs_std":    float(np.nanstd(ihs_std)),
                    "time_s":     ihs_time,
                },
                f"xpehh_{pop1}_{pop2}": {
                    "n_variants":  int(xp_valid.sum()),
                    "xpehh_mean":  float(np.nanmean(xp_raw)),
                    "xpehh_std":   float(np.nanstd(xp_raw)),
                    "time_s":      xp_time,
                },
                "load_time_s":  load_s,
                "total_time_s": time.time() - t0,
            }
            return key, result
        except Exception as e:
            log.error("  hap_stats %s: FAILED — %s", key, e, exc_info=True)
            return key, {"error": str(e), "total_time_s": time.time() - t0}

    n_done  = [0]
    n_total = len(tasks)
    with ThreadPoolExecutor(max_workers=GSCAN_WORKERS) as pool:
        futures = {pool.submit(run_one, *t): t for t in tasks}
        for fut in as_completed(futures):
            key, val = fut.result()
            with lock:
                existing[key] = val
                n_done[0] += 1
                log.info("  [%d/%d] %s: total=%.1fs",
                         n_done[0], n_total, key, val.get("total_time_s", 0))
                if n_done[0] % 10 == 0:
                    output["hap_stats"] = existing
                    save_output(output)

    output["hap_stats"] = existing
    save_output(output)
    log.info("hap_stats: done. %d results saved.", len(existing))


# ── Main ──────────────────────────────────────────────────────────────

SUBCOMMANDS = {
    "annotate": cmd_annotate,
    "pinsps": cmd_pinsps,
    "divergence": cmd_divergence,
    "gscan": cmd_gscan,
    "pbs": cmd_pbs,
    "tree": cmd_tree,
    "fay_wu": cmd_fay_wu,
    "usfs": cmd_usfs,
    "hwscan": cmd_hwscan,
    "roh_hmm": cmd_roh_hmm,
    "daf_enrichment": cmd_daf_enrichment,
    "report": cmd_report,
    "hap_stats": cmd_hap_stats,
}


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description="1000 Genomes Biological Interpretation")
    parser.add_argument(
        "command", choices=list(SUBCOMMANDS.keys()),
        help="Subcommand to run")
    args = parser.parse_args()
    SUBCOMMANDS[args.command](args)


if __name__ == "__main__":
    main()
