#!/usr/bin/env python3
"""
Benchmark: GraphPop nSL vs scikit-allel, GraphPop ROH vs PLINK 1.9.

Measures wall-clock time, peak memory, and numerical agreement.

Regions:
  - medium: chr22:16,000,000-16,500,000 (~5K variants)
  - large:  chr22:16,000,000-17,000,000 (~10K variants)
  - full:   full chr22 (~1.07M variants, nSL whole-chromosome)

Usage:
    conda run -n graphevo python scripts/benchmark-nsl-roh.py
    conda run -n graphevo python scripts/benchmark-nsl-roh.py --region full
    conda run -n graphevo python scripts/benchmark-nsl-roh.py --all-regions

Requires:
    - Neo4j running with GraphPop procedures deployed
    - 1000G chr22 VCF + panel file
    - conda env 'graphevo' with cyvcf2, scikit-allel, numpy
    - plink (1.9) for ROH comparison
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import tracemalloc
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent
VCF = str(ROOT / "data/raw/1000g/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz")
PANEL = str(ROOT / "data/raw/1000g/integrated_call_samples_v3.20130502.ALL.panel")
RESULTS_DIR = ROOT / "benchmarks" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CHR = "chr22"
POP = "EUR"

NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"

REGIONS = {
    "medium": (16_000_000, 16_500_000),
    "large":  (16_000_000, 17_000_000),
    "full":   (10_510_000, 50_810_000),
}

# ---------------------------------------------------------------------------
# Utility helpers (mirrored from benchmark-vs-tools.py)
# ---------------------------------------------------------------------------

def measure_python_call(func, *args, **kwargs):
    """Run a Python function, return (result, elapsed_s, peak_mem_mb)."""
    tracemalloc.start()
    t0 = time.perf_counter()
    result = func(*args, **kwargs)
    elapsed = time.perf_counter() - t0
    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result, elapsed, peak_bytes / (1024 * 1024)


def measure_subprocess(cmd, timeout=600, cwd=None):
    """Run a subprocess via /usr/bin/time, return (stdout, elapsed_s, peak_rss_mb)."""
    full_cmd = ["/usr/bin/time", "-v"] + cmd
    t0 = time.perf_counter()
    proc = subprocess.run(full_cmd, capture_output=True, text=True,
                          timeout=timeout, cwd=cwd)
    elapsed = time.perf_counter() - t0

    peak_rss_mb = 0.0
    for line in proc.stderr.split("\n"):
        m = re.search(r"Maximum resident set size.*?:\s*(\d+)", line)
        if m:
            peak_rss_mb = int(m.group(1)) / 1024.0
            break

    if proc.returncode != 0:
        actual_stderr = "\n".join(
            l for l in proc.stderr.split("\n")
            if not any(kw in l for kw in [
                "Maximum resident", "Command being timed", "wall clock",
                "Elapsed", "Minor", "Major", "Voluntary", "Involuntary",
                "File system", "Socket", "Signals", "Swaps", "Page size",
                "Exit status", "Percent of CPU", "User time", "System time",
                "Average",
            ])
        ).strip()
        raise RuntimeError(f"Command failed (rc={proc.returncode}): {actual_stderr[:500]}")

    return proc.stdout, elapsed, peak_rss_mb


def run_cypher(query, timeout=600):
    """Execute Cypher via cypher-shell, return (output, elapsed_s, peak_rss_mb)."""
    cmd = ["cypher-shell", "-u", NEO4J_USER, "-p", NEO4J_PASS,
           "--format", "plain", query]
    return measure_subprocess(cmd, timeout=timeout)


def parse_cypher_row(output):
    """Parse single-row cypher-shell plain output into dict."""
    lines = [l for l in output.strip().split("\n") if l.strip()]
    if len(lines) < 2:
        return {}
    keys = [k.strip().strip('"') for k in lines[0].split(",")]
    vals = [v.strip().strip('"') for v in lines[1].split(",")]
    result = {}
    for k, v in zip(keys, vals):
        try:
            result[k] = float(v)
        except ValueError:
            result[k] = v
    return result


def parse_cypher_rows(output):
    """Parse multi-row cypher-shell plain output into list of dicts."""
    lines = [l for l in output.strip().split("\n") if l.strip()]
    if len(lines) < 2:
        return []
    keys = [k.strip().strip('"') for k in lines[0].split(",")]
    rows = []
    for line in lines[1:]:
        vals = [v.strip().strip('"') for v in line.split(",")]
        row = {}
        for k, v in zip(keys, vals):
            try:
                row[k] = float(v)
            except ValueError:
                row[k] = v
        rows.append(row)
    return rows


_panel_cache = None

def load_panel():
    global _panel_cache
    if _panel_cache is not None:
        return _panel_cache
    sample_pop = {}
    with open(PANEL) as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                sample_pop[parts[0]] = parts[2]  # super_pop
    _panel_cache = sample_pop
    return sample_pop


def load_haplotype_data(start, end, pop):
    """Load haplotype data from VCF for a region and population."""
    import cyvcf2

    sample_pop = load_panel()
    vcf = cyvcf2.VCF(VCF)
    all_samples = list(vcf.samples)
    pop_indices = [i for i, s in enumerate(all_samples) if sample_pop.get(s) == pop]

    positions = []
    haplotypes_0 = []
    haplotypes_1 = []

    region = f"{CHR}:{start}-{end}"
    for v in vcf(region):
        if not v.is_snp or len(v.ALT) != 1:
            continue
        gt = v.genotypes
        h0 = []
        h1 = []
        for idx in pop_indices:
            g = gt[idx]
            h0.append(g[0])
            h1.append(g[1])
        h0 = np.array(h0, dtype=np.int8)
        h1 = np.array(h1, dtype=np.int8)
        # Skip monomorphic
        combined = np.concatenate([h0, h1])
        if combined.min() == combined.max():
            continue
        positions.append(v.POS)
        haplotypes_0.append(h0)
        haplotypes_1.append(h1)

    vcf.close()

    h0_arr = np.array(haplotypes_0)
    h1_arr = np.array(haplotypes_1)
    h = np.concatenate([h0_arr, h1_arr], axis=1)

    return {
        "positions": np.array(positions),
        "haplotypes": h,
        "n_samples": len(pop_indices),
        "n_variants": len(positions),
    }


def write_pop_keep_file(tmpdir, pop):
    """Write a keep file with one sample ID per line."""
    sample_pop = load_panel()
    path = os.path.join(tmpdir, f"{pop}_samples.txt")
    with open(path, "w") as f:
        for s, p in sample_pop.items():
            if p == pop:
                f.write(s + "\n")
    return path


def write_pop_keep_file_fam(tmpdir, pop):
    """Write a keep file in PLINK format (FID IID, one per line)."""
    sample_pop = load_panel()
    path = os.path.join(tmpdir, f"{pop}_keep.txt")
    with open(path, "w") as f:
        for s, p in sample_pop.items():
            if p == pop:
                f.write(f"{s} {s}\n")
    return path


def fmt_time(t):
    if t is None:
        return "—"
    return f"{t:.3f}s"


def fmt_mem(m):
    if m is None:
        return "—"
    return f"{m:.0f}MB"


def fmt_val(v, decimals=4):
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.{decimals}f}"
    return str(v)


def print_group_table(title, headers, rows):
    """Print a formatted comparison table."""
    print(f"\n{'='*90}")
    print(f"  {title}")
    print(f"{'='*90}")
    all_rows = [headers] + rows
    widths = [max(len(str(row[i])) for row in all_rows) for i in range(len(headers))]
    fmt = "  ".join(f"{{:<{w}}}" for w in widths)
    print(fmt.format(*headers))
    print("-" * sum(widths + [2 * (len(widths) - 1)]))
    for row in rows:
        print(fmt.format(*[str(x) for x in row]))


# ---------------------------------------------------------------------------
# Group F: nSL — scikit-allel vs GraphPop
# ---------------------------------------------------------------------------

def nsl_scikit_allel(start, end, pop):
    """Compute nSL using scikit-allel."""
    import allel

    def _compute():
        data = load_haplotype_data(start, end, pop)
        h = allel.HaplotypeArray(data["haplotypes"])

        ac = h.count_alleles()
        is_seg = ac.is_segregating()
        af = ac.to_frequencies()[:, 1]
        flt = is_seg & (af >= 0.05) & (af <= 0.95)

        h_flt = h[flt]
        pos_flt = data["positions"][flt]

        if len(pos_flt) < 10:
            return {"n_variants": 0, "nsl_mean": None}

        nsl_raw = allel.nsl(h_flt)

        # Standardize by allele count bins
        try:
            nsl_std, bins = allel.standardize_by_allele_count(
                nsl_raw, ac[flt][:, 1], diagnostics=False
            )
        except (ValueError, ZeroDivisionError):
            valid_raw = nsl_raw[~np.isnan(nsl_raw)]
            if len(valid_raw) > 1:
                nsl_std = (nsl_raw - np.nanmean(nsl_raw)) / np.nanstd(nsl_raw)
            else:
                nsl_std = nsl_raw

        valid = ~np.isnan(nsl_std)
        return {
            "n_variants": int(valid.sum()),
            "nsl_mean": float(np.nanmean(nsl_std)),
            "nsl_std": float(np.nanstd(nsl_std)),
        }

    result, elapsed, peak_mb = measure_python_call(_compute)
    return {**result, "time_s": elapsed, "peak_mem_mb": peak_mb}


def nsl_graphpop(start, end, pop):
    """Run graphpop.nsl (whole chromosome, then filter to region)."""
    # GraphPop nSL is whole-chromosome. We query and filter to region.
    query = (
        f"CALL graphpop.nsl('{CHR}', '{pop}', {{min_af: 0.05}}) "
        f"YIELD variantId, pos, af, nsl_unstd, nsl "
        f"WITH pos, nsl WHERE pos >= {start} AND pos <= {end} "
        f"RETURN count(*) AS n_variants, avg(nsl) AS nsl_mean, "
        f"stDev(nsl) AS nsl_std;"
    )
    out, elapsed, peak_mb = run_cypher(query, timeout=600)
    result = parse_cypher_row(out)
    return {
        "n_variants": result.get("n_variants"),
        "nsl_mean": result.get("nsl_mean"),
        "nsl_std": result.get("nsl_std"),
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


def nsl_graphpop_full(pop):
    """Run graphpop.nsl for whole chromosome (full region benchmark)."""
    query = (
        f"CALL graphpop.nsl('{CHR}', '{pop}', {{min_af: 0.05}}) "
        f"YIELD variantId, pos, af, nsl_unstd, nsl "
        f"RETURN count(*) AS n_variants, avg(nsl) AS nsl_mean, "
        f"stDev(nsl) AS nsl_std;"
    )
    out, elapsed, peak_mb = run_cypher(query, timeout=600)
    result = parse_cypher_row(out)
    return {
        "n_variants": result.get("n_variants"),
        "nsl_mean": result.get("nsl_mean"),
        "nsl_std": result.get("nsl_std"),
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


def nsl_scikit_allel_full(pop):
    """Compute nSL for full chromosome with scikit-allel.

    This is expensive: loads ~1M haplotypes into memory and computes
    pairwise EHH decay per variant. May take 10+ minutes.
    """
    import allel

    def _compute():
        start, end = REGIONS["full"]
        print("\n    Loading haplotype data...", end="", flush=True)
        data = load_haplotype_data(start, end, pop)
        print(f" {data['n_variants']} variants", flush=True)
        h = allel.HaplotypeArray(data["haplotypes"])

        ac = h.count_alleles()
        is_seg = ac.is_segregating()
        af = ac.to_frequencies()[:, 1]
        flt = is_seg & (af >= 0.05) & (af <= 0.95)

        h_flt = h[flt]
        print(f"    Computing nSL on {h_flt.shape[0]} variants...", flush=True)

        if h_flt.shape[0] < 10:
            return {"n_variants": 0, "nsl_mean": None}

        nsl_raw = allel.nsl(h_flt)

        try:
            nsl_std, bins = allel.standardize_by_allele_count(
                nsl_raw, ac[flt][:, 1], diagnostics=False
            )
        except (ValueError, ZeroDivisionError):
            valid_raw = nsl_raw[~np.isnan(nsl_raw)]
            if len(valid_raw) > 1:
                nsl_std = (nsl_raw - np.nanmean(nsl_raw)) / np.nanstd(nsl_raw)
            else:
                nsl_std = nsl_raw

        valid = ~np.isnan(nsl_std)
        return {
            "n_variants": int(valid.sum()),
            "nsl_mean": float(np.nanmean(nsl_std)),
            "nsl_std": float(np.nanstd(nsl_std)),
        }

    result, elapsed, peak_mb = measure_python_call(_compute)
    return {**result, "time_s": elapsed, "peak_mem_mb": peak_mb}


# ---------------------------------------------------------------------------
# Correlation helpers
# ---------------------------------------------------------------------------

def nsl_correlation(start, end, pop):
    """Compute Pearson correlation of nSL scores between scikit-allel and GraphPop.

    Matches variants by position, computes r on the unstandardized scores.

    Both scikit-allel and GraphPop now compute log(SL_1 / SL_0) using the
    original Ferrer-Admetlla (2014) pairwise SSL algorithm.
    """
    import allel

    # --- scikit-allel side ---
    data = load_haplotype_data(start, end, pop)
    h = allel.HaplotypeArray(data["haplotypes"])
    ac = h.count_alleles()
    af = ac.to_frequencies()[:, 1]
    flt = ac.is_segregating() & (af >= 0.05) & (af <= 0.95)

    h_flt = h[flt]
    pos_flt = data["positions"][flt]
    nsl_raw = allel.nsl(h_flt)

    sa_scores = {}
    for i, pos in enumerate(pos_flt):
        if not np.isnan(nsl_raw[i]):
            sa_scores[int(pos)] = float(nsl_raw[i])

    # --- GraphPop side: read previously computed nsl_unstd from Variant nodes ---
    query = (
        f"MATCH (v:Variant) "
        f"WHERE v.chr = '{CHR}' AND v.pos >= {start} AND v.pos <= {end} "
        f"AND v.nsl_unstd_{pop} IS NOT NULL "
        f"RETURN v.pos AS pos, v.nsl_unstd_{pop} AS nsl_unstd "
        f"ORDER BY pos;"
    )
    out, _, _ = run_cypher(query, timeout=120)
    gp_scores = {}
    lines = out.strip().split("\n")
    if len(lines) > 1:
        for line in lines[1:]:
            parts = line.strip().split(",")
            if len(parts) == 2:
                try:
                    gp_scores[int(float(parts[0].strip()))] = float(parts[1].strip())
                except ValueError:
                    continue

    # Match by position
    common_pos = sorted(set(sa_scores.keys()) & set(gp_scores.keys()))
    if len(common_pos) < 10:
        return {"n_common": len(common_pos), "pearson_r": None}

    sa_vals = np.array([sa_scores[p] for p in common_pos])
    gp_vals = np.array([gp_scores[p] for p in common_pos])

    r = np.corrcoef(sa_vals, gp_vals)[0, 1]
    return {
        "n_common": len(common_pos),
        "pearson_r": float(r),
    }


# ---------------------------------------------------------------------------
# Group G: ROH — PLINK 1.9 vs GraphPop
# ---------------------------------------------------------------------------

def roh_plink(pop, tmpdir, min_length_kb=500, min_snps=25, window_snps=50,
              max_het=1):
    """Run PLINK 1.9 --homozyg for ROH detection."""
    keep_file = write_pop_keep_file_fam(tmpdir, pop)

    cmd = [
        "plink",
        "--vcf", VCF,
        "--chr", "22",
        "--keep", keep_file,
        "--homozyg",
        "--homozyg-snp", str(min_snps),
        "--homozyg-kb", str(min_length_kb),
        "--homozyg-window-snp", str(window_snps),
        "--homozyg-window-het", str(max_het),
        "--homozyg-density", "50",      # max kb/SNP within ROH
        "--homozyg-gap", "1000",        # max gap kb between SNPs
        "--allow-extra-chr",
        "--out", os.path.join(tmpdir, "plink_roh"),
    ]

    try:
        _, elapsed, peak_mb = measure_subprocess(cmd, timeout=600, cwd=tmpdir)
    except RuntimeError as e:
        return {"error": str(e), "time_s": 0, "peak_mem_mb": 0}

    # Parse PLINK .hom.indiv output (per-individual summary)
    indiv_file = os.path.join(tmpdir, "plink_roh.hom.indiv")
    results = {}
    if os.path.exists(indiv_file):
        with open(indiv_file) as f:
            header = f.readline().split()
            for line in f:
                parts = line.split()
                if len(parts) >= 5:
                    sample_id = parts[1]  # IID
                    n_roh = int(parts[3])  # NSEG
                    total_kb = float(parts[4])  # KB
                    results[sample_id] = {
                        "n_roh": n_roh,
                        "total_length_kb": total_kb,
                    }

    # Parse detailed .hom output for segment-level data
    hom_file = os.path.join(tmpdir, "plink_roh.hom")
    n_segments = 0
    if os.path.exists(hom_file):
        with open(hom_file) as f:
            f.readline()  # skip header
            for line in f:
                if line.strip():
                    n_segments += 1

    n_with_roh = sum(1 for v in results.values() if v["n_roh"] > 0)
    total_samples = len(results)
    total_roh = sum(v["n_roh"] for v in results.values())
    mean_froh = 0.0
    if results:
        chr22_len_kb = 50818.468  # chr22 in kb
        frohs = [v["total_length_kb"] / chr22_len_kb for v in results.values()]
        mean_froh = float(np.mean(frohs))

    # Top 5 by total ROH length
    top5 = sorted(results.items(), key=lambda x: x[1]["total_length_kb"], reverse=True)[:5]

    return {
        "n_samples": total_samples,
        "n_with_roh": n_with_roh,
        "total_segments": total_roh,
        "n_segments_detailed": n_segments,
        "mean_froh": mean_froh,
        "top5": {s: v for s, v in top5},
        "all_results": results,
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


def roh_graphpop(pop, min_length=500000, min_snps=25, window_snps=50,
                 max_het=1):
    """Run graphpop.roh."""
    query = (
        f"CALL graphpop.roh('{CHR}', '{pop}', "
        f"{{min_length: {min_length}, min_snps: {min_snps}, "
        f"window_snps: {window_snps}, max_het: {max_het}}}) "
        f"YIELD sampleId, n_roh, total_length, froh, mean_length, max_length "
        f"RETURN sampleId, n_roh, total_length, froh, mean_length, max_length;"
    )
    out, elapsed, peak_mb = run_cypher(query, timeout=600)
    rows = parse_cypher_rows(out)

    results = {}
    for row in rows:
        results[row["sampleId"]] = {
            "n_roh": int(row["n_roh"]),
            "total_length_kb": row["total_length"] / 1000.0,
            "froh": row["froh"],
            "mean_length": row["mean_length"],
            "max_length": row.get("max_length", 0),
        }

    n_with_roh = sum(1 for v in results.values() if v["n_roh"] > 0)
    total_roh = sum(v["n_roh"] for v in results.values())
    mean_froh = float(np.mean([v["froh"] for v in results.values()])) if results else 0.0

    top5 = sorted(results.items(), key=lambda x: x[1]["total_length_kb"], reverse=True)[:5]

    return {
        "n_samples": len(results),
        "n_with_roh": n_with_roh,
        "total_segments": total_roh,
        "mean_froh": mean_froh,
        "top5": {s: v for s, v in top5},
        "all_results": results,
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


def roh_correlation(plink_results, graphpop_results):
    """Compute per-sample correlation of ROH metrics between PLINK and GraphPop."""
    pr = plink_results.get("all_results", {})
    gr = graphpop_results.get("all_results", {})

    common = sorted(set(pr.keys()) & set(gr.keys()))
    if len(common) < 10:
        return {"n_common": len(common), "pearson_r_total_kb": None,
                "pearson_r_n_roh": None}

    plink_total = np.array([pr[s]["total_length_kb"] for s in common])
    gp_total = np.array([gr[s]["total_length_kb"] for s in common])

    plink_nroh = np.array([pr[s]["n_roh"] for s in common], dtype=float)
    gp_nroh = np.array([gr[s]["n_roh"] for s in common], dtype=float)

    r_total = np.corrcoef(plink_total, gp_total)[0, 1]
    r_nroh = np.corrcoef(plink_nroh, gp_nroh)[0, 1]

    # Mean absolute difference in FROH
    chr22_len_kb = 50818.468
    plink_froh = plink_total / chr22_len_kb
    gp_froh = np.array([gr[s].get("froh", gr[s]["total_length_kb"] / chr22_len_kb)
                         for s in common])
    mae_froh = float(np.mean(np.abs(plink_froh - gp_froh)))

    return {
        "n_common": len(common),
        "pearson_r_total_kb": float(r_total),
        "pearson_r_n_roh": float(r_nroh),
        "mae_froh": mae_froh,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Benchmark GraphPop nSL and ROH vs external tools")
    parser.add_argument("--region", default="large",
                        choices=REGIONS.keys(),
                        help="Region size for nSL benchmark (default: large)")
    parser.add_argument("--all-regions", action="store_true",
                        help="Run all region sizes")
    parser.add_argument("--output", default=None,
                        help="Output JSONL file path")
    parser.add_argument("--nsl-only", action="store_true",
                        help="Run only nSL benchmarks")
    parser.add_argument("--roh-only", action="store_true",
                        help="Run only ROH benchmarks")
    args = parser.parse_args()

    if args.all_regions:
        regions_to_run = list(REGIONS.keys())
    else:
        regions_to_run = [args.region]

    has_plink = shutil.which("plink") is not None

    print("=" * 90)
    print("  GraphPop nSL + ROH Benchmark")
    print("=" * 90)
    print(f"  VCF: {VCF}")
    print(f"  Population: {POP}")
    print(f"  Regions: {', '.join(regions_to_run)}")
    print(f"  Tools: scikit-allel (nSL), PLINK 1.9 (ROH): {'available' if has_plink else 'NOT FOUND'}")
    print()

    all_results = []

    with tempfile.TemporaryDirectory(prefix="graphpop_bench_") as tmpdir:
        for region_name in regions_to_run:
            start, end = REGIONS[region_name]
            region_record = {
                "region": region_name, "chr": CHR,
                "start": start, "end": end, "pop": POP,
            }

            # ================================================================
            # nSL Benchmark
            # ================================================================
            if not args.roh_only:
                print(f"\n{'#'*90}")
                print(f"  nSL — {region_name} ({start:,}-{end:,})")
                print(f"{'#'*90}")

                nsl_results = {}

                # scikit-allel
                print("  scikit-allel...", end="", flush=True)
                if region_name == "full":
                    nsl_results["scikit-allel"] = nsl_scikit_allel_full(POP)
                else:
                    nsl_results["scikit-allel"] = nsl_scikit_allel(start, end, POP)
                sa = nsl_results["scikit-allel"]
                print(f" {fmt_time(sa['time_s'])}  ({sa.get('n_variants', 0)} variants)")

                # GraphPop
                print("  GraphPop...", end="", flush=True)
                if region_name == "full":
                    nsl_results["GraphPop"] = nsl_graphpop_full(POP)
                else:
                    nsl_results["GraphPop"] = nsl_graphpop(start, end, POP)
                gp = nsl_results["GraphPop"]
                print(f" {fmt_time(gp['time_s'])}  ({gp.get('n_variants', 0)} variants)")

                # Correlation
                print("  Computing nSL correlation...", end="", flush=True)
                if region_name == "full":
                    corr = nsl_correlation(REGIONS["full"][0], REGIONS["full"][1], POP)
                else:
                    corr = nsl_correlation(start, end, POP)
                nsl_results["correlation"] = corr
                print(f" r = {fmt_val(corr.get('pearson_r'), 4)} "
                      f"(n={corr.get('n_common', 0)} common variants)")

                # Speedup
                sa_t = sa.get("time_s", 0)
                gp_t = gp.get("time_s", 0)
                speedup = sa_t / gp_t if gp_t > 0 else float("inf")

                rows = [
                    ["scikit-allel",
                     fmt_val(sa.get("n_variants"), 0),
                     fmt_val(sa.get("nsl_mean"), 4),
                     fmt_val(sa.get("nsl_std"), 4),
                     fmt_time(sa.get("time_s")),
                     fmt_mem(sa.get("peak_mem_mb"))],
                    ["GraphPop",
                     fmt_val(gp.get("n_variants"), 0),
                     fmt_val(gp.get("nsl_mean"), 4),
                     fmt_val(gp.get("nsl_std"), 4),
                     fmt_time(gp.get("time_s")),
                     fmt_mem(gp.get("peak_mem_mb"))],
                ]
                print_group_table(
                    f"nSL — {region_name}  |  Speedup: {speedup:.1f}x  |  "
                    f"Correlation: r={fmt_val(corr.get('pearson_r'), 4)}",
                    ["Tool", "Variants", "Mean nSL", "Std nSL", "Time", "Peak Mem"],
                    rows,
                )

                region_record["nsl"] = nsl_results

            # ================================================================
            # ROH Benchmark (always whole-chromosome, run once for "large")
            # ================================================================
            if not args.nsl_only and region_name in ("large", "full"):
                print(f"\n{'#'*90}")
                print(f"  ROH — whole chr22  (min_length=500kb, min_snps=25, "
                      f"window_snps=50, max_het=1)")
                print(f"{'#'*90}")

                roh_results = {}

                # PLINK 1.9
                if has_plink:
                    print("  PLINK 1.9...", end="", flush=True)
                    roh_results["PLINK"] = roh_plink(
                        POP, tmpdir, min_length_kb=500, min_snps=25,
                        window_snps=50, max_het=1)
                    pk = roh_results["PLINK"]
                    print(f" {fmt_time(pk['time_s'])}  "
                          f"({pk.get('n_with_roh', 0)}/{pk.get('n_samples', 0)} "
                          f"samples with ROH)")
                else:
                    print("  PLINK 1.9... SKIPPED (not installed)")

                # GraphPop
                print("  GraphPop...", end="", flush=True)
                roh_results["GraphPop"] = roh_graphpop(
                    POP, min_length=500000, min_snps=25,
                    window_snps=50, max_het=1)
                gr = roh_results["GraphPop"]
                print(f" {fmt_time(gr['time_s'])}  "
                      f"({gr.get('n_with_roh', 0)}/{gr.get('n_samples', 0)} "
                      f"samples with ROH)")

                # Correlation
                if has_plink and "PLINK" in roh_results:
                    print("  Computing ROH correlation...", end="", flush=True)
                    corr = roh_correlation(roh_results["PLINK"],
                                          roh_results["GraphPop"])
                    roh_results["correlation"] = corr
                    print(f" r(total_kb)={fmt_val(corr.get('pearson_r_total_kb'), 4)}, "
                          f"r(n_roh)={fmt_val(corr.get('pearson_r_n_roh'), 4)}, "
                          f"MAE(FROH)={fmt_val(corr.get('mae_froh'), 6)}")

                # Speedup
                pk_t = roh_results.get("PLINK", {}).get("time_s", 0)
                gp_t = gr.get("time_s", 0)
                speedup_str = ""
                if pk_t > 0 and gp_t > 0:
                    speedup = pk_t / gp_t
                    speedup_str = f"  |  Speedup: {speedup:.1f}x"

                # Summary table
                rows = []
                if has_plink and "PLINK" in roh_results:
                    pk = roh_results["PLINK"]
                    rows.append([
                        "PLINK 1.9",
                        fmt_val(pk.get("n_samples"), 0),
                        fmt_val(pk.get("n_with_roh"), 0),
                        fmt_val(pk.get("total_segments"), 0),
                        fmt_val(pk.get("mean_froh"), 6),
                        fmt_time(pk.get("time_s")),
                        fmt_mem(pk.get("peak_mem_mb")),
                    ])
                rows.append([
                    "GraphPop",
                    fmt_val(gr.get("n_samples"), 0),
                    fmt_val(gr.get("n_with_roh"), 0),
                    fmt_val(gr.get("total_segments"), 0),
                    fmt_val(gr.get("mean_froh"), 6),
                    fmt_time(gr.get("time_s")),
                    fmt_mem(gr.get("peak_mem_mb")),
                ])

                corr_str = ""
                if "correlation" in roh_results:
                    c = roh_results["correlation"]
                    corr_str = (f"  |  r(total_kb)={fmt_val(c.get('pearson_r_total_kb'), 4)}"
                                f"  MAE(FROH)={fmt_val(c.get('mae_froh'), 6)}")

                print_group_table(
                    f"ROH — chr22 whole-chromosome{speedup_str}{corr_str}",
                    ["Tool", "Samples", "With ROH", "Segments",
                     "Mean FROH", "Time", "Peak Mem"],
                    rows,
                )

                # Top-5 comparison
                if has_plink and "PLINK" in roh_results:
                    print("\n  Top 5 samples by total ROH length:")
                    print(f"  {'Sample':<12s} {'PLINK (kb)':>12s} {'GraphPop (kb)':>14s} {'Diff (kb)':>10s}")
                    print(f"  {'-'*12} {'-'*12} {'-'*14} {'-'*10}")
                    pk_all = roh_results["PLINK"].get("all_results", {})
                    gp_all = gr.get("all_results", {})
                    # Merge top5 from both
                    all_samples = set(
                        list(roh_results["PLINK"].get("top5", {}).keys()) +
                        list(gr.get("top5", {}).keys())
                    )
                    combined = []
                    for s in all_samples:
                        pk_kb = pk_all.get(s, {}).get("total_length_kb", 0)
                        gp_kb = gp_all.get(s, {}).get("total_length_kb", 0)
                        combined.append((s, pk_kb, gp_kb))
                    combined.sort(key=lambda x: max(x[1], x[2]), reverse=True)
                    for s, pk_kb, gp_kb in combined[:5]:
                        diff = gp_kb - pk_kb
                        print(f"  {s:<12s} {pk_kb:>12.1f} {gp_kb:>14.1f} {diff:>+10.1f}")

                # Strip all_results before serializing (too large for JSONL)
                for tool_key in ["PLINK", "GraphPop"]:
                    if tool_key in roh_results:
                        roh_results[tool_key].pop("all_results", None)
                        roh_results[tool_key].pop("top5", None)

                region_record["roh"] = roh_results

            all_results.append(region_record)

    # Write JSONL
    outfile = args.output or str(RESULTS_DIR / "nsl_roh_benchmark.jsonl")
    with open(outfile, "w") as f:
        for r in all_results:
            f.write(json.dumps(r, default=str) + "\n")
    print(f"\nResults written to {outfile}")

    # Final summary
    print(f"\n{'='*90}")
    print("  FINAL SUMMARY")
    print(f"{'='*90}")
    for record in all_results:
        rname = record["region"]
        print(f"\n  Region: {rname}")
        if "nsl" in record:
            nsl = record["nsl"]
            sa_t = nsl.get("scikit-allel", {}).get("time_s")
            gp_t = nsl.get("GraphPop", {}).get("time_s")
            speedup = sa_t / gp_t if (sa_t and gp_t and gp_t > 0) else None
            corr_r = nsl.get("correlation", {}).get("pearson_r")
            print(f"    nSL:  scikit-allel={fmt_time(sa_t)}  "
                  f"GraphPop={fmt_time(gp_t)}  "
                  f"Speedup={fmt_val(speedup, 1) + 'x' if speedup else '—'}  "
                  f"r={fmt_val(corr_r, 4)}")
        if "roh" in record:
            roh = record["roh"]
            pk_t = roh.get("PLINK", {}).get("time_s")
            gp_t = roh.get("GraphPop", {}).get("time_s")
            speedup = pk_t / gp_t if (pk_t and gp_t and gp_t > 0) else None
            corr_r = roh.get("correlation", {}).get("pearson_r_total_kb")
            print(f"    ROH:  PLINK 1.9={fmt_time(pk_t)}  "
                  f"GraphPop={fmt_time(gp_t)}  "
                  f"Speedup={fmt_val(speedup, 1) + 'x' if speedup else '—'}  "
                  f"r(total_kb)={fmt_val(corr_r, 4)}")


if __name__ == "__main__":
    main()
