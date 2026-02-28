#!/usr/bin/env python3
"""
Head-to-head benchmark: GraphPop vs scikit-allel, vcftools, pixy, PLINK2, selscan.

Measures three dimensions per test:
  1. Wall-clock time (seconds)
  2. Peak memory (MB) — RSS for subprocesses, tracemalloc for in-process
  3. Accuracy (numerical agreement vs reference value)

Benchmark groups:
  A. Diversity   (π, θ_W, Tajima's D)   — scikit-allel, vcftools, GraphPop
  B. Differentiation (Fst, Dxy)          — scikit-allel, vcftools, pixy, PLINK2, GraphPop
  C. SFS                                 — scikit-allel, GraphPop
  D. LD (r²)                             — scikit-allel, PLINK2, GraphPop
  E. Selection (iHS, XP-EHH)            — scikit-allel, selscan, GraphPop

Usage:
    conda run -n graphevo python scripts/benchmark-vs-tools.py [--region small|medium|large]

Requires:
    - Neo4j running with GraphPop procedures deployed
    - 1000G chr22 VCF + panel file
    - conda env 'graphevo' with cyvcf2, scikit-allel, numpy, scipy
    - Optional CLI: vcftools, pixy, plink2, selscan (skipped if absent)
"""

import argparse
import json
import os
import re
import resource
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
POP1 = "AFR"
POP2 = "EUR"

NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"

REGIONS = {
    "small":  (16_000_000, 16_100_000),
    "medium": (16_000_000, 16_500_000),
    "large":  (16_000_000, 17_000_000),
}

# ---------------------------------------------------------------------------
# Utility: tool detection
# ---------------------------------------------------------------------------

def which(cmd):
    return shutil.which(cmd)


def tool_available(name):
    """Check if a tool is available."""
    if name == "scikit-allel":
        try:
            import allel
            return True
        except ImportError:
            return False
    return which(name) is not None


TOOLS = {}
for t in ["vcftools", "pixy", "plink2", "selscan"]:
    TOOLS[t] = tool_available(t)
TOOLS["scikit-allel"] = tool_available("scikit-allel")
TOOLS["graphpop"] = True  # always available if Neo4j is up


def check_neo4j():
    try:
        run_cypher("RETURN 1;", timeout=10)
        return True
    except Exception:
        return False

# ---------------------------------------------------------------------------
# Utility: measurement helpers
# ---------------------------------------------------------------------------

class Timer:
    """Context manager for wall-clock timing."""
    def __init__(self):
        self.elapsed = 0.0
    def __enter__(self):
        self._t0 = time.perf_counter()
        return self
    def __exit__(self, *args):
        self.elapsed = time.perf_counter() - self._t0


def measure_python_call(func, *args, **kwargs):
    """Run a Python function, return (result, elapsed_s, peak_mem_mb)."""
    tracemalloc.start()
    # Reset RSS baseline
    t0 = time.perf_counter()
    result = func(*args, **kwargs)
    elapsed = time.perf_counter() - t0
    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    peak_mb = peak_bytes / (1024 * 1024)
    return result, elapsed, peak_mb


def measure_subprocess(cmd, timeout=600, cwd=None):
    """Run a subprocess via /usr/bin/time, return (stdout, elapsed_s, peak_rss_mb)."""
    # Use GNU time for peak RSS
    time_cmd = "/usr/bin/time"
    full_cmd = [time_cmd, "-v"] + cmd

    t0 = time.perf_counter()
    proc = subprocess.run(
        full_cmd, capture_output=True, text=True,
        timeout=timeout, cwd=cwd,
    )
    elapsed = time.perf_counter() - t0

    # Parse peak RSS from /usr/bin/time stderr
    peak_rss_mb = 0.0
    for line in proc.stderr.split("\n"):
        m = re.search(r"Maximum resident set size.*?:\s*(\d+)", line)
        if m:
            peak_rss_mb = int(m.group(1)) / 1024.0  # kbytes → MB
            break

    if proc.returncode != 0:
        # Separate /usr/bin/time output from actual stderr
        actual_stderr = "\n".join(
            l for l in proc.stderr.split("\n")
            if "Maximum resident" not in l and "Command being timed" not in l
            and "wall clock" not in l and "Elapsed" not in l
            and "Minor" not in l and "Major" not in l
            and "Voluntary" not in l and "Involuntary" not in l
            and "File system" not in l and "Socket" not in l
            and "Signals" not in l and "Swaps" not in l
            and "Page size" not in l and "Exit status" not in l
            and "Percent of CPU" not in l and "User time" not in l
            and "System time" not in l and "Average" not in l
        ).strip()
        raise RuntimeError(f"Command failed (rc={proc.returncode}): {actual_stderr[:500]}")

    return proc.stdout, elapsed, peak_rss_mb


def run_cypher(query, timeout=600):
    """Execute Cypher via cypher-shell, return (output, elapsed_s, peak_rss_mb)."""
    cmd = [
        "cypher-shell",
        "-u", NEO4J_USER,
        "-p", NEO4J_PASS,
        "--format", "plain",
        query,
    ]
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


# ---------------------------------------------------------------------------
# Data loading helpers (shared across Python-based tools)
# ---------------------------------------------------------------------------

_panel_cache = None

def load_panel():
    global _panel_cache
    if _panel_cache is not None:
        return _panel_cache
    sample_pop = {}
    with open(PANEL) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                sample_pop[parts[0]] = parts[2]  # super_pop
    _panel_cache = sample_pop
    return sample_pop


def load_genotype_data(start, end, pop):
    """Load genotype data from VCF for a region and population.
    Returns (geno_012 array shape (n_variants, n_samples),
             positions array, variant_ids list, ac array, an_total).
    """
    import cyvcf2
    sample_pop = load_panel()

    vcf = cyvcf2.VCF(VCF)
    all_samples = list(vcf.samples)
    pop_indices = [i for i, s in enumerate(all_samples) if sample_pop.get(s) == pop]

    positions = []
    variant_ids = []
    genotypes_012 = []  # dosage: 0, 1, 2
    haplotypes_0 = []   # phase 0
    haplotypes_1 = []   # phase 1

    region = f"{CHR}:{start}-{end}"
    for v in vcf(region):
        if not v.is_snp or len(v.ALT) != 1:
            continue
        gt = v.genotypes
        dos = []
        h0 = []
        h1 = []
        for idx in pop_indices:
            g = gt[idx]
            dos.append(g[0] + g[1])
            h0.append(g[0])
            h1.append(g[1])
        dos = np.array(dos, dtype=np.int8)
        # skip monomorphic
        if dos.min() == dos.max():
            continue
        positions.append(v.POS)
        variant_ids.append(f"{v.CHROM}:{v.POS}:{v.REF}:{v.ALT[0]}")
        genotypes_012.append(dos)
        haplotypes_0.append(np.array(h0, dtype=np.int8))
        haplotypes_1.append(np.array(h1, dtype=np.int8))

    vcf.close()

    return {
        "positions": np.array(positions),
        "variant_ids": variant_ids,
        "genotypes": np.array(genotypes_012),
        "haplotypes_0": np.array(haplotypes_0),
        "haplotypes_1": np.array(haplotypes_1),
        "n_samples": len(pop_indices),
        "n_variants": len(positions),
    }


def load_genotype_data_pair(start, end, pop1, pop2):
    """Load genotype data for TWO populations from the same VCF pass.
    Keeps all biallelic SNPs (no per-pop monomorphic filtering) so variant
    sets are aligned between populations.
    """
    import cyvcf2
    sample_pop = load_panel()

    vcf = cyvcf2.VCF(VCF)
    all_samples = list(vcf.samples)
    idx1 = [i for i, s in enumerate(all_samples) if sample_pop.get(s) == pop1]
    idx2 = [i for i, s in enumerate(all_samples) if sample_pop.get(s) == pop2]

    positions = []
    variant_ids = []
    geno1 = []
    geno2 = []
    hap1_0, hap1_1 = [], []
    hap2_0, hap2_1 = [], []

    region = f"{CHR}:{start}-{end}"
    for v in vcf(region):
        if not v.is_snp or len(v.ALT) != 1:
            continue
        gt = v.genotypes
        d1, d2 = [], []
        h10, h11, h20, h21 = [], [], [], []
        for i in idx1:
            g = gt[i]; d1.append(g[0]+g[1]); h10.append(g[0]); h11.append(g[1])
        for i in idx2:
            g = gt[i]; d2.append(g[0]+g[1]); h20.append(g[0]); h21.append(g[1])
        d1a = np.array(d1, dtype=np.int8)
        d2a = np.array(d2, dtype=np.int8)
        # Keep if segregating in at least one population
        if d1a.min() == d1a.max() and d2a.min() == d2a.max():
            continue
        positions.append(v.POS)
        variant_ids.append(f"{v.CHROM}:{v.POS}:{v.REF}:{v.ALT[0]}")
        geno1.append(d1a); geno2.append(d2a)
        hap1_0.append(np.array(h10, dtype=np.int8))
        hap1_1.append(np.array(h11, dtype=np.int8))
        hap2_0.append(np.array(h20, dtype=np.int8))
        hap2_1.append(np.array(h21, dtype=np.int8))

    vcf.close()
    pos = np.array(positions)
    return {
        "positions": pos,
        "variant_ids": variant_ids,
        "n_variants": len(positions),
        "pop1": {
            "genotypes": np.array(geno1), "n_samples": len(idx1),
            "haplotypes_0": np.array(hap1_0), "haplotypes_1": np.array(hap1_1),
        },
        "pop2": {
            "genotypes": np.array(geno2), "n_samples": len(idx2),
            "haplotypes_0": np.array(hap2_0), "haplotypes_1": np.array(hap2_1),
        },
    }


def write_pop_keep_file(tmpdir, pop):
    """Write a keep file for vcftools (one sample per line)."""
    sample_pop = load_panel()
    path = os.path.join(tmpdir, f"{pop}_samples.txt")
    with open(path, "w") as f:
        for s, p in sample_pop.items():
            if p == pop:
                f.write(s + "\n")
    return path


def write_pixy_popfile(tmpdir):
    """Write pixy population file (tab-delimited: sample_id, population)."""
    sample_pop = load_panel()
    path = os.path.join(tmpdir, "populations.txt")
    with open(path, "w") as f:
        for s, p in sample_pop.items():
            if p in (POP1, POP2):
                f.write(f"{s}\t{p}\n")
    return path


# ---------------------------------------------------------------------------
# Group A: Diversity (π, θ_W, Tajima's D)
# ---------------------------------------------------------------------------

def diversity_scikit_allel(start, end, pop):
    """Compute π, θ_W, Tajima's D using scikit-allel.

    Uses ALL biallelic SNP variants (including pop-monomorphic) so the
    denominator matches GraphPop's nVariants-based normalization.
    Reports π = sum(pi_i) / n_variants (per-variant-site average).
    """
    import allel

    def _compute():
        # Load ALL biallelic SNPs (no monomorphic filtering)
        import cyvcf2
        sample_pop = load_panel()
        vcf = cyvcf2.VCF(VCF)
        all_samples = list(vcf.samples)
        pop_idx = [i for i, s in enumerate(all_samples) if sample_pop.get(s) == pop]
        n_samples = len(pop_idx)
        n = 2 * n_samples

        positions = []
        ac_alt_list = []
        region = f"{CHR}:{start}-{end}"
        for v in vcf(region):
            if not v.is_snp or len(v.ALT) != 1:
                continue
            gt = v.genotypes
            ac = sum(gt[i][0] + gt[i][1] for i in pop_idx)
            positions.append(v.POS)
            ac_alt_list.append(ac)
        vcf.close()

        pos = np.array(positions)
        ac_alt = np.array(ac_alt_list)
        ac_ref = n - ac_alt
        ac = np.column_stack([ac_ref, ac_alt])
        n_total = len(positions)
        n_seg = int(np.sum((ac_alt > 0) & (ac_alt < n)))

        pi_arr = allel.mean_pairwise_difference(ac)
        pi_sum = float(np.nansum(pi_arr))
        pi_per_variant = pi_sum / n_total if n_total > 0 else 0.0

        # θ_W and Tajima's D — allel computes over segregating sites
        seg_mask = (ac_alt > 0) & (ac_alt < n)
        ac_seg = ac[seg_mask]
        pos_seg = pos[seg_mask]

        # θ_W = S / a_n, normalized by n_total
        from scipy.special import digamma
        a_n = sum(1.0/k for k in range(1, n))
        theta_w = n_seg / a_n / n_total if n_total > 0 and a_n > 0 else 0.0

        tajima_d = allel.tajima_d(ac_seg, pos=pos_seg) if n_seg > 3 else float('nan')

        return {
            "pi": float(pi_per_variant),
            "theta_w": float(theta_w),
            "tajima_d": float(tajima_d),
            "n_variants": n_total,
            "n_segregating": n_seg,
        }

    result, elapsed, peak_mb = measure_python_call(_compute)
    return {**result, "time_s": elapsed, "peak_mem_mb": peak_mb}


def diversity_vcftools(start, end, pop, tmpdir):
    """Compute π, Tajima's D using vcftools."""
    keep_file = write_pop_keep_file(tmpdir, pop)
    outprefix = os.path.join(tmpdir, "vcftools_div")

    # π (per-site)
    cmd_pi = [
        "vcftools", "--gzvcf", VCF,
        "--chr", CHR, "--from-bp", str(start), "--to-bp", str(end),
        "--keep", keep_file,
        "--site-pi", "--out", outprefix + "_pi",
    ]
    try:
        _, time_pi, mem_pi = measure_subprocess(cmd_pi)
    except RuntimeError:
        time_pi, mem_pi = 0.0, 0.0

    # Parse π output
    pi_file = outprefix + "_pi.sites.pi"
    pi_vals = []
    if os.path.exists(pi_file):
        with open(pi_file) as f:
            f.readline()  # header
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    pi_vals.append(float(parts[2]))

    # Tajima's D (windowed — use large window to get single value)
    window = end - start
    cmd_td = [
        "vcftools", "--gzvcf", VCF,
        "--chr", CHR, "--from-bp", str(start), "--to-bp", str(end),
        "--keep", keep_file,
        "--TajimaD", str(window), "--out", outprefix + "_td",
    ]
    try:
        _, time_td, mem_td = measure_subprocess(cmd_td)
    except RuntimeError:
        time_td, mem_td = 0.0, 0.0

    td_file = outprefix + "_td.Tajima.D"
    tajima_d = None
    if os.path.exists(td_file):
        with open(td_file) as f:
            f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 4 and parts[3] != "nan":
                    tajima_d = float(parts[3])

    return {
        "pi": float(np.mean(pi_vals)) if pi_vals else None,
        "tajima_d": tajima_d,
        "n_variants": len(pi_vals),
        "time_s": time_pi + time_td,
        "peak_mem_mb": max(mem_pi, mem_td),
    }


def diversity_graphpop(start, end, pop):
    """Query graphpop.diversity."""
    query = (
        f"CALL graphpop.diversity('{CHR}', {start}, {end}, '{pop}') "
        "YIELD pi, theta_w, tajima_d, n_variants, n_segregating "
        "RETURN pi, theta_w, tajima_d, n_variants, n_segregating;"
    )
    out, elapsed, peak_mb = run_cypher(query)
    result = parse_cypher_row(out)
    return {
        "pi": result.get("pi"),
        "theta_w": result.get("theta_w"),
        "tajima_d": result.get("tajima_d"),
        "n_variants": result.get("n_variants"),
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


# ---------------------------------------------------------------------------
# Group B: Differentiation (Fst, Dxy)
# ---------------------------------------------------------------------------

def fst_scikit_allel(start, end, pop1, pop2):
    """Compute Hudson's Fst using scikit-allel."""
    import allel

    def _compute():
        data = load_genotype_data_pair(start, end, pop1, pop2)
        d1 = data["pop1"]
        d2 = data["pop2"]

        n1 = 2 * d1["n_samples"]
        n2 = 2 * d2["n_samples"]

        ac1_alt = d1["genotypes"].sum(axis=1)
        ac1_ref = n1 - ac1_alt
        ac1 = np.column_stack([ac1_ref, ac1_alt])

        ac2_alt = d2["genotypes"].sum(axis=1)
        ac2_ref = n2 - ac2_alt
        ac2 = np.column_stack([ac2_ref, ac2_alt])

        num, den = allel.hudson_fst(ac1, ac2)
        fst = np.sum(num) / np.sum(den)

        return {"fst": float(fst), "n_variants": data["n_variants"]}

    result, elapsed, peak_mb = measure_python_call(_compute)
    return {**result, "time_s": elapsed, "peak_mem_mb": peak_mb}


def fst_vcftools(start, end, pop1, pop2, tmpdir):
    """Compute Weir & Cockerham Fst using vcftools.

    NOTE: vcftools uses Weir & Cockerham Fst, not Hudson's. Values will
    differ from GraphPop/scikit-allel (both use Hudson's estimator).
    """
    keep1 = write_pop_keep_file(tmpdir, pop1)
    keep2 = write_pop_keep_file(tmpdir, pop2)
    outprefix = os.path.join(tmpdir, "vcftools_fst")

    # vcftools prints weighted Fst to stderr — use direct subprocess
    # to capture stderr separately from /usr/bin/time
    cmd = [
        "/usr/bin/time", "-v",
        "vcftools", "--gzvcf", VCF,
        "--chr", CHR, "--from-bp", str(start), "--to-bp", str(end),
        "--weir-fst-pop", keep1,
        "--weir-fst-pop", keep2,
        "--out", outprefix,
    ]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    elapsed = time.perf_counter() - t0

    # Parse peak RSS from /usr/bin/time in stderr
    peak_rss_mb = 0.0
    fst_val = None
    for line in proc.stderr.split("\n"):
        m = re.search(r"Maximum resident set size.*?:\s*(\d+)", line)
        if m:
            peak_rss_mb = int(m.group(1)) / 1024.0
        m = re.search(r"[Ww]eighted Fst estimate:\s*([-\d.eE+]+)", line)
        if m:
            fst_val = float(m.group(1))

    return {
        "fst_wc": fst_val,
        "time_s": elapsed,
        "peak_mem_mb": peak_rss_mb,
    }


def fst_pixy(start, end, pop1, pop2, tmpdir):
    """Compute π, Fst, Dxy using pixy."""
    popfile = write_pixy_popfile(tmpdir)
    window = end - start  # single window

    cmd = [
        "pixy",
        "--stats", "pi", "fst", "dxy",
        "--vcf", VCF,
        "--populations", popfile,
        "--chromosomes", CHR,
        "--window_size", str(window),
        "--interval_start", str(start),
        "--interval_end", str(end),
        "--output_folder", tmpdir,
        "--output_prefix", "pixy_bench",
    ]
    _, elapsed, peak_mb = measure_subprocess(cmd, timeout=600)

    results = {"time_s": elapsed, "peak_mem_mb": peak_mb}

    # Parse pixy output files
    for stat in ["pi", "fst", "dxy"]:
        outfile = os.path.join(tmpdir, f"pixy_bench_{stat}.txt")
        if os.path.exists(outfile):
            with open(outfile) as f:
                f.readline()  # header
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 5 and parts[-1] != "NA":
                        key = f"{stat}_{parts[0]}_{parts[1]}" if len(parts) > 5 else stat
                        # pixy output varies: pop, chrom, window_pos_1, window_pos_2, value
                        # For fst/dxy: pop1, pop2, chrom, ..., value
                        results[f"pixy_{stat}"] = float(parts[-1])

    return results


def fst_plink2(start, end, pop1, pop2, tmpdir):
    """Compute Hudson's Fst using PLINK2."""
    # PLINK2 --fst needs a phenotype file: #IID POP
    # (PLINK2 sets FID=0 when loading from VCF, so match on IID only)
    pheno_file = os.path.join(tmpdir, "plink2_pop.txt")
    sample_pop = load_panel()
    with open(pheno_file, "w") as f:
        f.write("#IID\tPOP\n")
        for s, p in sample_pop.items():
            if p in (pop1, pop2):
                f.write(f"{s}\t{p}\n")

    cmd = [
        "plink2",
        "--vcf", VCF,
        "--chr", CHR,
        "--from-bp", str(start), "--to-bp", str(end),
        "--pheno", pheno_file,
        "--fst", "POP", "method=hudson",
        "--allow-extra-chr",
        "--out", os.path.join(tmpdir, "plink2_fst"),
    ]
    try:
        _, elapsed, peak_mb = measure_subprocess(cmd, timeout=600)
    except RuntimeError as e:
        return {"error": str(e), "time_s": 0, "peak_mem_mb": 0}

    # Parse Fst from .fst.summary (#POP1 POP2 HUDSON_FST)
    fst_val = None
    summary = os.path.join(tmpdir, "plink2_fst.fst.summary")
    if os.path.exists(summary):
        with open(summary) as f:
            header = f.readline().strip().split("\t")
            data_line = f.readline()
            if data_line:
                parts = data_line.strip().split("\t")
                # Last column is the Fst value
                try:
                    fst_val = float(parts[-1])
                except (ValueError, IndexError):
                    pass

    return {
        "fst_hudson": fst_val,
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


def fst_graphpop(start, end, pop1, pop2):
    """Query graphpop.divergence for Fst and Dxy."""
    query = (
        f"CALL graphpop.divergence('{CHR}', {start}, {end}, '{pop1}', '{pop2}') "
        "YIELD fst_hudson, dxy, da, n_variants "
        "RETURN fst_hudson, dxy, da, n_variants;"
    )
    out, elapsed, peak_mb = run_cypher(query)
    result = parse_cypher_row(out)
    return {
        "fst": result.get("fst_hudson"),
        "dxy": result.get("dxy"),
        "n_variants": result.get("n_variants"),
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


# ---------------------------------------------------------------------------
# Group C: SFS
# ---------------------------------------------------------------------------

def sfs_scikit_allel(start, end, pop):
    """Compute SFS using scikit-allel."""
    import allel

    def _compute():
        data = load_genotype_data(start, end, pop)
        n = 2 * data["n_samples"]
        ac_alt = data["genotypes"].sum(axis=1)
        ac_ref = n - ac_alt
        ac = np.column_stack([ac_ref, ac_alt])
        sfs = allel.sfs(ac[:, 1], n=n)
        return {"sfs_len": len(sfs), "sfs_sum": int(sfs.sum()), "n_variants": data["n_variants"]}

    result, elapsed, peak_mb = measure_python_call(_compute)
    return {**result, "time_s": elapsed, "peak_mem_mb": peak_mb}


def sfs_graphpop(start, end, pop):
    """Query graphpop.sfs."""
    query = (
        f"CALL graphpop.sfs('{CHR}', {start}, {end}, '{pop}') "
        "YIELD n_variants RETURN count(*) AS n_bins, sum(n_variants) AS total;"
    )
    out, elapsed, peak_mb = run_cypher(query)
    result = parse_cypher_row(out)
    return {
        "n_bins": result.get("n_bins"),
        "total": result.get("total"),
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


# ---------------------------------------------------------------------------
# Group D: LD (r²)
# ---------------------------------------------------------------------------

def ld_scikit_allel(start, end, pop):
    """Compute pairwise r² using scikit-allel."""
    import allel

    def _compute():
        data = load_genotype_data(start, end, pop)
        g = data["genotypes"].astype(float)
        r = allel.rogers_huff_r(g)
        from scipy.spatial.distance import squareform
        r2 = squareform(r) ** 2
        n_pairs = int(np.sum(r2 >= 0.2) / 2)
        return {"n_pairs_r2_02": n_pairs, "n_variants": data["n_variants"]}

    result, elapsed, peak_mb = measure_python_call(_compute)
    return {**result, "time_s": elapsed, "peak_mem_mb": peak_mb}


def ld_plink2(start, end, pop, tmpdir):
    """Compute LD using PLINK2."""
    keep_file = write_pop_keep_file(tmpdir, pop)

    # PLINK2 --keep format: #IID (match on IID since FID=0 from VCF)
    keep_plink = os.path.join(tmpdir, f"{pop}_keep_plink.txt")
    with open(keep_file) as fin, open(keep_plink, "w") as fout:
        fout.write("#IID\n")
        for line in fin:
            s = line.strip()
            if s:
                fout.write(f"{s}\n")

    cmd = [
        "plink2",
        "--vcf", VCF,
        "--chr", CHR,
        "--from-bp", str(start), "--to-bp", str(end),
        "--keep", keep_plink,
        "--r2",
        "--ld-window-r2", "0.2",
        "--allow-extra-chr",
        "--out", os.path.join(tmpdir, "plink2_ld"),
    ]
    try:
        _, elapsed, peak_mb = measure_subprocess(cmd, timeout=600)
    except RuntimeError as e:
        return {"error": str(e), "time_s": 0, "peak_mem_mb": 0}

    # Count LD pairs from output
    n_pairs = 0
    for ext in [".vcor", ".ld"]:
        ld_file = os.path.join(tmpdir, f"plink2_ld{ext}")
        if os.path.exists(ld_file):
            with open(ld_file) as f:
                for _ in f:
                    n_pairs += 1
            n_pairs -= 1  # subtract header
            break

    return {
        "n_pairs_r2_02": n_pairs,
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


def ld_graphpop(start, end, pop):
    """Run graphpop.ld and count pairs."""
    # First clean any existing LD edges
    try:
        run_cypher("MATCH ()-[r:LD]->() DELETE r;", timeout=60)
    except Exception:
        pass

    query = (
        f"CALL graphpop.ld('{CHR}', {start}, {end}, '{pop}', 500000, 0.2) "
        "YIELD variant1 RETURN count(*) AS n_pairs;"
    )
    out, elapsed, peak_mb = run_cypher(query, timeout=600)
    result = parse_cypher_row(out)

    # Clean up
    try:
        run_cypher("MATCH ()-[r:LD]->() DELETE r;", timeout=60)
    except Exception:
        pass

    return {
        "n_pairs_r2_02": result.get("n_pairs"),
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


# ---------------------------------------------------------------------------
# Group E: Selection (iHS, XP-EHH)
# ---------------------------------------------------------------------------

def ihs_scikit_allel(start, end, pop):
    """Compute iHS using scikit-allel."""
    import allel

    def _compute():
        data = load_genotype_data(start, end, pop)

        # Build haplotype array: (n_variants, n_haplotypes)
        h = np.concatenate([data["haplotypes_0"], data["haplotypes_1"]], axis=1)
        hap = allel.HaplotypeArray(h)

        # Allele counts for filtering
        ac = hap.count_alleles()

        # Filter: biallelic, common (AF >= 0.05)
        is_seg = ac.is_segregating()
        af = ac.to_frequencies()[:, 1]
        flt = is_seg & (af >= 0.05) & (af <= 0.95)

        hap_flt = hap[flt]
        pos_flt = data["positions"][flt]

        if len(pos_flt) < 10:
            return {"n_variants": 0, "ihs_mean": None}

        ihs_raw = allel.ihs(hap_flt, pos_flt, include_edges=True)
        # Standardize — may fail on small regions with too few variants
        try:
            ihs_std, bins = allel.standardize_by_allele_count(
                ihs_raw, ac[flt][:, 1], diagnostics=False
            )
        except (ValueError, ZeroDivisionError):
            # Fall back to simple z-score standardization
            valid_raw = ihs_raw[~np.isnan(ihs_raw)]
            if len(valid_raw) > 1:
                ihs_std = (ihs_raw - np.nanmean(ihs_raw)) / np.nanstd(ihs_raw)
            else:
                ihs_std = ihs_raw
        valid = ~np.isnan(ihs_std)
        return {
            "n_variants": int(valid.sum()),
            "ihs_mean": float(np.nanmean(ihs_std)),
            "ihs_std": float(np.nanstd(ihs_std)),
        }

    result, elapsed, peak_mb = measure_python_call(_compute)
    return {**result, "time_s": elapsed, "peak_mem_mb": peak_mb}


def ihs_selscan(start, end, pop, tmpdir):
    """Compute iHS using selscan."""
    # selscan needs: --vcf, --map, --out
    # Write a genetic map (physical positions as genetic positions, 1cM/Mb)
    data = load_genotype_data(start, end, pop)
    map_file = os.path.join(tmpdir, "selscan.map")
    with open(map_file, "w") as f:
        for pos in data["positions"]:
            f.write(f"{CHR}\t.\t{pos/1e6:.6f}\t{pos}\n")

    # Write a minimal VCF with only the target population
    keep_file = write_pop_keep_file(tmpdir, pop)

    # Use bcftools to subset (if available) or write manually
    pop_vcf = os.path.join(tmpdir, "selscan_pop.vcf.gz")
    bcf = which("bcftools")
    if bcf:
        subprocess.run([
            bcf, "view", "-S", keep_file,
            "-r", f"{CHR}:{start}-{end}",
            VCF, "-Oz", "-o", pop_vcf
        ], check=True, capture_output=True)
    else:
        return {"error": "bcftools required for VCF subsetting", "time_s": 0, "peak_mem_mb": 0}

    cmd = [
        "selscan", "--ihs",
        "--vcf", pop_vcf,
        "--map", map_file,
        "--out", os.path.join(tmpdir, "selscan_ihs"),
    ]
    try:
        _, elapsed, peak_mb = measure_subprocess(cmd, timeout=600)
    except RuntimeError as e:
        return {"error": str(e), "time_s": 0, "peak_mem_mb": 0}

    # Count results
    out_file = os.path.join(tmpdir, "selscan_ihs.ihs.out")
    n_variants = 0
    if os.path.exists(out_file):
        with open(out_file) as f:
            for _ in f:
                n_variants += 1

    return {
        "n_variants": n_variants,
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


def ihs_graphpop(start, end, pop):
    """Run graphpop.ihs."""
    query = (
        f"CALL graphpop.ihs('{CHR}', '{pop}', "
        f"{{start: {start}, end: {end}, min_af: 0.05}}) "
        "YIELD variantId, ihs "
        "RETURN count(*) AS n_variants, avg(ihs) AS ihs_mean, stDev(ihs) AS ihs_std;"
    )
    out, elapsed, peak_mb = run_cypher(query, timeout=600)
    result = parse_cypher_row(out)
    return {
        "n_variants": result.get("n_variants"),
        "ihs_mean": result.get("ihs_mean"),
        "ihs_std": result.get("ihs_std"),
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


def xpehh_scikit_allel(start, end, pop1, pop2):
    """Compute XP-EHH using scikit-allel."""
    import allel

    def _compute():
        data = load_genotype_data_pair(start, end, pop1, pop2)
        d1 = data["pop1"]
        d2 = data["pop2"]

        h1 = np.concatenate([d1["haplotypes_0"], d1["haplotypes_1"]], axis=1)
        h2 = np.concatenate([d2["haplotypes_0"], d2["haplotypes_1"]], axis=1)

        hap1 = allel.HaplotypeArray(h1)
        hap2 = allel.HaplotypeArray(h2)

        pos = data["positions"]
        if len(pos) < 10:
            return {"n_variants": 0}

        xp = allel.xpehh(hap1, hap2, pos, include_edges=True)
        valid = ~np.isnan(xp)
        return {
            "n_variants": int(valid.sum()),
            "xpehh_mean": float(np.nanmean(xp)),
            "xpehh_std": float(np.nanstd(xp)),
        }

    result, elapsed, peak_mb = measure_python_call(_compute)
    return {**result, "time_s": elapsed, "peak_mem_mb": peak_mb}


def xpehh_selscan(start, end, pop1, pop2, tmpdir):
    """Compute XP-EHH using selscan."""
    data = load_genotype_data(start, end, pop1)
    map_file = os.path.join(tmpdir, "selscan_xp.map")
    with open(map_file, "w") as f:
        for pos in data["positions"]:
            f.write(f"{CHR}\t.\t{pos/1e6:.6f}\t{pos}\n")

    keep1 = write_pop_keep_file(tmpdir, pop1)
    keep2 = write_pop_keep_file(tmpdir, pop2)

    bcf = which("bcftools")
    if not bcf:
        return {"error": "bcftools required", "time_s": 0, "peak_mem_mb": 0}

    vcf1 = os.path.join(tmpdir, "xp_pop1.vcf.gz")
    vcf2 = os.path.join(tmpdir, "xp_pop2.vcf.gz")
    region = f"{CHR}:{start}-{end}"

    subprocess.run([bcf, "view", "-S", keep1, "-r", region, VCF, "-Oz", "-o", vcf1],
                   check=True, capture_output=True)
    subprocess.run([bcf, "view", "-S", keep2, "-r", region, VCF, "-Oz", "-o", vcf2],
                   check=True, capture_output=True)

    cmd = [
        "selscan", "--xpehh",
        "--vcf", vcf1, "--vcf-ref", vcf2,
        "--map", map_file,
        "--out", os.path.join(tmpdir, "selscan_xpehh"),
    ]
    try:
        _, elapsed, peak_mb = measure_subprocess(cmd, timeout=600)
    except RuntimeError as e:
        return {"error": str(e), "time_s": 0, "peak_mem_mb": 0}

    out_file = os.path.join(tmpdir, "selscan_xpehh.xpehh.out")
    n_variants = 0
    if os.path.exists(out_file):
        with open(out_file) as f:
            for _ in f:
                n_variants += 1

    return {
        "n_variants": n_variants,
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


def xpehh_graphpop(start, end, pop1, pop2):
    """Run graphpop.xpehh."""
    query = (
        f"CALL graphpop.xpehh('{CHR}', '{pop1}', '{pop2}', "
        f"{{start: {start}, end: {end}}}) "
        "YIELD variantId, xpehh "
        "RETURN count(*) AS n_variants, avg(xpehh) AS xpehh_mean, stDev(xpehh) AS xpehh_std;"
    )
    out, elapsed, peak_mb = run_cypher(query, timeout=600)
    result = parse_cypher_row(out)
    return {
        "n_variants": result.get("n_variants"),
        "xpehh_mean": result.get("xpehh_mean"),
        "xpehh_std": result.get("xpehh_std"),
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def fmt_time(t):
    if t is None:
        return "—"
    return f"{t:.3f}s"


def fmt_mem(m):
    if m is None:
        return "—"
    return f"{m:.0f}MB"


def fmt_val(v, decimals=6):
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.{decimals}f}"
    return str(v)


def print_group_table(title, headers, rows):
    """Print a formatted comparison table."""
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}")

    # Compute column widths
    all_rows = [headers] + rows
    widths = [max(len(str(row[i])) for row in all_rows) for i in range(len(headers))]

    fmt = "  ".join(f"{{:<{w}}}" for w in widths)
    print(fmt.format(*headers))
    print("-" * sum(widths + [2 * (len(widths) - 1)]))
    for row in rows:
        print(fmt.format(*[str(x) for x in row]))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="GraphPop vs external tools benchmark")
    parser.add_argument("--region", default="small", choices=REGIONS.keys(),
                        help="Region size to benchmark")
    parser.add_argument("--all-regions", action="store_true",
                        help="Run all region sizes")
    parser.add_argument("--output", default=None,
                        help="Output JSONL file path")
    args = parser.parse_args()

    regions_to_run = list(REGIONS.keys()) if args.all_regions else [args.region]

    print("=" * 80)
    print("  GraphPop Head-to-Head Benchmark")
    print("=" * 80)
    print()
    print("  Tool availability:")
    for tool, avail in TOOLS.items():
        status = "FOUND" if avail else "NOT FOUND (will skip)"
        print(f"    {tool:16s}  {status}")
    print()

    if not check_neo4j():
        print("ERROR: Neo4j is not responding. Start Neo4j first.", file=sys.stderr)
        sys.exit(1)

    all_results = []

    for region_name in regions_to_run:
        start, end = REGIONS[region_name]
        print(f"\n{'#'*80}")
        print(f"  Region: {region_name} — {CHR}:{start:,}-{end:,}")
        print(f"{'#'*80}")

        with tempfile.TemporaryDirectory(prefix="graphpop_bench_") as tmpdir:

            # ---- Group A: Diversity ----
            print(f"\n--- Group A: Diversity (π, θ_W, Tajima's D) — pop={POP1} ---")
            div_results = {}

            if TOOLS["scikit-allel"]:
                print("  scikit-allel...", end="", flush=True)
                div_results["scikit-allel"] = diversity_scikit_allel(start, end, POP1)
                print(f" {fmt_time(div_results['scikit-allel']['time_s'])}")

            if TOOLS["vcftools"]:
                print("  vcftools...", end="", flush=True)
                div_results["vcftools"] = diversity_vcftools(start, end, POP1, tmpdir)
                print(f" {fmt_time(div_results['vcftools']['time_s'])}")

            print("  GraphPop...", end="", flush=True)
            div_results["GraphPop"] = diversity_graphpop(start, end, POP1)
            print(f" {fmt_time(div_results['GraphPop']['time_s'])}")

            rows = []
            for tool, r in div_results.items():
                rows.append([
                    tool,
                    fmt_val(r.get("pi")), fmt_val(r.get("theta_w")),
                    fmt_val(r.get("tajima_d")),
                    fmt_time(r.get("time_s")), fmt_mem(r.get("peak_mem_mb")),
                ])
            print_group_table(
                f"Diversity — {region_name} ({CHR}:{start:,}-{end:,})",
                ["Tool", "π", "θ_W", "Tajima's D", "Time", "Peak Mem"],
                rows,
            )

            # ---- Group B: Differentiation ----
            print(f"\n--- Group B: Differentiation (Fst, Dxy) — {POP1} vs {POP2} ---")
            fst_results = {}

            if TOOLS["scikit-allel"]:
                print("  scikit-allel...", end="", flush=True)
                fst_results["scikit-allel"] = fst_scikit_allel(start, end, POP1, POP2)
                print(f" {fmt_time(fst_results['scikit-allel']['time_s'])}")

            if TOOLS["vcftools"]:
                print("  vcftools...", end="", flush=True)
                fst_results["vcftools"] = fst_vcftools(start, end, POP1, POP2, tmpdir)
                print(f" {fmt_time(fst_results['vcftools']['time_s'])}")

            if TOOLS["pixy"]:
                print("  pixy...", end="", flush=True)
                fst_results["pixy"] = fst_pixy(start, end, POP1, POP2, tmpdir)
                print(f" {fmt_time(fst_results['pixy']['time_s'])}")

            if TOOLS["plink2"]:
                print("  PLINK2...", end="", flush=True)
                fst_results["PLINK2"] = fst_plink2(start, end, POP1, POP2, tmpdir)
                print(f" {fmt_time(fst_results['PLINK2']['time_s'])}")

            print("  GraphPop...", end="", flush=True)
            fst_results["GraphPop"] = fst_graphpop(start, end, POP1, POP2)
            print(f" {fmt_time(fst_results['GraphPop']['time_s'])}")

            rows = []
            for tool, r in fst_results.items():
                fst_val = r.get("fst") or r.get("fst_wc") or r.get("fst_hudson") or r.get("pixy_fst")
                dxy_val = r.get("dxy") or r.get("pixy_dxy")
                rows.append([
                    tool, fmt_val(fst_val), fmt_val(dxy_val),
                    fmt_time(r.get("time_s")), fmt_mem(r.get("peak_mem_mb")),
                ])
            print_group_table(
                f"Differentiation — {region_name} ({POP1} vs {POP2})",
                ["Tool", "Fst", "Dxy", "Time", "Peak Mem"],
                rows,
            )

            # ---- Group C: SFS ----
            print(f"\n--- Group C: SFS — pop={POP2} ---")
            sfs_results = {}

            if TOOLS["scikit-allel"]:
                print("  scikit-allel...", end="", flush=True)
                sfs_results["scikit-allel"] = sfs_scikit_allel(start, end, POP2)
                print(f" {fmt_time(sfs_results['scikit-allel']['time_s'])}")

            print("  GraphPop...", end="", flush=True)
            sfs_results["GraphPop"] = sfs_graphpop(start, end, POP2)
            print(f" {fmt_time(sfs_results['GraphPop']['time_s'])}")

            rows = []
            for tool, r in sfs_results.items():
                rows.append([
                    tool, fmt_val(r.get("sfs_sum") or r.get("total"), 0),
                    fmt_time(r.get("time_s")), fmt_mem(r.get("peak_mem_mb")),
                ])
            print_group_table(
                f"SFS — {region_name}",
                ["Tool", "Total Variants", "Time", "Peak Mem"],
                rows,
            )

            # ---- Group D: LD ----
            # Only run LD on small/medium (large is too expensive for pairwise)
            if region_name != "large":
                print(f"\n--- Group D: LD (r² ≥ 0.2) — pop={POP2} ---")
                ld_results = {}

                if TOOLS["scikit-allel"]:
                    print("  scikit-allel...", end="", flush=True)
                    ld_results["scikit-allel"] = ld_scikit_allel(start, end, POP2)
                    print(f" {fmt_time(ld_results['scikit-allel']['time_s'])}")

                if TOOLS["plink2"]:
                    print("  PLINK2...", end="", flush=True)
                    ld_results["PLINK2"] = ld_plink2(start, end, POP2, tmpdir)
                    print(f" {fmt_time(ld_results['PLINK2']['time_s'])}")

                print("  GraphPop...", end="", flush=True)
                ld_results["GraphPop"] = ld_graphpop(start, end, POP2)
                print(f" {fmt_time(ld_results['GraphPop']['time_s'])}")

                rows = []
                for tool, r in ld_results.items():
                    rows.append([
                        tool, fmt_val(r.get("n_pairs_r2_02"), 0),
                        fmt_time(r.get("time_s")), fmt_mem(r.get("peak_mem_mb")),
                    ])
                print_group_table(
                    f"LD — {region_name}",
                    ["Tool", "Pairs (r²≥0.2)", "Time", "Peak Mem"],
                    rows,
                )

            # ---- Group E: Selection ----
            # Only on small/medium (expensive)
            if region_name != "large":
                print(f"\n--- Group E: iHS — pop={POP2} ---")
                ihs_results = {}

                if TOOLS["scikit-allel"]:
                    print("  scikit-allel...", end="", flush=True)
                    ihs_results["scikit-allel"] = ihs_scikit_allel(start, end, POP2)
                    print(f" {fmt_time(ihs_results['scikit-allel']['time_s'])}")

                if TOOLS["selscan"]:
                    print("  selscan...", end="", flush=True)
                    ihs_results["selscan"] = ihs_selscan(start, end, POP2, tmpdir)
                    print(f" {fmt_time(ihs_results['selscan']['time_s'])}")

                print("  GraphPop...", end="", flush=True)
                ihs_results["GraphPop"] = ihs_graphpop(start, end, POP2)
                print(f" {fmt_time(ihs_results['GraphPop']['time_s'])}")

                rows = []
                for tool, r in ihs_results.items():
                    rows.append([
                        tool, fmt_val(r.get("n_variants"), 0),
                        fmt_val(r.get("ihs_mean"), 4), fmt_val(r.get("ihs_std"), 4),
                        fmt_time(r.get("time_s")), fmt_mem(r.get("peak_mem_mb")),
                    ])
                print_group_table(
                    f"iHS — {region_name}",
                    ["Tool", "Variants", "Mean iHS", "Std iHS", "Time", "Peak Mem"],
                    rows,
                )

                print(f"\n--- Group E: XP-EHH — {POP2} vs {POP1} ---")
                xp_results = {}

                if TOOLS["scikit-allel"]:
                    print("  scikit-allel...", end="", flush=True)
                    xp_results["scikit-allel"] = xpehh_scikit_allel(start, end, POP2, POP1)
                    print(f" {fmt_time(xp_results['scikit-allel']['time_s'])}")

                if TOOLS["selscan"]:
                    print("  selscan...", end="", flush=True)
                    xp_results["selscan"] = xpehh_selscan(start, end, POP2, POP1, tmpdir)
                    print(f" {fmt_time(xp_results['selscan']['time_s'])}")

                print("  GraphPop...", end="", flush=True)
                xp_results["GraphPop"] = xpehh_graphpop(start, end, POP2, POP1)
                print(f" {fmt_time(xp_results['GraphPop']['time_s'])}")

                rows = []
                for tool, r in xp_results.items():
                    rows.append([
                        tool, fmt_val(r.get("n_variants"), 0),
                        fmt_val(r.get("xpehh_mean"), 4), fmt_val(r.get("xpehh_std"), 4),
                        fmt_time(r.get("time_s")), fmt_mem(r.get("peak_mem_mb")),
                    ])
                print_group_table(
                    f"XP-EHH — {region_name}",
                    ["Tool", "Variants", "Mean XP-EHH", "Std XP-EHH", "Time", "Peak Mem"],
                    rows,
                )

            # Collect all results for JSONL output
            region_record = {
                "region": region_name, "chr": CHR, "start": start, "end": end,
                "diversity": {k: v for k, v in div_results.items()},
                "differentiation": {k: v for k, v in fst_results.items()},
                "sfs": {k: v for k, v in sfs_results.items()},
            }
            if region_name != "large":
                region_record["ld"] = {k: v for k, v in ld_results.items()}
                region_record["ihs"] = {k: v for k, v in ihs_results.items()}
                region_record["xpehh"] = {k: v for k, v in xp_results.items()}
            all_results.append(region_record)

    # Write JSONL
    outfile = args.output or str(RESULTS_DIR / "vs_tools.jsonl")
    with open(outfile, "w") as f:
        for r in all_results:
            f.write(json.dumps(r, default=str) + "\n")
    print(f"\nResults written to {outfile}")

    # Final summary table
    print(f"\n{'='*80}")
    print("  FINAL SUMMARY — Wall-clock time (seconds)")
    print(f"{'='*80}")
    for record in all_results:
        rname = record["region"]
        print(f"\n  Region: {rname}")
        print(f"  {'Statistic':<20s}", end="")
        tool_names = set()
        for group in ["diversity", "differentiation", "sfs", "ld", "ihs", "xpehh"]:
            if group in record:
                tool_names.update(record[group].keys())
        tool_names = sorted(tool_names)
        for t in tool_names:
            print(f"  {t:>14s}", end="")
        print()
        print(f"  {'-'*20}", end="")
        for _ in tool_names:
            print(f"  {'-'*14}", end="")
        print()

        stat_groups = [
            ("π/θ_W/D", "diversity"),
            ("Fst/Dxy", "differentiation"),
            ("SFS", "sfs"),
            ("LD (r²)", "ld"),
            ("iHS", "ihs"),
            ("XP-EHH", "xpehh"),
        ]
        for label, group in stat_groups:
            if group not in record:
                continue
            print(f"  {label:<20s}", end="")
            for t in tool_names:
                if t in record[group]:
                    v = record[group][t].get("time_s")
                    print(f"  {fmt_time(v):>14s}", end="")
                else:
                    print(f"  {'—':>14s}", end="")
            print()


if __name__ == "__main__":
    main()
