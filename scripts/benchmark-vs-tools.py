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
  D. LD scenarios:
     D1. Regional Pairwise (100kb, r²≥0.1, MAF≥0.05) — scikit-allel, vcftools, PLINK2, GraphPop
     D2. LD Decay (500kb, r²≥0, MAF≥0.05)            — scikit-allel, vcftools, PLINK2, GraphPop
     D3. LD Pruning (1000kb/step50, r²=0.8)           — PLINK2 only
  E. Selection (iHS, XP-EHH)            — scikit-allel, selscan, GraphPop

Usage:
    conda run -n graphevo python scripts/benchmark-vs-tools.py [--region small|medium|large]
    conda run -n graphevo python scripts/benchmark-vs-tools.py --ld-scenarios

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
PGEN_PREFIX = str(ROOT / "data/raw/1000g/chr22_zst")  # .pgen/.pvar.zst/.psam
PANEL = str(ROOT / "data/raw/1000g/integrated_call_samples_v3.20130502.ALL.panel")
RESULTS_DIR = ROOT / "benchmarks" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CHR = "chr22"
POP1 = "YRI"   # sub-pop (Yoruba), used for diversity / single-pop stats
POP2 = "CEU"   # sub-pop (CEPH European), used for selection / pairwise stats
HAP_CACHE_DIR = str(ROOT / "data/hap_cache")

NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"

REGIONS = {
    "small":  (16_000_000, 16_100_000),
    "medium": (16_000_000, 16_500_000),
    "large":  (16_000_000, 17_000_000),
    "full":   (10_510_000, 50_810_000),   # full chr22 (~1.07M variants)
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
TOOLS["plink2-pgen"] = (TOOLS["plink2"]
                         and os.path.exists(PGEN_PREFIX + ".pgen")
                         and os.path.exists(PGEN_PREFIX + ".pvar.zst"))


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
            if len(parts) >= 2:
                sample_pop[parts[0]] = parts[1]  # sub-pop (ACB, ASW, CEU, YRI, …)
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
        _, time_pi, mem_pi = measure_subprocess(cmd_pi, timeout=3600)
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
        _, time_td, mem_td = measure_subprocess(cmd_td, timeout=3600)
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
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
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


def _write_plink2_pheno(tmpdir, pop1, pop2):
    """Write #IID/POP pheno file for PLINK2 (cached per tmpdir)."""
    path = os.path.join(tmpdir, "plink2_pop_pheno.txt")
    if os.path.exists(path):
        return path
    sample_pop = load_panel()
    with open(path, "w") as f:
        f.write("#IID\tPOP\n")
        for s, p in sample_pop.items():
            if p in (pop1, pop2):
                f.write(f"{s}\t{p}\n")
    return path


def _write_plink2_keep(tmpdir, pop):
    """Write #IID keep file for PLINK2 (cached per tmpdir)."""
    path = os.path.join(tmpdir, f"plink2_keep_{pop}.txt")
    if os.path.exists(path):
        return path
    sample_pop = load_panel()
    with open(path, "w") as f:
        f.write("#IID\n")
        for s, p in sample_pop.items():
            if p == pop:
                f.write(f"{s}\n")
    return path


def fst_plink2_pgen(start, end, pop1, pop2, tmpdir):
    """Compute Hudson's Fst using PLINK2 from native .pgen files."""
    pheno_file = _write_plink2_pheno(tmpdir, pop1, pop2)
    cmd = [
        "plink2",
        "--pfile", PGEN_PREFIX, "vzs",
        "--chr", CHR,
        "--from-bp", str(start), "--to-bp", str(end),
        "--pheno", pheno_file,
        "--fst", "POP", "method=hudson",
        "--allow-extra-chr",
        "--out", os.path.join(tmpdir, "plink2_pgen_fst"),
    ]
    try:
        _, elapsed, peak_mb = measure_subprocess(cmd, timeout=120)
    except RuntimeError as e:
        return {"error": str(e), "time_s": 0, "peak_mem_mb": 0}

    fst_val = None
    summary = os.path.join(tmpdir, "plink2_pgen_fst.fst.summary")
    if os.path.exists(summary):
        with open(summary) as f:
            header = f.readline()
            data_line = f.readline()
            if data_line:
                parts = data_line.strip().split("\t")
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
# Group D: LD scenarios
#   D1 — Regional Pairwise LD (100kb window, r²≥0.1, MAF≥0.05)
#   D2 — LD Decay            (500kb window, r²≥0,   MAF≥0.05)
#   D3 — LD Pruning          (1000kb/step50, r²=0.8, MAF≥0.05) — PLINK2 only
# ---------------------------------------------------------------------------

def _count_file_lines(path, skip_header=True):
    """Count data lines in a file (minus header if requested)."""
    n = 0
    with open(path) as f:
        for _ in f:
            n += 1
    return max(0, n - 1) if skip_header else n


def ld_scikit_allel_windowed(start, end, pop, max_dist_bp, r2_threshold, min_maf):
    """Compute windowed pairwise r² using scikit-allel (Rogers-Huff estimator).

    Only feasible for ≤1Mb regions (~3K variants after MAF filter).
    """
    import allel

    def _compute():
        data = load_genotype_data(start, end, pop)
        positions = data["positions"]
        genotypes = data["genotypes"].astype(float)
        n_samples = data["n_samples"]
        n_total = 2 * n_samples

        # MAF filter
        ac = genotypes.sum(axis=1)
        af = ac / n_total
        maf = np.minimum(af, 1.0 - af)
        mask = maf >= min_maf
        positions = positions[mask]
        genotypes = genotypes[mask]
        n_variants = len(positions)

        if n_variants < 2:
            return {"n_pairs": 0, "n_variants": n_variants}

        r = allel.rogers_huff_r(genotypes)
        from scipy.spatial.distance import squareform
        r2_mat = squareform(r) ** 2

        # Apply distance window and r² threshold
        n_pairs = 0
        for i in range(n_variants):
            for j in range(i + 1, n_variants):
                dist = positions[j] - positions[i]
                if dist > max_dist_bp:
                    break
                if r2_mat[i, j] >= r2_threshold:
                    n_pairs += 1

        return {"n_pairs": n_pairs, "n_variants": n_variants}

    result, elapsed, peak_mb = measure_python_call(_compute)
    return {**result, "time_s": elapsed, "peak_mem_mb": peak_mb}


def ld_vcftools(start, end, pop, tmpdir, max_dist_bp, r2_threshold, min_maf):
    """Compute windowed LD using vcftools --geno-r2.

    vcftools uses an EM algorithm for unphased genotype r².
    Output: CHR POS1 POS2 N_INDV R^2
    """
    keep_file = write_pop_keep_file(tmpdir, pop)
    tag = f"vcftools_ld_{start}_{end}_{max_dist_bp}"
    outprefix = os.path.join(tmpdir, tag)

    cmd = [
        "vcftools", "--gzvcf", VCF,
        "--chr", CHR, "--from-bp", str(start), "--to-bp", str(end),
        "--keep", keep_file,
        "--geno-r2",
        "--ld-window-bp", str(max_dist_bp),
        "--min-r2", str(r2_threshold),
        "--maf", str(min_maf),
        "--out", outprefix,
    ]
    try:
        _, elapsed, peak_mb = measure_subprocess(cmd, timeout=3600)
    except RuntimeError as e:
        return {"error": str(e), "time_s": 0, "peak_mem_mb": 0}

    # Count pairs from .geno.ld output
    ld_file = outprefix + ".geno.ld"
    n_pairs = 0
    if os.path.exists(ld_file):
        n_pairs = _count_file_lines(ld_file, skip_header=True)

    return {
        "n_pairs": n_pairs,
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


def ld_plink2_windowed(start, end, pop, tmpdir, max_dist_kb, r2_threshold, min_maf,
                        use_pgen=False, tag_suffix=""):
    """Compute windowed LD using PLINK2 --r2-unphased.

    Parameterized for VCF vs pgen input via use_pgen flag.
    """
    keep_file = _write_plink2_keep(tmpdir, pop)
    tag = f"plink2_ld{'_pgen' if use_pgen else ''}_{start}_{end}{tag_suffix}"
    outprefix = os.path.join(tmpdir, tag)

    if use_pgen:
        input_args = ["--pfile", PGEN_PREFIX, "vzs"]
    else:
        input_args = ["--vcf", VCF]

    cmd = [
        "plink2",
        *input_args,
        "--chr", CHR,
        "--from-bp", str(start), "--to-bp", str(end),
        "--keep", keep_file,
        "--r2-unphased",
        "--ld-window-kb", str(max_dist_kb),
        "--ld-window-r2", str(r2_threshold),
        "--maf", str(min_maf),
        "--allow-extra-chr",
        "--out", outprefix,
    ]
    try:
        _, elapsed, peak_mb = measure_subprocess(cmd, timeout=3600)
    except RuntimeError as e:
        return {"error": str(e), "time_s": 0, "peak_mem_mb": 0}

    # Count pairs from output (.vcor or .ld)
    n_pairs = 0
    for ext in [".vcor", ".ld"]:
        ld_file = outprefix + ext
        if os.path.exists(ld_file):
            n_pairs = _count_file_lines(ld_file, skip_header=True)
            break

    return {
        "n_pairs": n_pairs,
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


def ld_graphpop_windowed(start, end, pop, max_dist, r2_threshold, min_af):
    """Run graphpop.ld with MAF filter and write_edges=false (read-only benchmark)."""
    opts = f"{{min_af: {min_af}, write_edges: false}}"
    query = (
        f"CALL graphpop.ld('{CHR}', {start}, {end}, '{pop}', {max_dist}, "
        f"{r2_threshold}, {opts}) "
        "YIELD variant1 RETURN count(*) AS n_pairs;"
    )
    out, elapsed, peak_mb = run_cypher(query, timeout=3600)
    result = parse_cypher_row(out)

    return {
        "n_pairs": result.get("n_pairs"),
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


def ld_graphpop_numpy(start, end, pop, max_dist_bp, r2_threshold, min_maf):
    """Compute windowed pairwise r² from the numpy haplotype cache (bypasses Neo4j).

    Uses phased haplotypes: r²(i,j) = cov(h_i, h_j)² / (var_i × var_j)
    where h_i, h_j are 0/1 vectors across 2*n_pop haplotypes.

    For M ≤ 5000 variants: computes full M×M matrix at once.
    For M > 5000: streams in row-blocks of 512 to bound peak memory.
    """
    def _compute():
        pos, hap = load_hap_cache(pop, start, end)
        if len(pos) < 2:
            return {"n_pairs": 0, "n_variants": 0}

        hap_f = hap.astype(np.float32)   # (M, 2K)
        n_hap = float(hap_f.shape[1])

        # MAF filter
        freq = hap_f.mean(axis=1)        # (M,)
        maf  = np.minimum(freq, 1.0 - freq)
        mask = maf >= min_maf
        hap_f = hap_f[mask]
        pos   = pos[mask]
        freq  = freq[mask]
        M = len(pos)
        if M < 2:
            return {"n_pairs": 0, "n_variants": M}

        # Clamp to prevent divide-by-zero at fixed sites
        var = (freq * (1.0 - freq)).clip(min=1e-9).astype(np.float32)  # (M,)

        FULL_THRESH = 5000
        n_pairs = 0

        if M <= FULL_THRESH:
            # Full matrix in one shot
            cov = hap_f @ hap_f.T / n_hap - np.outer(freq, freq)  # (M, M)
            r2  = (cov ** 2) / np.outer(var, var)
            i_idx, j_idx = np.triu_indices(M, k=1)
            in_window = (pos[j_idx] - pos[i_idx]) <= max_dist_bp
            n_pairs = int((r2[i_idx[in_window], j_idx[in_window]] >= r2_threshold).sum())
        else:
            # Streaming block approach: process BLOCK rows at a time
            BLOCK = 512
            # For each variant i, precompute the last j within max_dist_bp
            j_max = np.searchsorted(pos, pos + max_dist_bp, side="right")  # (M,)

            for i0 in range(0, M, BLOCK):
                i1    = min(i0 + BLOCK, M)
                # j_end is the max col index reachable from any row in this block
                j_end = int(j_max[i1 - 1])
                b     = i1 - i0
                J     = j_end - i0

                hi = hap_f[i0:i1]       # (b, 2K)
                hj = hap_f[i0:j_end]    # (J, 2K) — aligned to global index i0

                cov   = hi @ hj.T / n_hap - np.outer(freq[i0:i1], freq[i0:j_end])
                denom = np.outer(var[i0:i1], var[i0:j_end])
                r2    = cov ** 2 / denom  # (b, J)

                # Validity mask: strict upper triangle and within distance window
                gi    = np.arange(i0, i1)[:, None]     # (b, 1) global row
                gj    = np.arange(i0, j_end)[None, :]  # (1, J) global col
                valid = (gj > gi) & (
                    (pos[i0:j_end][None, :] - pos[i0:i1, None]) <= max_dist_bp)
                n_pairs += int((r2[valid] >= r2_threshold).sum())

        return {"n_pairs": n_pairs, "n_variants": M}

    result, elapsed, peak_mb = measure_python_call(_compute)
    return {**result, "time_s": elapsed, "peak_mem_mb": peak_mb}


def ld_plink2_pruning(start, end, pop, tmpdir, window_kb, step, r2, min_maf,
                       use_pgen=False):
    """Run PLINK2 --indep-pairwise for LD pruning.

    Reports retained variant count and time.
    """
    keep_file = _write_plink2_keep(tmpdir, pop)
    tag = f"plink2_prune{'_pgen' if use_pgen else ''}_{start}_{end}"
    outprefix = os.path.join(tmpdir, tag)

    if use_pgen:
        input_args = ["--pfile", PGEN_PREFIX, "vzs"]
    else:
        input_args = ["--vcf", VCF]

    cmd = [
        "plink2",
        *input_args,
        "--chr", CHR,
        "--from-bp", str(start), "--to-bp", str(end),
        "--keep", keep_file,
        "--indep-pairwise", f"{window_kb}kb", "1", str(r2),
        "--maf", str(min_maf),
        "--allow-extra-chr",
        "--out", outprefix,
    ]
    try:
        _, elapsed, peak_mb = measure_subprocess(cmd, timeout=3600)
    except RuntimeError as e:
        return {"error": str(e), "time_s": 0, "peak_mem_mb": 0}

    # Count retained variants from .prune.in
    prune_file = outprefix + ".prune.in"
    n_retained = 0
    if os.path.exists(prune_file):
        n_retained = _count_file_lines(prune_file, skip_header=False)

    return {
        "n_retained": n_retained,
        "time_s": elapsed,
        "peak_mem_mb": peak_mb,
    }


# ---------------------------------------------------------------------------
# Haplotype cache helper (for numpy bypasses in Group E)
# ---------------------------------------------------------------------------

def load_hap_cache(pop, start=None, end=None, chrom=None):
    """Load population haplotypes from the numpy haplotype cache.

    Returns (pos, hap) where:
      pos: int64 array (M,) — variant positions
      hap: int8 array (M, 2*n_pop) — haplotype matrix for the population
    """
    import json
    chrom = chrom or CHR
    meta_path = os.path.join(HAP_CACHE_DIR, "sample_meta.json")
    npz_path  = os.path.join(HAP_CACHE_DIR, f"{chrom}.npz")
    if not os.path.exists(meta_path) or not os.path.exists(npz_path):
        raise FileNotFoundError(
            f"Hap cache missing. Run: python scripts/export_hap_cache.py --chr {chrom}")

    with open(meta_path) as f:
        meta = json.load(f)
    pop_ids   = meta["pop_ids"]
    n_samples = meta["n_samples"]

    sample_idx = np.array([i for i, p in enumerate(pop_ids) if p == pop], dtype=np.int32)
    if len(sample_idx) == 0:
        raise ValueError(f"Population '{pop}' not in hap cache. "
                         f"Available: {sorted(set(pop_ids))[:10]}")

    d          = np.load(npz_path, mmap_mode="r")
    pos        = d["pos"]
    hap_packed = d["hap"]

    if start is not None and end is not None:
        mask       = (pos >= start) & (pos <= end)
        pos        = np.array(pos[mask])
        hap_packed = np.array(hap_packed[mask])
    else:
        pos        = np.array(pos)
        hap_packed = np.array(hap_packed)

    if len(pos) == 0:
        return pos, np.empty((0, 2 * len(sample_idx)), dtype=np.int8)

    hap_full = np.unpackbits(hap_packed, axis=1, bitorder="big")[:, :2 * n_samples]

    hap_idx = np.empty(2 * len(sample_idx), dtype=np.int32)
    hap_idx[0::2] = 2 * sample_idx
    hap_idx[1::2] = 2 * sample_idx + 1

    return pos, hap_full[:, hap_idx].astype(np.int8)


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


def ihs_graphpop_numpy(start, end, pop):
    """Compute iHS using scikit-allel on the numpy haplotype cache (bypasses Neo4j)."""
    import allel

    def _compute():
        pos, hap = load_hap_cache(pop, start, end)
        if len(pos) == 0:
            return {"n_variants": 0, "ihs_mean": None}

        h  = allel.HaplotypeArray(hap)
        ac = h.count_alleles()
        af = ac.to_frequencies()[:, 1]
        flt = ac.is_segregating() & (af >= 0.05) & (af <= 0.95)

        h_flt   = h[flt]
        pos_flt = pos[flt]
        if len(pos_flt) < 10:
            return {"n_variants": 0, "ihs_mean": None}

        ihs_raw = allel.ihs(h_flt, pos_flt, include_edges=True)
        try:
            ihs_std, _ = allel.standardize_by_allele_count(
                ihs_raw, ac[flt][:, 1], diagnostics=False)
        except (ValueError, ZeroDivisionError):
            valid_raw = ihs_raw[~np.isnan(ihs_raw)]
            ihs_std = (ihs_raw - np.nanmean(ihs_raw)) / np.nanstd(ihs_raw) \
                      if len(valid_raw) > 1 else ihs_raw

        valid = ~np.isnan(ihs_std)
        return {
            "n_variants": int(valid.sum()),
            "ihs_mean": float(np.nanmean(ihs_std)),
            "ihs_std": float(np.nanstd(ihs_std)),
        }

    result, elapsed, peak_mb = measure_python_call(_compute)
    return {**result, "time_s": elapsed, "peak_mem_mb": peak_mb}


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


def xpehh_graphpop_numpy(start, end, pop1, pop2):
    """Compute XP-EHH using scikit-allel on the numpy haplotype cache (bypasses Neo4j)."""
    import allel

    def _compute():
        pos1, hap1 = load_hap_cache(pop1, start, end)
        pos2, hap2 = load_hap_cache(pop2, start, end)

        # Use common positions (identical since same cache, but guard against edge cases)
        pos = pos1
        if len(pos) < 10:
            return {"n_variants": 0}

        xp = allel.xpehh(
            allel.HaplotypeArray(hap1),
            allel.HaplotypeArray(hap2),
            pos, include_edges=True,
        )
        valid = ~np.isnan(xp)
        return {
            "n_variants": int(valid.sum()),
            "xpehh_mean": float(np.nanmean(xp)),
            "xpehh_std": float(np.nanstd(xp)),
        }

    result, elapsed, peak_mb = measure_python_call(_compute)
    return {**result, "time_s": elapsed, "peak_mem_mb": peak_mb}


def sel_stats_numpy_joint(start, end, pop1, pop2, chrom=None):
    """
    Load haplotype cache ONCE for pop1+pop2, then compute nSL, iHS, XP-EHH
    sharing the same unpacked arrays.
    Returns sequential (shared load) timing and parallel timing.
    """
    import allel
    from concurrent.futures import ThreadPoolExecutor as _TPE

    tracemalloc.start()

    # ── Phase 1: Load once ────────────────────────────────────────────────────
    t_load0 = time.perf_counter()
    pos1, hap1 = load_hap_cache(pop1, start, end, chrom=chrom or CHR)
    pos2, hap2 = load_hap_cache(pop2, start, end, chrom=chrom or CHR)
    load_time = time.perf_counter() - t_load0

    if len(pos1) == 0:
        tracemalloc.stop()
        return {"load_time_s": load_time, "nsl": None, "ihs": None, "xpehh": None}

    # ── Phase 2: Shared filter ────────────────────────────────────────────────
    h1 = allel.HaplotypeArray(hap1)
    h2 = allel.HaplotypeArray(hap2)
    ac1 = h1.count_alleles()
    af1 = ac1.to_frequencies()[:, 1]
    flt = ac1.is_segregating() & (af1 >= 0.05) & (af1 <= 0.95)

    h1f = h1[flt]; h2f = h2[flt]; pos_flt = pos1[flt]; ac1f = ac1[flt]

    # ── Phase 3: Compute functions ────────────────────────────────────────────
    def _nsl():
        t0 = time.perf_counter()
        raw = allel.nsl(h1f)
        try:
            std, _ = allel.standardize_by_allele_count(raw, ac1f[:, 1], diagnostics=False)
        except Exception:
            std = (raw - np.nanmean(raw)) / np.nanstd(raw)
        valid = ~np.isnan(std)
        return {"n_variants": int(valid.sum()),
                "nsl_mean": float(np.nanmean(std)),
                "time_s": time.perf_counter() - t0}

    def _ihs():
        t0 = time.perf_counter()
        raw = allel.ihs(h1f, pos_flt, include_edges=True)
        try:
            std, _ = allel.standardize_by_allele_count(raw, ac1f[:, 1], diagnostics=False)
        except Exception:
            std = (raw - np.nanmean(raw)) / np.nanstd(raw)
        valid = ~np.isnan(std)
        return {"n_variants": int(valid.sum()),
                "ihs_mean": float(np.nanmean(std)),
                "time_s": time.perf_counter() - t0}

    def _xpehh():
        t0 = time.perf_counter()
        xp = allel.xpehh(h1f, h2f, pos_flt, include_edges=True)
        valid = ~np.isnan(xp)
        return {"n_variants": int(valid.sum()),
                "xpehh_mean": float(np.nanmean(xp)),
                "time_s": time.perf_counter() - t0}

    # Sequential (shared load, sequential compute)
    nsl_r  = _nsl()
    ihs_r  = _ihs()
    xp_r   = _xpehh()
    seq_compute = nsl_r["time_s"] + ihs_r["time_s"] + xp_r["time_s"]
    seq_total   = load_time + seq_compute

    # Parallel (shared load, parallel compute with 3 threads)
    t_par0 = time.perf_counter()
    with _TPE(max_workers=3) as pool:
        f_nsl  = pool.submit(_nsl)
        f_ihs  = pool.submit(_ihs)
        f_xp   = pool.submit(_xpehh)
        nsl_rp  = f_nsl.result()
        ihs_rp  = f_ihs.result()
        xp_rp   = f_xp.result()
    par_compute = time.perf_counter() - t_par0
    par_total   = load_time + par_compute

    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return {
        "load_time_s":   load_time,
        "nsl":           nsl_r,
        "ihs":           ihs_r,
        "xpehh":         xp_r,
        "seq_compute_s": seq_compute,
        "seq_total_s":   seq_total,
        "par_compute_s": par_compute,
        "par_total_s":   par_total,
        "peak_mem_mb":   peak_bytes / (1024 * 1024),
        "pop1": pop1, "pop2": pop2,
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
    parser.add_argument("--ld-scenarios", action="store_true",
                        help="Run only LD scenario benchmarks (D1/D2/D3)")
    args = parser.parse_args()

    if args.ld_scenarios:
        regions_to_run = ["large", "full"]
    elif args.all_regions:
        regions_to_run = list(REGIONS.keys())
    else:
        regions_to_run = [args.region]

    print("=" * 80)
    print("  GraphPop Head-to-Head Benchmark")
    if args.ld_scenarios:
        print("  (LD scenarios only)")
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

            div_results = {}
            fst_results = {}
            sfs_results = {}
            ihs_results = {}
            xp_results = {}
            joint_result = None
            ld_only = args.ld_scenarios

            if not ld_only:
                # ---- Group A: Diversity ----
                print(f"\n--- Group A: Diversity (π, θ_W, Tajima's D) — pop={POP1} ---")

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
                    print("  PLINK2 (VCF)...", end="", flush=True)
                    fst_results["PLINK2-VCF"] = fst_plink2(start, end, POP1, POP2, tmpdir)
                    print(f" {fmt_time(fst_results['PLINK2-VCF']['time_s'])}")

                if TOOLS["plink2-pgen"]:
                    print("  PLINK2 (.pgen)...", end="", flush=True)
                    fst_results["PLINK2-pgen"] = fst_plink2_pgen(start, end, POP1, POP2, tmpdir)
                    print(f" {fmt_time(fst_results['PLINK2-pgen']['time_s'])}")

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

            # ---- Group D: LD Scenarios ----
            ld_scenarios = {}
            _hap_cache_ok = (
                os.path.exists(os.path.join(HAP_CACHE_DIR, "sample_meta.json")) and
                os.path.exists(os.path.join(HAP_CACHE_DIR, f"{CHR}.npz"))
            )

            # D1: Regional Pairwise LD (100kb window, r²≥0.1, MAF≥0.05)
            # Run on large (1Mb) with all tools; on full (chr22) skip scikit-allel
            if region_name in ("large", "full"):
                print(f"\n--- D1: Regional Pairwise LD (100kb, r²≥0.1, MAF≥0.05) — pop={POP2} ---")
                d1 = {}

                if TOOLS["scikit-allel"] and region_name == "large":
                    print("  scikit-allel...", end="", flush=True)
                    d1["scikit-allel"] = ld_scikit_allel_windowed(
                        start, end, POP2, max_dist_bp=100_000, r2_threshold=0.1, min_maf=0.05)
                    print(f" {fmt_time(d1['scikit-allel']['time_s'])}")

                # Skip vcftools on full chr22 — too slow (100M+ output pairs)
                if TOOLS["vcftools"] and region_name != "full":
                    print("  vcftools...", end="", flush=True)
                    d1["vcftools"] = ld_vcftools(
                        start, end, POP2, tmpdir, max_dist_bp=100_000,
                        r2_threshold=0.1, min_maf=0.05)
                    print(f" {fmt_time(d1['vcftools']['time_s'])}")

                if TOOLS["plink2"]:
                    print("  PLINK2 (VCF)...", end="", flush=True)
                    d1["PLINK2-VCF"] = ld_plink2_windowed(
                        start, end, POP2, tmpdir, max_dist_kb=100,
                        r2_threshold=0.1, min_maf=0.05, use_pgen=False,
                        tag_suffix="_d1")
                    print(f" {fmt_time(d1['PLINK2-VCF']['time_s'])}")

                if TOOLS["plink2-pgen"]:
                    print("  PLINK2 (.pgen)...", end="", flush=True)
                    d1["PLINK2-pgen"] = ld_plink2_windowed(
                        start, end, POP2, tmpdir, max_dist_kb=100,
                        r2_threshold=0.1, min_maf=0.05, use_pgen=True,
                        tag_suffix="_d1")
                    print(f" {fmt_time(d1['PLINK2-pgen']['time_s'])}")

                print("  GraphPop (Neo4j)...", end="", flush=True)
                d1["GraphPop"] = ld_graphpop_windowed(
                    start, end, POP2, max_dist=100_000, r2_threshold=0.1, min_af=0.05)
                print(f" {fmt_time(d1['GraphPop']['time_s'])}")

                if _hap_cache_ok:
                    print("  GraphPop (numpy)...", end="", flush=True)
                    d1["GraphPop-numpy"] = ld_graphpop_numpy(
                        start, end, POP2, max_dist_bp=100_000, r2_threshold=0.1, min_maf=0.05)
                    print(f" {fmt_time(d1['GraphPop-numpy']['time_s'])}")
                else:
                    print("  GraphPop (numpy)... SKIPPED (hap cache not ready)")

                rows = []
                for tool, r in d1.items():
                    rows.append([
                        tool, fmt_val(r.get("n_pairs"), 0),
                        fmt_time(r.get("time_s")), fmt_mem(r.get("peak_mem_mb")),
                    ])
                print_group_table(
                    f"D1: Regional Pairwise LD — {region_name}",
                    ["Tool", "Pairs (r²≥0.1)", "Time", "Peak Mem"],
                    rows,
                )
                ld_scenarios["d1"] = d1

            # D2: LD Decay (500kb window, r²≥0 = all pairs, MAF≥0.05)
            # Only feasible on large (1Mb) — full chr22 produces too many pairs
            if region_name == "large":
                print(f"\n--- D2: LD Decay (500kb, r²≥0, MAF≥0.05) — pop={POP2} ---")
                d2 = {}

                if TOOLS["scikit-allel"]:
                    print("  scikit-allel...", end="", flush=True)
                    d2["scikit-allel"] = ld_scikit_allel_windowed(
                        start, end, POP2, max_dist_bp=500_000, r2_threshold=0.0, min_maf=0.05)
                    print(f" {fmt_time(d2['scikit-allel']['time_s'])}")

                if TOOLS["vcftools"]:
                    print("  vcftools...", end="", flush=True)
                    d2["vcftools"] = ld_vcftools(
                        start, end, POP2, tmpdir, max_dist_bp=500_000,
                        r2_threshold=0.0, min_maf=0.05)
                    print(f" {fmt_time(d2['vcftools']['time_s'])}")

                if TOOLS["plink2"]:
                    print("  PLINK2 (VCF)...", end="", flush=True)
                    d2["PLINK2-VCF"] = ld_plink2_windowed(
                        start, end, POP2, tmpdir, max_dist_kb=500,
                        r2_threshold=0.0, min_maf=0.05, use_pgen=False,
                        tag_suffix="_d2")
                    print(f" {fmt_time(d2['PLINK2-VCF']['time_s'])}")

                if TOOLS["plink2-pgen"]:
                    print("  PLINK2 (.pgen)...", end="", flush=True)
                    d2["PLINK2-pgen"] = ld_plink2_windowed(
                        start, end, POP2, tmpdir, max_dist_kb=500,
                        r2_threshold=0.0, min_maf=0.05, use_pgen=True,
                        tag_suffix="_d2")
                    print(f" {fmt_time(d2['PLINK2-pgen']['time_s'])}")

                print("  GraphPop (Neo4j)...", end="", flush=True)
                d2["GraphPop"] = ld_graphpop_windowed(
                    start, end, POP2, max_dist=500_000, r2_threshold=0.0, min_af=0.05)
                print(f" {fmt_time(d2['GraphPop']['time_s'])}")

                if _hap_cache_ok:
                    print("  GraphPop (numpy)...", end="", flush=True)
                    d2["GraphPop-numpy"] = ld_graphpop_numpy(
                        start, end, POP2, max_dist_bp=500_000, r2_threshold=0.0, min_maf=0.05)
                    print(f" {fmt_time(d2['GraphPop-numpy']['time_s'])}")
                else:
                    print("  GraphPop (numpy)... SKIPPED (hap cache not ready)")

                rows = []
                for tool, r in d2.items():
                    rows.append([
                        tool, fmt_val(r.get("n_pairs"), 0),
                        fmt_time(r.get("time_s")), fmt_mem(r.get("peak_mem_mb")),
                    ])
                print_group_table(
                    f"D2: LD Decay — {region_name}",
                    ["Tool", "Pairs (all)", "Time", "Peak Mem"],
                    rows,
                )
                ld_scenarios["d2"] = d2

            # D3: LD Pruning (1000kb window, step 50, r²=0.8, MAF≥0.05)
            # PLINK2 only (GraphPop has no pruning proc). Full chr22 only.
            if region_name == "full" and TOOLS["plink2"]:
                print(f"\n--- D3: LD Pruning (1000kb/step50, r²=0.8, MAF≥0.05) — pop={POP2} ---")
                d3 = {}

                print("  PLINK2 (VCF)...", end="", flush=True)
                d3["PLINK2-VCF"] = ld_plink2_pruning(
                    start, end, POP2, tmpdir, window_kb=1000, step=50,
                    r2=0.8, min_maf=0.05, use_pgen=False)
                print(f" {fmt_time(d3['PLINK2-VCF']['time_s'])}")

                if TOOLS["plink2-pgen"]:
                    print("  PLINK2 (.pgen)...", end="", flush=True)
                    d3["PLINK2-pgen"] = ld_plink2_pruning(
                        start, end, POP2, tmpdir, window_kb=1000, step=50,
                        r2=0.8, min_maf=0.05, use_pgen=True)
                    print(f" {fmt_time(d3['PLINK2-pgen']['time_s'])}")

                rows = []
                for tool, r in d3.items():
                    rows.append([
                        tool, fmt_val(r.get("n_retained"), 0),
                        fmt_time(r.get("time_s")), fmt_mem(r.get("peak_mem_mb")),
                    ])
                print_group_table(
                    f"D3: LD Pruning — {region_name}",
                    ["Tool", "Retained Variants", "Time", "Peak Mem"],
                    rows,
                )
                ld_scenarios["d3"] = d3

            # ---- Group E: Selection ----
            # small/medium for quick smoke tests; full for paper benchmark
            if region_name in ("small", "medium", "full") and not ld_only:
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

                print("  GraphPop (Neo4j)...", end="", flush=True)
                ihs_results["GraphPop"] = ihs_graphpop(start, end, POP2)
                print(f" {fmt_time(ihs_results['GraphPop']['time_s'])}")

                if _hap_cache_ok:
                    print("  GraphPop (numpy)...", end="", flush=True)
                    ihs_results["GraphPop-numpy"] = ihs_graphpop_numpy(start, end, POP2)
                    print(f" {fmt_time(ihs_results['GraphPop-numpy']['time_s'])}")
                else:
                    print("  GraphPop (numpy)... SKIPPED (hap cache not ready)")

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

                print("  GraphPop (Neo4j)...", end="", flush=True)
                xp_results["GraphPop"] = xpehh_graphpop(start, end, POP2, POP1)
                print(f" {fmt_time(xp_results['GraphPop']['time_s'])}")

                if _hap_cache_ok:
                    print("  GraphPop (numpy)...", end="", flush=True)
                    xp_results["GraphPop-numpy"] = xpehh_graphpop_numpy(start, end, POP2, POP1)
                    print(f" {fmt_time(xp_results['GraphPop-numpy']['time_s'])}")
                else:
                    print("  GraphPop (numpy)... SKIPPED (hap cache not ready)")

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

                # ---- Group F: Joint Selection Stats ----
                joint_result = None
                if _hap_cache_ok:
                    print(f"\n--- Group F: Joint nSL+iHS+XP-EHH — {POP2} vs {POP1} ---")
                    print("  Loading cache + computing jointly...", end="", flush=True)
                    joint_result = sel_stats_numpy_joint(start, end, POP2, POP1)
                    print(f"  load={joint_result['load_time_s']:.1f}s  "
                          f"seq={joint_result['seq_total_s']:.1f}s  "
                          f"par={joint_result['par_total_s']:.1f}s")

                    indep_ihs_t  = ihs_results.get("GraphPop-numpy", {}).get("time_s") or 0
                    indep_xp_t   = xp_results.get("GraphPop-numpy",  {}).get("time_s") or 0
                    indep_sum    = indep_ihs_t + indep_xp_t

                    def fmt(v):
                        return f"{v:.4f}" if v is not None else "—"

                    rows_f = [
                        ["Load (shared)",    "—", "—",
                         f"{joint_result['load_time_s']:.2f}s", "—"],
                        ["nSL (shared)",     "—", fmt(joint_result["nsl"]["nsl_mean"]),
                         f"{joint_result['nsl']['time_s']:.2f}s", "—"],
                        ["iHS (shared)",     "—", fmt(joint_result["ihs"]["ihs_mean"]),
                         f"{joint_result['ihs']['time_s']:.2f}s", "—"],
                        ["XP-EHH (shared)", "—", fmt(joint_result["xpehh"]["xpehh_mean"]),
                         f"{joint_result['xpehh']['time_s']:.2f}s", "—"],
                        ["TOTAL sequential", "—", "—",
                         f"{joint_result['seq_total_s']:.2f}s", "—"],
                        ["TOTAL parallel",   "—", "—",
                         f"{joint_result['par_total_s']:.2f}s",
                         f"{indep_sum/joint_result['par_total_s']:.1f}×"
                         if indep_sum > 0 and joint_result['par_total_s'] > 0 else "—"],
                    ]
                    print_group_table(
                        f"Group F: Joint nSL+iHS+XP-EHH (GraphPop-numpy) — {region_name}",
                        ["Stat", "Variants", "Mean", "Time", "Speedup vs indep"],
                        rows_f,
                    )

            # Collect all results for JSONL output
            region_record = {
                "region": region_name, "chr": CHR, "start": start, "end": end,
                "diversity": {k: v for k, v in div_results.items()},
                "differentiation": {k: v for k, v in fst_results.items()},
                "sfs": {k: v for k, v in sfs_results.items()},
            }
            if ld_scenarios:
                region_record["ld_scenarios"] = ld_scenarios
            if region_name in ("small", "medium", "full"):
                region_record["ihs"] = {k: v for k, v in ihs_results.items()}
                region_record["xpehh"] = {k: v for k, v in xp_results.items()}
                if joint_result is not None:
                    region_record["joint"] = joint_result
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
        for group in ["diversity", "differentiation", "sfs", "ihs", "xpehh"]:
            if group in record:
                tool_names.update(record[group].keys())
        if "ld_scenarios" in record:
            for scenario in record["ld_scenarios"].values():
                tool_names.update(scenario.keys())
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

        # LD scenarios
        if "ld_scenarios" in record:
            scenario_labels = {"d1": "D1: Pairwise LD", "d2": "D2: LD Decay", "d3": "D3: LD Pruning"}
            for key, label in scenario_labels.items():
                if key not in record["ld_scenarios"]:
                    continue
                scenario = record["ld_scenarios"][key]
                print(f"  {label:<20s}", end="")
                for t in tool_names:
                    if t in scenario:
                        v = scenario[t].get("time_s")
                        print(f"  {fmt_time(v):>14s}", end="")
                    else:
                        print(f"  {'—':>14s}", end="")
                print()


if __name__ == "__main__":
    main()
