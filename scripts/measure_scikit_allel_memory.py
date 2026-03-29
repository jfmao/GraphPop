#!/usr/bin/env python3
"""Measure scikit-allel peak memory for all statistics on full chr22 (CEU).

Reports per-statistic tracemalloc peak and overall worst-case.
Uses the same data loading + computation functions as benchmark-vs-tools.py.
"""

import json
import sys
import time
import tracemalloc
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
VCF = str(ROOT / "data" / "raw" / "1000g" / "CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz")
CHR = "chr22"
START = 10_510_000
END = 50_810_000
POP1 = "CEU"
POP2 = "YRI"


PANEL = str(ROOT / "data" / "raw" / "1000g" / "integrated_call_samples_v3.20130502.ALL.panel")


def load_panel():
    mapping = {}
    with open(PANEL) as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping


def measure(label, func):
    """Run func(), return (result, time_s, peak_mb)."""
    tracemalloc.start()
    t0 = time.perf_counter()
    result = func()
    elapsed = time.perf_counter() - t0
    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    peak_mb = peak_bytes / (1024 * 1024)
    print(f"  {label}: {elapsed:.1f}s, peak_mem={peak_mb:.1f} MB")
    return {"time_s": elapsed, "peak_mem_mb": peak_mb, **result}


def run_diversity():
    """Full chr22 diversity (pi, theta_w, Tajima's D) — CEU."""
    import allel
    import cyvcf2
    from scipy.special import digamma

    def _compute():
        sample_pop = load_panel()
        vcf = cyvcf2.VCF(VCF)
        all_samples = list(vcf.samples)
        pop_idx = [i for i, s in enumerate(all_samples) if sample_pop.get(s) == POP1]
        n = 2 * len(pop_idx)

        positions = []
        ac_alt_list = []
        region = f"{CHR}:{START}-{END}"
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
        pi = pi_sum / n_total if n_total > 0 else 0.0

        seg_mask = (ac_alt > 0) & (ac_alt < n)
        ac_seg = ac[seg_mask]
        pos_seg = pos[seg_mask]
        a_n = sum(1.0 / k for k in range(1, n))
        theta_w = n_seg / a_n / n_total if n_total > 0 else 0.0
        tajima_d = allel.tajima_d(ac_seg, pos=pos_seg) if n_seg > 3 else float('nan')

        return {"stat": "diversity", "pi": float(pi), "theta_w": float(theta_w),
                "tajima_d": float(tajima_d), "n_variants": n_total}

    return measure("Diversity", _compute)


def run_fst():
    """Full chr22 Fst (Hudson) — CEU vs YRI."""
    import allel
    import cyvcf2

    def _compute():
        sample_pop = load_panel()
        vcf = cyvcf2.VCF(VCF)
        all_samples = list(vcf.samples)
        idx1 = [i for i, s in enumerate(all_samples) if sample_pop.get(s) == POP1]
        idx2 = [i for i, s in enumerate(all_samples) if sample_pop.get(s) == POP2]
        n1, n2 = 2 * len(idx1), 2 * len(idx2)

        geno1, geno2 = [], []
        region = f"{CHR}:{START}-{END}"
        for v in vcf(region):
            if not v.is_snp or len(v.ALT) != 1:
                continue
            gt = v.genotypes
            g1 = np.array([gt[i][0] + gt[i][1] for i in idx1], dtype=np.int8)
            g2 = np.array([gt[i][0] + gt[i][1] for i in idx2], dtype=np.int8)
            geno1.append(g1)
            geno2.append(g2)
        vcf.close()

        g1_arr = np.array(geno1)
        g2_arr = np.array(geno2)
        ac1_alt = g1_arr.sum(axis=1)
        ac1 = np.column_stack([n1 - ac1_alt, ac1_alt])
        ac2_alt = g2_arr.sum(axis=1)
        ac2 = np.column_stack([n2 - ac2_alt, ac2_alt])
        num, den = allel.hudson_fst(ac1, ac2)
        fst = np.sum(num) / np.sum(den)
        return {"stat": "fst", "fst": float(fst), "n_variants": len(geno1)}

    return measure("Fst/Dxy", _compute)


def run_sfs():
    """Full chr22 SFS — CEU."""
    import allel
    import cyvcf2

    def _compute():
        sample_pop = load_panel()
        vcf = cyvcf2.VCF(VCF)
        all_samples = list(vcf.samples)
        pop_idx = [i for i, s in enumerate(all_samples) if sample_pop.get(s) == POP1]
        n = 2 * len(pop_idx)

        ac_alt_list = []
        region = f"{CHR}:{START}-{END}"
        for v in vcf(region):
            if not v.is_snp or len(v.ALT) != 1:
                continue
            gt = v.genotypes
            ac = sum(gt[i][0] + gt[i][1] for i in pop_idx)
            ac_alt_list.append(ac)
        vcf.close()

        ac_alt = np.array(ac_alt_list)
        ac_ref = n - ac_alt
        ac = np.column_stack([ac_ref, ac_alt])
        sfs = allel.sfs(ac_alt, n=n)
        return {"stat": "sfs", "sfs_sum": int(sfs.sum()), "n_variants": len(ac_alt_list)}

    return measure("SFS", _compute)


def run_ihs():
    """Full chr22 iHS — CEU."""
    import allel
    import cyvcf2

    def _compute():
        sample_pop = load_panel()
        vcf = cyvcf2.VCF(VCF)
        all_samples = list(vcf.samples)
        pop_idx = [i for i, s in enumerate(all_samples) if sample_pop.get(s) == POP1]

        positions = []
        hap_list = []
        region = f"{CHR}:{START}-{END}"
        for v in vcf(region):
            if not v.is_snp or len(v.ALT) != 1:
                continue
            gt = v.genotypes
            h = []
            for i in pop_idx:
                h.append(gt[i][0])
                h.append(gt[i][1])
            hap_list.append(np.array(h, dtype=np.int8))
            positions.append(v.POS)
        vcf.close()

        hap = np.array(hap_list)
        pos = np.array(positions)

        # Filter to segregating
        ac = hap.sum(axis=1)
        n_hap = hap.shape[1]
        seg = (ac > 0) & (ac < n_hap)
        hap_seg = hap[seg]
        pos_seg = pos[seg]

        ihs_raw = allel.ihs(hap_seg, pos_seg, include_edges=True)
        valid = ~np.isnan(ihs_raw)
        ihs_std = (ihs_raw[valid] - ihs_raw[valid].mean()) / ihs_raw[valid].std()
        return {"stat": "ihs", "n_variants": int(valid.sum())}

    return measure("iHS", _compute)


def run_xpehh():
    """Full chr22 XP-EHH — CEU vs YRI."""
    import allel
    import cyvcf2

    def _compute():
        sample_pop = load_panel()
        vcf = cyvcf2.VCF(VCF)
        all_samples = list(vcf.samples)
        idx1 = [i for i, s in enumerate(all_samples) if sample_pop.get(s) == POP1]
        idx2 = [i for i, s in enumerate(all_samples) if sample_pop.get(s) == POP2]

        positions = []
        hap1_list, hap2_list = [], []
        region = f"{CHR}:{START}-{END}"
        for v in vcf(region):
            if not v.is_snp or len(v.ALT) != 1:
                continue
            gt = v.genotypes
            h1, h2 = [], []
            for i in idx1:
                h1.append(gt[i][0])
                h1.append(gt[i][1])
            for i in idx2:
                h2.append(gt[i][0])
                h2.append(gt[i][1])
            hap1_list.append(np.array(h1, dtype=np.int8))
            hap2_list.append(np.array(h2, dtype=np.int8))
            positions.append(v.POS)
        vcf.close()

        hap1 = np.array(hap1_list)
        hap2 = np.array(hap2_list)
        pos = np.array(positions)

        xp = allel.xpehh(hap1, hap2, pos, include_edges=True)
        valid = ~np.isnan(xp)
        return {"stat": "xpehh", "n_variants": int(valid.sum())}

    return measure("XP-EHH", _compute)


def main():
    print(f"scikit-allel memory measurement — full chr22, {POP1} (vs {POP2} for pairwise)")
    print(f"Region: {CHR}:{START}-{END}\n")

    results = {}

    print("Running diversity...")
    results["diversity"] = run_diversity()

    print("Running Fst...")
    results["fst"] = run_fst()

    print("Running SFS...")
    results["sfs"] = run_sfs()

    print("Running iHS...")
    results["ihs"] = run_ihs()

    print("Running XP-EHH...")
    results["xpehh"] = run_xpehh()

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY — scikit-allel peak memory (tracemalloc)")
    print("=" * 60)
    worst = 0
    worst_stat = ""
    for stat, r in results.items():
        mem = r["peak_mem_mb"]
        t = r["time_s"]
        print(f"  {stat:12s}: {mem:8.1f} MB  ({t:.1f}s)")
        if mem > worst:
            worst = mem
            worst_stat = stat
    print(f"\n  WORST CASE:   {worst:.1f} MB  ({worst_stat})")

    out_path = ROOT / "benchmarks" / "results" / "scikit_allel_memory_ceu.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
