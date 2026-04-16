#!/usr/bin/env python3
"""Detailed comparison of iHS/XP-EHH: unstandardized vs standardized correlations.

Decomposes the correlation gap to show that the core EHH computation agrees
closely, and most deviation comes from standardization on a sub-chromosomal region.
"""

import os
import subprocess
import sys
import time

import allel
import numpy as np
from cyvcf2 import VCF

PANEL_PATH = "data/raw/1000g/20130606_g1k_3202_samples_ped_population.txt"
VCF_PATH = "data/raw/1000g/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz"
REGION = "chr22:16000000-17000000"


def load_panel():
    panel = {}
    with open(PANEL_PATH) as f:
        header = f.readline().strip().split()
        si, pi = header.index("SampleID"), header.index("Superpopulation")
        for line in f:
            parts = line.strip().split()
            panel[parts[si]] = parts[pi]
    return panel


def load_haplotypes(region, panel, pops):
    vcf = VCF(VCF_PATH)
    all_samples = list(vcf.samples)
    pop_idx = {p: [i for i, s in enumerate(all_samples) if panel.get(s) == p] for p in pops}

    positions = []
    hap_data = {p: [] for p in pops}

    for v in vcf(region):
        if v.FILTER is not None:
            continue
        if not v.ALT or v.ALT[0] in (".", "<*>", "*"):
            continue
        if len(v.ALT) > 1:
            continue
        if len(v.REF) != 1 or len(v.ALT[0]) != 1:
            continue

        positions.append(v.POS)
        genotypes = v.genotypes
        for p in pops:
            row = []
            for i in pop_idx[p]:
                a0, a1, _ = genotypes[i]
                row.extend([max(0, a0), max(0, a1)])
            hap_data[p].append(row)

    vcf.close()
    pos = np.array(positions)
    haps = {p: np.array(hap_data[p], dtype=np.int8) for p in pops}
    return pos, haps


def query_gp(cypher):
    result = subprocess.run(
        ["cypher-shell", "-u", "neo4j", "-p", os.environ.get("GRAPHPOP_PASSWORD", "graphpop"), "--format", "plain", cypher],
        capture_output=True, text=True, timeout=120,
    )
    return result.stdout.strip().split("\n")


def main():
    print("=" * 72)
    print("Detailed iHS / XP-EHH unstandardized vs standardized comparison")
    print(f"Region: {REGION}")
    print("=" * 72)
    print()

    panel = load_panel()
    t0 = time.time()
    pos, haps = load_haplotypes(REGION, panel, ["EUR", "EAS"])
    print(f"Loaded {len(pos)} SNPs in {time.time() - t0:.1f}s")
    print(f"  EUR: {haps['EUR'].shape[1]} haplotypes")
    print(f"  EAS: {haps['EAS'].shape[1]} haplotypes")

    # ---- scikit-allel iHS (returns unstandardized) ----
    ihs_unstd_allel = allel.ihs(
        allel.HaplotypeArray(haps["EUR"]), pos,
        min_maf=0.05, include_edges=False,
    )
    # Standardize manually
    ac_eur = allel.HaplotypeArray(haps["EUR"]).count_alleles()[:, 1]
    ihs_std_allel, _ = allel.standardize_by_allele_count(
        ihs_unstd_allel, ac_eur, n_bins=20, diagnostics=False,
    )

    # ---- scikit-allel XP-EHH (returns unstandardized) ----
    xpehh_unstd_allel = allel.xpehh(
        allel.HaplotypeArray(haps["EUR"]),
        allel.HaplotypeArray(haps["EAS"]),
        pos, include_edges=False,
    )
    # XP-EHH standardized: genome-wide z-score
    valid_xp = ~np.isnan(xpehh_unstd_allel)
    xpehh_std_allel = np.full_like(xpehh_unstd_allel, np.nan)
    if np.sum(valid_xp) > 1:
        m = np.nanmean(xpehh_unstd_allel[valid_xp])
        s = np.nanstd(xpehh_unstd_allel[valid_xp], ddof=1)
        xpehh_std_allel[valid_xp] = (xpehh_unstd_allel[valid_xp] - m) / s

    # ---- GraphPop iHS ----
    lines = query_gp(
        "CALL graphpop.ihs('chr22', 'EUR', {start: 16000000, end: 17000000}) "
        "YIELD pos, ihs_unstd, ihs RETURN pos, ihs_unstd, ihs;"
    )
    gp_ihs = {}
    for line in lines[1:]:
        parts = line.split(",")
        p = int(parts[0].strip())
        gp_ihs[p] = (float(parts[1].strip()), float(parts[2].strip()))

    # ---- GraphPop XP-EHH ----
    lines = query_gp(
        "CALL graphpop.xpehh('chr22', 'EUR', 'EAS', {start: 16000000, end: 17000000}) "
        "YIELD pos, xpehh_unstd, xpehh RETURN pos, xpehh_unstd, xpehh;"
    )
    gp_xpehh = {}
    for line in lines[1:]:
        parts = line.split(",")
        p = int(parts[0].strip())
        gp_xpehh[p] = (float(parts[1].strip()), float(parts[2].strip()))

    # ---- Match iHS ----
    valid_allel = ~np.isnan(ihs_unstd_allel)
    ihs_u_a, ihs_u_g, ihs_s_a, ihs_s_g = [], [], [], []
    for i, p in enumerate(pos):
        if valid_allel[i] and p in gp_ihs:
            ihs_u_a.append(ihs_unstd_allel[i])
            ihs_u_g.append(gp_ihs[p][0])
            if not np.isnan(ihs_std_allel[i]):
                ihs_s_a.append(ihs_std_allel[i])
                ihs_s_g.append(gp_ihs[p][1])

    ihs_u_a, ihs_u_g = np.array(ihs_u_a), np.array(ihs_u_g)
    ihs_s_a, ihs_s_g = np.array(ihs_s_a), np.array(ihs_s_g)

    # ---- Match XP-EHH ----
    xp_u_a, xp_u_g, xp_s_a, xp_s_g = [], [], [], []
    for i, p in enumerate(pos):
        if valid_xp[i] and p in gp_xpehh:
            xp_u_a.append(xpehh_unstd_allel[i])
            xp_u_g.append(gp_xpehh[p][0])
            if not np.isnan(xpehh_std_allel[i]):
                xp_s_a.append(xpehh_std_allel[i])
                xp_s_g.append(gp_xpehh[p][1])

    xp_u_a, xp_u_g = np.array(xp_u_a), np.array(xp_u_g)
    xp_s_a, xp_s_g = np.array(xp_s_a), np.array(xp_s_g)

    # ---- Variant set overlap ----
    gp_ihs_set = set(gp_ihs.keys())
    allel_ihs_set = set(pos[i] for i in range(len(pos)) if valid_allel[i])
    gp_xp_set = set(gp_xpehh.keys())
    allel_xp_set = set(pos[i] for i in range(len(pos)) if valid_xp[i])

    print()
    print("-" * 72)
    print("iHS")
    print("-" * 72)
    print(f"  Variant sets: allel={len(allel_ihs_set)}, GP={len(gp_ihs_set)}, "
          f"overlap={len(allel_ihs_set & gp_ihs_set)}")
    print(f"    allel-only: {len(allel_ihs_set - gp_ihs_set)}, "
          f"GP-only: {len(gp_ihs_set - allel_ihs_set)}")
    print(f"  Matched unstandardized: {len(ihs_u_a)}")
    print(f"  Matched standardized:   {len(ihs_s_a)}")
    r_u = np.corrcoef(ihs_u_a, ihs_u_g)[0, 1]
    r_s = np.corrcoef(ihs_s_a, ihs_s_g)[0, 1]
    mae_u = np.mean(np.abs(ihs_u_a - ihs_u_g))
    mae_s = np.mean(np.abs(ihs_s_a - ihs_s_g))
    print(f"  UNSTANDARDIZED  r = {r_u:.6f}  MAE = {mae_u:.4f}")
    print(f"  STANDARDIZED    r = {r_s:.6f}  MAE = {mae_s:.4f}")

    # Rank correlation (Spearman)
    from scipy.stats import spearmanr
    rho_u, _ = spearmanr(ihs_u_a, ihs_u_g)
    rho_s, _ = spearmanr(ihs_s_a, ihs_s_g)
    print(f"  UNSTANDARDIZED  rho = {rho_u:.6f}  (Spearman)")
    print(f"  STANDARDIZED    rho = {rho_s:.6f}  (Spearman)")

    print()
    print("-" * 72)
    print("XP-EHH")
    print("-" * 72)
    print(f"  Variant sets: allel={len(allel_xp_set)}, GP={len(gp_xp_set)}, "
          f"overlap={len(allel_xp_set & gp_xp_set)}")
    print(f"    allel-only: {len(allel_xp_set - gp_xp_set)}, "
          f"GP-only: {len(gp_xp_set - allel_xp_set)}")
    print(f"  Matched unstandardized: {len(xp_u_a)}")
    print(f"  Matched standardized:   {len(xp_s_a)}")
    r_u = np.corrcoef(xp_u_a, xp_u_g)[0, 1]
    r_s = np.corrcoef(xp_s_a, xp_s_g)[0, 1]
    mae_u = np.mean(np.abs(xp_u_a - xp_u_g))
    mae_s = np.mean(np.abs(xp_s_a - xp_s_g))
    print(f"  UNSTANDARDIZED  r = {r_u:.6f}  MAE = {mae_u:.4f}")
    print(f"  STANDARDIZED    r = {r_s:.6f}  MAE = {mae_s:.4f}")

    rho_u, _ = spearmanr(xp_u_a, xp_u_g)
    rho_s, _ = spearmanr(xp_s_a, xp_s_g)
    print(f"  UNSTANDARDIZED  rho = {rho_u:.6f}  (Spearman)")
    print(f"  STANDARDIZED    rho = {rho_s:.6f}  (Spearman)")

    print()
    print("=" * 72)


if __name__ == "__main__":
    main()
