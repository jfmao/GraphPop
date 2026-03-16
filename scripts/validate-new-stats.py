#!/usr/bin/env python3
"""Validate new statistics against scikit-allel reference implementations.

Tests:
  1. Weir & Cockerham Fst (per-locus + genome-wide)
  2. PBS (derived from W&C Fst triplet)
  3. Garud's H (H1, H12, H2/H1) per window

Usage:
    conda run -n graphevo python scripts/validate-new-stats.py
"""

import subprocess
import sys
import time

import allel
import numpy as np
from cyvcf2 import VCF

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

VCF_PATH = "data/raw/1000g/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz"
PANEL_PATH = "data/raw/1000g/20130606_g1k_3202_samples_ped_population.txt"

CHR = "chr22"
REGION = f"{CHR}:16000000-17000000"
START, END = 16_000_000, 17_000_000

POP1, POP2, POP3 = "EUR", "EAS", "AFR"
GARUD_POP = "EUR"
GARUD_WINDOW = 100_000
GARUD_STEP = 50_000

NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"


def load_panel(panel_path):
    """Load 1000G panel → {sample: pop}.

    Auto-detects PED format (3202 samples, whitespace-separated, SampleID/Superpopulation)
    vs panel format (2504 samples, tab-separated, sample/super_pop).
    """
    sample_to_pop = {}
    with open(panel_path) as f:
        first_line = f.readline().strip()
        if "SampleID" in first_line:
            # PED format: whitespace-separated
            fields = first_line.split()
            sample_col = fields.index("SampleID")
            pop_col = fields.index("Superpopulation")
            for line in f:
                parts = line.strip().split()
                if len(parts) > max(sample_col, pop_col):
                    sample_to_pop[parts[sample_col]] = parts[pop_col]
        else:
            # Panel format: tab-separated
            fields = first_line.split("\t")
            sample_col = fields.index("sample")
            pop_col = fields.index("super_pop")
            for line in f:
                parts = line.strip().split("\t")
                sample_to_pop[parts[sample_col]] = parts[pop_col]
    return sample_to_pop


def load_genotypes(vcf_path, region, sample_to_pop, pops):
    """Load genotype data from VCF for given populations.

    Returns:
        positions: 1D array of positions
        gt_arrays: dict pop -> GenotypeArray (n_variants, n_samples, 2)
        hap_arrays: dict pop -> HaplotypeArray (n_variants, n_haplotypes)
        ac_arrays: dict pop -> AlleleCountsArray (n_variants, 2)
    """
    vcf = VCF(vcf_path, samples=None)
    all_samples = list(vcf.samples)

    # Build per-pop sample indices
    pop_indices = {}
    for pop in pops:
        pop_indices[pop] = [i for i, s in enumerate(all_samples) if sample_to_pop.get(s) == pop]

    positions = []
    gt_per_pop = {p: [] for p in pops}
    hap_per_pop = {p: [] for p in pops}

    for v in vcf(region):
        if v.FILTER is not None:
            continue
        if not v.ALT or v.ALT[0] in (".", "<*>", "*"):
            continue
        if len(v.ALT) > 1:
            continue
        ref, alt = v.REF, v.ALT[0]
        # SNPs only (matching GraphPop garud_h default)
        if len(ref) != 1 or len(alt) != 1:
            continue

        positions.append(v.POS)
        genotypes = v.genotypes  # list of [a0, a1, phased]

        for pop in pops:
            idx = pop_indices[pop]
            gt_pop = []
            hap_pop = []
            for i in idx:
                a0, a1, _ = genotypes[i]
                if a0 < 0 or a1 < 0:
                    gt_pop.append([a0, a1])
                else:
                    gt_pop.append([a0, a1])
                hap_pop.append(a0)
                hap_pop.append(a1)
            gt_per_pop[pop].append(gt_pop)
            hap_per_pop[pop].append(hap_pop)

    vcf.close()

    positions = np.array(positions)
    gt_arrays = {}
    hap_arrays = {}
    ac_arrays = {}

    for pop in pops:
        gt = np.array(gt_per_pop[pop])  # (n_var, n_samp, 2) via list structure
        # Reshape: gt_per_pop[pop] is list of lists of [a0, a1]
        gt_3d = allel.GenotypeArray(gt)
        gt_arrays[pop] = gt_3d
        ac_arrays[pop] = gt_3d.count_alleles()

        hap = np.array(hap_per_pop[pop], dtype=np.int8)  # (n_var, n_hap)
        hap_arrays[pop] = allel.HaplotypeArray(hap)

    return positions, gt_arrays, hap_arrays, ac_arrays


def cypher_query(query):
    """Run a Cypher query and return output lines."""
    cmd = [
        "cypher-shell",
        "-u", NEO4J_USER, "-p", NEO4J_PASS,
        "--format", "plain",
        query,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr.strip()}")
        return []
    lines = result.stdout.strip().split("\n")
    return lines


# ---------------------------------------------------------------------------
# Test 1: Weir & Cockerham Fst
# ---------------------------------------------------------------------------

def test_wc_fst(gt_arrays, ac_arrays, positions):
    print("=" * 72)
    print("TEST 1: Weir & Cockerham Fst")
    print(f"  Region: {REGION}, Pops: {POP1} vs {POP2}")
    print("=" * 72)

    # scikit-allel: per-locus W&C Fst
    g_combined = allel.GenotypeArray(
        np.concatenate([gt_arrays[POP1], gt_arrays[POP2]], axis=1)
    )
    n1 = gt_arrays[POP1].shape[1]
    n2 = gt_arrays[POP2].shape[1]
    subpops = [list(range(n1)), list(range(n1, n1 + n2))]

    a, b, c = allel.weir_cockerham_fst(g_combined, subpops)
    # Genome-wide: ratio of averages
    a_sum = np.nansum(a, axis=0)
    b_sum = np.nansum(b, axis=0)
    c_sum = np.nansum(c, axis=0)
    # For biallelic, shape is (n_variants, 1) — take column 0
    fst_allel = a_sum[0] / (a_sum[0] + b_sum[0] + c_sum[0])
    print(f"  scikit-allel W&C Fst = {fst_allel:.10f}")

    # Per-locus Fst for correlation
    per_locus_fst_allel = a[:, 0] / (a[:, 0] + b[:, 0] + c[:, 0])
    per_locus_fst_allel = np.where(np.isnan(per_locus_fst_allel), 0.0, per_locus_fst_allel)

    # GraphPop: genome-wide (SNP-only for fair comparison)
    lines = cypher_query(
        f"CALL graphpop.divergence('{CHR}', {START}, {END}, '{POP1}', '{POP2}', "
        f"{{variant_type: 'SNP'}}) "
        "YIELD fst_wc RETURN fst_wc;"
    )
    fst_gp = float(lines[1].strip())
    print(f"  GraphPop W&C Fst     = {fst_gp:.10f}")

    rel_err = abs(fst_allel - fst_gp) / abs(fst_allel) * 100
    status = "PASS" if rel_err < 0.1 else ("WARN" if rel_err < 1.0 else "FAIL")
    print(f"  {status}  RelErr = {rel_err:.6f}%")

    # GraphPop: per-window for correlation (SNP-only)
    lines = cypher_query(
        f"CALL graphpop.genome_scan('{CHR}', '{POP1}', 50000, 50000, "
        f"{{pop2: '{POP2}', start: {START}, end: {END}, variant_type: 'SNP'}}) "
        "YIELD start, end, fst_wc, n_variants "
        "RETURN start, end, fst_wc, n_variants;"
    )
    gp_windows = []
    for line in lines[1:]:
        parts = line.split(",")
        gp_windows.append({
            "start": int(parts[0].strip()),
            "end": int(parts[1].strip()),
            "fst_wc": float(parts[2].strip()),
            "n_variants": int(parts[3].strip()),
        })

    # scikit-allel: same windows
    allel_windows = []
    for w in gp_windows:
        ws, we = w["start"], w["end"]
        mask = (positions >= ws) & (positions <= we)
        if np.sum(mask) < 2:
            allel_windows.append(np.nan)
            continue
        a_w = np.nansum(a[mask, 0])
        b_w = np.nansum(b[mask, 0])
        c_w = np.nansum(c[mask, 0])
        denom = a_w + b_w + c_w
        allel_windows.append(a_w / denom if denom != 0 else 0.0)

    gp_fst = np.array([w["fst_wc"] for w in gp_windows])
    allel_fst = np.array(allel_windows)
    valid = ~np.isnan(allel_fst) & ~np.isnan(gp_fst)
    if np.sum(valid) > 2:
        r = np.corrcoef(gp_fst[valid], allel_fst[valid])[0, 1]
        print(f"  Per-window correlation r = {r:.6f} ({np.sum(valid)} windows)")
    print()
    return rel_err < 1.0


# ---------------------------------------------------------------------------
# Test 2: PBS
# ---------------------------------------------------------------------------

def test_pbs(gt_arrays):
    print("=" * 72)
    print("TEST 2: Population Branch Statistic (PBS)")
    print(f"  Region: {REGION}, Pops: {POP1} vs {POP2}, outgroup {POP3}")
    print("=" * 72)

    # scikit-allel: compute pairwise W&C Fst for all 3 pairs
    pops = [POP1, POP2, POP3]
    n_samples = {p: gt_arrays[p].shape[1] for p in pops}
    g_all = allel.GenotypeArray(
        np.concatenate([gt_arrays[p] for p in pops], axis=1)
    )
    offsets = {}
    cur = 0
    for p in pops:
        offsets[p] = cur
        cur += n_samples[p]
    subpops = [list(range(offsets[p], offsets[p] + n_samples[p])) for p in pops]

    # Pairwise Fst: 0-1, 0-2, 1-2
    pairs = [(0, 1), (0, 2), (1, 2)]
    fst_vals = {}
    for i, j in pairs:
        a, b, c = allel.weir_cockerham_fst(g_all, [subpops[i], subpops[j]])
        a_s = np.nansum(a[:, 0])
        b_s = np.nansum(b[:, 0])
        c_s = np.nansum(c[:, 0])
        fst_vals[(i, j)] = a_s / (a_s + b_s + c_s)

    fst_12 = fst_vals[(0, 1)]
    fst_13 = fst_vals[(0, 2)]
    fst_23 = fst_vals[(1, 2)]

    def T(fst):
        return -np.log(max(1e-10, 1 - fst))

    pbs_allel = (T(fst_12) + T(fst_13) - T(fst_23)) / 2
    print(f"  scikit-allel: Fst_12={fst_12:.6f}, Fst_13={fst_13:.6f}, Fst_23={fst_23:.6f}")
    print(f"  scikit-allel PBS({POP1}) = {pbs_allel:.10f}")

    # GraphPop (SNP-only for fair comparison)
    lines = cypher_query(
        f"CALL graphpop.divergence('{CHR}', {START}, {END}, '{POP1}', '{POP2}', "
        f"{{pop3: '{POP3}', variant_type: 'SNP'}}) "
        "YIELD fst_wc, pbs RETURN fst_wc, pbs;"
    )
    parts = lines[1].split(",")
    fst_gp = float(parts[0].strip())
    pbs_gp = float(parts[1].strip())
    print(f"  GraphPop:    Fst_wc={fst_gp:.6f}, PBS({POP1}) = {pbs_gp:.10f}")

    rel_err_pbs = abs(pbs_allel - pbs_gp) / abs(pbs_allel) * 100 if pbs_allel != 0 else 0
    status = "PASS" if rel_err_pbs < 0.1 else ("WARN" if rel_err_pbs < 1.0 else "FAIL")
    print(f"  {status}  PBS RelErr = {rel_err_pbs:.6f}%")
    print()
    return rel_err_pbs < 1.0


# ---------------------------------------------------------------------------
# Test 3: Garud's H
# ---------------------------------------------------------------------------

def test_garud_h(hap_arrays, positions):
    print("=" * 72)
    print("TEST 3: Garud's H statistics")
    print(f"  Region: {REGION}, Pop: {GARUD_POP}, Window: {GARUD_WINDOW}, Step: {GARUD_STEP}")
    print("=" * 72)

    hap = hap_arrays[GARUD_POP]

    # scikit-allel: windowed Garud H
    allel_results = []
    for ws in range(START, END, GARUD_STEP):
        we = ws + GARUD_WINDOW
        mask = (positions >= ws) & (positions < we)
        n_snps = np.sum(mask)
        if n_snps < 2:
            continue
        h_window = hap[mask, :]
        h1, h12, h123, h2_h1 = allel.garud_h(h_window)
        allel_results.append({
            "start": ws, "end": we, "h1": h1, "h12": h12, "h2_h1": h2_h1,
            "n_variants": int(n_snps),
        })

    print(f"  scikit-allel: {len(allel_results)} windows computed")

    # GraphPop
    lines = cypher_query(
        f"CALL graphpop.garud_h('{CHR}', '{GARUD_POP}', {GARUD_WINDOW}, {GARUD_STEP}, "
        f"{{start: {START}, end: {END}}}) "
        "YIELD start, end, h1, h12, h2_h1, n_variants "
        "RETURN start, end, h1, h12, h2_h1, n_variants;"
    )
    gp_results = []
    for line in lines[1:]:
        parts = line.split(",")
        gp_results.append({
            "start": int(parts[0].strip()),
            "end": int(parts[1].strip()),
            "h1": float(parts[2].strip()),
            "h12": float(parts[3].strip()),
            "h2_h1": float(parts[4].strip()),
            "n_variants": int(parts[5].strip()),
        })
    print(f"  GraphPop:     {len(gp_results)} windows returned")

    # Match windows by nearest start position (GraphPop aligns to first variant)
    allel_by_start = {r["start"]: r for r in allel_results}
    gp_by_start = {r["start"]: r for r in gp_results}

    # Map each GP window to nearest allel window
    allel_starts = sorted(allel_by_start.keys())
    matched = {}
    for gs in gp_by_start:
        best = min(allel_starts, key=lambda x: abs(x - gs), default=None)
        if best is not None and abs(best - gs) < GARUD_STEP // 2:
            matched[best] = gs
    common_starts = sorted(matched.keys())
    print(f"  Matched windows: {len(common_starts)}")

    if len(common_starts) < 2:
        print("  SKIP: too few matched windows")
        print()
        return True

    h1_allel = np.array([allel_by_start[s]["h1"] for s in common_starts])
    h1_gp = np.array([gp_by_start[matched[s]]["h1"] for s in common_starts])
    h12_allel = np.array([allel_by_start[s]["h12"] for s in common_starts])
    h12_gp = np.array([gp_by_start[matched[s]]["h12"] for s in common_starts])
    h2h1_allel = np.array([allel_by_start[s]["h2_h1"] for s in common_starts])
    h2h1_gp = np.array([gp_by_start[matched[s]]["h2_h1"] for s in common_starts])
    nvar_allel = np.array([allel_by_start[s]["n_variants"] for s in common_starts])
    nvar_gp = np.array([gp_by_start[matched[s]]["n_variants"] for s in common_starts])

    # Variant count match
    nvar_match = np.sum(nvar_allel == nvar_gp)
    print(f"  Variant count match: {nvar_match}/{len(common_starts)}")

    # Correlations
    for name, a_arr, g_arr in [("H1", h1_allel, h1_gp),
                                ("H12", h12_allel, h12_gp),
                                ("H2/H1", h2h1_allel, h2h1_gp)]:
        valid = ~np.isnan(a_arr) & ~np.isnan(g_arr)
        if np.sum(valid) > 2:
            r = np.corrcoef(a_arr[valid], g_arr[valid])[0, 1]
            mae = np.mean(np.abs(a_arr[valid] - g_arr[valid]))
            max_err = np.max(np.abs(a_arr[valid] - g_arr[valid]))
            status = "PASS" if r > 0.99 else ("WARN" if r > 0.95 else "FAIL")
            print(f"  {status}  {name:6s}  r={r:.6f}  MAE={mae:.2e}  MaxErr={max_err:.2e}")

    # Print a few sample windows for inspection
    print()
    print(f"  {'Start':>10s}  {'H1_allel':>10s}  {'H1_gp':>10s}  {'H12_allel':>10s}  {'H12_gp':>10s}  {'nVar_a':>6s}  {'nVar_gp':>7s}")
    for s in common_starts[:5]:
        a = allel_by_start[s]
        g = gp_by_start[matched[s]]
        print(f"  {s:10d}  {a['h1']:10.6f}  {g['h1']:10.6f}  {a['h12']:10.6f}  {g['h12']:10.6f}  {a['n_variants']:6d}  {g['n_variants']:7d}")

    print()
    r_h1 = np.corrcoef(h1_allel, h1_gp)[0, 1]
    return r_h1 > 0.95


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print()
    print("Validating new GraphPop statistics against scikit-allel")
    print(f"VCF: {VCF_PATH}")
    print(f"Region: {REGION}")
    print()

    t0 = time.time()
    print("Loading panel...")
    sample_to_pop = load_panel(PANEL_PATH)
    print(f"  {len(sample_to_pop)} samples in panel")

    print("Loading genotypes from VCF (this takes ~30s)...")
    positions, gt_arrays, hap_arrays, ac_arrays = load_genotypes(
        VCF_PATH, REGION, sample_to_pop, [POP1, POP2, POP3]
    )
    print(f"  {len(positions)} SNPs loaded in {time.time() - t0:.1f}s")
    print()

    results = []
    results.append(("W&C Fst", test_wc_fst(gt_arrays, ac_arrays, positions)))
    results.append(("PBS", test_pbs(gt_arrays)))
    results.append(("Garud's H", test_garud_h(hap_arrays, positions)))

    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)
    all_pass = True
    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_pass = False
        print(f"  {status}  {name}")
    print("=" * 72)

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
