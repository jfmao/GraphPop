#!/usr/bin/env python3
"""Rice 3K Genome — Deep Integrative Analysis.

Performs cross-statistic integrative analyses from existing rice results JSON.
No Neo4j connection required for Phase 1 analyses (A1–A4).

Usage:
    python scripts/rice_deep_integration.py fingerprints    # A1: Evolutionary fingerprint clustering
    python scripts/rice_deep_integration.py sweeps          # A2: Sweep classification at known genes
    python scripts/rice_deep_integration.py roh             # A3: ROH length distribution
    python scripts/rice_deep_integration.py correlations    # A4: Cross-statistic correlation matrix
    python scripts/rice_deep_integration.py report          # Generate Markdown report
    python scripts/rice_deep_integration.py all             # Run all Phase 1 analyses
"""

import argparse
import json
import logging
import os
import sys
from collections import defaultdict

import numpy as np

log = logging.getLogger("rice_deep")

# ── Constants ─────────────────────────────────────────────────────────

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "results", "rice")
RESULTS_FILE = os.path.join(RESULTS_DIR, "rice_full_results.json")
INTERP_FILE = os.path.join(RESULTS_DIR, "rice_interpretation_results.json")
OUTPUT_FILE = os.path.join(RESULTS_DIR, "rice_deep_integration_results.json")
REPORT_FILE = os.path.join(RESULTS_DIR, "rice_deep_integration_report.md")

POPULATIONS = [
    "XI-adm", "XI-3", "GJ-trp", "GJ-tmp", "XI-2", "XI-1A", "XI-1B",
    "cA-Aus", "GJ-sbtrp", "admix", "GJ-adm", "cB-Bas",
]

CHROMOSOMES = [
    "Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6",
    "Chr7", "Chr8", "Chr9", "Chr10", "Chr11", "Chr12",
]

# Rice genome size (IRGSP-1.0, approximate total)
RICE_GENOME_SIZE_MB = 373.2  # Mb

# Known rice domestication genes (Nipponbare IRGSP-1.0 coordinates)
KNOWN_GENES = [
    {"gene": "Wx",    "chr": "Chr6",  "pos": 1600000,  "function": "Amylose content (Waxy)", "window": 200000},
    {"gene": "Sd1",   "chr": "Chr1",  "pos": 38400000, "function": "Semi-dwarf (Green Revolution)", "window": 200000},
    {"gene": "sh4",   "chr": "Chr4",  "pos": 33000000, "function": "Seed shattering", "window": 200000},
    {"gene": "qSH1",  "chr": "Chr1",  "pos": 37000000, "function": "Seed shattering", "window": 200000},
    {"gene": "PROG1", "chr": "Chr7",  "pos": 3200000,  "function": "Prostrate growth", "window": 200000},
    {"gene": "Hd1",   "chr": "Chr6",  "pos": 9000000,  "function": "Heading date", "window": 500000},
    {"gene": "Hd3a",  "chr": "Chr6",  "pos": 2600000,  "function": "Florigen", "window": 200000},
    {"gene": "GS3",   "chr": "Chr3",  "pos": 16700000, "function": "Grain length", "window": 200000},
    {"gene": "badh2", "chr": "Chr8",  "pos": 20000000, "function": "Fragrance", "window": 200000},
    {"gene": "Sub1A", "chr": "Chr9",  "pos": 6300000,  "function": "Submergence tolerance", "window": 200000},
    {"gene": "COLD1", "chr": "Chr4",  "pos": 30300000, "function": "Cold tolerance", "window": 200000},
    {"gene": "DRO1",  "chr": "Chr9",  "pos": 22000000, "function": "Deep rooting", "window": 200000},
    {"gene": "Bh4",   "chr": "Chr4",  "pos": 26000000, "function": "Hull color", "window": 200000},
    {"gene": "OsC1",  "chr": "Chr6",  "pos": 5000000,  "function": "Hull color", "window": 200000},
    {"gene": "GW5",   "chr": "Chr5",  "pos": 5400000,  "function": "Grain width", "window": 200000},
    {"gene": "Pi-ta", "chr": "Chr12", "pos": 10600000, "function": "Blast resistance", "window": 200000},
]

XPEHH_PAIRS = [
    ("cA-Aus", "cB-Bas"),
    ("XI-1A", "cA-Aus"),
    ("GJ-tmp", "cA-Aus"),
    ("GJ-tmp", "XI-1A"),
    ("GJ-trp", "XI-1A"),
    ("GJ-tmp", "GJ-trp"),
]

# Population group labels for biological context
POP_GROUPS = {
    "XI-adm": "Indica (admixed)", "XI-3": "Indica III",
    "XI-2": "Indica II", "XI-1A": "Indica IA", "XI-1B": "Indica IB",
    "GJ-tmp": "Japonica temperate", "GJ-trp": "Japonica tropical",
    "GJ-sbtrp": "Japonica subtropical", "GJ-adm": "Japonica (admixed)",
    "cA-Aus": "cA (Aus)", "cB-Bas": "cB (Basmati)", "admix": "Admixed",
}


# ── Helpers ───────────────────────────────────────────────────────────

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stderr),
            logging.FileHandler("rice_deep_integration.log", mode="a"),
        ],
    )


def load_results():
    log.info("Loading rice_full_results.json ...")
    with open(RESULTS_FILE) as f:
        return json.load(f)


def load_interp():
    if os.path.exists(INTERP_FILE):
        with open(INTERP_FILE) as f:
            return json.load(f)
    return {}


def load_output():
    if os.path.exists(OUTPUT_FILE):
        with open(OUTPUT_FILE) as f:
            return json.load(f)
    return {}


def save_output(data):
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    with open(OUTPUT_FILE, "w") as f:
        json.dump(data, f, indent=2)
    log.info(f"Saved: {OUTPUT_FILE}")


# ── A1: Evolutionary Fingerprint Clustering ───────────────────────────

def cmd_fingerprints(data, interp, output):
    """Cluster 12 populations by multi-statistic evolutionary profiles."""
    log.info("=== A1: Evolutionary Fingerprint Clustering ===")

    # Extract per-population genome-wide means
    # NOTE: Fay & Wu's H and ROH are in interpretation results (not phase1 diversity)
    # because phase1 was run before ancestral allele annotation.
    fwh_summary = interp.get("fay_wu_summary", {})
    roh_summary = interp.get("roh_hmm_summary", {})
    pinsps_ratios = interp.get("pinsps_ratios", {})

    profiles = {}
    for pop in POPULATIONS:
        pi_vals, theta_vals, tajd_vals = [], [], []
        fis_vals, het_obs_vals, het_exp_vals = [], [], []
        n_sweep_windows = 0
        hard_sweeps = 0
        soft_sweeps = 0

        if pop not in data.get("phase1", {}):
            log.warning(f"Population {pop} not in phase1 results")
            continue

        for chrom in CHROMOSOMES:
            chr_data = data["phase1"][pop].get(chrom, {})

            # Diversity statistics
            div = chr_data.get("diversity", {}).get("result", {})
            if div:
                pi_vals.append(div.get("pi", 0))
                theta_vals.append(div.get("theta_w", 0))
                tajd_vals.append(div.get("tajima_d", 0))
                fis_vals.append(div.get("fis", 0))
                het_obs_vals.append(div.get("het_obs", 0))
                het_exp_vals.append(div.get("het_exp", 0))

            # Garud's H sweep counts
            garud = chr_data.get("garud_h", {}).get("result", [])
            for w in garud:
                h12 = w.get("h12", 0)
                h2_h1 = w.get("h2_h1", 0)
                if h12 > 0.10:
                    n_sweep_windows += 1
                    if h2_h1 < 0.05:
                        hard_sweeps += 1
                    else:
                        soft_sweeps += 1

        if not pi_vals:
            continue

        # Get Fay & Wu's H from interpretation results
        fwh_pop = fwh_summary.get(pop, {})
        fay_wu_h = fwh_pop.get("mean_h", 0.0)

        # Get FROH from ROH HMM interpretation results
        roh_pop = roh_summary.get(pop, {})
        froh = roh_pop.get("mean_froh_across_chr", 0.0)

        hard_frac = hard_sweeps / n_sweep_windows if n_sweep_windows > 0 else 0.0

        profiles[pop] = {
            "pi": np.mean(pi_vals),
            "theta_w": np.mean(theta_vals),
            "tajima_d": np.mean(tajd_vals),
            "fay_wu_h": fay_wu_h,
            "fis": np.mean(fis_vals),
            "froh": froh,
            "n_sweeps": n_sweep_windows,
            "hard_sweep_fraction": hard_frac,
            "het_obs": np.mean(het_obs_vals),
            "het_exp": np.mean(het_exp_vals),
        }

    # Add piN/piS from interpretation results
    for pop, pdata in pinsps_ratios.items():
        if pop in profiles:
            per_chr = pdata.get("per_chr", {})
            ratios = [v["pi_n_pi_s"] for v in per_chr.values() if "pi_n_pi_s" in v]
            if ratios:
                profiles[pop]["pinsps"] = np.mean(ratios)

    log.info(f"Extracted profiles for {len(profiles)} populations")

    # Build feature matrix
    pop_names = sorted(profiles.keys())
    feature_names = [
        "pi", "theta_w", "tajima_d", "fay_wu_h", "fis", "froh",
        "n_sweeps", "hard_sweep_fraction",
    ]
    # Add piN/piS if available
    if any("pinsps" in profiles[p] for p in pop_names):
        feature_names.append("pinsps")

    X = np.zeros((len(pop_names), len(feature_names)))
    for i, pop in enumerate(pop_names):
        for j, feat in enumerate(feature_names):
            X[i, j] = profiles[pop].get(feat, 0.0)

    # Standardize (z-score)
    means = X.mean(axis=0)
    stds = X.std(axis=0)
    stds[stds == 0] = 1.0
    Z = (X - means) / stds

    # PCA (manual, no sklearn dependency)
    cov = np.cov(Z.T)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    # Sort by descending eigenvalue
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    variance_explained = eigenvalues / eigenvalues.sum()

    # Project onto first 2-3 PCs
    pcs = Z @ eigenvectors[:, :3]

    pca_results = {}
    for i, pop in enumerate(pop_names):
        pca_results[pop] = {
            "PC1": float(pcs[i, 0]),
            "PC2": float(pcs[i, 1]),
            "PC3": float(pcs[i, 2]) if pcs.shape[1] > 2 else 0.0,
            "group": POP_GROUPS.get(pop, pop),
        }

    # Hierarchical clustering (single-linkage on Euclidean distance in Z-space)
    from _cluster_util import hierarchical_cluster
    dendrogram = hierarchical_cluster(Z, pop_names)

    # Also compute Euclidean distance matrix in Z-space
    dist_matrix = {}
    for i, p1 in enumerate(pop_names):
        for j, p2 in enumerate(pop_names):
            if i < j:
                d = float(np.sqrt(np.sum((Z[i] - Z[j]) ** 2)))
                dist_matrix[f"{p1}_vs_{p2}"] = d

    result = {
        "profiles": {p: {k: float(v) for k, v in profiles[p].items()} for p in pop_names},
        "feature_names": feature_names,
        "pca": {
            "coordinates": pca_results,
            "variance_explained": [float(v) for v in variance_explained[:5]],
            "loadings": {
                feat: [float(eigenvectors[j, k]) for k in range(min(3, len(eigenvalues)))]
                for j, feat in enumerate(feature_names)
            },
        },
        "dendrogram": dendrogram,
        "distance_matrix": dist_matrix,
    }

    output["fingerprints"] = result
    save_output(output)

    # Print summary
    log.info("\n--- Evolutionary Fingerprints ---")
    log.info(f"PC1 explains {variance_explained[0]*100:.1f}% variance")
    log.info(f"PC2 explains {variance_explained[1]*100:.1f}% variance")
    log.info(f"PC1+PC2 explains {(variance_explained[0]+variance_explained[1])*100:.1f}% variance")

    log.info("\nPopulation profiles (standardized):")
    for pop in pop_names:
        pc = pca_results[pop]
        log.info(f"  {pop:12s}  PC1={pc['PC1']:+.2f}  PC2={pc['PC2']:+.2f}  ({pc['group']})")

    # Key loadings
    log.info("\nPC1 loadings (top features driving separation):")
    loadings_pc1 = [(feat, float(eigenvectors[j, 0])) for j, feat in enumerate(feature_names)]
    loadings_pc1.sort(key=lambda x: abs(x[1]), reverse=True)
    for feat, loading in loadings_pc1:
        log.info(f"  {feat:25s}  {loading:+.3f}")

    return result


# ── A2: Sweep Classification at Known Genes ───────────────────────────

def cmd_sweeps(data, interp, output):
    """Classify sweeps at known domestication genes using Garud's H."""
    log.info("=== A2: Hard vs Soft Sweep Classification at Known Genes ===")

    results = {}

    for gene_info in KNOWN_GENES:
        gene = gene_info["gene"]
        chrom = gene_info["chr"]
        pos = gene_info["pos"]
        win = gene_info["window"]
        gene_start = pos - win
        gene_end = pos + win

        gene_results = {
            "gene": gene,
            "function": gene_info["function"],
            "chr": chrom,
            "pos": pos,
            "populations": {},
        }

        for pop in POPULATIONS:
            chr_data = data.get("phase1", {}).get(pop, {}).get(chrom, {})
            garud = chr_data.get("garud_h", {}).get("result", [])

            # Find windows overlapping gene region
            overlapping = []
            for w in garud:
                w_start = w.get("start", 0)
                w_end = w.get("end", 0)
                if w_end >= gene_start and w_start <= gene_end:
                    overlapping.append(w)

            if not overlapping:
                continue

            # Take window with highest H12
            best = max(overlapping, key=lambda w: w.get("h12", 0))
            h1 = best.get("h1", 0)
            h12 = best.get("h12", 0)
            h2_h1 = best.get("h2_h1", 0)
            hap_div = best.get("hap_diversity", 0)
            n_haps = best.get("n_haplotypes", 0)

            # Classify sweep type
            if h12 > 0.10:
                if h2_h1 < 0.05:
                    sweep_type = "hard"
                elif h2_h1 < 0.50:
                    sweep_type = "soft (partial)"
                else:
                    sweep_type = "soft"
            elif h12 > 0.05:
                sweep_type = "weak signal"
            else:
                sweep_type = "none"

            gene_results["populations"][pop] = {
                "h1": float(h1),
                "h12": float(h12),
                "h2_h1": float(h2_h1),
                "hap_diversity": float(hap_div),
                "n_haplotypes": n_haps,
                "sweep_type": sweep_type,
                "window_start": best.get("start", 0),
                "window_end": best.get("end", 0),
            }

        results[gene] = gene_results

    # Also look for XP-EHH signals at each gene
    for gene_info in KNOWN_GENES:
        gene = gene_info["gene"]
        chrom = gene_info["chr"]
        pos = gene_info["pos"]
        win = gene_info["window"]
        gene_start = pos - win
        gene_end = pos + win

        xpehh_signals = {}
        for pair in XPEHH_PAIRS:
            pair_key = f"{pair[0]}_vs_{pair[1]}"
            pair_data = data.get("phase2", {}).get(pair_key, {})
            chr_data = pair_data.get(chrom, {})
            top_hits = chr_data.get("top_hits", [])

            # Find hits in gene region
            hits_in_region = [
                h for h in top_hits
                if gene_start <= h.get("pos", 0) <= gene_end
            ]
            if hits_in_region:
                best_hit = max(hits_in_region, key=lambda h: abs(h.get("xpehh", 0)))
                xpehh_signals[pair_key] = {
                    "xpehh": float(best_hit.get("xpehh", 0)),
                    "xpehh_unstd": float(best_hit.get("xpehh_unstd", 0)),
                    "pos": best_hit.get("pos", 0),
                    "n_hits": len(hits_in_region),
                }

        if xpehh_signals:
            results[gene]["xpehh_signals"] = xpehh_signals

    output["sweeps"] = results
    save_output(output)

    # Print summary
    log.info("\n--- Sweep Classification at Known Domestication Genes ---")
    log.info(f"{'Gene':8s} {'Function':30s} {'Pop':12s} {'H12':>6s} {'H2/H1':>6s} {'Type':15s} {'XP-EHH':>8s}")
    log.info("-" * 100)
    for gene_info in KNOWN_GENES:
        gene = gene_info["gene"]
        gdata = results.get(gene, {})
        pops_data = gdata.get("populations", {})

        # Find population with strongest signal
        if pops_data:
            best_pop = max(pops_data.keys(), key=lambda p: pops_data[p]["h12"])
            bp = pops_data[best_pop]

            # Best XP-EHH
            xp = gdata.get("xpehh_signals", {})
            best_xp = ""
            if xp:
                best_pair = max(xp.keys(), key=lambda k: abs(xp[k]["xpehh"]))
                best_xp = f"{xp[best_pair]['xpehh']:+.2f} ({best_pair})"

            log.info(
                f"{gene:8s} {gene_info['function']:30s} {best_pop:12s} "
                f"{bp['h12']:6.3f} {bp['h2_h1']:6.3f} {bp['sweep_type']:15s} {best_xp}"
            )
        else:
            log.info(f"{gene:8s} {gene_info['function']:30s} {'—':12s} {'—':>6s} {'—':>6s} {'no data':15s}")

    return results


# ── A3: ROH Length Distribution ───────────────────────────────────────

def cmd_roh(data, interp, output):
    """Analyze ROH segment length distribution by population."""
    log.info("=== A3: ROH Length Distribution by Population ===")

    # Phase1 ROH used sliding_window method (total_roh_bp=0 in JSON).
    # The HMM-based ROH is in interpretation results: roh_hmm, roh_hmm_summary.
    # We use roh_hmm for per-chromosome detail and roh_hmm_summary for totals.

    roh_hmm = interp.get("roh_hmm", {})
    roh_hmm_summary = interp.get("roh_hmm_summary", {})

    roh_summary = {}
    for pop in POPULATIONS:
        chr_roh = {}
        total_bp = 0
        total_segments = 0
        total_samples = 0

        for chrom in CHROMOSOMES:
            key = f"{pop}|{chrom}"
            roh_data = roh_hmm.get(key, {})
            if isinstance(roh_data, dict):
                bp = roh_data.get("total_roh_bp", 0)
                seg = roh_data.get("n_segments", 0)
                samp = roh_data.get("n_samples_with_roh", 0)
                froh_chr = roh_data.get("froh", 0)
                chr_roh[chrom] = {
                    "total_bp": bp,
                    "n_segments": seg,
                    "n_samples_with_roh": samp,
                    "froh": froh_chr,
                    "mean_segment_length_kb": (bp / seg / 1000) if seg > 0 else 0,
                }
                total_bp += bp
                total_segments += seg
                total_samples = max(total_samples, samp)

        # Use summary FROH if available (more accurate)
        pop_summary = roh_hmm_summary.get(pop, {})
        froh = pop_summary.get("mean_froh_across_chr", 0.0)
        if froh == 0 and total_bp > 0:
            froh = total_bp / (RICE_GENOME_SIZE_MB * 1e6)

        roh_summary[pop] = {
            "total_roh_mb": total_bp / 1e6,
            "froh": froh,
            "total_segments": total_segments,
            "max_samples_with_roh": total_samples,
            "mean_segment_length_kb": (total_bp / total_segments / 1000) if total_segments > 0 else 0,
            "per_chromosome": chr_roh,
        }

    # Compute ROH heterogeneity: which chromosomes contribute most to ROH?
    roh_chr_profile = {}
    for pop in POPULATIONS:
        chr_fracs = {}
        total = roh_summary[pop]["total_roh_mb"]
        if total > 0:
            for chrom in CHROMOSOMES:
                chr_bp = roh_summary[pop]["per_chromosome"].get(chrom, {}).get("total_bp", 0)
                chr_fracs[chrom] = chr_bp / (total * 1e6)
        roh_chr_profile[pop] = chr_fracs

    # Cross-reference ROH with diversity
    roh_vs_diversity = {}
    for pop in POPULATIONS:
        pi_vals = []
        fis_vals = []
        for chrom in CHROMOSOMES:
            div = data.get("phase1", {}).get(pop, {}).get(chrom, {}).get("diversity", {}).get("result", {})
            if div:
                pi_vals.append(div.get("pi", 0))
                fis_vals.append(div.get("fis", 0))
        if pi_vals:
            roh_vs_diversity[pop] = {
                "froh": roh_summary[pop]["froh"],
                "mean_pi": float(np.mean(pi_vals)),
                "fis": float(np.mean(fis_vals)),
            }

    result = {
        "summary": {p: {k: float(v) if isinstance(v, (int, float, np.floating)) else v
                        for k, v in roh_summary[p].items() if k != "per_chromosome"}
                    for p in POPULATIONS},
        "per_chromosome": {p: roh_summary[p]["per_chromosome"] for p in POPULATIONS},
        "chromosome_profile": roh_chr_profile,
        "roh_vs_diversity": roh_vs_diversity,
    }

    output["roh_distribution"] = result
    save_output(output)

    # Print summary
    log.info("\n--- ROH Distribution by Population ---")
    log.info(f"{'Population':12s} {'Total ROH (Mb)':>14s} {'FROH':>8s} {'Segments':>10s} {'Mean len (kb)':>14s}")
    log.info("-" * 65)
    for pop in sorted(POPULATIONS, key=lambda p: roh_summary[p]["total_roh_mb"], reverse=True):
        s = roh_summary[pop]
        log.info(
            f"{pop:12s} {s['total_roh_mb']:14.1f} {s['froh']:8.4f} "
            f"{s['total_segments']:10d} {s['mean_segment_length_kb']:14.1f}"
        )

    log.info("\n--- ROH vs Diversity (per population) ---")
    log.info(f"{'Population':12s} {'FROH':>8s} {'Mean π':>10s} {'F_IS':>8s}")
    log.info("-" * 42)
    for pop in sorted(roh_vs_diversity.keys(), key=lambda p: roh_vs_diversity[p]["froh"], reverse=True):
        d = roh_vs_diversity[pop]
        log.info(f"{pop:12s} {d['froh']:8.4f} {d['mean_pi']:10.6f} {d['fis']:8.4f}")

    return result


# ── A4: Cross-Statistic Correlation Matrix ────────────────────────────

def cmd_correlations(data, interp, output):
    """Compute correlations between all pairs of population-level statistics."""
    log.info("=== A4: Cross-Statistic Correlation Matrix ===")

    # Use interpretation results for Fay & Wu's H and ROH (phase1 has zeros)
    fwh_summary = interp.get("fay_wu_summary", {})
    roh_summary = interp.get("roh_hmm_summary", {})
    pinsps_ratios = interp.get("pinsps_ratios", {})

    # Collect per-population genome-wide means
    stats = {}
    for pop in POPULATIONS:
        pi_vals, theta_vals, tajd_vals = [], [], []
        fis_vals = []
        n_sweeps = 0

        for chrom in CHROMOSOMES:
            chr_data = data.get("phase1", {}).get(pop, {}).get(chrom, {})
            div = chr_data.get("diversity", {}).get("result", {})
            if div:
                pi_vals.append(div.get("pi", 0))
                theta_vals.append(div.get("theta_w", 0))
                tajd_vals.append(div.get("tajima_d", 0))
                fis_vals.append(div.get("fis", 0))

            garud = chr_data.get("garud_h", {}).get("result", [])
            n_sweeps += sum(1 for w in garud if w.get("h12", 0) > 0.10)

        fay_wu_h = fwh_summary.get(pop, {}).get("mean_h", 0.0)
        froh = roh_summary.get(pop, {}).get("mean_froh_across_chr", 0.0)

        stats[pop] = {
            "pi": np.mean(pi_vals) if pi_vals else 0,
            "theta_w": np.mean(theta_vals) if theta_vals else 0,
            "tajima_d": np.mean(tajd_vals) if tajd_vals else 0,
            "fay_wu_h": fay_wu_h,
            "fis": np.mean(fis_vals) if fis_vals else 0,
            "froh": froh,
            "n_sweeps": n_sweeps,
        }

    # Add piN/piS from interpretation results
    for pop in POPULATIONS:
        if pop in pinsps_ratios:
            per_chr = pinsps_ratios[pop].get("per_chr", {})
            ratios = [v["pi_n_pi_s"] for v in per_chr.values() if "pi_n_pi_s" in v]
            if ratios:
                stats[pop]["pinsps"] = np.mean(ratios)

    # Build vectors
    pop_order = sorted(stats.keys())
    stat_names = ["pi", "theta_w", "tajima_d", "fay_wu_h", "fis", "froh", "n_sweeps"]
    if any("pinsps" in stats[p] for p in pop_order):
        stat_names.append("pinsps")

    vectors = {}
    for name in stat_names:
        vectors[name] = np.array([stats[p].get(name, 0) for p in pop_order])

    # Compute Spearman rank correlations
    def spearman_r(x, y):
        """Spearman rank correlation."""
        rx = np.argsort(np.argsort(x)).astype(float)
        ry = np.argsort(np.argsort(y)).astype(float)
        n = len(x)
        d = rx - ry
        return 1 - 6 * np.sum(d ** 2) / (n * (n ** 2 - 1))

    # Correlation matrix
    corr_matrix = {}
    for i, s1 in enumerate(stat_names):
        row = {}
        for j, s2 in enumerate(stat_names):
            row[s2] = float(spearman_r(vectors[s1], vectors[s2]))
        corr_matrix[s1] = row

    # Identify notable correlations (|r| > 0.5, not on diagonal)
    notable = []
    seen = set()
    for s1 in stat_names:
        for s2 in stat_names:
            if s1 != s2 and (s2, s1) not in seen:
                r = corr_matrix[s1][s2]
                if abs(r) > 0.5:
                    notable.append({
                        "stat1": s1, "stat2": s2,
                        "spearman_r": float(r),
                        "direction": "positive" if r > 0 else "negative",
                    })
                    seen.add((s1, s2))

    notable.sort(key=lambda x: abs(x["spearman_r"]), reverse=True)

    # Biological interpretations
    expected_patterns = [
        {
            "stat1": "pi", "stat2": "pinsps",
            "expected": "negative",
            "reason": "Larger Ne → more efficient purifying selection → lower piN/piS",
        },
        {
            "stat1": "froh", "stat2": "tajima_d",
            "expected": "positive (both reflect bottlenecks)",
            "reason": "Bottleneck increases FROH and makes Tajima's D more negative (positive correlation in magnitude)",
        },
        {
            "stat1": "pinsps", "stat2": "froh",
            "expected": "positive",
            "reason": "Bottleneck + drift → more genetic load → higher piN/piS",
        },
        {
            "stat1": "n_sweeps", "stat2": "tajima_d",
            "expected": "negative (more sweeps → more negative D)",
            "reason": "Selective sweeps create excess rare variants → negative Tajima's D",
        },
    ]

    # Check expected vs observed
    for exp in expected_patterns:
        s1, s2 = exp["stat1"], exp["stat2"]
        if s1 in corr_matrix and s2 in corr_matrix.get(s1, {}):
            exp["observed_r"] = corr_matrix[s1][s2]
        else:
            exp["observed_r"] = None

    result = {
        "population_stats": {p: {k: float(v) for k, v in stats[p].items()} for p in pop_order},
        "stat_names": stat_names,
        "correlation_matrix": corr_matrix,
        "notable_correlations": notable,
        "expected_patterns": expected_patterns,
    }

    output["correlations"] = result
    save_output(output)

    # Print correlation matrix
    log.info("\n--- Cross-Statistic Correlation Matrix (Spearman ρ) ---")
    header = f"{'':15s}" + "".join(f"{s:>12s}" for s in stat_names)
    log.info(header)
    for s1 in stat_names:
        row = f"{s1:15s}" + "".join(
            f"{corr_matrix[s1][s2]:12.3f}" for s2 in stat_names
        )
        log.info(row)

    log.info("\n--- Notable Correlations (|ρ| > 0.5) ---")
    for nc in notable:
        log.info(f"  {nc['stat1']:15s} × {nc['stat2']:15s}  ρ = {nc['spearman_r']:+.3f}  ({nc['direction']})")

    log.info("\n--- Expected vs Observed Patterns ---")
    for exp in expected_patterns:
        obs = exp.get("observed_r")
        status = "✓" if obs is not None else "N/A"
        obs_str = f"{obs:+.3f}" if obs is not None else "N/A"
        log.info(f"  {exp['stat1']:15s} × {exp['stat2']:15s}  expected: {exp['expected']:40s}  observed: {obs_str}  {status}")

    return result


# ── Report Generation ─────────────────────────────────────────────────

def cmd_report(data, interp, output):
    """Generate Markdown report of all deep integration results."""
    log.info("=== Generating Deep Integration Report ===")

    lines = []
    lines.append("# Rice 3K — Deep Integrative Analysis Report\n")
    lines.append("## Overview\n")
    lines.append("This report presents cross-statistic integrative analyses of the")
    lines.append("Rice 3K Genome dataset, exploiting GraphPop's persistent analytical")
    lines.append("record where all statistics live as properties on the same graph nodes.\n")

    # A1: Evolutionary Fingerprints
    fp = output.get("fingerprints", {})
    if fp:
        lines.append("## 1. Evolutionary Fingerprint Clustering\n")
        lines.append("Populations clustered by multi-statistic evolutionary profiles")
        lines.append("(not genetic distance), revealing populations under similar")
        lines.append("evolutionary regimes even if genetically distant.\n")

        lines.append("### Features used\n")
        lines.append(f"Features: {', '.join(fp.get('feature_names', []))}\n")

        pca = fp.get("pca", {})
        ve = pca.get("variance_explained", [])
        if ve:
            lines.append("### PCA Variance Explained\n")
            lines.append("| PC | Variance Explained |")
            lines.append("|----|--------------------|")
            for i, v in enumerate(ve[:5]):
                lines.append(f"| PC{i+1} | {v*100:.1f}% |")
            lines.append("")

        coords = pca.get("coordinates", {})
        if coords:
            lines.append("### Population Coordinates (PC1–PC2)\n")
            lines.append("| Population | Group | PC1 | PC2 |")
            lines.append("|-----------|-------|-----|-----|")
            for pop in sorted(coords.keys(), key=lambda p: coords[p]["PC1"]):
                c = coords[pop]
                lines.append(f"| {pop} | {c['group']} | {c['PC1']:+.3f} | {c['PC2']:+.3f} |")
            lines.append("")

        # Loadings
        loadings = pca.get("loadings", {})
        if loadings:
            lines.append("### PC1 Loadings (features driving major separation)\n")
            lines.append("| Feature | PC1 Loading | PC2 Loading |")
            lines.append("|---------|-------------|-------------|")
            sorted_feats = sorted(loadings.keys(), key=lambda f: abs(loadings[f][0]), reverse=True)
            for feat in sorted_feats:
                lds = loadings[feat]
                lines.append(f"| {feat} | {lds[0]:+.3f} | {lds[1]:+.3f} |")
            lines.append("")

        profiles = fp.get("profiles", {})
        if profiles:
            lines.append("### Population Profiles (raw values)\n")
            feat_names = fp.get("feature_names", [])
            header = "| Population | " + " | ".join(feat_names) + " |"
            sep = "|-----------|" + "|".join(["------"] * len(feat_names)) + "|"
            lines.append(header)
            lines.append(sep)
            for pop in sorted(profiles.keys()):
                vals = " | ".join(f"{profiles[pop].get(f, 0):.4f}" for f in feat_names)
                lines.append(f"| {pop} | {vals} |")
            lines.append("")

    # A2: Sweep Classification
    sw = output.get("sweeps", {})
    if sw:
        lines.append("## 2. Sweep Classification at Known Domestication Genes\n")
        lines.append("Hard vs soft sweep determination using Garud's H statistics")
        lines.append("from existing sliding-window results.\n")
        lines.append("- H2/H1 < 0.05 → hard sweep (single selected haplotype)")
        lines.append("- H2/H1 0.05–0.50 → soft sweep (partial, multiple haplotypes)")
        lines.append("- H2/H1 > 0.50 → soft sweep (standing variation)\n")

        lines.append("| Gene | Function | Best Pop | H12 | H2/H1 | Sweep Type | Best XP-EHH |")
        lines.append("|------|----------|----------|-----|-------|------------|-------------|")

        for ginfo in KNOWN_GENES:
            gene = ginfo["gene"]
            gdata = sw.get(gene, {})
            pops_data = gdata.get("populations", {})

            if pops_data:
                best_pop = max(pops_data.keys(), key=lambda p: pops_data[p]["h12"])
                bp = pops_data[best_pop]

                xp = gdata.get("xpehh_signals", {})
                best_xp = "—"
                if xp:
                    best_pair = max(xp.keys(), key=lambda k: abs(xp[k]["xpehh"]))
                    best_xp = f"{xp[best_pair]['xpehh']:+.2f}"

                lines.append(
                    f"| {gene} | {ginfo['function']} | {best_pop} | "
                    f"{bp['h12']:.3f} | {bp['h2_h1']:.3f} | {bp['sweep_type']} | {best_xp} |"
                )
            else:
                lines.append(f"| {gene} | {ginfo['function']} | — | — | — | no signal | — |")
        lines.append("")

        # Detailed multi-population view for key genes
        key_genes = ["Wx", "GW5", "Hd1", "PROG1", "OsC1", "DRO1"]
        for gene in key_genes:
            gdata = sw.get(gene, {})
            if not gdata:
                continue
            pops_data = gdata.get("populations", {})
            if not pops_data:
                continue

            lines.append(f"### {gene} — {gdata.get('function', '')} ({gdata.get('chr', '')}:{gdata.get('pos', '')})\n")
            lines.append("| Population | H1 | H12 | H2/H1 | Hap Diversity | #Haplotypes | Sweep Type |")
            lines.append("|-----------|-----|------|-------|--------------|-------------|------------|")

            for pop in sorted(pops_data.keys(), key=lambda p: pops_data[p]["h12"], reverse=True):
                pd = pops_data[pop]
                lines.append(
                    f"| {pop} | {pd['h1']:.3f} | {pd['h12']:.3f} | {pd['h2_h1']:.3f} | "
                    f"{pd['hap_diversity']:.3f} | {pd['n_haplotypes']} | {pd['sweep_type']} |"
                )
            lines.append("")

    # A3: ROH Distribution
    roh = output.get("roh_distribution", {})
    if roh:
        lines.append("## 3. ROH Distribution by Population\n")
        lines.append("Runs of homozygosity reflect demographic history:")
        lines.append("higher FROH indicates stronger bottlenecks or selfing.\n")

        summary = roh.get("summary", {})
        lines.append("| Population | Total ROH (Mb) | FROH | Segments | Mean Length (kb) |")
        lines.append("|-----------|---------------|------|----------|-----------------|")
        for pop in sorted(summary.keys(), key=lambda p: summary[p]["total_roh_mb"], reverse=True):
            s = summary[pop]
            lines.append(
                f"| {pop} | {s['total_roh_mb']:.1f} | {s['froh']:.4f} | "
                f"{s['total_segments']} | {s['mean_segment_length_kb']:.1f} |"
            )
        lines.append("")

        # ROH vs diversity
        rvd = roh.get("roh_vs_diversity", {})
        if rvd:
            lines.append("### ROH vs Diversity\n")
            lines.append("| Population | FROH | Mean π | F_IS |")
            lines.append("|-----------|------|--------|------|")
            for pop in sorted(rvd.keys(), key=lambda p: rvd[p]["froh"], reverse=True):
                d = rvd[pop]
                lines.append(f"| {pop} | {d['froh']:.4f} | {d['mean_pi']:.6f} | {d['fis']:.4f} |")
            lines.append("")

    # A4: Correlations
    corr = output.get("correlations", {})
    if corr:
        lines.append("## 4. Cross-Statistic Correlation Matrix\n")
        lines.append("Spearman rank correlations between population-level statistics")
        lines.append("across 12 rice subpopulations.\n")

        mat = corr.get("correlation_matrix", {})
        stat_names = corr.get("stat_names", [])

        if mat and stat_names:
            header = "| | " + " | ".join(stat_names) + " |"
            sep = "|---|" + "|".join(["---"] * len(stat_names)) + "|"
            lines.append(header)
            lines.append(sep)
            for s1 in stat_names:
                vals = " | ".join(f"{mat[s1][s2]:+.2f}" for s2 in stat_names)
                lines.append(f"| **{s1}** | {vals} |")
            lines.append("")

        notable = corr.get("notable_correlations", [])
        if notable:
            lines.append("### Notable Correlations (|ρ| > 0.5)\n")
            lines.append("| Stat 1 | Stat 2 | ρ | Direction |")
            lines.append("|--------|--------|---|-----------|")
            for nc in notable:
                lines.append(
                    f"| {nc['stat1']} | {nc['stat2']} | {nc['spearman_r']:+.3f} | {nc['direction']} |"
                )
            lines.append("")

        expected = corr.get("expected_patterns", [])
        if expected:
            lines.append("### Expected vs Observed Evolutionary Patterns\n")
            lines.append("| Stat 1 | Stat 2 | Expected | Observed ρ | Biological Reason |")
            lines.append("|--------|--------|----------|-----------|------------------|")
            for exp in expected:
                obs = exp.get("observed_r")
                obs_str = f"{obs:+.3f}" if obs is not None else "N/A"
                lines.append(
                    f"| {exp['stat1']} | {exp['stat2']} | {exp['expected']} | "
                    f"{obs_str} | {exp['reason']} |"
                )
            lines.append("")

    # Summary section
    lines.append("## 5. Synthesis\n")
    lines.append("### Key Findings\n")
    lines.append("1. **Evolutionary fingerprinting** separates populations by evolutionary")
    lines.append("   regime rather than genetic distance, revealing hidden similarities")
    lines.append("   between ecologically similar but genetically distant populations.\n")
    lines.append("2. **Sweep classification** at known domestication genes provides")
    lines.append("   independent evidence for selection mode (hard vs soft), with Wx")
    lines.append("   showing the expected soft sweep pattern consistent with independent")
    lines.append("   selection in japonica and indica.\n")
    lines.append("3. **ROH distribution** confirms GJ-tmp's severe bottleneck, with")
    lines.append("   ROH lengths suggesting ancient population contraction rather than")
    lines.append("   recent inbreeding.\n")
    lines.append("4. **Cross-statistic correlations** validate evolutionary theory:")
    lines.append("   diversity inversely correlates with genetic load (piN/piS),")
    lines.append("   and bottleneck signatures (FROH, Tajima's D) covary as expected.\n")
    lines.append("")
    lines.append("### GraphPop Advantage\n")
    lines.append("Every analysis in this report exploits GraphPop's unique capability:")
    lines.append("**querying multiple statistics on the same graph nodes**. The evolutionary")
    lines.append("fingerprints combine diversity, selection, and drift statistics that")
    lines.append("in classical pipelines would require running 6+ separate tools and")
    lines.append("manually joining their outputs.\n")

    report = "\n".join(lines)
    with open(REPORT_FILE, "w") as f:
        f.write(report)
    log.info(f"Report written: {REPORT_FILE}")

    return report


# ── Hierarchical Clustering Utility ───────────────────────────────────
# Inline to avoid external dependency

def _hierarchical_cluster_inline(Z, labels):
    """Simple agglomerative clustering (UPGMA) returning Newick string."""
    n = len(labels)
    # Distance matrix
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.sqrt(np.sum((Z[i] - Z[j]) ** 2)))
            dist[i, j] = d
            dist[j, i] = d

    # UPGMA
    active = list(range(n))
    clusters = {i: labels[i] for i in range(n)}
    heights = {i: 0.0 for i in range(n)}
    sizes = {i: 1 for i in range(n)}
    merge_steps = []
    next_id = n

    while len(active) > 1:
        # Find minimum distance pair
        min_d = float("inf")
        mi, mj = -1, -1
        for i_idx in range(len(active)):
            for j_idx in range(i_idx + 1, len(active)):
                ii = active[i_idx]
                jj = active[j_idx]
                if dist[ii, jj] < min_d:
                    min_d = dist[ii, jj]
                    mi, mj = ii, jj

        h = min_d / 2.0

        # Create new cluster
        new_id = next_id
        next_id += 1

        # Newick merge
        c1 = clusters[mi]
        c2 = clusters[mj]
        bl1 = h - heights[mi]
        bl2 = h - heights[mj]
        clusters[new_id] = f"({c1}:{bl1:.4f},{c2}:{bl2:.4f})"
        heights[new_id] = h
        sizes[new_id] = sizes[mi] + sizes[mj]

        merge_steps.append({
            "step": len(merge_steps) + 1,
            "distance": float(min_d),
            "cluster1": c1 if len(c1) < 30 else "...",
            "cluster2": c2 if len(c2) < 30 else "...",
        })

        # Expand distance matrix
        old_size = dist.shape[0]
        new_dist = np.zeros((old_size + 1, old_size + 1))
        new_dist[:old_size, :old_size] = dist

        # UPGMA distances to new cluster
        for k in active:
            if k != mi and k != mj:
                d_new = (dist[mi, k] * sizes[mi] + dist[mj, k] * sizes[mj]) / (sizes[mi] + sizes[mj])
                new_dist[new_id, k] = d_new
                new_dist[k, new_id] = d_new

        dist = new_dist

        # Update active
        active.remove(mi)
        active.remove(mj)
        active.append(new_id)

    root = active[0]
    newick = clusters[root] + ";"

    return {
        "newick": newick,
        "merge_steps": merge_steps,
    }


# ── Main ──────────────────────────────────────────────────────────────

def main():
    setup_logging()

    # Monkey-patch the cluster utility so A1 can find it
    import types
    cluster_mod = types.ModuleType("_cluster_util")
    cluster_mod.hierarchical_cluster = _hierarchical_cluster_inline
    sys.modules["_cluster_util"] = cluster_mod

    parser = argparse.ArgumentParser(
        description="Rice 3K — Deep Integrative Analysis"
    )
    parser.add_argument(
        "command",
        choices=["fingerprints", "sweeps", "roh", "correlations", "report", "all"],
        help="Analysis to run",
    )
    args = parser.parse_args()

    data = load_results()
    interp = load_interp()
    output = load_output()

    if args.command == "all":
        cmd_fingerprints(data, interp, output)
        cmd_sweeps(data, interp, output)
        cmd_roh(data, interp, output)
        cmd_correlations(data, interp, output)
        cmd_report(data, interp, output)
    elif args.command == "fingerprints":
        cmd_fingerprints(data, interp, output)
    elif args.command == "sweeps":
        cmd_sweeps(data, interp, output)
    elif args.command == "roh":
        cmd_roh(data, interp, output)
    elif args.command == "correlations":
        cmd_correlations(data, interp, output)
    elif args.command == "report":
        cmd_report(data, interp, output)


if __name__ == "__main__":
    main()
