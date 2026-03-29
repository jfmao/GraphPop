#!/usr/bin/env python3
"""Rice Inv.R11: Multi-Evidence Synthesis.

Integrates results from all rice investigations (R01-R10) into a unified
population x evidence heatmap with indica vs japonica comparison.

Data source: data/results/rice/rice_inv*.json (all investigation outputs)
Output:      data/results/rice/rice_inv11_synthesis.{tsv,json}
Figure:      data/results/rice/figures/rice_fig11_synthesis.png
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "data" / "results" / "rice"
FIG_DIR = RESULTS_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

GROUP_COLORS = {
    "Indica": "#e41a1c",
    "Japonica": "#377eb8",
    "Aus/Basmati": "#4daf4a",
    "Admixed": "#984ea3",
}

POP_GROUPS = {
    "XI-1A": "Indica", "XI-1B": "Indica", "XI-2": "Indica",
    "XI-3": "Indica", "XI-adm": "Indica",
    "GJ-tmp": "Japonica", "GJ-trp": "Japonica",
    "GJ-sbtrp": "Japonica", "GJ-adm": "Japonica",
    "cA-Aus": "Aus/Basmati", "cB-Bas": "Aus/Basmati",
    "admix": "Admixed",
}

ALL_POPS = sorted(POP_GROUPS.keys())


def load_json_safe(path):
    """Load a JSON file, returning empty dict if not found."""
    try:
        with open(path) as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"  WARNING: {path.name} not found, skipping")
        return {}


def main():
    # Load all investigation results
    inv01 = load_json_safe(RESULTS_DIR / "rice_inv01_evo_fingerprint.json")
    inv02 = load_json_safe(RESULTS_DIR / "rice_inv02_pinsps.json")
    inv03 = load_json_safe(RESULTS_DIR / "rice_inv03_fst_landscape.json")
    inv04 = load_json_safe(RESULTS_DIR / "rice_inv04_conditioned_fst.json")
    inv05 = load_json_safe(RESULTS_DIR / "rice_inv05_domestication_genes.json")
    inv06 = load_json_safe(RESULTS_DIR / "rice_inv06_roh_landscape.json")
    inv07 = load_json_safe(RESULTS_DIR / "rice_inv07_pbs_sweeps.json")
    inv08 = load_json_safe(RESULTS_DIR / "rice_inv08_sweep_types.json")
    inv09 = load_json_safe(RESULTS_DIR / "rice_inv09_stat_correlations.json")
    inv10 = load_json_safe(RESULTS_DIR / "rice_inv10_fay_wu_daf.json")

    # Also load deep integration for per-pop stats
    di_path = ROOT / "results" / "rice" / "rice_deep_integration_results.json"
    with open(di_path) as f:
        di = json.load(f)

    pop_stats = di["correlations"]["population_stats"]
    stat_names = di["correlations"]["stat_names"]

    # Build evidence matrix: population x evidence dimensions
    # Evidence dimensions (z-scored across populations):
    # 1. pi (diversity)
    # 2. theta_W
    # 3. Tajima's D
    # 4. Fay & Wu H
    # 5. FIS (inbreeding)
    # 6. FROH
    # 7. N sweeps
    # 8. piN/piS (genetic load)
    evidence_names = stat_names  # ['pi', 'theta_w', 'tajima_d', 'fay_wu_h', 'fis', 'froh', 'n_sweeps', 'pinsps']

    # Build matrix for populations present in pop_stats
    available_pops = [p for p in ALL_POPS if p in pop_stats]
    n_pops = len(available_pops)
    n_evidence = len(evidence_names)

    raw_matrix = np.zeros((n_pops, n_evidence))
    for i, pop in enumerate(available_pops):
        for j, stat in enumerate(evidence_names):
            raw_matrix[i, j] = pop_stats[pop].get(stat, float("nan"))

    # Z-score across populations
    means = np.nanmean(raw_matrix, axis=0)
    stds = np.nanstd(raw_matrix, axis=0)
    stds[stds == 0] = 1
    z_matrix = (raw_matrix - means) / stds

    # Save TSV
    header = ["population", "group"] + [f"{s}_raw" for s in evidence_names] + \
             [f"{s}_z" for s in evidence_names]
    rows = []
    for i, pop in enumerate(available_pops):
        row = {"population": pop, "group": POP_GROUPS.get(pop, "Unknown")}
        for j, stat in enumerate(evidence_names):
            row[f"{stat}_raw"] = raw_matrix[i, j]
            row[f"{stat}_z"] = z_matrix[i, j]
        rows.append(row)

    with open(RESULTS_DIR / "rice_inv11_synthesis.tsv", "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r[h]) for h in header) + "\n")

    # Indica vs Japonica comparison
    indica_pops = [p for p in available_pops if POP_GROUPS.get(p) == "Indica"]
    japonica_pops = [p for p in available_pops if POP_GROUPS.get(p) == "Japonica"]
    indica_idx = [available_pops.index(p) for p in indica_pops]
    japonica_idx = [available_pops.index(p) for p in japonica_pops]

    indica_means = np.nanmean(raw_matrix[indica_idx], axis=0) if indica_idx else np.zeros(n_evidence)
    japonica_means = np.nanmean(raw_matrix[japonica_idx], axis=0) if japonica_idx else np.zeros(n_evidence)

    comparison = {}
    for j, stat in enumerate(evidence_names):
        comparison[stat] = {
            "indica_mean": float(indica_means[j]),
            "japonica_mean": float(japonica_means[j]),
            "difference": float(indica_means[j] - japonica_means[j]),
        }

    # Count investigations completed
    inv_status = {
        "R01_evo_fingerprint": bool(inv01),
        "R02_pinsps": bool(inv02),
        "R03_fst_landscape": bool(inv03),
        "R04_conditioned_fst": bool(inv04),
        "R05_domestication_genes": bool(inv05),
        "R06_roh_landscape": bool(inv06),
        "R07_pbs_sweeps": bool(inv07),
        "R08_sweep_types": bool(inv08),
        "R09_stat_correlations": bool(inv09),
        "R10_fay_wu_daf": bool(inv10),
    }

    # Key findings from each investigation
    key_findings = []
    if inv01:
        key_findings.append(f"R01: PCA separates indica-japonica axis "
                            f"({inv01.get('n_features', '?')} features)")
    if inv02:
        key_findings.append(f"R02: piN/piS analysis complete")
    if inv03:
        key_findings.append(f"R03: FST landscape mapped")
    if inv05:
        n_genes = inv05.get("n_genes", 0)
        n_swept = inv05.get("genes_with_sweeps", 0)
        key_findings.append(f"R05: {n_swept}/{n_genes} domestication genes show sweeps")
    if inv07:
        n_comp = inv07.get("n_comparisons", 0)
        key_findings.append(f"R07: {n_comp} PBS comparisons analyzed")
    if inv08:
        tc = inv08.get("type_counts", {})
        key_findings.append(f"R08: {tc.get('hard', 0)} hard, {tc.get('soft', 0)} soft sweeps")
    if inv09:
        sp = inv09.get("strongest_positive", {})
        sn = inv09.get("strongest_negative", {})
        key_findings.append(
            f"R09: Strongest correlation: {sp.get('stat1', '')} x {sp.get('stat2', '')} "
            f"(r={sp.get('spearman_r', 0):.2f})")
    if inv10:
        key_findings.append(
            f"R10: PBS DAF enrichment={inv10.get('mean_daf_enrichment_pbs', 0):.3f}, "
            f"XP-EHH={inv10.get('mean_daf_enrichment_xpehh', 0):.3f}")

    summary = {
        "n_populations": n_pops,
        "n_evidence_dimensions": n_evidence,
        "evidence_names": evidence_names,
        "indica_vs_japonica": comparison,
        "investigations_completed": sum(inv_status.values()),
        "investigation_status": inv_status,
        "key_findings": key_findings,
    }
    with open(RESULTS_DIR / "rice_inv11_synthesis.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ── Figure: 2 panels ──────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Friendly labels for display
    display_names = {
        "pi": "pi", "theta_w": "theta_W", "tajima_d": "Tajima's D",
        "fay_wu_h": "Fay&Wu H", "fis": "FIS", "froh": "FROH",
        "n_sweeps": "N sweeps", "pinsps": "piN/piS",
    }

    # Panel (a): Population x Evidence z-score heatmap
    ax = axes[0]
    # Sort populations by group for visual clustering
    sort_order = []
    for grp in ["Indica", "Japonica", "Aus/Basmati", "Admixed"]:
        sort_order.extend([i for i, p in enumerate(available_pops)
                           if POP_GROUPS.get(p) == grp])
    sorted_pops = [available_pops[i] for i in sort_order]
    sorted_z = z_matrix[sort_order]

    im = ax.imshow(sorted_z, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(n_evidence))
    ax.set_xticklabels([display_names.get(s, s) for s in evidence_names],
                       rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(n_pops))
    # Color y-tick labels by group
    ylabels = ax.set_yticklabels(sorted_pops, fontsize=7)
    for lbl, pop in zip(ylabels, sorted_pops):
        grp = POP_GROUPS.get(pop, "Unknown")
        lbl.set_color(GROUP_COLORS.get(grp, "black"))
    plt.colorbar(im, ax=ax, label="Z-score", shrink=0.8)
    ax.set_title("a  Population x evidence z-scores", fontweight="bold", loc="left")

    # Panel (b): Indica vs Japonica comparison (grouped bar)
    ax = axes[1]
    x = np.arange(n_evidence)
    width = 0.35
    ax.bar(x - width / 2, indica_means, width, label="Indica", color="#e41a1c",
           edgecolor="black", linewidth=0.5, alpha=0.8)
    ax.bar(x + width / 2, japonica_means, width, label="Japonica", color="#377eb8",
           edgecolor="black", linewidth=0.5, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([display_names.get(s, s) for s in evidence_names],
                       rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Mean statistic value")
    ax.axhline(0, color="black", linewidth=0.5)
    ax.legend(fontsize=8)
    ax.set_title("b  Indica vs Japonica mean statistics", fontweight="bold", loc="left")

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig11_synthesis.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Print summary
    print("=== Rice Inv.R11: Multi-Evidence Synthesis ===")
    print(f"Populations: {n_pops}")
    print(f"Evidence dimensions: {n_evidence}")
    print(f"Investigations completed: {sum(inv_status.values())}/{len(inv_status)}")
    print(f"\nIndica vs Japonica comparison:")
    for stat in evidence_names:
        c = comparison[stat]
        direction = ">" if c["difference"] > 0 else "<"
        print(f"  {display_names.get(stat, stat):12s}: "
              f"indica={c['indica_mean']:.4f} {direction} japonica={c['japonica_mean']:.4f} "
              f"(diff={c['difference']:+.4f})")
    print(f"\nKey findings:")
    for kf in key_findings:
        print(f"  {kf}")
    print(f"\nSaved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
