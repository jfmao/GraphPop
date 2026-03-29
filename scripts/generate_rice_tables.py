#!/usr/bin/env python3
"""Generate rice supplementary tables for the GraphPop paper.

Reads rice investigation JSON/TSV results and creates properly formatted
supplementary table TSV files in paper/supplementary_tables/.

Tables produced:
  ST14: Rice population summary (12 pops)
  ST15: Rice piN/piS
  ST16: Rice Fst matrix (12x12)
  ST17: Rice conditioned Fst (missense vs synonymous)
  ST18: Rice domestication genes (16 genes)
  ST19: Rice consequence bias (HIGH/MODERATE/LOW)
  ST20: Rice gene-level Fst top 50
  ST21: Rice rare variant burden per population
  ST22: Rice computational timing
"""

import json
import csv
import sys
from pathlib import Path

ROOT = Path("/mnt/data/GraphPop")
RESULTS_RICE = ROOT / "data" / "results" / "rice"
INTERP_RICE = ROOT / "results" / "rice"
OUT = ROOT / "paper" / "supplementary_tables"
OUT.mkdir(parents=True, exist_ok=True)


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    """Write a TSV file with header and rows."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
    print(f"  wrote {path.name} ({len(rows)} rows)")


def fmt(v, decimals=6):
    """Format a numeric value."""
    if v is None:
        return "NA"
    if isinstance(v, float):
        return f"{v:.{decimals}f}"
    return str(v)


# ── ST14: Population summary ──────────────────────────────────────────────
def make_st14():
    deep = json.loads((INTERP_RICE / "rice_deep_integration_results.json").read_text())
    profiles = deep["fingerprints"]["profiles"]
    header = ["population", "group", "pi", "theta_w", "tajima_d", "fis",
              "froh", "n_sweeps", "pinsps"]
    rows = []
    # Get group info from the fingerprint TSV
    pop_group = {}
    with open(RESULTS_RICE / "rice_inv01_evo_fingerprint.tsv") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            pop_group[r["population"]] = r["group"]

    for pop in sorted(profiles.keys()):
        p = profiles[pop]
        group = pop_group.get(pop, "Unknown")
        rows.append([
            pop, group,
            fmt(p.get("pi")), fmt(p.get("theta_w")), fmt(p.get("tajima_d")),
            fmt(p.get("fis")), fmt(p.get("froh")),
            str(int(p.get("n_sweeps", 0))),
            fmt(p.get("pinsps"), 4),
        ])
    write_tsv(OUT / "SupplementaryTable14_rice_population_summary.tsv", header, rows)


# ── ST15: piN/piS ─────────────────────────────────────────────────────────
def make_st15():
    # Direct copy from the investigation TSV (already well-formatted)
    src = RESULTS_RICE / "rice_inv02_pinsps.tsv"
    with open(src) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        rows = [r for r in reader]
    write_tsv(OUT / "SupplementaryTable15_rice_pinsps.tsv", header, rows)


# ── ST16: Fst matrix ──────────────────────────────────────────────────────
def make_st16():
    interp = json.loads((INTERP_RICE / "rice_interpretation_results.json").read_text())
    fst = interp["fst_matrix"]
    # Extract unique populations
    pops = set()
    for key in fst:
        p1, p2 = key.split("_vs_")
        pops.add(p1)
        pops.add(p2)
    pops = sorted(pops)

    # Build matrix
    matrix = {}
    for key, val in fst.items():
        p1, p2 = key.split("_vs_")
        matrix[(p1, p2)] = val["mean_fst"]
        matrix[(p2, p1)] = val["mean_fst"]

    header = ["population"] + pops
    rows = []
    for p1 in pops:
        row = [p1]
        for p2 in pops:
            if p1 == p2:
                row.append("0.000000")
            else:
                v = matrix.get((p1, p2))
                row.append(fmt(v) if v is not None else "NA")
        rows.append(row)
    write_tsv(OUT / "SupplementaryTable16_rice_fst_matrix.tsv", header, rows)


# ── ST17: Conditioned Fst ─────────────────────────────────────────────────
def make_st17():
    data = json.loads((RESULTS_RICE / "rice_inv04_conditioned_fst.json").read_text())
    pairs = data["pairs"]
    header = ["pair", "n_windows_missense", "n_windows_synonymous",
              "mean_fst_missense", "mean_fst_synonymous", "ratio"]
    rows = []
    for pair_name in sorted(pairs.keys()):
        p = pairs[pair_name]
        rows.append([
            pair_name,
            str(p["n_windows_missense"]),
            str(p["n_windows_synonymous"]),
            fmt(p["mean_fst_missense"]),
            fmt(p["mean_fst_synonymous"]),
            fmt(p["ratio"], 4),
        ])
    write_tsv(OUT / "SupplementaryTable17_rice_conditioned_fst.tsv", header, rows)


# ── ST18: Domestication genes ─────────────────────────────────────────────
def make_st18():
    data = json.loads((RESULTS_RICE / "rice_inv05_domestication_genes.json").read_text())
    genes = data["genes"]
    header = ["gene", "function", "chromosome", "hard_sweeps", "soft_sweeps",
              "total_sweeps"]
    rows = []
    for gene in sorted(genes.keys()):
        g = genes[gene]
        rows.append([
            gene, g["function"], g["chr"],
            str(g["hard"]), str(g["soft"]), str(g["total"]),
        ])
    write_tsv(OUT / "SupplementaryTable18_rice_domestication_genes.tsv", header, rows)


# ── ST19: Consequence bias ────────────────────────────────────────────────
def make_st19():
    data = json.loads((RESULTS_RICE / "rice_inv14_consequence_bias.json").read_text())
    impacts = data["impacts"]
    mw = data.get("mann_whitney", {})
    pair = data.get("pair", "GJ-tmp_vs_XI-1A")

    header = ["pair", "impact", "n_variants", "mean_fst", "median_fst", "std_fst"]
    rows = []
    for impact in ["HIGH", "MODERATE", "LOW"]:
        if impact in impacts:
            imp = impacts[impact]
            rows.append([
                pair, impact,
                str(imp["n_variants"]),
                fmt(imp["mean_fst"]),
                fmt(imp["median_fst"]),
                fmt(imp["std_fst"]),
            ])

    # Add Mann-Whitney results as additional rows
    if mw:
        header_mw = ["comparison", "U_statistic", "p_value"]
        rows_mw = []
        for comp, vals in sorted(mw.items()):
            rows_mw.append([comp, fmt(vals["U"], 1), f"{vals['p']:.2e}"])

        # Write both sections: impacts then MannWhitney
        with open(OUT / "SupplementaryTable19_rice_consequence_bias.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["# FST by functional impact category"])
            w.writerow(header)
            w.writerows(rows)
            w.writerow([])
            w.writerow(["# Mann-Whitney U tests"])
            w.writerow(header_mw)
            w.writerows(rows_mw)
        print(f"  wrote SupplementaryTable19_rice_consequence_bias.tsv ({len(rows)} impact rows, {len(rows_mw)} test rows)")
        return

    write_tsv(OUT / "SupplementaryTable19_rice_consequence_bias.tsv", header, rows)


# ── ST20: Gene Fst top 50 ────────────────────────────────────────────────
def make_st20():
    src = RESULTS_RICE / "rice_inv13_gene_fst.tsv"
    with open(src) as f:
        reader = csv.DictReader(f, delimiter="\t")
        all_rows = list(reader)

    # Sort by fst descending
    all_rows.sort(key=lambda r: float(r.get("fst", 0)), reverse=True)
    top50 = all_rows[:50]

    header = ["rank", "pair", "geneId", "symbol", "n_variants", "n_used", "fst"]
    rows = []
    for i, r in enumerate(top50, 1):
        rows.append([
            str(i), r["pair"], r["geneId"], r["symbol"],
            r["n_variants"], r["n_used"], r["fst"],
        ])
    write_tsv(OUT / "SupplementaryTable20_rice_gene_fst_top50.tsv", header, rows)


# ── ST21: Rare burden ─────────────────────────────────────────────────────
def make_st21():
    src = RESULTS_RICE / "rice_inv15_rare_burden.tsv"
    with open(src) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        rows = [r for r in reader]
    write_tsv(OUT / "SupplementaryTable21_rice_rare_burden.tsv", header, rows)


# ── ST22: Computational timing ────────────────────────────────────────────
def make_st22():
    header = ["stage", "detail", "time_minutes", "items"]
    rows = [
        ["CSV generation", "29,635,224 variants, 8 workers", "274.5", "29635224"],
        ["Annotation loading", "55,328 genes + 7,004,802 HAS_CONSEQUENCE edges", "4.4", "7060130"],
        ["Phase 1 (summary stats)", "12 pops x 12 chr x 7 procedures", "48.0", "1008"],
        ["Phase 2 (XP-EHH)", "6 pairs x 12 chr", "78.0", "72"],
        ["Interpretation pipeline", "All automated analyses", "60.0", "—"],
        ["Pass 2 investigations", "R13 gene FST + R14 consequence + R15 rare burden", "70.0", "3"],
        ["  R13 gene-level FST", "Per-gene Fst for all pairs", "15.0", "—"],
        ["  R14 consequence bias", "HIGH/MODERATE/LOW FST stratification", "0.3", "—"],
        ["  R15 rare burden", "Population-specific rare variant counts", "52.0", "—"],
        ["TOTAL", "End-to-end pipeline", "534.9", "—"],
    ]
    write_tsv(OUT / "SupplementaryTable22_rice_computational_timing.tsv", header, rows)


# ── Main ──────────────────────────────────────────────────────────────────
def main():
    print("Generating rice supplementary tables...")
    make_st14()
    make_st15()
    make_st16()
    make_st17()
    make_st18()
    make_st19()
    make_st20()
    make_st21()
    make_st22()
    print("\nDone. All tables written to paper/supplementary_tables/")

    # Verify all files exist
    expected = [f"SupplementaryTable{i}_*.tsv" for i in range(14, 23)]
    for i in range(14, 23):
        matches = list(OUT.glob(f"SupplementaryTable{i}_*.tsv"))
        if matches:
            for m in matches:
                size = m.stat().st_size
                print(f"  OK  {m.name} ({size:,} bytes)")
        else:
            print(f"  MISSING  SupplementaryTable{i}_*.tsv")
            sys.exit(1)


if __name__ == "__main__":
    main()
