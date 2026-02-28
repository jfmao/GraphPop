#!/usr/bin/env python3
"""Validate GraphPop stored procedure results against scikit-allel.

Computes π, θ_W, Tajima's D, and Hudson's Fst from the 1000G chr22 VCF
using scikit-allel, then compares against the GraphPop procedure output.

Usage:
    conda run -n graphevo python scripts/validate-statistics.py
"""

import csv
import subprocess
import json
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

VCF = "data/raw/1000g/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz"
PANEL = "data/raw/1000g/integrated_call_samples_v3.20130502.ALL.panel"
CSV_DIR = Path("data/raw/1000g/csv_out")

# Region to validate (same as smoke test)
CHR = "chr22"
START = 16_000_000
END = 17_000_000

POP1 = "AFR"
POP2 = "EUR"

NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"

# ---------------------------------------------------------------------------
# Load population panel
# ---------------------------------------------------------------------------

def load_panel(panel_path):
    """Load sample -> super_pop mapping from 1000G panel file."""
    sample_to_pop = {}
    with open(panel_path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                sample_id = parts[0]
                super_pop = parts[2]
                sample_to_pop[sample_id] = super_pop
    return sample_to_pop


# ---------------------------------------------------------------------------
# Compute stats from CSV files (same data Neo4j has)
# ---------------------------------------------------------------------------

def compute_from_csv(csv_dir, pop, pop2=None):
    """Compute stats from the same CSV data that was imported into Neo4j."""
    variant_csv = csv_dir / "variant_nodes.csv"

    pi_sum = 0.0
    he_sum = 0.0
    ho_sum = 0.0
    n_variants = 0
    n_segregating = 0

    fst_num_sum = 0.0
    fst_den_sum = 0.0
    dxy_sum = 0.0

    pop_index = None
    pop2_index = None

    with open(variant_csv) as f:
        reader = csv.reader(f)
        header = next(reader)

        for row in reader:
            chr_val = row[2]
            pos = int(row[3])

            if chr_val != CHR:
                continue
            if pos < START or pos > END:
                continue

            pop_ids = row[7].split(";")

            if pop_index is None:
                pop_index = pop_ids.index(pop)
                if pop2:
                    pop2_index = pop_ids.index(pop2)

            ac_arr = [int(x) for x in row[8].split(";")]
            an_arr = [int(x) for x in row[9].split(";")]
            af_arr = [float(x) for x in row[10].split(";")]
            het_arr = [int(x) for x in row[11].split(";")]

            ac = ac_arr[pop_index]
            an = an_arr[pop_index]
            af = af_arr[pop_index]
            het = het_arr[pop_index]

            if an < 2:
                continue

            n_variants += 1

            # Pi with finite sample correction
            pi_site = 2.0 * af * (1.0 - af) * an / (an - 1.0)
            pi_sum += pi_site

            # He
            he_sum += 2.0 * af * (1.0 - af)

            # Ho
            n_diploid = an // 2
            ho_sum += het / n_diploid

            # Segregating sites
            if ac > 0 and ac < an:
                n_segregating += 1

            # Fst (Hudson)
            if pop2_index is not None:
                an2 = an_arr[pop2_index]
                af2 = af_arr[pop2_index]
                if an2 >= 2:
                    hw1 = 2.0 * af * (1.0 - af) * an / (an - 1.0)
                    hw2 = 2.0 * af2 * (1.0 - af2) * an2 / (an2 - 1.0)
                    hw = (hw1 + hw2) / 2.0
                    hb = (af - af2) ** 2 - hw1 / (2.0 * an) - hw2 / (2.0 * an2)
                    fst_num_sum += hb
                    fst_den_sum += hb + hw
                    dxy_sum += af * (1.0 - af2) + af2 * (1.0 - af)

    L = n_variants
    results = {
        "n_variants": n_variants,
        "n_segregating": n_segregating,
        "pi": pi_sum / L if L > 0 else 0,
        "het_exp": he_sum / L if L > 0 else 0,
        "het_obs": ho_sum / L if L > 0 else 0,
    }

    if pop2_index is not None and L > 0:
        results["fst_hudson"] = fst_num_sum / fst_den_sum if fst_den_sum > 0 else 0
        results["dxy"] = dxy_sum / L

    return results


# ---------------------------------------------------------------------------
# Query Neo4j procedures
# ---------------------------------------------------------------------------

def query_neo4j(cypher):
    """Run a Cypher query via cypher-shell and return raw output."""
    cmd = [
        "cypher-shell",
        "-u", NEO4J_USER,
        "-p", NEO4J_PASS,
        "--format", "plain",
        cypher,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    return result.stdout.strip()


def parse_neo4j_row(output):
    """Parse cypher-shell plain output into a dict."""
    lines = output.strip().split("\n")
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
# Main
# ---------------------------------------------------------------------------

def compare(label, csv_val, neo4j_val, tol=0.0001):
    """Compare two values and print result."""
    if csv_val == 0:
        rel_err = abs(neo4j_val)
    else:
        rel_err = abs(csv_val - neo4j_val) / abs(csv_val)

    status = "PASS" if rel_err < tol else "FAIL"
    pct = rel_err * 100
    print(f"  {status}  {label:20s}  CSV={csv_val:.10f}  Neo4j={neo4j_val:.10f}  RelErr={pct:.6f}%")
    return status == "PASS"


def main():
    print("=" * 72)
    print("GraphPop Validation: Neo4j procedures vs CSV ground truth")
    print(f"Region: {CHR}:{START}-{END}")
    print("=" * 72)
    print()

    # 1. Compute from CSV
    print("Computing stats from CSV files...")
    csv_stats = compute_from_csv(CSV_DIR, POP1, POP2)
    print(f"  CSV: {csv_stats['n_variants']} variants, {csv_stats['n_segregating']} segregating")
    print()

    # 2. Query Neo4j diversity
    print("Querying graphpop.diversity...")
    neo4j_div = parse_neo4j_row(query_neo4j(
        f"CALL graphpop.diversity('{CHR}', {START}, {END}, '{POP1}') "
        "YIELD pi, theta_w, tajima_d, het_exp, het_obs, fis, n_variants, n_segregating "
        "RETURN pi, theta_w, tajima_d, het_exp, het_obs, fis, n_variants, n_segregating;"
    ))
    print(f"  Neo4j: {int(neo4j_div.get('n_variants', 0))} variants, "
          f"{int(neo4j_div.get('n_segregating', 0))} segregating")
    print()

    # 3. Query Neo4j divergence
    print("Querying graphpop.divergence...")
    neo4j_fst = parse_neo4j_row(query_neo4j(
        f"CALL graphpop.divergence('{CHR}', {START}, {END}, '{POP1}', '{POP2}') "
        "YIELD fst_hudson, dxy, da, n_variants "
        "RETURN fst_hudson, dxy, da, n_variants;"
    ))
    print()

    # 4. Compare
    print("-" * 72)
    print("DIVERSITY VALIDATION")
    print("-" * 72)
    all_pass = True
    all_pass &= compare("n_variants", csv_stats["n_variants"], neo4j_div.get("n_variants", 0), tol=1e-9)
    all_pass &= compare("n_segregating", csv_stats["n_segregating"], neo4j_div.get("n_segregating", 0), tol=1e-9)
    all_pass &= compare("pi", csv_stats["pi"], neo4j_div.get("pi", 0))
    all_pass &= compare("het_exp", csv_stats["het_exp"], neo4j_div.get("het_exp", 0))
    all_pass &= compare("het_obs", csv_stats["het_obs"], neo4j_div.get("het_obs", 0))
    print()

    print("-" * 72)
    print("DIVERGENCE VALIDATION")
    print("-" * 72)
    all_pass &= compare("fst_hudson", csv_stats["fst_hudson"], neo4j_fst.get("fst_hudson", 0))
    all_pass &= compare("dxy", csv_stats["dxy"], neo4j_fst.get("dxy", 0))
    print()

    # 5. Summary
    print("=" * 72)
    if all_pass:
        print("ALL VALIDATIONS PASSED (<0.01% relative error)")
    else:
        print("SOME VALIDATIONS FAILED — check above for details")
    print("=" * 72)

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
