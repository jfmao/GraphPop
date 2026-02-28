#!/usr/bin/env python3
"""Validate Tajima's D: scikit-allel vs GraphPop.

Reads allele counts from the CSV files (same data Neo4j has) and
computes Tajima's D using scikit-allel's implementation, then compares
against the Neo4j graphpop.diversity procedure.

Usage:
    conda run -n graphevo python scripts/validate-tajimad.py
"""

import csv
import subprocess
import sys
from pathlib import Path

import numpy as np
import allel

CSV_DIR = Path("data/raw/1000g/csv_out")
CHR = "chr22"
START = 16_000_000
END = 17_000_000
POP = "AFR"

NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"


def load_ac_from_csv(csv_dir, pop, chr_filter, start, end):
    """Load allele counts for a population from variant_nodes.csv."""
    variant_csv = csv_dir / "variant_nodes.csv"
    pop_index = None
    ac_list = []
    an_list = []

    with open(variant_csv) as f:
        reader = csv.reader(f)
        next(reader)  # skip header

        for row in reader:
            chr_val = row[2]
            pos = int(row[3])
            if chr_val != chr_filter or pos < start or pos > end:
                continue

            pop_ids = row[7].split(";")
            if pop_index is None:
                pop_index = pop_ids.index(pop)

            ac_arr = [int(x) for x in row[8].split(";")]
            an_arr = [int(x) for x in row[9].split(";")]

            ac = ac_arr[pop_index]
            an = an_arr[pop_index]
            if an < 2:
                continue

            ac_list.append(ac)
            an_list.append(an)

    return np.array(ac_list), np.array(an_list)


def tajima_d_scikit_allel(ac_alt, an):
    """Compute Tajima's D using scikit-allel."""
    # scikit-allel expects a 2D allele counts array: shape (n_variants, 2)
    # Column 0 = ref count, Column 1 = alt count
    ac_ref = an - ac_alt
    ac_2d = np.column_stack([ac_ref, ac_alt])
    ac_obj = allel.AlleleCountsArray(ac_2d)

    # Tajima's D
    D = allel.tajima_d(ac_obj)
    return D


def query_neo4j_tajimad():
    """Get Tajima's D from Neo4j."""
    cmd = [
        "cypher-shell",
        "-u", NEO4J_USER, "-p", NEO4J_PASS,
        "--format", "plain",
        f"CALL graphpop.diversity('{CHR}', {START}, {END}, '{POP}') "
        "YIELD tajima_d RETURN tajima_d;",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    lines = result.stdout.strip().split("\n")
    return float(lines[1].strip())


def main():
    print("=" * 72)
    print("Tajima's D Validation: scikit-allel vs GraphPop")
    print(f"Region: {CHR}:{START}-{END}, Population: {POP}")
    print("=" * 72)
    print()

    # Load data
    print("Loading allele counts from CSV...")
    ac_alt, an = load_ac_from_csv(CSV_DIR, POP, CHR, START, END)
    print(f"  {len(ac_alt)} variants loaded")
    print()

    # scikit-allel
    print("Computing Tajima's D with scikit-allel...")
    D_allel = tajima_d_scikit_allel(ac_alt, an)
    print(f"  scikit-allel Tajima's D = {D_allel:.10f}")
    print()

    # Neo4j
    print("Querying graphpop.diversity for Tajima's D...")
    D_neo4j = query_neo4j_tajimad()
    print(f"  Neo4j Tajima's D       = {D_neo4j:.10f}")
    print()

    # Compare
    rel_err = abs(D_allel - D_neo4j) / abs(D_allel) if D_allel != 0 else abs(D_neo4j)
    pct = rel_err * 100

    print("-" * 72)
    status = "PASS" if pct < 0.01 else ("WARN" if pct < 1.0 else "FAIL")
    print(f"  {status}  Tajima's D  allel={D_allel:.10f}  neo4j={D_neo4j:.10f}  RelErr={pct:.6f}%")

    if pct >= 0.01:
        print()
        print("  Note: small discrepancies are expected due to:")
        print("   - scikit-allel uses a single n for the whole region (modal AN)")
        print("   - GraphPop uses per-population n_samples from the Population node")
        print("   - Both are valid estimators of Tajima's D")

    print()
    print("=" * 72)

    return 0 if pct < 1.0 else 1


if __name__ == "__main__":
    sys.exit(main())
