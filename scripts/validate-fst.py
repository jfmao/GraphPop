#!/usr/bin/env python3
"""Validate Hudson's Fst: scikit-allel vs GraphPop.

Usage:
    conda run -n graphevo python scripts/validate-fst.py
"""

import csv
import subprocess
import sys
from pathlib import Path

import numpy as np
import allel
import os

CSV_DIR = Path("data/raw/1000g/csv_out")
CHR = "chr22"
START = 16_000_000
END = 17_000_000
POP1 = "AFR"
POP2 = "EUR"

NEO4J_USER = os.environ.get("GRAPHPOP_USER", "neo4j")
NEO4J_PASS = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")


def load_ac_pair(csv_dir, pop1, pop2, chr_filter, start, end):
    """Load allele count arrays for two populations."""
    variant_csv = csv_dir / "variant_nodes.csv"
    p1_idx = None
    p2_idx = None
    ac1_list, an1_list = [], []
    ac2_list, an2_list = [], []

    with open(variant_csv) as f:
        reader = csv.reader(f)
        next(reader)

        for row in reader:
            chr_val = row[2]
            pos = int(row[3])
            if chr_val != chr_filter or pos < start or pos > end:
                continue

            pop_ids = row[7].split(";")
            if p1_idx is None:
                p1_idx = pop_ids.index(pop1)
                p2_idx = pop_ids.index(pop2)

            ac = [int(x) for x in row[8].split(";")]
            an = [int(x) for x in row[9].split(";")]

            a1, n1 = ac[p1_idx], an[p1_idx]
            a2, n2 = ac[p2_idx], an[p2_idx]

            if n1 < 2 or n2 < 2:
                continue

            ac1_list.append(a1)
            an1_list.append(n1)
            ac2_list.append(a2)
            an2_list.append(n2)

    return (np.array(ac1_list), np.array(an1_list),
            np.array(ac2_list), np.array(an2_list))


def hudson_fst_scikit_allel(ac1_alt, an1, ac2_alt, an2):
    """Compute Hudson's Fst using scikit-allel."""
    ac1_ref = an1 - ac1_alt
    ac2_ref = an2 - ac2_alt

    ac1_2d = np.column_stack([ac1_ref, ac1_alt])
    ac2_2d = np.column_stack([ac2_ref, ac2_alt])

    num, den = allel.hudson_fst(
        allel.AlleleCountsArray(ac1_2d),
        allel.AlleleCountsArray(ac2_2d),
    )
    fst = np.sum(num) / np.sum(den)
    return fst


def query_neo4j_fst():
    """Get Fst from Neo4j."""
    cmd = [
        "cypher-shell",
        "-u", NEO4J_USER, "-p", NEO4J_PASS,
        "--format", "plain",
        f"CALL graphpop.divergence('{CHR}', {START}, {END}, '{POP1}', '{POP2}') "
        "YIELD fst_hudson RETURN fst_hudson;",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    lines = result.stdout.strip().split("\n")
    return float(lines[1].strip())


def main():
    print("=" * 72)
    print(f"Hudson's Fst Validation: scikit-allel vs GraphPop")
    print(f"Region: {CHR}:{START}-{END}, Pops: {POP1} vs {POP2}")
    print("=" * 72)
    print()

    print("Loading allele counts from CSV...")
    ac1, an1, ac2, an2 = load_ac_pair(CSV_DIR, POP1, POP2, CHR, START, END)
    print(f"  {len(ac1)} variants loaded")
    print()

    print("Computing Hudson's Fst with scikit-allel...")
    fst_allel = hudson_fst_scikit_allel(ac1, an1, ac2, an2)
    print(f"  scikit-allel Fst = {fst_allel:.10f}")
    print()

    print("Querying graphpop.divergence for Fst...")
    fst_neo4j = query_neo4j_fst()
    print(f"  Neo4j Fst        = {fst_neo4j:.10f}")
    print()

    rel_err = abs(fst_allel - fst_neo4j) / abs(fst_allel) * 100
    status = "PASS" if rel_err < 0.01 else ("WARN" if rel_err < 1.0 else "FAIL")

    print("-" * 72)
    print(f"  {status}  Hudson's Fst  allel={fst_allel:.10f}  neo4j={fst_neo4j:.10f}  RelErr={rel_err:.6f}%")
    print("=" * 72)

    return 0 if rel_err < 1.0 else 1


if __name__ == "__main__":
    sys.exit(main())
