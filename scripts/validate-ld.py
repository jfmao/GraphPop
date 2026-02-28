#!/usr/bin/env python3
"""Validate GraphPop LD (r²) against scikit-allel's rogers_huff_r.

Computes pairwise r² for a small region from the 1000G chr22 VCF using
scikit-allel, then compares against GraphPop graphpop.ld procedure output.

Usage:
    conda run -n graphevo python scripts/validate-ld.py
"""

import subprocess
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

VCF = "data/raw/1000g/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz"
PANEL = "data/raw/1000g/integrated_call_samples_v3.20130502.ALL.panel"

CHR = "chr22"
START = 16_000_000
END = 16_100_000  # Small region to keep validation fast

POP = "EUR"

NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"

MAX_DIST = 500_000
R2_THRESHOLD = 0.2

# ---------------------------------------------------------------------------
# Load population panel
# ---------------------------------------------------------------------------

def load_panel(panel_path):
    """Load sample→super_pop mapping."""
    sample_pop = {}
    with open(panel_path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                sample_pop[parts[0]] = parts[2]  # super_pop
    return sample_pop

# ---------------------------------------------------------------------------
# Compute LD with scikit-allel
# ---------------------------------------------------------------------------

def compute_ld_scikit_allel():
    """Compute pairwise r² using scikit-allel."""
    import allel
    import cyvcf2

    sample_pop = load_panel(PANEL)

    vcf = cyvcf2.VCF(VCF)
    all_samples = list(vcf.samples)

    # Find indices for target population
    pop_indices = [i for i, s in enumerate(all_samples) if sample_pop.get(s) == POP]
    print(f"  {POP}: {len(pop_indices)} samples")

    # Extract genotypes for the region
    positions = []
    variant_ids = []
    genotypes = []

    region = f"chr22:{START}-{END}"
    for v in vcf(region):
        if v.is_snp and len(v.ALT) == 1:
            gt = v.genotypes  # list of [a1, a2, phased]
            # Extract dosages for pop samples
            dosage = []
            for idx in pop_indices:
                g = gt[idx]
                dosage.append(g[0] + g[1])
            dosage = np.array(dosage, dtype=float)

            # Skip monomorphic
            if dosage.min() == dosage.max():
                continue

            positions.append(v.POS)
            variant_ids.append(f"{v.CHROM}:{v.POS}:{v.REF}:{v.ALT[0]}")
            genotypes.append(dosage)

    vcf.close()

    if len(genotypes) < 2:
        print("  Not enough polymorphic variants for LD computation")
        return {}

    genotypes = np.array(genotypes)
    positions = np.array(positions)
    print(f"  {len(positions)} polymorphic variants in {CHR}:{START}-{END}")

    # Compute pairwise r² using rogers_huff_r from scikit-allel
    # allel.rogers_huff_r expects a 2D array of shape (n_variants, n_samples)
    r_matrix = allel.rogers_huff_r(genotypes)

    # r_matrix is a condensed distance matrix (upper triangle)
    # Convert to square form
    from scipy.spatial.distance import squareform
    r_square = squareform(r_matrix)
    r2_matrix = r_square ** 2

    # Build dict of (var1, var2) → r² for pairs above threshold
    ld_pairs = {}
    n = len(variant_ids)
    for i in range(n):
        for j in range(i + 1, n):
            dist = abs(positions[j] - positions[i])
            if dist > MAX_DIST:
                break
            r2_val = r2_matrix[i, j]
            if r2_val >= R2_THRESHOLD:
                key = (variant_ids[i], variant_ids[j])
                ld_pairs[key] = r2_val

    print(f"  {len(ld_pairs)} pairs with r² >= {R2_THRESHOLD}")
    return ld_pairs

# ---------------------------------------------------------------------------
# Query GraphPop LD
# ---------------------------------------------------------------------------

def query_graphpop_ld():
    """Query LD edges from Neo4j."""
    cypher = (
        f"MATCH (v1:Variant)-[r:LD {{population: '{POP}'}}]->(v2:Variant) "
        f"WHERE v1.chr = '{CHR}' AND v1.pos >= {START} AND v1.pos <= {END} "
        f"AND v2.pos >= {START} AND v2.pos <= {END} "
        f"RETURN v1.variantId AS v1, v2.variantId AS v2, r.r2 AS r2 "
        f"ORDER BY v1.pos, v2.pos"
    )

    result = subprocess.run(
        ["cypher-shell", "-u", NEO4J_USER, "-p", NEO4J_PASS,
         "--format", "plain", cypher],
        capture_output=True, text=True
    )

    if result.returncode != 0:
        print(f"  cypher-shell error: {result.stderr}")
        return {}

    ld_pairs = {}
    for line in result.stdout.strip().split("\n"):
        line = line.strip()
        if not line or line.startswith("v1") or line.startswith("+"):
            continue
        parts = [p.strip().strip('"') for p in line.split(",")]
        if len(parts) >= 3:
            try:
                v1, v2, r2 = parts[0], parts[1], float(parts[2])
                ld_pairs[(v1, v2)] = r2
            except (ValueError, IndexError):
                continue

    print(f"  {len(ld_pairs)} LD edges from GraphPop")
    return ld_pairs

# ---------------------------------------------------------------------------
# Compare
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("Validating graphpop.ld against scikit-allel rogers_huff_r")
    print("=" * 60)
    print()

    print("Computing LD with scikit-allel...")
    allel_ld = compute_ld_scikit_allel()
    print()

    print("Querying GraphPop LD edges...")
    gp_ld = query_graphpop_ld()
    print()

    if not allel_ld or not gp_ld:
        print("WARNING: Not enough data for comparison.")
        print("  Run graphpop.ld first: CALL graphpop.ld('chr22', 16000000, 16100000, 'EUR', 500000, 0.2)")
        return

    # Find common pairs (accounting for direction — GraphPop may store A→B, allel has A,B)
    matched = 0
    max_err = 0.0
    sum_abs_err = 0.0
    errors = []

    for pair, allel_r2 in allel_ld.items():
        gp_r2 = gp_ld.get(pair) or gp_ld.get((pair[1], pair[0]))
        if gp_r2 is not None:
            matched += 1
            err = abs(allel_r2 - gp_r2)
            max_err = max(max_err, err)
            sum_abs_err += err
            errors.append(err)

    print(f"Matched pairs: {matched}")
    print(f"scikit-allel pairs (r² >= {R2_THRESHOLD}): {len(allel_ld)}")
    print(f"GraphPop pairs: {len(gp_ld)}")

    if matched > 0:
        mean_err = sum_abs_err / matched
        print(f"Mean absolute error: {mean_err:.6f}")
        print(f"Max absolute error:  {max_err:.6f}")

        if max_err < 0.001:
            print("\nPASS: r² values match scikit-allel within 0.001 absolute error")
        elif max_err < 0.01:
            print("\nPASS (loose): r² values match within 0.01 absolute error")
        else:
            print(f"\nFAIL: max error {max_err:.6f} exceeds tolerance")
    else:
        print("\nWARNING: No overlapping pairs found. Check variant IDs match.")

if __name__ == "__main__":
    main()
