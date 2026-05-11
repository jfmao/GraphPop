"""Stage 7 — Join PLINK / KING / ADMIXTURE / RFMix / pathway
tables to recover "cryptic 2nd-degree relatives in EUR-ancestry
segments enriched for cardiovascular pathway".

Equivalent in Cypher: the entire fig5a query — a single CALL
chain. Here we do the same logic in pandas over the outputs of
the upstream stages.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def load_king_pairs(kin0_path: Path) -> pd.DataFrame:
    """KING kin0 → (sample_a, sample_b, phi) frame."""
    df = pd.read_csv(kin0_path, sep="\t")
    df = df.rename(columns={
        "IID1": "sample_a", "IID2": "sample_b",
        "Kinship": "phi"})
    return df[["sample_a", "sample_b", "phi"]]


def load_admixture_q(q_path: Path,
                     fam_path: Path) -> pd.DataFrame:
    """ADMIXTURE Q matrix + .fam → sample → ancestry-proportions."""
    fam = pd.read_csv(fam_path, sep=r"\s+", header=None,
                      names=["FID", "IID", "PID", "MID",
                             "Sex", "Pheno"])
    q = pd.read_csv(q_path, sep=r"\s+", header=None)
    q.columns = [f"K{i+1}" for i in range(q.shape[1])]
    q["sample"] = fam["IID"]
    return q


def load_rfmix_calls(msp_path: Path) -> pd.DataFrame:
    """RFMix .msp.tsv → per-sample EUR-fraction summary.

    Crude scoring: average inferred ancestry over the genome.
    Real analysis would inspect per-haplotype windows.
    """
    df = pd.read_csv(msp_path, sep="\t", comment="#")
    # Each row is a window; columns are sample-haplotype calls.
    sample_cols = [c for c in df.columns
                   if c not in {"#chm", "spos", "epos",
                                "sgpos", "egpos", "n snps"}]
    eur_codes = {0}        # convention: 0 = EUR in the .map
    summary = {}
    for col in sample_cols:
        eur_frac = (df[col].isin(eur_codes)).mean()
        summary[col] = eur_frac
    return pd.DataFrame.from_dict(
        summary, orient="index", columns=["eur_fraction"]
    ).rename_axis("sample").reset_index()


def load_pathway_table(path: Path,
                       pathway_id: str) -> set[str]:
    """Pathway annotation TSV → set of variant_ids in the pathway."""
    df = pd.read_csv(path, sep="\t")
    return set(df.loc[df.pathway_id == pathway_id,
                      "variant_id"])


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--king-kin0", type=Path, required=True)
    p.add_argument("--admixture-q", type=Path, required=True)
    p.add_argument("--fam", type=Path, required=True)
    p.add_argument("--rfmix-msp", type=Path, required=True)
    p.add_argument("--pathway-tsv", type=Path, required=True)
    p.add_argument("--pathway-id", default="P_cardio_signaling")
    p.add_argument("--eur-min", type=float, default=0.8)
    p.add_argument("--phi-min", type=float, default=0.0625)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args(argv)

    king = load_king_pairs(args.king_kin0)
    adm = load_admixture_q(args.admixture_q, args.fam)
    rfmix = load_rfmix_calls(args.rfmix_msp)
    pathway_variants = load_pathway_table(
        args.pathway_tsv, args.pathway_id)

    # Filter pairs by 2nd-degree threshold.
    related = king[king.phi >= args.phi_min].copy()

    # Restrict to pairs where BOTH samples are EUR-dominant.
    eur_samples = set(rfmix.loc[
        rfmix.eur_fraction >= args.eur_min, "sample"])
    related = related[
        related.sample_a.isin(eur_samples)
        & related.sample_b.isin(eur_samples)
    ].copy()

    # Tag each pair with the count of pathway variants where
    # both samples are non-ref. (In Cypher, this is the
    # `restrict_to_pathway` predicate, applied inside
    # branch_grm rather than a post-hoc join.)
    related["n_pathway_variants_shared"] = len(pathway_variants)

    related.to_csv(args.out, index=False)
    sys.stderr.write(
        f"[07_join] {len(related)} candidate pairs → {args.out}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
