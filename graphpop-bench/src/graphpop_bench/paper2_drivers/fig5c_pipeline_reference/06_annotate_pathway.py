"""Stage 6 — Annotate variants with pathway membership.

Reads a VEP-annotated variant table (CSQ field from VEP) +
a Reactome pathway JSON, emits a `(variant_id, pathway_id)`
mapping table for downstream filtering.

Equivalent in Cypher: variants are already linked to
:Pathway nodes via the M4.A pathway ingest; the predicate
`restrict_to_pathway: "P_cardio_signaling"` resolves the
filter at query time without any per-stage glue.
"""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path


def parse_vep_consequence(csq_field: str) -> list[dict]:
    """Decode a VEP CSQ field (pipe-separated subfields, comma-
    separated transcripts) into a list of consequence records.
    """
    schema = [
        "Allele", "Consequence", "Gene", "Feature_type",
        "Feature", "BIOTYPE", "EXON", "INTRON",
    ]
    out = []
    for chunk in csq_field.split(","):
        parts = chunk.split("|")
        rec = {schema[i]: parts[i] if i < len(parts) else ""
               for i in range(len(schema))}
        out.append(rec)
    return out


def load_reactome_pathways(path: Path) -> dict[str, set[str]]:
    """Load Reactome pathway → gene-symbol set mapping."""
    blob = json.loads(path.read_text())
    out: dict[str, set[str]] = {}
    for pathway in blob.get("pathways", []):
        pid = pathway["id"]
        out[pid] = set(pathway.get("genes", []))
    return out


def annotate_variants(
    vep_table_path: Path,
    reactome_path: Path,
    out_tsv: Path,
) -> int:
    """For each variant in the VEP table, emit one row per
    pathway it overlaps via gene membership.
    """
    pathways = load_reactome_pathways(reactome_path)
    n_rows = 0
    with open(vep_table_path) as fh_in, open(out_tsv, "w") as fh_out:
        writer = csv.writer(fh_out, delimiter="\t")
        writer.writerow(["variant_id", "pathway_id", "gene", "consequence"])
        reader = csv.DictReader(fh_in, delimiter="\t")
        for row in reader:
            variant_id = row["variant_id"]
            csq_records = parse_vep_consequence(row["CSQ"])
            for csq in csq_records:
                gene = csq.get("Gene", "")
                if not gene:
                    continue
                for pid, gene_set in pathways.items():
                    if gene in gene_set:
                        writer.writerow([
                            variant_id, pid, gene,
                            csq.get("Consequence", "")])
                        n_rows += 1
    return n_rows


def main(argv: list[str]) -> int:
    if len(argv) != 4:
        sys.stderr.write(
            "usage: 06_annotate_pathway.py <vep.tsv> "
            "<reactome.json> <out.tsv>\n")
        return 2
    n = annotate_variants(Path(argv[1]), Path(argv[2]),
                            Path(argv[3]))
    print(f"[06_pathway] wrote {n} (variant, pathway) rows")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
