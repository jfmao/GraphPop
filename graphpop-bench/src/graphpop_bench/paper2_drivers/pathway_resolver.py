"""Reactome pathway → variant resolver for 1000G real-data panels.

Used by the R2 Fig 3e real-data driver to translate a Reactome
pathway ID into the set of 1000G variant IDs that hit genes in
that pathway on a specific chromosome.

Inputs (pre-existing on disk):
- `data/raw/1000g/csv_out/in_pathway_edges.csv` —
  `:START_ID(Gene),:END_ID(Pathway),:TYPE,evidence`
- `data/raw/1000g/csv_out/pathway_nodes.csv` —
  `pathwayId:ID(Pathway),:LABEL,name,source`
- `data/raw/1000g/csv_out/has_consequence_edges.csv` —
  `:START_ID(Variant),:END_ID(Gene),:TYPE,consequence,impact,...`
  (9.4 GB; streamed, not loaded fully into memory)

Variant IDs use the convention `chr{N}:{pos}:{ref}:{alt}` (Paper-1
ingest pipeline schema).
"""
from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, Set


DEFAULT_PATHWAY_NODES = Path(
    "/mnt/data/GraphPop/data/raw/1000g/csv_out/pathway_nodes.csv")
DEFAULT_IN_PATHWAY_EDGES = Path(
    "/mnt/data/GraphPop/data/raw/1000g/csv_out/in_pathway_edges.csv")
DEFAULT_HAS_CONSEQUENCE_EDGES = Path(
    "/mnt/data/GraphPop/data/raw/1000g/csv_out/has_consequence_edges.csv")


# ---------------------------------------------------------------------------
# Pathway-gene mapping
# ---------------------------------------------------------------------------

def load_pathway_names(
    path: Path = DEFAULT_PATHWAY_NODES,
) -> Dict[str, str]:
    """`pathwayId → name` for all Reactome pathways in the catalog."""
    if not path.exists():
        raise FileNotFoundError(f"pathway nodes CSV: {path}")
    out: Dict[str, str] = {}
    with open(path) as fh:
        reader = csv.reader(fh)
        header = next(reader)
        for row in reader:
            if len(row) >= 3:
                out[row[0]] = row[2]
    return out


def load_pathway_to_genes(
    path: Path = DEFAULT_IN_PATHWAY_EDGES,
) -> Dict[str, Set[str]]:
    """`pathwayId → {geneId, ...}` reverse-mapping the Reactome edges.

    Reads the IN_PATHWAY edges (~1500 rows on the existing
    1000G ingest); memory cost is trivial.
    """
    if not path.exists():
        raise FileNotFoundError(f"in_pathway_edges CSV: {path}")
    out: Dict[str, Set[str]] = defaultdict(set)
    with open(path) as fh:
        reader = csv.reader(fh)
        header = next(reader)
        for row in reader:
            if len(row) >= 2:
                gene_id, pathway_id = row[0], row[1]
                out[pathway_id].add(gene_id)
    return dict(out)


# ---------------------------------------------------------------------------
# chr-restricted gene → variant stream
# ---------------------------------------------------------------------------

@dataclass
class VariantRow:
    """One streamed row of the has_consequence catalog."""

    variant_id: str
    gene_id: str
    consequence: str


def stream_has_consequence(
    path: Path = DEFAULT_HAS_CONSEQUENCE_EDGES,
    chr_filter: str | None = None,
) -> Iterator[VariantRow]:
    """Stream rows from the 9.4 GB has_consequence_edges CSV.

    chr_filter : if set (e.g. "22"), only yield rows whose
                 variant_id starts with `chr{chr_filter}:`.
    """
    if not path.exists():
        raise FileNotFoundError(f"has_consequence CSV: {path}")
    prefix = f"chr{chr_filter}:" if chr_filter else None
    with open(path) as fh:
        # We DON'T use csv.reader for the bulk stream — string-prefix
        # match on the raw line is dramatically faster than per-row
        # CSV parse on a 9.4 GB file.
        header_line = fh.readline()
        for line in fh:
            if prefix is not None and not line.startswith(prefix):
                continue
            # parts[0]=variant_id, parts[1]=gene_id, parts[2]=:TYPE,
            # parts[3]=consequence
            parts = line.rstrip("\n").split(",", 4)
            if len(parts) < 4:
                continue
            yield VariantRow(
                variant_id=parts[0],
                gene_id=parts[1],
                consequence=parts[3],
            )


def chr_pathway_variants(
    has_consequence_csv: Path,
    gene_set: Set[str],
    chr_filter: str = "22",
) -> Set[str]:
    """Return the set of `variant_id`s on `chr_filter` whose target
    gene is in `gene_set`.

    Stream-filters the has_consequence CSV; single-pass.
    """
    keep_genes = set(gene_set)
    if not keep_genes:
        return set()
    out: Set[str] = set()
    for row in stream_has_consequence(
            has_consequence_csv, chr_filter=chr_filter):
        if row.gene_id in keep_genes:
            out.add(row.variant_id)
    return out


# ---------------------------------------------------------------------------
# Convenience top-level
# ---------------------------------------------------------------------------

def chr_pathway_variants_by_id(
    pathway_id: str,
    chr_filter: str = "22",
    pathway_edges_csv: Path = DEFAULT_IN_PATHWAY_EDGES,
    has_consequence_csv: Path = DEFAULT_HAS_CONSEQUENCE_EDGES,
) -> tuple[Set[str], Set[str]]:
    """One-shot helper: (gene_set, variant_set) for a given pathway.

    `variant_set` is restricted to `chr_filter`; `gene_set` is the
    full Reactome gene set for the pathway (not chr-restricted).
    """
    pathway_to_genes = load_pathway_to_genes(pathway_edges_csv)
    gene_set = pathway_to_genes.get(pathway_id, set())
    variant_set = chr_pathway_variants(
        has_consequence_csv, gene_set, chr_filter=chr_filter)
    return gene_set, variant_set
