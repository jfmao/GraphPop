"""Download Reactome pathway and GO term annotations and emit CSVs for Neo4j LOAD CSV."""

from __future__ import annotations

import csv
import logging
import urllib.request
from pathlib import Path

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# CSV headers (LOAD CSV compatible)
# ---------------------------------------------------------------------------

PATHWAY_HEADER = ["pathwayId:ID(Pathway)", ":LABEL", "name", "source"]
GOTERM_HEADER = ["goId:ID(GOTerm)", ":LABEL", "name", "aspect"]
IN_PATHWAY_HEADER = [":START_ID(Gene)", ":END_ID(Pathway)", ":TYPE", "evidence"]
HAS_GO_TERM_HEADER = [":START_ID(Gene)", ":END_ID(GOTerm)", ":TYPE"]

# ---------------------------------------------------------------------------
# URLs
# ---------------------------------------------------------------------------

REACTOME_URL = "https://reactome.org/download/current/Ensembl2Reactome.txt"

# BioMart REST query for human gene → GO term mapping
_BIOMART_XML = """\
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1"
       uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="go_id"/>
    <Attribute name="name_1006"/>
    <Attribute name="namespace_1003"/>
  </Dataset>
</Query>"""

BIOMART_URL = (
    "https://www.ensembl.org/biomart/martservice?query="
    + urllib.request.quote(_BIOMART_XML.strip())
)

# Map BioMart namespace values to shorter aspect names
_ASPECT_MAP = {
    "biological_process": "biological_process",
    "molecular_function": "molecular_function",
    "cellular_component": "cellular_component",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load_gene_id_set(gene_csv_path: str | Path) -> set[str]:
    """Read geneId column from gene_nodes.csv."""
    ids: set[str] = set()
    with open(gene_csv_path) as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            if row:
                ids.add(row[0])
    logger.info("Loaded %d gene IDs from %s", len(ids), gene_csv_path)
    return ids


# ---------------------------------------------------------------------------
# PathwayParser
# ---------------------------------------------------------------------------


class PathwayParser:
    """Download Reactome and GO annotations and emit CSVs for LOAD CSV.

    Filters output to genes already in the graph (reads gene_nodes.csv).

    Usage::

        parser = PathwayParser("output/csv", gene_csv_path="output/csv/gene_nodes.csv")
        parser.run()
    """

    def __init__(
        self,
        out_dir: str | Path,
        *,
        gene_id_set: set[str] | None = None,
        gene_csv_path: str | Path | None = None,
        cache_dir: str | Path | None = None,
    ) -> None:
        self._out_dir = Path(out_dir)
        self._out_dir.mkdir(parents=True, exist_ok=True)

        if gene_id_set is not None:
            self._gene_ids = gene_id_set
        elif gene_csv_path is not None:
            self._gene_ids = _load_gene_id_set(gene_csv_path)
        else:
            self._gene_ids = None  # no filtering

        if cache_dir is not None:
            self._cache_dir = Path(cache_dir)
        else:
            self._cache_dir = self._out_dir / "_cache"
        self._cache_dir.mkdir(parents=True, exist_ok=True)

        # Stats
        self._n_pathways = 0
        self._n_goterms = 0
        self._n_in_pathway = 0
        self._n_has_go_term = 0

    # -- Properties ---------------------------------------------------------

    @property
    def n_pathways(self) -> int:
        return self._n_pathways

    @property
    def n_goterms(self) -> int:
        return self._n_goterms

    @property
    def n_in_pathway(self) -> int:
        return self._n_in_pathway

    @property
    def n_has_go_term(self) -> int:
        return self._n_has_go_term

    # -- Downloads ----------------------------------------------------------

    def download_reactome(self) -> Path:
        """Download Ensembl2Reactome.txt, skip if cached."""
        dest = self._cache_dir / "Ensembl2Reactome.txt"
        if dest.exists():
            logger.info("Reactome file cached: %s", dest)
            return dest
        logger.info("Downloading Reactome from %s ...", REACTOME_URL)
        urllib.request.urlretrieve(REACTOME_URL, dest)
        logger.info("Reactome downloaded: %s", dest)
        return dest

    def download_biomart_go(self) -> Path:
        """Download GO annotations from BioMart, skip if cached."""
        dest = self._cache_dir / "biomart_go.tsv"
        if dest.exists():
            logger.info("BioMart GO file cached: %s", dest)
            return dest
        logger.info("Downloading GO terms from BioMart ...")
        urllib.request.urlretrieve(BIOMART_URL, dest)
        logger.info("BioMart GO downloaded: %s", dest)
        return dest

    # -- Parsing ------------------------------------------------------------

    def parse_reactome(self, path: Path) -> None:
        """Parse Reactome TSV, filter for Homo sapiens + known genes, emit CSVs."""
        pathways: dict[str, str] = {}  # pathwayId → name
        edges: set[tuple[str, str, str]] = set()  # (geneId, pathwayId, evidence)

        with open(path) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 6:
                    continue

                gene_id, reactome_id, _url, pathway_name, evidence, species = parts[:6]

                if species != "Homo sapiens":
                    continue
                if self._gene_ids is not None and gene_id not in self._gene_ids:
                    continue

                pathways[reactome_id] = pathway_name
                edges.add((gene_id, reactome_id, evidence))

        # Write pathway_nodes.csv
        node_path = self._out_dir / "pathway_nodes.csv"
        with open(node_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(PATHWAY_HEADER)
            for pid in sorted(pathways):
                w.writerow([pid, "Pathway", pathways[pid], "Reactome"])
        self._n_pathways = len(pathways)

        # Write in_pathway_edges.csv
        edge_path = self._out_dir / "in_pathway_edges.csv"
        with open(edge_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(IN_PATHWAY_HEADER)
            for gene_id, pid, evidence in sorted(edges):
                w.writerow([gene_id, pid, "IN_PATHWAY", evidence])
        self._n_in_pathway = len(edges)

        logger.info(
            "Reactome: %d pathways, %d IN_PATHWAY edges",
            self._n_pathways, self._n_in_pathway,
        )

    def parse_go(self, path: Path) -> None:
        """Parse BioMart GO TSV, filter for known genes, emit CSVs."""
        goterms: dict[str, tuple[str, str]] = {}  # goId → (name, aspect)
        edges: set[tuple[str, str]] = set()  # (geneId, goId)

        with open(path) as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, None)
            if header is None:
                logger.warning("BioMart GO file is empty: %s", path)
                return

            for row in reader:
                if len(row) < 4:
                    continue

                gene_id, go_id, go_name, go_namespace = row[:4]

                if not go_id:
                    continue
                if self._gene_ids is not None and gene_id not in self._gene_ids:
                    continue

                aspect = _ASPECT_MAP.get(go_namespace, go_namespace)
                goterms[go_id] = (go_name, aspect)
                edges.add((gene_id, go_id))

        # Write goterm_nodes.csv
        node_path = self._out_dir / "goterm_nodes.csv"
        with open(node_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(GOTERM_HEADER)
            for gid in sorted(goterms):
                name, aspect = goterms[gid]
                w.writerow([gid, "GOTerm", name, aspect])
        self._n_goterms = len(goterms)

        # Write has_go_term_edges.csv
        edge_path = self._out_dir / "has_go_term_edges.csv"
        with open(edge_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(HAS_GO_TERM_HEADER)
            for gene_id, gid in sorted(edges):
                w.writerow([gene_id, gid, "HAS_GO_TERM"])
        self._n_has_go_term = len(edges)

        logger.info(
            "GO: %d terms, %d HAS_GO_TERM edges",
            self._n_goterms, self._n_has_go_term,
        )

    # -- Orchestration ------------------------------------------------------

    def run(self) -> None:
        """Download and parse all sources. BioMart failure does not block Reactome."""
        # Reactome
        reactome_path = self.download_reactome()
        self.parse_reactome(reactome_path)

        # BioMart GO (can be flaky)
        try:
            go_path = self.download_biomart_go()
            self.parse_go(go_path)
        except Exception as exc:
            logger.warning(
                "BioMart GO download/parse failed (Reactome data still available): %s",
                exc,
            )

        logger.info(
            "PathwayParser complete: %d pathways, %d GO terms, "
            "%d IN_PATHWAY edges, %d HAS_GO_TERM edges",
            self._n_pathways, self._n_goterms,
            self._n_in_pathway, self._n_has_go_term,
        )
