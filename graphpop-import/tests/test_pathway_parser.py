"""Tests for pathway_parser.py — Reactome pathway and GO term parsing."""

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from graphpop_import.pathway_parser import (
    GOTERM_HEADER,
    HAS_GO_TERM_HEADER,
    IN_PATHWAY_HEADER,
    PATHWAY_HEADER,
    PathwayParser,
    _load_gene_id_set,
)


# ---------------------------------------------------------------------------
# _load_gene_id_set tests
# ---------------------------------------------------------------------------


class TestLoadGeneIdSet:
    def test_loads_gene_ids(self, tmp_path):
        """Reads geneId from the first column (skips header)."""
        gene_csv = tmp_path / "gene_nodes.csv"
        gene_csv.write_text(
            "geneId:ID(Gene),:LABEL,symbol,biotype\n"
            "ENSG00000001,Gene,BRCA1,protein_coding\n"
            "ENSG00000002,Gene,TP53,protein_coding\n"
        )
        ids = _load_gene_id_set(gene_csv)
        assert ids == {"ENSG00000001", "ENSG00000002"}

    def test_empty_file(self, tmp_path):
        """Header-only file returns empty set."""
        gene_csv = tmp_path / "gene_nodes.csv"
        gene_csv.write_text("geneId:ID(Gene),:LABEL,symbol,biotype\n")
        ids = _load_gene_id_set(gene_csv)
        assert ids == set()


# ---------------------------------------------------------------------------
# Reactome parsing tests
# ---------------------------------------------------------------------------

# Reactome TSV format: gene_id \t reactome_id \t url \t pathway_name \t evidence \t species
_REACTOME_CONTENT = """\
ENSG00000001\tR-HSA-001\thttp://reactome.org\tApoptosis\tIEA\tHomo sapiens
ENSG00000002\tR-HSA-001\thttp://reactome.org\tApoptosis\tTAS\tHomo sapiens
ENSG00000002\tR-HSA-002\thttp://reactome.org\tCell cycle\tIEA\tHomo sapiens
ENSG00000003\tR-HSA-003\thttp://reactome.org\tRice pathway\tIEA\tOryza sativa
ENSG00000099\tR-HSA-004\thttp://reactome.org\tSignaling\tIEA\tHomo sapiens
"""


class TestParseReactome:
    def test_basic_parsing(self, tmp_path):
        """Parses Reactome TSV, filters by species and gene set."""
        reactome_file = tmp_path / "Ensembl2Reactome.txt"
        reactome_file.write_text(_REACTOME_CONTENT)

        parser = PathwayParser(
            tmp_path / "output",
            gene_id_set={"ENSG00000001", "ENSG00000002"},
        )
        parser.parse_reactome(reactome_file)

        assert parser.n_pathways == 2  # R-HSA-001, R-HSA-002
        assert parser.n_in_pathway == 3  # 3 gene-pathway edges (deduplicated by set)

    def test_pathway_nodes_csv(self, tmp_path):
        """pathway_nodes.csv has correct header and content."""
        reactome_file = tmp_path / "Ensembl2Reactome.txt"
        reactome_file.write_text(_REACTOME_CONTENT)

        out_dir = tmp_path / "output"
        parser = PathwayParser(out_dir, gene_id_set={"ENSG00000001", "ENSG00000002"})
        parser.parse_reactome(reactome_file)

        rows = list(csv.reader(open(out_dir / "pathway_nodes.csv")))
        assert rows[0] == PATHWAY_HEADER
        # 2 pathways
        assert len(rows) == 3
        pathway_ids = {r[0] for r in rows[1:]}
        assert pathway_ids == {"R-HSA-001", "R-HSA-002"}
        # Check source column
        assert all(r[3] == "Reactome" for r in rows[1:])

    def test_in_pathway_edges_csv(self, tmp_path):
        """in_pathway_edges.csv has correct edges with evidence codes."""
        reactome_file = tmp_path / "Ensembl2Reactome.txt"
        reactome_file.write_text(_REACTOME_CONTENT)

        out_dir = tmp_path / "output"
        parser = PathwayParser(out_dir, gene_id_set={"ENSG00000001", "ENSG00000002"})
        parser.parse_reactome(reactome_file)

        rows = list(csv.reader(open(out_dir / "in_pathway_edges.csv")))
        assert rows[0] == IN_PATHWAY_HEADER
        assert len(rows) == 4  # header + 3 edges
        # All edges should have IN_PATHWAY as type
        assert all(r[2] == "IN_PATHWAY" for r in rows[1:])

    def test_species_filter(self, tmp_path):
        """Non-Homo sapiens entries are excluded."""
        reactome_file = tmp_path / "Ensembl2Reactome.txt"
        reactome_file.write_text(_REACTOME_CONTENT)

        out_dir = tmp_path / "output"
        # Include ENSG00000003 in gene set, but it's Oryza sativa
        parser = PathwayParser(
            out_dir,
            gene_id_set={"ENSG00000003"},
        )
        parser.parse_reactome(reactome_file)

        assert parser.n_pathways == 0
        assert parser.n_in_pathway == 0

    def test_gene_filter(self, tmp_path):
        """Only genes in gene_id_set are included."""
        reactome_file = tmp_path / "Ensembl2Reactome.txt"
        reactome_file.write_text(_REACTOME_CONTENT)

        out_dir = tmp_path / "output"
        parser = PathwayParser(out_dir, gene_id_set={"ENSG00000001"})
        parser.parse_reactome(reactome_file)

        assert parser.n_pathways == 1  # only R-HSA-001
        assert parser.n_in_pathway == 1

    def test_no_gene_filter(self, tmp_path):
        """Without gene_id_set, all Homo sapiens entries are included."""
        reactome_file = tmp_path / "Ensembl2Reactome.txt"
        reactome_file.write_text(_REACTOME_CONTENT)

        out_dir = tmp_path / "output"
        parser = PathwayParser(out_dir)  # no gene filtering
        parser.parse_reactome(reactome_file)

        # ENSG00000001→R-HSA-001, ENSG00000002→R-HSA-001, ENSG00000002→R-HSA-002, ENSG00000099→R-HSA-004
        assert parser.n_pathways == 3  # R-HSA-001, R-HSA-002, R-HSA-004
        assert parser.n_in_pathway == 4


# ---------------------------------------------------------------------------
# GO term parsing tests
# ---------------------------------------------------------------------------

_GO_TSV_CONTENT = """\
Gene stable ID\tGO term accession\tGO term name\tGO domain
ENSG00000001\tGO:0006915\tapoptotic process\tbiological_process
ENSG00000001\tGO:0005737\tcytoplasm\tcellular_component
ENSG00000002\tGO:0006915\tapoptotic process\tbiological_process
ENSG00000002\tGO:0003677\tDNA binding\tmolecular_function
ENSG00000099\tGO:0007049\tcell cycle\tbiological_process
"""


class TestParseGO:
    def test_basic_parsing(self, tmp_path):
        """Parses BioMart GO TSV, filters by gene set."""
        go_file = tmp_path / "biomart_go.tsv"
        go_file.write_text(_GO_TSV_CONTENT)

        out_dir = tmp_path / "output"
        parser = PathwayParser(
            out_dir,
            gene_id_set={"ENSG00000001", "ENSG00000002"},
        )
        parser.parse_go(go_file)

        assert parser.n_goterms == 3  # GO:0006915, GO:0005737, GO:0003677
        assert parser.n_has_go_term == 4  # 4 gene-GO edges

    def test_goterm_nodes_csv(self, tmp_path):
        """goterm_nodes.csv has correct header and content."""
        go_file = tmp_path / "biomart_go.tsv"
        go_file.write_text(_GO_TSV_CONTENT)

        out_dir = tmp_path / "output"
        parser = PathwayParser(
            out_dir,
            gene_id_set={"ENSG00000001", "ENSG00000002"},
        )
        parser.parse_go(go_file)

        rows = list(csv.reader(open(out_dir / "goterm_nodes.csv")))
        assert rows[0] == GOTERM_HEADER
        assert len(rows) == 4  # header + 3 GO terms
        go_ids = {r[0] for r in rows[1:]}
        assert go_ids == {"GO:0006915", "GO:0005737", "GO:0003677"}

    def test_aspect_mapping(self, tmp_path):
        """GO aspects are preserved correctly."""
        go_file = tmp_path / "biomart_go.tsv"
        go_file.write_text(_GO_TSV_CONTENT)

        out_dir = tmp_path / "output"
        parser = PathwayParser(out_dir, gene_id_set={"ENSG00000001"})
        parser.parse_go(go_file)

        rows = list(csv.reader(open(out_dir / "goterm_nodes.csv")))
        aspect_by_id = {r[0]: r[3] for r in rows[1:]}
        assert aspect_by_id["GO:0006915"] == "biological_process"
        assert aspect_by_id["GO:0005737"] == "cellular_component"

    def test_has_go_term_edges_csv(self, tmp_path):
        """has_go_term_edges.csv has correct structure."""
        go_file = tmp_path / "biomart_go.tsv"
        go_file.write_text(_GO_TSV_CONTENT)

        out_dir = tmp_path / "output"
        parser = PathwayParser(
            out_dir,
            gene_id_set={"ENSG00000001", "ENSG00000002"},
        )
        parser.parse_go(go_file)

        rows = list(csv.reader(open(out_dir / "has_go_term_edges.csv")))
        assert rows[0] == HAS_GO_TERM_HEADER
        assert len(rows) == 5  # header + 4 edges
        assert all(r[2] == "HAS_GO_TERM" for r in rows[1:])

    def test_gene_filter_go(self, tmp_path):
        """Only genes in gene_id_set get GO term edges."""
        go_file = tmp_path / "biomart_go.tsv"
        go_file.write_text(_GO_TSV_CONTENT)

        out_dir = tmp_path / "output"
        parser = PathwayParser(out_dir, gene_id_set={"ENSG00000001"})
        parser.parse_go(go_file)

        assert parser.n_goterms == 2  # GO:0006915, GO:0005737
        assert parser.n_has_go_term == 2

    def test_empty_go_file(self, tmp_path):
        """Empty GO file (header only) produces 0 terms."""
        go_file = tmp_path / "biomart_go.tsv"
        go_file.write_text(
            "Gene stable ID\tGO term accession\tGO term name\tGO domain\n"
        )

        out_dir = tmp_path / "output"
        parser = PathwayParser(out_dir)
        parser.parse_go(go_file)

        assert parser.n_goterms == 0
        assert parser.n_has_go_term == 0

    def test_skip_rows_without_go_id(self, tmp_path):
        """Rows with empty GO ID are skipped."""
        go_file = tmp_path / "biomart_go.tsv"
        go_file.write_text(
            "Gene stable ID\tGO term accession\tGO term name\tGO domain\n"
            "ENSG00000001\t\t\t\n"
            "ENSG00000001\tGO:0006915\tapoptotic process\tbiological_process\n"
        )

        out_dir = tmp_path / "output"
        parser = PathwayParser(out_dir)
        parser.parse_go(go_file)

        assert parser.n_goterms == 1
        assert parser.n_has_go_term == 1


# ---------------------------------------------------------------------------
# PathwayParser initialization tests
# ---------------------------------------------------------------------------


class TestPathwayParserInit:
    def test_gene_id_set_from_csv(self, tmp_path):
        """Parser loads gene IDs from gene_csv_path."""
        gene_csv = tmp_path / "gene_nodes.csv"
        gene_csv.write_text(
            "geneId:ID(Gene),:LABEL,symbol,biotype\n"
            "ENSG00000001,Gene,BRCA1,protein_coding\n"
        )

        parser = PathwayParser(tmp_path / "output", gene_csv_path=gene_csv)
        assert parser._gene_ids == {"ENSG00000001"}

    def test_gene_id_set_direct(self, tmp_path):
        """Parser accepts direct gene_id_set."""
        parser = PathwayParser(
            tmp_path / "output",
            gene_id_set={"GENE1", "GENE2"},
        )
        assert parser._gene_ids == {"GENE1", "GENE2"}

    def test_no_gene_filter(self, tmp_path):
        """Without gene_id_set or gene_csv_path, no filtering."""
        parser = PathwayParser(tmp_path / "output")
        assert parser._gene_ids is None

    def test_cache_dir(self, tmp_path):
        """Custom cache_dir is used."""
        cache = tmp_path / "my_cache"
        parser = PathwayParser(tmp_path / "output", cache_dir=cache)
        assert parser._cache_dir == cache
        assert cache.exists()

    def test_default_cache_dir(self, tmp_path):
        """Default cache is _cache under out_dir."""
        out_dir = tmp_path / "output"
        parser = PathwayParser(out_dir)
        assert parser._cache_dir == out_dir / "_cache"
