"""Unit tests for the R2 pathway-resolver + Fig 3e real driver."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from graphpop_bench.paper2_drivers import (
    fig3e_real_panels as f3e_real,
    pathway_resolver as pr,
)


# ---------------------------------------------------------------------------
# Pathway resolver — uses real 1000G catalog when available
# ---------------------------------------------------------------------------

def _pathway_csvs_available() -> bool:
    return (
        pr.DEFAULT_PATHWAY_NODES.exists()
        and pr.DEFAULT_IN_PATHWAY_EDGES.exists()
    )


@pytest.mark.skipif(
    not _pathway_csvs_available(),
    reason="1000G pathway CSVs not present on this host",
)
def test_load_pathway_names_returns_reactome_pathways():
    names = pr.load_pathway_names()
    assert len(names) > 500    # 1000G ingest has ~663 pathways
    # Spot-check a known pathway.
    assert "R-HSA-202733" in names
    assert "vascular wall" in names["R-HSA-202733"].lower()


@pytest.mark.skipif(
    not _pathway_csvs_available(),
    reason="1000G pathway CSVs not present on this host",
)
def test_load_pathway_to_genes_finds_pathway():
    p2g = pr.load_pathway_to_genes()
    genes = p2g.get("R-HSA-202733", set())
    assert len(genes) >= 10
    # All entries are ENSG-style IDs.
    assert all(g.startswith("ENSG") for g in genes)


# ---------------------------------------------------------------------------
# Synthetic-CSV unit tests (no real 1000G needed)
# ---------------------------------------------------------------------------

def test_load_pathway_to_genes_synthetic(tmp_path):
    p = tmp_path / "edges.csv"
    p.write_text(
        ":START_ID(Gene),:END_ID(Pathway),:TYPE,evidence\n"
        "ENSG_A,P1,IN_PATHWAY,TAS\n"
        "ENSG_B,P1,IN_PATHWAY,TAS\n"
        "ENSG_C,P2,IN_PATHWAY,TAS\n"
    )
    out = pr.load_pathway_to_genes(p)
    assert out == {
        "P1": {"ENSG_A", "ENSG_B"},
        "P2": {"ENSG_C"},
    }


def test_load_pathway_names_synthetic(tmp_path):
    p = tmp_path / "nodes.csv"
    p.write_text(
        "pathwayId:ID(Pathway),:LABEL,name,source\n"
        "P1,Pathway,Cell signalling,Reactome\n"
        "P2,Pathway,DNA repair,Reactome\n"
    )
    out = pr.load_pathway_names(p)
    assert out == {"P1": "Cell signalling", "P2": "DNA repair"}


def test_stream_has_consequence_chr_filter(tmp_path):
    p = tmp_path / "hc.csv"
    p.write_text(
        ":START_ID(Variant),:END_ID(Gene),:TYPE,consequence\n"
        "chr22:100:A:G,ENSG_A,HAS_CONSEQUENCE,missense\n"
        "chr1:200:T:C,ENSG_B,HAS_CONSEQUENCE,synonymous\n"
        "chr22:300:G:A,ENSG_C,HAS_CONSEQUENCE,intron\n"
    )
    rows = list(pr.stream_has_consequence(p, chr_filter="22"))
    assert len(rows) == 2
    assert {r.variant_id for r in rows} == {
        "chr22:100:A:G", "chr22:300:G:A"}


def test_chr_pathway_variants_synthetic(tmp_path):
    p = tmp_path / "hc.csv"
    p.write_text(
        ":START_ID(Variant),:END_ID(Gene),:TYPE,consequence\n"
        "chr22:100:A:G,ENSG_A,HAS_CONSEQUENCE,missense\n"
        "chr22:200:T:C,ENSG_X,HAS_CONSEQUENCE,synonymous\n"
        "chr22:300:G:A,ENSG_B,HAS_CONSEQUENCE,intron\n"
        "chr1:400:T:C,ENSG_A,HAS_CONSEQUENCE,intron\n"
    )
    out = pr.chr_pathway_variants(p, {"ENSG_A", "ENSG_B"},
                                    chr_filter="22")
    assert out == {"chr22:100:A:G", "chr22:300:G:A"}


# ---------------------------------------------------------------------------
# parse_variant_position
# ---------------------------------------------------------------------------

def test_parse_variant_position_canonical():
    assert f3e_real.parse_variant_position("chr22:1234567:A:G") == 1234567


def test_parse_variant_position_malformed():
    assert f3e_real.parse_variant_position("not-a-variant") is None
    assert f3e_real.parse_variant_position("chr22:") is None
    assert f3e_real.parse_variant_position("chr22:NaN:A:G") is None


def test_variant_set_to_positions_with_region():
    vids = {
        "chr22:100:A:G",
        "chr22:5000:A:G",
        "chr22:10000:A:G",
        "chr22:99999:A:G",
    }
    out = f3e_real.variant_set_to_positions(
        vids, in_region=(1000, 50000))
    assert out == {5000, 10000}


# ---------------------------------------------------------------------------
# pathway_child_nodes_from_ts (gated on msprime for the fixture)
# ---------------------------------------------------------------------------

def _msprime_available() -> bool:
    try:
        import msprime  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.mark.skipif(
    not _msprime_available(),
    reason="msprime required for synthetic .trees fixture",
)
def test_pathway_child_nodes_from_ts_filters_by_position(tmp_path):
    import msprime
    import tskit
    ts = msprime.sim_ancestry(
        samples=5, sequence_length=10_000,
        recombination_rate=1e-5, random_seed=42)
    ts = msprime.sim_mutations(
        ts, rate=1e-3, random_seed=42,
        model=msprime.BinaryMutationModel())
    # Pick the first 3 site positions as "pathway".
    site_positions = [int(s.position) for s in ts.sites()][:3]
    children = f3e_real.pathway_child_nodes_from_ts(
        ts, set(site_positions))
    assert len(children) > 0
    # Empty pathway set → empty children.
    empty = f3e_real.pathway_child_nodes_from_ts(ts, set())
    assert empty == set()
