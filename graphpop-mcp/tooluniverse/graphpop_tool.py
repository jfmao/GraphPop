"""ToolUniverse BaseTool wrapper for GraphPop.

Provides 12 population genomics tools that map to GraphPop's Neo4j stored procedures:
  - GraphPop_Diversity — pi, theta_W, Tajima's D, Fay & Wu's H, heterozygosity, F_IS
  - GraphPop_Divergence — Hudson Fst, W&C Fst, Dxy, Da, PBS
  - GraphPop_SFS — site frequency spectrum (folded/unfolded)
  - GraphPop_JointSFS — 2D joint SFS
  - GraphPop_GenomeScan — sliding-window genome scan
  - GraphPop_LD — pairwise r2, D'
  - GraphPop_iHS — integrated haplotype score
  - GraphPop_XPEHH — cross-population EHH
  - GraphPop_nSL — number of segregating sites by length
  - GraphPop_ROH — runs of homozygosity (HMM or sliding-window)
  - GraphPop_GarudH — Garud's H1/H12/H2_H1 (sweep detection)
  - GraphPop_PopSummary — whole-chromosome population summary

Each tool connects directly to a Neo4j database containing a GraphPop graph.
All tools support annotation-conditioned queries via consequence/pathway parameters.

Environment variables:
    GRAPHPOP_NEO4J_URI      (default: bolt://localhost:7687)
    GRAPHPOP_NEO4J_USER     (default: neo4j)
    GRAPHPOP_NEO4J_PASSWORD (default: graphpop)
"""

from __future__ import annotations

import json
import os

_driver = None


def _get_driver():
    """Lazy-initialize a shared Neo4j driver."""
    global _driver
    if _driver is None:
        from neo4j import GraphDatabase

        _driver = GraphDatabase.driver(
            os.environ.get("GRAPHPOP_NEO4J_URI", "bolt://localhost:7687"),
            auth=(
                os.environ.get("GRAPHPOP_NEO4J_USER", "neo4j"),
                os.environ.get("GRAPHPOP_NEO4J_PASSWORD", "graphpop"),
            ),
        )
    return _driver


def _run(cypher: str, params: dict | None = None) -> list[dict]:
    """Execute a Cypher statement and return records as dicts."""
    with _get_driver().session() as session:
        result = session.run(cypher, params or {})
        return [r.data() for r in result]


def _ok(result):
    return {"result": result, "success": True}


def _err(msg: str):
    return {"result": msg, "success": False}


def _build_options(args: dict, keys: list[str]) -> dict:
    """Extract option keys from arguments dict, skipping None values."""
    return {k: args[k] for k in keys if args.get(k) is not None}


# ---------------------------------------------------------------------------
# Tool implementations
# ---------------------------------------------------------------------------


def diversity(arguments: dict) -> dict:
    """GraphPop_Diversity — diversity statistics with optional annotation conditioning."""
    try:
        opts = _build_options(arguments, [
            "min_af", "max_af", "variant_type", "consequence", "pathway", "samples",
        ])
        results = _run(
            "CALL graphpop.diversity($chr, $start, $end, $pop, $options)",
            {"chr": arguments["chr"], "start": arguments["start"],
             "end": arguments["end"], "pop": arguments["pop"], "options": opts},
        )
        return _ok(results)
    except Exception as e:
        return _err(str(e))


def divergence(arguments: dict) -> dict:
    """GraphPop_Divergence — Fst, Dxy, Da, PBS."""
    try:
        opts = _build_options(arguments, [
            "min_af", "max_af", "variant_type", "consequence", "pathway",
            "samples1", "samples2", "pop3",
        ])
        results = _run(
            "CALL graphpop.divergence($chr, $start, $end, $pop1, $pop2, $options)",
            {"chr": arguments["chr"], "start": arguments["start"],
             "end": arguments["end"], "pop1": arguments["pop1"],
             "pop2": arguments["pop2"], "options": opts},
        )
        return _ok(results)
    except Exception as e:
        return _err(str(e))


def sfs(arguments: dict) -> dict:
    """GraphPop_SFS — site frequency spectrum."""
    try:
        opts = _build_options(arguments, [
            "min_af", "max_af", "variant_type", "consequence", "pathway", "samples",
        ])
        results = _run(
            "CALL graphpop.sfs($chr, $start, $end, $pop, $folded, $options)",
            {"chr": arguments["chr"], "start": arguments["start"],
             "end": arguments["end"], "pop": arguments["pop"],
             "folded": arguments.get("folded", False), "options": opts},
        )
        return _ok(results)
    except Exception as e:
        return _err(str(e))


def joint_sfs(arguments: dict) -> dict:
    """GraphPop_JointSFS — 2D joint SFS."""
    try:
        opts = _build_options(arguments, [
            "min_af", "max_af", "variant_type", "consequence", "pathway",
            "samples1", "samples2",
        ])
        results = _run(
            "CALL graphpop.joint_sfs($chr, $start, $end, $pop1, $pop2, $folded, $options)",
            {"chr": arguments["chr"], "start": arguments["start"],
             "end": arguments["end"], "pop1": arguments["pop1"],
             "pop2": arguments["pop2"],
             "folded": arguments.get("folded", False), "options": opts},
        )
        return _ok(results)
    except Exception as e:
        return _err(str(e))


def genome_scan(arguments: dict) -> dict:
    """GraphPop_GenomeScan — sliding-window genome scan."""
    try:
        opts = _build_options(arguments, [
            "min_af", "max_af", "variant_type", "consequence", "pathway",
            "pop2", "pop3", "run_id",
        ])
        results = _run(
            "CALL graphpop.genome_scan($chr, $pop, $window, $step, $options)",
            {"chr": arguments["chr"], "pop": arguments["pop"],
             "window": arguments.get("window", 100000),
             "step": arguments.get("step", 50000), "options": opts},
        )
        return _ok(results)
    except Exception as e:
        return _err(str(e))


def ld(arguments: dict) -> dict:
    """GraphPop_LD — pairwise linkage disequilibrium."""
    try:
        opts = _build_options(arguments, ["min_af", "variant_type", "samples", "write_edges"])
        results = _run(
            "CALL graphpop.ld($chr, $start, $end, $pop, $max_dist, $r2_threshold, $options)",
            {"chr": arguments["chr"], "start": arguments["start"],
             "end": arguments["end"], "pop": arguments["pop"],
             "max_dist": arguments.get("max_dist", 500000),
             "r2_threshold": arguments.get("r2_threshold", 0.2), "options": opts},
        )
        return _ok(results)
    except Exception as e:
        return _err(str(e))


def ihs(arguments: dict) -> dict:
    """GraphPop_iHS — integrated haplotype score."""
    try:
        opts = _build_options(arguments, [
            "min_af", "max_af", "variant_type", "min_ehh", "max_gap", "n_af_bins",
        ])
        results = _run(
            "CALL graphpop.ihs($chr, $pop, $options)",
            {"chr": arguments["chr"], "pop": arguments["pop"], "options": opts},
        )
        return _ok(results)
    except Exception as e:
        return _err(str(e))


def xpehh(arguments: dict) -> dict:
    """GraphPop_XPEHH — cross-population EHH."""
    try:
        opts = _build_options(arguments, [
            "min_af", "variant_type", "min_ehh", "max_gap", "chunk_size", "ehh_margin",
        ])
        results = _run(
            "CALL graphpop.xpehh($chr, $pop1, $pop2, $options)",
            {"chr": arguments["chr"], "pop1": arguments["pop1"],
             "pop2": arguments["pop2"], "options": opts},
        )
        return _ok(results)
    except Exception as e:
        return _err(str(e))


def nsl(arguments: dict) -> dict:
    """GraphPop_nSL — number of segregating sites by length."""
    try:
        opts = _build_options(arguments, ["min_af", "max_af", "variant_type", "n_af_bins"])
        results = _run(
            "CALL graphpop.nsl($chr, $pop, $options)",
            {"chr": arguments["chr"], "pop": arguments["pop"], "options": opts},
        )
        return _ok(results)
    except Exception as e:
        return _err(str(e))


def roh(arguments: dict) -> dict:
    """GraphPop_ROH — runs of homozygosity."""
    try:
        opts = _build_options(arguments, [
            "method", "min_length", "min_snps", "het_error_rate",
            "hmm_min_af", "variant_type",
        ])
        results = _run(
            "CALL graphpop.roh($chr, $pop, $options)",
            {"chr": arguments["chr"], "pop": arguments["pop"], "options": opts},
        )
        return _ok(results)
    except Exception as e:
        return _err(str(e))


def garud_h(arguments: dict) -> dict:
    """GraphPop_GarudH — Garud's H sweep statistics."""
    try:
        opts = _build_options(arguments, ["min_af", "variant_type"])
        results = _run(
            "CALL graphpop.garud_h($chr, $pop, $window, $step, $options)",
            {"chr": arguments["chr"], "pop": arguments["pop"],
             "window": arguments.get("window", 100000),
             "step": arguments.get("step", 50000), "options": opts},
        )
        return _ok(results)
    except Exception as e:
        return _err(str(e))


def pop_summary(arguments: dict) -> dict:
    """GraphPop_PopSummary — whole-chromosome population summary."""
    try:
        opts = _build_options(arguments, [
            "min_af", "max_af", "variant_type", "consequence", "pathway", "samples",
        ])
        results = _run(
            "CALL graphpop.pop_summary($chr, $pop, $options)",
            {"chr": arguments["chr"], "pop": arguments["pop"], "options": opts},
        )
        return _ok(results)
    except Exception as e:
        return _err(str(e))


# ---------------------------------------------------------------------------
# Tool registry
# ---------------------------------------------------------------------------

TOOL_REGISTRY = {
    "GraphPop_Diversity": diversity,
    "GraphPop_Divergence": divergence,
    "GraphPop_SFS": sfs,
    "GraphPop_JointSFS": joint_sfs,
    "GraphPop_GenomeScan": genome_scan,
    "GraphPop_LD": ld,
    "GraphPop_iHS": ihs,
    "GraphPop_XPEHH": xpehh,
    "GraphPop_nSL": nsl,
    "GraphPop_ROH": roh,
    "GraphPop_GarudH": garud_h,
    "GraphPop_PopSummary": pop_summary,
}


def run(tool_name: str, arguments: dict) -> dict:
    """ToolUniverse entry point — dispatch to the appropriate tool function."""
    handler = TOOL_REGISTRY.get(tool_name)
    if handler is None:
        return _err(f"Unknown tool: {tool_name}")
    return handler(arguments)
