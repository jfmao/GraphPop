#!/usr/bin/env python3
"""Rice 3K Genome Analysis — run all 12 GraphPop procedures with resource tracking.

Demonstrates GraphPop on the 3K Rice Genome Project (3kRG):
  - 29.6M biallelic SNPs, 3,024 accessions, 12 chromosomes
  - 13 fine-grained populations (GJ-tmp, XI-1A, XI-3, etc.)
  - SnpEff functional annotations (HAS_CONSEQUENCE edges)

Usage:
    python scripts/rice_analysis.py [--chr Chr1] [--pop1 GJ-tmp] [--pop2 XI-1A]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass, field

from neo4j import GraphDatabase

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"
NEO4J_DB = "neo4j"


@dataclass
class RunResult:
    procedure: str
    parameters: str
    wall_sec: float = 0.0
    heap_used_mb: float = 0.0
    heap_max_mb: float = 0.0
    result_summary: str = ""
    error: str = ""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def get_heap_stats(session) -> tuple[float, float]:
    """Call graphpop.heap_stats() and return (used_mb, max_mb)."""
    try:
        rec = session.run("CALL graphpop.heap_stats()").single()
        if rec:
            return (rec["usedMB"], rec["maxMB"])
    except Exception:
        pass
    return (0.0, 0.0)


def run_procedure(session, cypher: str, name: str, params_desc: str) -> RunResult:
    """Execute a Cypher procedure call with timing and heap tracking."""
    result = RunResult(procedure=name, parameters=params_desc)
    try:
        t0 = time.time()
        records = list(session.run(cypher))
        result.wall_sec = time.time() - t0

        used, mx = get_heap_stats(session)
        result.heap_used_mb = used
        result.heap_max_mb = mx

        # Summarize result
        if records:
            r0 = records[0]
            keys = r0.keys()
            if len(records) == 1:
                summary_parts = []
                for k in keys[:6]:
                    v = r0[k]
                    if isinstance(v, float):
                        summary_parts.append(f"{k}={v:.6f}")
                    else:
                        summary_parts.append(f"{k}={v}")
                result.result_summary = ", ".join(summary_parts)
            else:
                result.result_summary = f"{len(records)} rows, cols={list(keys)}"
        else:
            result.result_summary = "(no rows)"
    except Exception as e:
        result.wall_sec = time.time() - t0 if 't0' in dir() else 0.0
        result.error = str(e)
    return result


# ---------------------------------------------------------------------------
# Procedure runs
# ---------------------------------------------------------------------------


def _get_chr_length(session, chrom: str) -> int:
    """Get chromosome length from the database."""
    rec = session.run(
        "MATCH (c:Chromosome {chromosomeId: $chr}) RETURN c.length AS len",
        chr=chrom,
    ).single()
    return rec["len"] if rec else 50_000_000


def run_all_procedures(
    session, chrom: str, pop1: str, pop2: str, pop3: str,
) -> list[RunResult]:
    results: list[RunResult] = []
    chr_len = _get_chr_length(session, chrom)

    # 1. Diversity (chr, start, end, pop)
    results.append(run_procedure(
        session,
        f"CALL graphpop.diversity('{chrom}', 1, {chr_len}, '{pop1}')",
        "graphpop.diversity",
        f"chr={chrom}, pop={pop1}",
    ))

    # 2. Divergence (chr, start, end, pop1, pop2)
    results.append(run_procedure(
        session,
        f"CALL graphpop.divergence('{chrom}', 1, {chr_len}, '{pop1}', '{pop2}')",
        "graphpop.divergence",
        f"chr={chrom}, pop1={pop1}, pop2={pop2}",
    ))

    # 3. SFS (chr, start, end, pop)
    results.append(run_procedure(
        session,
        f"CALL graphpop.sfs('{chrom}', 1, {chr_len}, '{pop3}')",
        "graphpop.sfs",
        f"chr={chrom}, pop={pop3}",
    ))

    # 4. Joint SFS (chr, start, end, pop1, pop2)
    results.append(run_procedure(
        session,
        f"CALL graphpop.joint_sfs('{chrom}', 1, {chr_len}, '{pop1}', '{pop2}')",
        "graphpop.joint_sfs",
        f"chr={chrom}, pop1={pop1}, pop2={pop2}",
    ))

    # 5. Genome scan (chr, pop, window, step)
    results.append(run_procedure(
        session,
        f"CALL graphpop.genome_scan('{chrom}', '{pop1}', 100000, 50000)",
        "graphpop.genome_scan",
        f"chr={chrom}, pop={pop1}, window=100kb",
    ))

    # 6. LD (chr, start, end, pop)
    results.append(run_procedure(
        session,
        f"CALL graphpop.ld('{chrom}', 1, 500000, '{pop1}')",
        "graphpop.ld",
        f"chr={chrom}, pop={pop1}, region=1-500kb",
    ))

    # 7. iHS (chr, pop)
    results.append(run_procedure(
        session,
        f"CALL graphpop.ihs('{chrom}', '{pop3}')",
        "graphpop.ihs",
        f"chr={chrom}, pop={pop3}",
    ))

    # 8. XP-EHH (chr, pop1, pop2)
    results.append(run_procedure(
        session,
        f"CALL graphpop.xpehh('{chrom}', '{pop1}', '{pop2}')",
        "graphpop.xpehh",
        f"chr={chrom}, pop1={pop1}, pop2={pop2}",
    ))

    # 9. nSL (chr, pop)
    results.append(run_procedure(
        session,
        f"CALL graphpop.nsl('{chrom}', '{pop3}')",
        "graphpop.nsl",
        f"chr={chrom}, pop={pop3}",
    ))

    # 10. ROH (chr, pop)
    results.append(run_procedure(
        session,
        f"CALL graphpop.roh('{chrom}', '{pop1}')",
        "graphpop.roh",
        f"chr={chrom}, pop={pop1}",
    ))

    # 11. Pop summary (chr, pop)
    results.append(run_procedure(
        session,
        f"CALL graphpop.pop_summary('{chrom}', '{pop1}')",
        "graphpop.pop_summary",
        f"chr={chrom}, pop={pop1}",
    ))

    # 12. Garud's H (chr, pop, window, step)
    results.append(run_procedure(
        session,
        f"CALL graphpop.garud_h('{chrom}', '{pop1}', 100000, 50000)",
        "graphpop.garud_h",
        f"chr={chrom}, pop={pop1}, window=100kb",
    ))

    return results


# ---------------------------------------------------------------------------
# Functional annotation demos
# ---------------------------------------------------------------------------


def run_annotation_demos(
    session, chrom: str, pop1: str, pop2: str,
) -> list[RunResult]:
    results: list[RunResult] = []
    chr_len = _get_chr_length(session, chrom)

    # Demo 1: Diversity at missense variants
    results.append(run_procedure(
        session,
        f"CALL graphpop.diversity('{chrom}', 1, {chr_len}, '{pop1}', "
        f"{{consequence: 'missense_variant'}})",
        "diversity@missense",
        f"chr={chrom}, pop={pop1}, consequence=missense_variant",
    ))

    # Demo 2: Diversity at synonymous variants
    results.append(run_procedure(
        session,
        f"CALL graphpop.diversity('{chrom}', 1, {chr_len}, '{pop1}', "
        f"{{consequence: 'synonymous_variant'}})",
        "diversity@synonymous",
        f"chr={chrom}, pop={pop1}, consequence=synonymous_variant",
    ))

    # Demo 3: Genome scan at missense variants
    results.append(run_procedure(
        session,
        f"CALL graphpop.genome_scan('{chrom}', '{pop1}', 100000, 50000, "
        f"{{consequence: 'missense_variant'}})",
        "genome_scan@missense",
        f"chr={chrom}, pop={pop1}, window=100kb, consequence=missense_variant",
    ))

    # Demo 4: Divergence at missense variants
    results.append(run_procedure(
        session,
        f"CALL graphpop.divergence('{chrom}', 1, {chr_len}, '{pop1}', '{pop2}', "
        f"{{consequence: 'missense_variant'}})",
        "divergence@missense",
        f"chr={chrom}, pop1={pop1}, pop2={pop2}, consequence=missense_variant",
    ))

    # Demo 5: Gene-level query (top differentiated genes)
    gene_cypher = """
    MATCH (v:Variant)-[c:HAS_CONSEQUENCE]->(g:Gene)
    WHERE v.chr = $chr AND c.impact IN ['HIGH', 'MODERATE']
    RETURN g.symbol AS gene, g.geneId AS geneId,
           count(v) AS n_variants, collect(DISTINCT c.consequence) AS consequences
    ORDER BY n_variants DESC LIMIT 20
    """
    results.append(run_procedure(
        session,
        gene_cypher.replace("$chr", f"'{chrom}'"),
        "gene_query@HIGH+MODERATE",
        f"chr={chrom}, impact=HIGH/MODERATE, top 20 genes",
    ))

    return results


# ---------------------------------------------------------------------------
# Database overview
# ---------------------------------------------------------------------------


def db_overview(session) -> dict:
    """Get node/edge counts and population list."""
    info: dict = {}

    # Node counts
    rec = session.run(
        "MATCH (n) RETURN labels(n)[0] AS label, count(n) AS cnt ORDER BY label"
    )
    info["node_counts"] = {r["label"]: r["cnt"] for r in rec}

    # Edge counts
    rec = session.run(
        "MATCH ()-[r]->() RETURN type(r) AS rtype, count(r) AS cnt ORDER BY rtype"
    )
    info["edge_counts"] = {r["rtype"]: r["cnt"] for r in rec}

    # Populations
    rec = session.run(
        "MATCH (p:Population) RETURN p.name AS name, p.n_samples AS n ORDER BY name"
    )
    info["populations"] = {r["name"]: r["n"] for r in rec}

    # Chromosomes
    rec = session.run(
        "MATCH (c:Chromosome) RETURN c.chromosomeId AS chr, c.length AS len ORDER BY chr"
    )
    info["chromosomes"] = {r["chr"]: r["len"] for r in rec}

    return info


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def print_markdown(
    overview: dict,
    proc_results: list[RunResult],
    demo_results: list[RunResult],
) -> None:
    print("# Rice 3K Genome — GraphPop Analysis Results\n")

    # Database overview
    print("## Database Overview\n")
    print("| Label | Count |")
    print("|-------|------:|")
    for label, cnt in sorted(overview.get("node_counts", {}).items()):
        print(f"| {label} | {cnt:,} |")
    print()
    print("| Relationship | Count |")
    print("|-------------|------:|")
    for rtype, cnt in sorted(overview.get("edge_counts", {}).items()):
        print(f"| {rtype} | {cnt:,} |")
    print()

    # Populations
    print("## Populations\n")
    print("| Population | Samples |")
    print("|-----------|--------:|")
    for name, n in sorted(overview.get("populations", {}).items()):
        print(f"| {name} | {n} |")
    print()

    # Procedure results
    print("## Procedure Results\n")
    print("| # | Procedure | Parameters | Wall (s) | Heap Used (MB) | Result |")
    print("|---|-----------|-----------|--------:|---------------:|--------|")
    for i, r in enumerate(proc_results, 1):
        status = r.error if r.error else r.result_summary
        if len(status) > 80:
            status = status[:77] + "..."
        print(f"| {i} | `{r.procedure}` | {r.parameters} | "
              f"{r.wall_sec:.1f} | {r.heap_used_mb:.0f} | {status} |")
    print()

    total_time = sum(r.wall_sec for r in proc_results)
    print(f"**Total procedure time: {total_time:.1f}s**\n")

    # Functional annotation demos
    if demo_results:
        print("## Functional Annotation Demos\n")
        print("| Demo | Parameters | Wall (s) | Result |")
        print("|------|-----------|--------:|--------|")
        for r in demo_results:
            status = r.error if r.error else r.result_summary
            if len(status) > 80:
                status = status[:77] + "..."
            print(f"| `{r.procedure}` | {r.parameters} | "
                  f"{r.wall_sec:.1f} | {status} |")
        print()

    # Biological validation
    print("## Biological Validation\n")
    div_result = next((r for r in proc_results if r.procedure == "graphpop.divergence"), None)
    if div_result and not div_result.error:
        print(f"- **Indica–Japonica divergence**: {div_result.result_summary}")
        print("  - Expected Fst ~0.3–0.5 for major subspecies split")
    div_result2 = next((r for r in proc_results if r.procedure == "graphpop.diversity"), None)
    if div_result2 and not div_result2.error:
        print(f"- **Diversity**: {div_result2.result_summary}")
        print("  - Rice is self-fertilizing → expect high F_IS (~0.9+)")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description="Rice 3K Genome GraphPop Analysis")
    parser.add_argument("--chr", default="Chr1", help="Chromosome (default: Chr1)")
    parser.add_argument("--pop1", default="GJ-tmp", help="Population 1 (default: GJ-tmp)")
    parser.add_argument("--pop2", default="XI-1A", help="Population 2 (default: XI-1A)")
    parser.add_argument("--pop3", default="XI-3", help="Population 3 for SFS/iHS/nSL (default: XI-3)")
    parser.add_argument("--json", action="store_true", help="Also output JSON results")
    parser.add_argument("--skip-procs", action="store_true", help="Skip procedure runs")
    parser.add_argument("--skip-demos", action="store_true", help="Skip annotation demos")
    args = parser.parse_args()

    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))

    with driver.session(database=NEO4J_DB) as session:
        # Database overview
        print("Collecting database overview...", file=sys.stderr)
        overview = db_overview(session)

        # Procedure runs
        proc_results: list[RunResult] = []
        if not args.skip_procs:
            print(f"Running 12 procedures on {args.chr}...", file=sys.stderr)
            proc_results = run_all_procedures(
                session, args.chr, args.pop1, args.pop2, args.pop3,
            )
            for r in proc_results:
                status = "OK" if not r.error else f"ERROR: {r.error[:60]}"
                print(f"  {r.procedure}: {r.wall_sec:.1f}s — {status}", file=sys.stderr)

        # Annotation demos
        demo_results: list[RunResult] = []
        if not args.skip_demos:
            print("Running functional annotation demos...", file=sys.stderr)
            demo_results = run_annotation_demos(
                session, args.chr, args.pop1, args.pop2,
            )
            for r in demo_results:
                status = "OK" if not r.error else f"ERROR: {r.error[:60]}"
                print(f"  {r.procedure}: {r.wall_sec:.1f}s — {status}", file=sys.stderr)

    driver.close()

    # Output
    print_markdown(overview, proc_results, demo_results)

    if args.json:
        all_results = [
            {
                "procedure": r.procedure,
                "parameters": r.parameters,
                "wall_sec": r.wall_sec,
                "heap_used_mb": r.heap_used_mb,
                "result_summary": r.result_summary,
                "error": r.error,
            }
            for r in proc_results + demo_results
        ]
        with open("rice_analysis_results.json", "w") as f:
            json.dump({"overview": overview, "results": all_results}, f, indent=2)
        print("JSON results written to rice_analysis_results.json", file=sys.stderr)


if __name__ == "__main__":
    main()
