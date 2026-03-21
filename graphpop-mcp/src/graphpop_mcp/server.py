"""MCP server wrapping GraphPop Neo4j stored procedures for AI agent access.

Each tool maps to one of GraphPop's 12 population genomics procedures.
Connection parameters come from environment variables.
"""

from __future__ import annotations

import json
import os

from mcp.server.fastmcp import FastMCP

mcp = FastMCP(
    "GraphPop",
    instructions=(
        "Graph-native population genomics analytical engine. Runs 12 stored procedures "
        "inside a Neo4j graph database: diversity, divergence, SFS, joint SFS, genome scan, "
        "LD, iHS, XP-EHH, nSL, ROH, Garud's H, and population summary. "
        "Supports annotation-conditioned queries (consequence, pathway, gene) on any procedure."
    ),
)

_driver = None


def _get_driver():
    """Lazy-initialize a shared Neo4j driver."""
    global _driver
    if _driver is None:
        from neo4j import GraphDatabase

        uri = os.environ.get("GRAPHPOP_NEO4J_URI", "bolt://localhost:7687")
        user = os.environ.get("GRAPHPOP_NEO4J_USER", "neo4j")
        password = os.environ.get("GRAPHPOP_NEO4J_PASSWORD", "graphpop")
        _driver = GraphDatabase.driver(uri, auth=(user, password))
    return _driver


def _run_procedure(cypher: str, params: dict | None = None) -> list[dict]:
    """Execute a Cypher procedure call and return results as a list of dicts."""
    driver = _get_driver()
    with driver.session() as session:
        result = session.run(cypher, params or {})
        return [record.data() for record in result]


def _build_options(
    *,
    min_af: float | None = None,
    max_af: float | None = None,
    min_call_rate: float | None = None,
    hwe_pvalue: float | None = None,
    variant_type: str | None = None,
    consequence: str | None = None,
    pathway: str | None = None,
    samples: list[str] | None = None,
    samples1: list[str] | None = None,
    samples2: list[str] | None = None,
    pop2: str | None = None,
    pop3: str | None = None,
    extra: dict | None = None,
) -> dict:
    """Build an options map from keyword arguments, omitting None values."""
    opts = {}
    if min_af is not None:
        opts["min_af"] = min_af
    if max_af is not None:
        opts["max_af"] = max_af
    if min_call_rate is not None:
        opts["min_call_rate"] = min_call_rate
    if hwe_pvalue is not None:
        opts["hwe_pvalue"] = hwe_pvalue
    if variant_type is not None:
        opts["variant_type"] = variant_type
    if consequence is not None:
        opts["consequence"] = consequence
    if pathway is not None:
        opts["pathway"] = pathway
    if samples is not None:
        opts["samples"] = samples
    if samples1 is not None:
        opts["samples1"] = samples1
    if samples2 is not None:
        opts["samples2"] = samples2
    if pop2 is not None:
        opts["pop2"] = pop2
    if pop3 is not None:
        opts["pop3"] = pop3
    if extra:
        opts.update(extra)
    return opts


# ---------------------------------------------------------------------------
# Tools — 12 population genomics procedures
# ---------------------------------------------------------------------------


@mcp.tool()
def graphpop_diversity(
    chr: str,
    start: int,
    end: int,
    pop: str,
    min_af: float | None = None,
    max_af: float | None = None,
    variant_type: str | None = None,
    consequence: str | None = None,
    pathway: str | None = None,
    samples: list[str] | None = None,
) -> str:
    """Compute diversity statistics (pi, theta_W, Tajima's D, Fay & Wu's H, heterozygosity, F_IS).

    FAST PATH: O(V*K) using allele-count arrays on Variant nodes.

    Args:
        chr: Chromosome ID (e.g. "22", "chr1").
        start: Start position (inclusive).
        end: End position (inclusive).
        pop: Population name (e.g. "AFR", "EUR", "GJ-tmp").
        min_af: Minimum allele frequency filter.
        max_af: Maximum allele frequency filter.
        variant_type: Filter by variant type ("SNP", "INDEL").
        consequence: VEP consequence filter (e.g. "missense_variant", "synonymous_variant").
        pathway: Reactome pathway filter.
        samples: Custom sample list (overrides population).

    Returns JSON with: pi, theta_w, tajima_d, fay_wu_h, fay_wu_h_norm,
    het_exp, het_obs, fis, n_variants, n_segregating, n_polarized.
    """
    opts = _build_options(
        min_af=min_af, max_af=max_af, variant_type=variant_type,
        consequence=consequence, pathway=pathway, samples=samples,
    )
    cypher = "CALL graphpop.diversity($chr, $start, $end, $pop, $options)"
    results = _run_procedure(cypher, {
        "chr": chr, "start": start, "end": end, "pop": pop, "options": opts,
    })
    return json.dumps(results)


@mcp.tool()
def graphpop_divergence(
    chr: str,
    start: int,
    end: int,
    pop1: str,
    pop2: str,
    pop3: str | None = None,
    min_af: float | None = None,
    max_af: float | None = None,
    variant_type: str | None = None,
    consequence: str | None = None,
    pathway: str | None = None,
    samples1: list[str] | None = None,
    samples2: list[str] | None = None,
) -> str:
    """Compute divergence statistics (Hudson Fst, W&C Fst, Dxy, Da, PBS) between two populations.

    FAST PATH: O(V*K) using allele-count arrays. PBS requires pop3.

    Args:
        chr: Chromosome ID.
        start: Start position (inclusive).
        end: End position (inclusive).
        pop1: First population name.
        pop2: Second population name.
        pop3: Third population for PBS (optional).
        min_af: Minimum allele frequency filter.
        max_af: Maximum allele frequency filter.
        variant_type: Filter by variant type.
        consequence: VEP consequence filter (annotation-conditioned).
        pathway: Reactome pathway filter.
        samples1: Custom sample list for pop1.
        samples2: Custom sample list for pop2.

    Returns JSON with: fst_hudson, fst_wc, dxy, da, pbs, n_variants.
    """
    opts = _build_options(
        min_af=min_af, max_af=max_af, variant_type=variant_type,
        consequence=consequence, pathway=pathway,
        samples1=samples1, samples2=samples2, pop3=pop3,
    )
    cypher = "CALL graphpop.divergence($chr, $start, $end, $pop1, $pop2, $options)"
    results = _run_procedure(cypher, {
        "chr": chr, "start": start, "end": end,
        "pop1": pop1, "pop2": pop2, "options": opts,
    })
    return json.dumps(results)


@mcp.tool()
def graphpop_sfs(
    chr: str,
    start: int,
    end: int,
    pop: str,
    folded: bool = False,
    min_af: float | None = None,
    max_af: float | None = None,
    variant_type: str | None = None,
    consequence: str | None = None,
    pathway: str | None = None,
    samples: list[str] | None = None,
) -> str:
    """Compute the site frequency spectrum (SFS) for a population.

    FAST PATH. Supports folded and unfolded (polarized via ancestral allele) SFS.

    Args:
        chr: Chromosome ID.
        start: Start position (inclusive).
        end: End position (inclusive).
        pop: Population name.
        folded: If True, return folded SFS; if False, unfolded (default).
        min_af: Minimum allele frequency filter.
        max_af: Maximum allele frequency filter.
        variant_type: Filter by variant type.
        consequence: VEP consequence filter.
        pathway: Reactome pathway filter.
        samples: Custom sample list.

    Returns JSON with: sfs (array of counts), n_variants, max_ac, n_polarized.
    """
    opts = _build_options(
        min_af=min_af, max_af=max_af, variant_type=variant_type,
        consequence=consequence, pathway=pathway, samples=samples,
    )
    cypher = "CALL graphpop.sfs($chr, $start, $end, $pop, $folded, $options)"
    results = _run_procedure(cypher, {
        "chr": chr, "start": start, "end": end,
        "pop": pop, "folded": folded, "options": opts,
    })
    return json.dumps(results)


@mcp.tool()
def graphpop_joint_sfs(
    chr: str,
    start: int,
    end: int,
    pop1: str,
    pop2: str,
    folded: bool = False,
    min_af: float | None = None,
    max_af: float | None = None,
    variant_type: str | None = None,
    consequence: str | None = None,
    pathway: str | None = None,
    samples1: list[str] | None = None,
    samples2: list[str] | None = None,
) -> str:
    """Compute the 2D joint site frequency spectrum between two populations.

    FAST PATH. Useful for demographic inference (dadi, fastsimcoal2).

    Args:
        chr: Chromosome ID.
        start: Start position (inclusive).
        end: End position (inclusive).
        pop1: First population name.
        pop2: Second population name.
        folded: If True, return folded joint SFS.
        min_af: Minimum allele frequency filter.
        max_af: Maximum allele frequency filter.
        variant_type: Filter by variant type.
        consequence: VEP consequence filter.
        pathway: Reactome pathway filter.
        samples1: Custom sample list for pop1.
        samples2: Custom sample list for pop2.

    Returns JSON with: joint_sfs (flattened 2D array), n_variants, max_ac1, max_ac2, dim1, dim2.
    """
    opts = _build_options(
        min_af=min_af, max_af=max_af, variant_type=variant_type,
        consequence=consequence, pathway=pathway,
        samples1=samples1, samples2=samples2,
    )
    cypher = "CALL graphpop.joint_sfs($chr, $start, $end, $pop1, $pop2, $folded, $options)"
    results = _run_procedure(cypher, {
        "chr": chr, "start": start, "end": end,
        "pop1": pop1, "pop2": pop2, "folded": folded, "options": opts,
    })
    return json.dumps(results)


@mcp.tool()
def graphpop_genome_scan(
    chr: str,
    pop: str,
    window: int = 100000,
    step: int = 50000,
    pop2: str | None = None,
    pop3: str | None = None,
    min_af: float | None = None,
    max_af: float | None = None,
    variant_type: str | None = None,
    consequence: str | None = None,
    pathway: str | None = None,
    run_id: str | None = None,
) -> str:
    """Sliding-window genome scan with per-window statistics.

    Materializes GenomicWindow nodes. Computes pi, theta_W, Tajima's D, Fay & Wu's H,
    heterozygosity, F_IS per window. If pop2 provided, also Fst (Hudson + W&C), Dxy, Da.
    If pop3 provided, also PBS.

    Args:
        chr: Chromosome ID.
        pop: Population name.
        window: Window size in bp (default 100,000).
        step: Step size in bp (default 50,000).
        pop2: Second population for Fst/Dxy (optional).
        pop3: Third population for PBS (optional).
        min_af: Minimum allele frequency filter.
        max_af: Maximum allele frequency filter.
        variant_type: Filter by variant type.
        consequence: VEP consequence filter (annotation-conditioned genome scan).
        pathway: Reactome pathway filter.
        run_id: Custom run identifier.

    Returns JSON array of window results.
    """
    opts = _build_options(
        min_af=min_af, max_af=max_af, variant_type=variant_type,
        consequence=consequence, pathway=pathway,
        pop2=pop2, pop3=pop3,
        extra={"run_id": run_id} if run_id else None,
    )
    cypher = "CALL graphpop.genome_scan($chr, $pop, $window, $step, $options)"
    results = _run_procedure(cypher, {
        "chr": chr, "pop": pop, "window": window, "step": step, "options": opts,
    })
    return json.dumps(results)


@mcp.tool()
def graphpop_ld(
    chr: str,
    start: int,
    end: int,
    pop: str,
    max_dist: int = 500000,
    r2_threshold: float = 0.2,
    variant_type: str | None = None,
    min_af: float | None = None,
    samples: list[str] | None = None,
    write_edges: bool = True,
) -> str:
    """Compute pairwise linkage disequilibrium (r2, D') within a genomic region.

    FULL PATH: uses bit-packed HaplotypeMatrix. Writes LD edges between Variant nodes.

    Args:
        chr: Chromosome ID.
        start: Start position (inclusive).
        end: End position (inclusive).
        pop: Population name.
        max_dist: Maximum distance between variant pairs in bp (default 500,000).
        r2_threshold: Minimum r2 to report/write (default 0.2).
        variant_type: Filter by variant type.
        min_af: Minimum allele frequency filter.
        samples: Custom sample list.
        write_edges: Whether to write LD edges to the graph (default True).

    Returns JSON array with: variant1, variant2, r2, dprime, distance.
    """
    opts = _build_options(
        min_af=min_af, variant_type=variant_type, samples=samples,
        extra={"write_edges": write_edges},
    )
    cypher = (
        "CALL graphpop.ld($chr, $start, $end, $pop, $max_dist, $r2_threshold, $options)"
    )
    results = _run_procedure(cypher, {
        "chr": chr, "start": start, "end": end, "pop": pop,
        "max_dist": max_dist, "r2_threshold": r2_threshold, "options": opts,
    })
    return json.dumps(results)


@mcp.tool()
def graphpop_ihs(
    chr: str,
    pop: str,
    min_af: float | None = None,
    max_af: float | None = None,
    variant_type: str | None = None,
    min_ehh: float | None = None,
    max_gap: int | None = None,
    n_af_bins: int | None = None,
) -> str:
    """Compute integrated haplotype score (iHS) for a population across a chromosome.

    FULL PATH: bit-packed HaplotypeMatrix, EHH decay, AF-bin standardization.
    Results written as ihs_{pop} properties on Variant nodes.

    Args:
        chr: Chromosome ID.
        pop: Population name.
        min_af: Minimum allele frequency (default 0.05).
        max_af: Maximum allele frequency.
        variant_type: Filter by variant type (default "SNP").
        min_ehh: EHH truncation threshold (default 0.05).
        max_gap: Maximum gap in bp before aborting EHH walk (default 200,000).
        n_af_bins: Number of allele frequency bins for standardization (default 20).

    Returns JSON array with per-variant: variantId, pos, af, ihs, ihs_unstd.
    """
    opts: dict = {}
    if min_af is not None:
        opts["min_af"] = min_af
    if max_af is not None:
        opts["max_af"] = max_af
    if variant_type is not None:
        opts["variant_type"] = variant_type
    if min_ehh is not None:
        opts["min_ehh"] = min_ehh
    if max_gap is not None:
        opts["max_gap"] = max_gap
    if n_af_bins is not None:
        opts["n_af_bins"] = n_af_bins
    cypher = "CALL graphpop.ihs($chr, $pop, $options)"
    results = _run_procedure(cypher, {"chr": chr, "pop": pop, "options": opts})
    return json.dumps(results)


@mcp.tool()
def graphpop_xpehh(
    chr: str,
    pop1: str,
    pop2: str,
    min_af: float | None = None,
    variant_type: str | None = None,
    min_ehh: float | None = None,
    max_gap: int | None = None,
    chunk_size: int | None = None,
    ehh_margin: int | None = None,
) -> str:
    """Compute cross-population extended haplotype homozygosity (XP-EHH).

    FULL PATH: loads two populations' HaplotypeMatrices, genome-wide standardized.
    Results written as xpehh_{pop1}_{pop2} on Variant nodes. Chunked for memory.

    Args:
        chr: Chromosome ID.
        pop1: First population (focal).
        pop2: Second population (reference).
        min_af: Minimum allele frequency.
        variant_type: Filter by variant type (default "SNP").
        min_ehh: EHH truncation threshold.
        max_gap: Maximum gap in bp.
        chunk_size: Core window size for chunking (default 5,000,000).
        ehh_margin: EHH margin on each side (default 2,000,000).

    Returns JSON array with per-variant: variantId, pos, xpehh, xpehh_unstd.
    """
    opts: dict = {}
    if min_af is not None:
        opts["min_af"] = min_af
    if variant_type is not None:
        opts["variant_type"] = variant_type
    if min_ehh is not None:
        opts["min_ehh"] = min_ehh
    if max_gap is not None:
        opts["max_gap"] = max_gap
    if chunk_size is not None:
        opts["chunk_size"] = chunk_size
    if ehh_margin is not None:
        opts["ehh_margin"] = ehh_margin
    cypher = "CALL graphpop.xpehh($chr, $pop1, $pop2, $options)"
    results = _run_procedure(cypher, {
        "chr": chr, "pop1": pop1, "pop2": pop2, "options": opts,
    })
    return json.dumps(results)


@mcp.tool()
def graphpop_nsl(
    chr: str,
    pop: str,
    min_af: float | None = None,
    max_af: float | None = None,
    variant_type: str | None = None,
    n_af_bins: int | None = None,
) -> str:
    """Compute nSL (number of segregating sites by length) for a population.

    FULL PATH: pairwise SSL algorithm (Ferrer-Admetlla 2014), AF-bin standardized.
    Results written as nsl_{pop} on Variant nodes.

    Args:
        chr: Chromosome ID.
        pop: Population name.
        min_af: Minimum allele frequency.
        max_af: Maximum allele frequency.
        variant_type: Filter by variant type (default "SNP").
        n_af_bins: Number of AF bins for standardization (default 20).

    Returns JSON array with per-variant: variantId, pos, af, nsl, nsl_unstd.
    """
    opts: dict = {}
    if min_af is not None:
        opts["min_af"] = min_af
    if max_af is not None:
        opts["max_af"] = max_af
    if variant_type is not None:
        opts["variant_type"] = variant_type
    if n_af_bins is not None:
        opts["n_af_bins"] = n_af_bins
    cypher = "CALL graphpop.nsl($chr, $pop, $options)"
    results = _run_procedure(cypher, {"chr": chr, "pop": pop, "options": opts})
    return json.dumps(results)


@mcp.tool()
def graphpop_roh(
    chr: str,
    pop: str,
    method: str | None = None,
    min_length: int | None = None,
    min_snps: int | None = None,
    het_error_rate: float | None = None,
    hmm_min_af: float | None = None,
    variant_type: str | None = None,
) -> str:
    """Detect runs of homozygosity (ROH) per sample in a population.

    FULL PATH: HMM (default, Narasimhan 2016 / bcftools roh) or sliding-window method.
    Parallel per-sample computation.

    Args:
        chr: Chromosome ID.
        pop: Population name.
        method: "hmm" (default) or "sliding_window".
        min_length: Minimum ROH length in bp (default 500,000).
        min_snps: Minimum SNPs in ROH (default 25).
        het_error_rate: Heterozygous error rate for HMM (default 0.001).
        hmm_min_af: Minimum AF for HMM variants (default 0.0).
        variant_type: Filter by variant type (default "SNP").

    Returns JSON array with per-sample ROH segments: sample, start, end, length_kb, n_snps.
    """
    opts: dict = {}
    if method is not None:
        opts["method"] = method
    if min_length is not None:
        opts["min_length"] = min_length
    if min_snps is not None:
        opts["min_snps"] = min_snps
    if het_error_rate is not None:
        opts["het_error_rate"] = het_error_rate
    if hmm_min_af is not None:
        opts["hmm_min_af"] = hmm_min_af
    if variant_type is not None:
        opts["variant_type"] = variant_type
    cypher = "CALL graphpop.roh($chr, $pop, $options)"
    results = _run_procedure(cypher, {"chr": chr, "pop": pop, "options": opts})
    return json.dumps(results)


@mcp.tool()
def graphpop_garud_h(
    chr: str,
    pop: str,
    window: int = 100000,
    step: int = 50000,
    min_af: float | None = None,
    variant_type: str | None = None,
) -> str:
    """Compute Garud's H statistics (H1, H12, H2/H1, haplotype diversity) in sliding windows.

    FULL PATH: haplotype hashing, detects hard and soft sweeps.

    Args:
        chr: Chromosome ID.
        pop: Population name.
        window: Window size in bp (default 100,000).
        step: Step size in bp (default 50,000).
        min_af: Minimum allele frequency filter.
        variant_type: Filter by variant type (default "SNP").

    Returns JSON array of windows with: chr, start, end, h1, h12, h2_h1, hap_diversity, n_variants.
    """
    opts: dict = {}
    if min_af is not None:
        opts["min_af"] = min_af
    if variant_type is not None:
        opts["variant_type"] = variant_type
    cypher = "CALL graphpop.garud_h($chr, $pop, $window, $step, $options)"
    results = _run_procedure(cypher, {
        "chr": chr, "pop": pop, "window": window, "step": step, "options": opts,
    })
    return json.dumps(results)


@mcp.tool()
def graphpop_pop_summary(
    chr: str,
    pop: str,
    min_af: float | None = None,
    max_af: float | None = None,
    variant_type: str | None = None,
    consequence: str | None = None,
    pathway: str | None = None,
    samples: list[str] | None = None,
) -> str:
    """Compute whole-chromosome population summary and write to Population node.

    FAST PATH: aggregates diversity across all variants on a chromosome.

    Args:
        chr: Chromosome ID.
        pop: Population name.
        min_af: Minimum allele frequency filter.
        max_af: Maximum allele frequency filter.
        variant_type: Filter by variant type.
        consequence: VEP consequence filter (annotation-conditioned summary).
        pathway: Reactome pathway filter.
        samples: Custom sample list.

    Returns JSON with: population, chr, pi, theta_w, tajima_d, het_exp, het_obs, fis,
    fay_wu_h, n_variants, n_segregating.
    """
    opts = _build_options(
        min_af=min_af, max_af=max_af, variant_type=variant_type,
        consequence=consequence, pathway=pathway, samples=samples,
    )
    cypher = "CALL graphpop.pop_summary($chr, $pop, $options)"
    results = _run_procedure(cypher, {"chr": chr, "pop": pop, "options": opts})
    return json.dumps(results)


# ---------------------------------------------------------------------------
# Utility tools
# ---------------------------------------------------------------------------


@mcp.tool()
def graphpop_status() -> str:
    """Return database summary: node counts per label and available populations.

    No parameters required.
    """
    cypher = """
    CALL db.labels() YIELD label
    CALL db.stats.retrieve('GRAPH COUNTS') YIELD data
    WITH label, data
    CALL {
        WITH label
        MATCH (n)
        WHERE label IN labels(n)
        RETURN count(n) AS cnt
    }
    RETURN label, cnt ORDER BY cnt DESC
    """
    # Simpler approach: count key labels
    counts_cypher = """
    OPTIONAL MATCH (v:Variant) WITH count(v) AS variants
    OPTIONAL MATCH (s:Sample) WITH variants, count(s) AS samples
    OPTIONAL MATCH (p:Population) WITH variants, samples, count(p) AS populations
    OPTIONAL MATCH (g:Gene) WITH variants, samples, populations, count(g) AS genes
    OPTIONAL MATCH (pw:Pathway) WITH variants, samples, populations, genes, count(pw) AS pathways
    OPTIONAL MATCH (c:Chromosome) WITH variants, samples, populations, genes, pathways, count(c) AS chromosomes
    RETURN variants, samples, populations, genes, pathways, chromosomes
    """
    results = _run_procedure(counts_cypher)

    pops_cypher = "MATCH (p:Population) RETURN p.name AS name ORDER BY p.name"
    pops = _run_procedure(pops_cypher)

    return json.dumps({
        "counts": results[0] if results else {},
        "populations": [p["name"] for p in pops],
    })


@mcp.tool()
def graphpop_query(cypher: str, params: str | None = None) -> str:
    """Run an arbitrary Cypher query against the GraphPop database.

    Args:
        cypher: Cypher query string. Use $param placeholders for parameters.
        params: Optional JSON string of query parameters.

    Returns JSON array of result records.
    """
    parsed_params = None
    if params:
        parsed_params = json.loads(params)
    results = _run_procedure(cypher, parsed_params)
    return json.dumps(results)


def main():
    """Entry point for the graphpop-mcp command."""
    mcp.run()


if __name__ == "__main__":
    main()
