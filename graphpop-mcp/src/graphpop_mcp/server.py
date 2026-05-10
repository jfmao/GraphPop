"""MCP server wrapping GraphPop for AI agent access.

Exposes 12 population genomics procedures, multi-statistic integration tools
(converge, rank-genes), annotation lookup, and graph query capabilities.
Connection parameters come from environment variables (shared with graphpop CLI).
"""

from __future__ import annotations

import json
import os

from mcp.server.fastmcp import FastMCP

mcp = FastMCP(
    "GraphPop",
    instructions=(
        "GraphPop is a graph-native population genomics platform with 12 stored "
        "procedures inside a Neo4j graph database: diversity, divergence, SFS, "
        "joint SFS, genome scan, LD, iHS, XP-EHH, nSL, ROH, Garud's H, and "
        "population summary. Beyond procedures, it provides multi-statistic "
        "convergence detection (converge), gene ranking (rank_genes), annotation "
        "lookup (lookup_gene, lookup_pathway, lookup_variant, lookup_region), "
        "post-hoc filtering of persisted results (filter), and database inventory. "
        "All FAST PATH procedures support annotation-conditioned queries "
        "(consequence, pathway, gene). FULL PATH results can be filtered post-hoc "
        "via the filter tool."
    ),
)

_driver = None
_database = None


def _get_driver():
    """Lazy-initialize a shared Neo4j driver."""
    global _driver, _database
    if _driver is None:
        from neo4j import GraphDatabase

        # Use same env var names as graphpop CLI for consistency
        uri = os.environ.get("GRAPHPOP_URI", "bolt://localhost:7687")
        user = os.environ.get("GRAPHPOP_USER", "neo4j")
        password = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")
        _database = os.environ.get("GRAPHPOP_DATABASE", "neo4j")
        _driver = GraphDatabase.driver(uri, auth=(user, password))
    return _driver


def _run_procedure(cypher: str, params: dict | None = None) -> list[dict]:
    """Execute a Cypher procedure call and return results as a list of dicts."""
    driver = _get_driver()
    with driver.session(database=_database) as session:
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

    Returns JSON array of window results with: window_id, chr, start, end, population,
    n_variants, n_segregating, pi, theta_w, tajima_d, fst, fst_wc, dxy, pbs, fay_wu_h.
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

    Returns JSON array with per-variant: variantId, pos, af, ihs_unstd, ihs.
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

    Returns JSON array with per-variant: variantId, pos, af_pop1, af_pop2, xpehh_unstd, xpehh.
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

    Returns JSON array with per-variant: variantId, pos, af, nsl_unstd, nsl.
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

    Returns JSON array with per-sample aggregates: sampleId, n_roh, total_length, froh, mean_length, max_length.
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

    Returns JSON array of windows with: chr, start, end, population, h1, h12, h2_h1, hap_diversity, n_haplotypes, n_variants.
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

    Returns JSON with: pi, theta_w, tajima_d, fay_wu_h, mean_he, mean_ho, mean_fis,
    n_variants, n_segregating, n_polarized.
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


# ---------------------------------------------------------------------------
# High-level analytical tools (beyond raw procedures)
# ---------------------------------------------------------------------------


@mcp.tool()
def graphpop_converge(
    pop: str,
    stats: str = "ihs,xpehh,h12",
    thresholds: str = "2.0,2.0,0.3",
    chr: str | None = None,
    pop2: str | None = None,
    limit: int = 100,
) -> str:
    """Find genomic regions where multiple selection statistics simultaneously exceed thresholds.

    This is GraphPop's signature capability: querying persisted results across
    independently computed statistics to identify convergent signals.

    Args:
        pop: Population name.
        stats: Comma-separated statistic names (ihs, xpehh, nsl, h12, fst, pi, tajima_d).
        thresholds: Comma-separated threshold values (matched positionally to stats).
        chr: Restrict to a chromosome (optional).
        pop2: Second population for xpehh (optional).
        limit: Maximum results to return (default 100).

    Returns JSON array of convergent regions with gene annotations.
    """
    stat_list = [s.strip() for s in stats.split(",")]
    thresh_list = [float(t.strip()) for t in thresholds.split(",")]

    variant_stats = {"ihs", "xpehh", "nsl"}
    window_stats = {"h12", "fst", "pi", "tajima_d"}

    where_parts = []
    if chr:
        where_parts.append(f"v.chr = '{chr}'")

    for stat, thresh in zip(stat_list, thresh_list):
        if stat in variant_stats:
            if stat == "xpehh" and pop2:
                prop = f"xpehh_{pop}_{pop2}"
            else:
                prop = f"{stat}_{pop}"
            where_parts.append(f"abs(v.{prop}) >= {thresh}")

    if not where_parts:
        return json.dumps({"error": "No variant-level filters specified"})

    cypher = (
        f"MATCH (v:Variant) WHERE {' AND '.join(where_parts)} "
        f"OPTIONAL MATCH (v)-[:HAS_CONSEQUENCE]->(g:Gene) "
        f"RETURN v.variantId AS variant, v.pos AS pos, v.chr AS chr, "
        f"g.symbol AS gene "
        f"ORDER BY v.chr, v.pos LIMIT {limit}"
    )
    results = _run_procedure(cypher)
    return json.dumps(results)


@mcp.tool()
def graphpop_rank_genes(
    pop: str,
    chr: str | None = None,
    top: int = 50,
    pop2: str | None = None,
) -> str:
    """Rank genes by composite selection evidence.

    Combines max |iHS|, max |XP-EHH|, mean Fst, and HIGH-impact variant count
    per gene into a composite score.

    Args:
        pop: Population name.
        chr: Restrict to a chromosome (optional).
        top: Number of top genes to return (default 50).
        pop2: Second population for XP-EHH (optional).

    Returns JSON array of ranked genes with per-statistic scores.
    """
    ihs_prop = f"ihs_{pop}"
    where_chr = f"AND v.chr = '{chr}'" if chr else ""

    cypher = (
        f"MATCH (v:Variant)-[hc:HAS_CONSEQUENCE]->(g:Gene) "
        f"WHERE v.{ihs_prop} IS NOT NULL {where_chr} "
        f"WITH g, "
        f"max(abs(v.{ihs_prop})) AS max_ihs, "
        f"count(CASE WHEN hc.impact = 'HIGH' THEN 1 END) AS n_high, "
        f"count(v) AS n_variants "
        f"RETURN g.symbol AS gene, g.chr AS chr, "
        f"max_ihs, n_high, n_variants "
        f"ORDER BY max_ihs DESC LIMIT {top}"
    )
    results = _run_procedure(cypher)
    return json.dumps(results)


@mcp.tool()
def graphpop_lookup_gene(gene: str) -> str:
    """Look up a gene: variant count, consequence types, pathways, and selection statistics.

    Args:
        gene: Gene symbol or ID (e.g. "KCNE1", "GW5", "LOC_Os05g09520").

    Returns JSON with gene details, variant summary, pathway memberships, and
    any persisted selection statistics.
    """
    cypher = (
        f"MATCH (g:Gene) WHERE g.symbol = '{gene}' OR g.geneId = '{gene}' "
        f"OPTIONAL MATCH (g)<-[hc:HAS_CONSEQUENCE]-(v:Variant) "
        f"WITH g, count(v) AS n_variants, "
        f"collect(DISTINCT hc.impact) AS impacts, "
        f"collect(DISTINCT hc.consequence) AS consequences "
        f"OPTIONAL MATCH (g)-[:IN_PATHWAY]->(p:Pathway) "
        f"RETURN g.symbol AS symbol, g.geneId AS geneId, "
        f"g.chr AS chr, g.start AS start, g.end AS end, "
        f"n_variants, impacts, consequences, "
        f"collect(p.name) AS pathways"
    )
    results = _run_procedure(cypher)
    return json.dumps(results)


@mcp.tool()
def graphpop_lookup_pathway(pathway: str) -> str:
    """Look up a pathway: member genes and variant counts.

    Args:
        pathway: Pathway name or partial name (case-insensitive search).

    Returns JSON with pathway details and member gene list.
    """
    cypher = (
        f"MATCH (p:Pathway) WHERE toLower(p.name) CONTAINS toLower('{pathway}') "
        f"OPTIONAL MATCH (p)<-[:IN_PATHWAY]-(g:Gene)<-[:HAS_CONSEQUENCE]-(v:Variant) "
        f"WITH p, g, count(v) AS n_variants "
        f"RETURN p.name AS pathway, p.pathwayId AS pathwayId, "
        f"collect({{gene: g.symbol, n_variants: n_variants}}) AS genes "
        f"ORDER BY p.name"
    )
    results = _run_procedure(cypher)
    return json.dumps(results)


@mcp.tool()
def graphpop_lookup_region(chr: str, start: int, end: int) -> str:
    """Look up genes and summary statistics in a genomic region.

    Args:
        chr: Chromosome ID.
        start: Start position.
        end: End position.

    Returns JSON with genes in the region and their variant/consequence summary.
    """
    cypher = (
        f"MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene) "
        f"WHERE v.chr = '{chr}' AND v.pos >= {start} AND v.pos <= {end} "
        f"WITH g, count(v) AS n_variants, "
        f"collect(DISTINCT v.variantId) AS sample_variants "
        f"RETURN g.symbol AS gene, g.chr AS chr, g.start AS start, "
        f"g.end AS end, n_variants "
        f"ORDER BY n_variants DESC"
    )
    results = _run_procedure(cypher)
    return json.dumps(results)


@mcp.tool()
def graphpop_filter(
    statistic: str,
    chr: str,
    pop: str,
    consequence: str | None = None,
    pathway: str | None = None,
    gene: str | None = None,
    min_score: float | None = None,
    pop2: str | None = None,
    limit: int = 1000,
) -> str:
    """Query persisted statistics (iHS, XP-EHH, nSL) filtered by annotation.

    Use this for conditioned analysis on FULL PATH statistics. These statistics
    must be computed genome-wide first (via graphpop_ihs etc.), then filtered
    by annotation post-hoc.

    Args:
        statistic: One of: ihs, xpehh, nsl.
        chr: Chromosome ID.
        pop: Population name.
        consequence: VEP consequence filter (e.g. "missense_variant").
        pathway: Pathway name filter.
        gene: Gene name filter.
        min_score: Minimum absolute score.
        pop2: Second population for xpehh.
        limit: Maximum results (default 1000).

    Returns JSON array of filtered variants with scores and annotations.
    """
    if statistic == "xpehh" and pop2:
        prop = f"xpehh_{pop}_{pop2}"
    else:
        prop = f"{statistic}_{pop}"

    match_clause = "MATCH (v:Variant)"
    where_parts = [f"v.chr = '{chr}'", f"v.{prop} IS NOT NULL"]

    if consequence:
        match_clause += "-[:HAS_CONSEQUENCE]->(hc)"
        where_parts.append(f"hc.consequence = '{consequence}'")
    if pathway:
        match_clause += "-[:HAS_CONSEQUENCE]->(:Gene)-[:IN_PATHWAY]->(pw:Pathway)"
        where_parts.append(f"pw.name CONTAINS '{pathway}'")
    if gene:
        match_clause += "-[:HAS_CONSEQUENCE]->(g:Gene)"
        where_parts.append(f"(g.symbol = '{gene}' OR g.geneId = '{gene}')")
    if min_score is not None:
        where_parts.append(f"abs(v.{prop}) >= {min_score}")

    return_cols = [f"v.variantId AS variant", "v.pos AS pos",
                   f"v.{prop} AS {statistic}"]
    if consequence:
        return_cols.append("hc.consequence AS consequence")
    if gene:
        return_cols.append("g.symbol AS gene")

    cypher = (
        f"{match_clause} WHERE {' AND '.join(where_parts)} "
        f"RETURN DISTINCT {', '.join(return_cols)} "
        f"ORDER BY v.pos LIMIT {limit}"
    )
    results = _run_procedure(cypher)
    return json.dumps(results)


@mcp.tool()
def graphpop_inventory() -> str:
    """Show comprehensive database inventory: what has been computed.

    Returns populations, chromosomes, annotation types loaded, and which
    selection statistics have been persisted on Variant/GenomicWindow nodes.
    No parameters required.
    """
    inventory = {}

    # Populations
    pops = _run_procedure(
        "MATCH (p:Population) RETURN p.populationId AS pop, "
        "p.n_samples AS n ORDER BY n DESC"
    )
    inventory["populations"] = pops

    # Chromosomes
    chrs = _run_procedure(
        "MATCH (c:Chromosome) RETURN c.chromosomeId AS chr, "
        "c.length AS length ORDER BY chr"
    )
    inventory["chromosomes"] = chrs

    # Node counts
    for label in ["Variant", "Sample", "Gene", "Pathway", "GOTerm", "GenomicWindow"]:
        result = _run_procedure(f"MATCH (n:{label}) RETURN count(n) AS cnt")
        inventory[f"n_{label.lower()}"] = result[0]["cnt"] if result else 0

    # Persisted statistics (check a sample variant)
    props = _run_procedure(
        "MATCH (v:Variant) WITH v LIMIT 1 RETURN keys(v) AS props"
    )
    if props:
        all_props = props[0]["props"]
        ihs_props = [p for p in all_props if p.startswith("ihs_")]
        xpehh_props = [p for p in all_props if p.startswith("xpehh_")]
        nsl_props = [p for p in all_props if p.startswith("nsl_")]
        inventory["persisted_ihs"] = ihs_props
        inventory["persisted_xpehh"] = xpehh_props
        inventory["persisted_nsl"] = nsl_props
        inventory["has_ancestral_allele"] = "ancestral_allele" in all_props

    return json.dumps(inventory, indent=2)


@mcp.tool()
def graphpop_kinship_king(
    chr: str,
    pop: str,
    start: int | None = None,
    end: int | None = None,
    min_snp: int | None = None,
    min_phi: float | None = None,
    include_self: bool = False,
    samples: list[str] | None = None,
) -> str:
    """KING-robust pairwise kinship (Manichaikul 2010) for samples in a population.

    Validation baseline of the GraphPop kinship suite — agrees with PLINK2
    --make-king to r >= 0.999. For graph-native, ARG-based kinship (the
    branch-GRM family that closes the open gaps in the literature) use
    graphpop_kinship_branch_grm once the M4.A ARG layer is in place.

    Args:
        chr: Chromosome ID.
        pop: Population name.
        start: Start position in bp (default 1).
        end: End position in bp (default end-of-chromosome).
        min_snp: Skip pairs with fewer informative variants (default 1000).
        min_phi: Skip pairs with kinship below this threshold (default -0.5).
        include_self: Emit (s, s) self-pair rows (default False).
        samples: Restrict to a subset of sample IDs (overrides ``pop``).

    Returns JSON array with one row per pair: sample_a, sample_b, phi, ibs0,
    het_het, n_snp, n_aa_min, method.
    """
    opts: dict = {}
    if start is not None:
        opts["start"] = start
    if end is not None:
        opts["end"] = end
    if min_snp is not None:
        opts["min_snp"] = min_snp
    if min_phi is not None:
        opts["min_phi"] = min_phi
    if include_self:
        opts["include_self"] = True
    if samples:
        opts["samples"] = samples
    cypher = "CALL graphpop.kinship.king($chr, $pop, $options)"
    results = _run_procedure(
        cypher, {"chr": chr, "pop": pop, "options": opts}
    )
    return json.dumps(results)


@mcp.tool()
def graphpop_arg_ingest(
    tree_sequence_path: str,
    run_id: str,
    source: str = "tsinfer",
    params: dict | None = None,
    sample_id_map: dict[int, str] | None = None,
) -> str:
    """Ingest a tskit TreeSequence into the GraphPop ARG layer (M4.A).

    Additive: leaves existing :Variant, :Sample, :Population nodes
    untouched. Creates one :ARGRun + many :TreeNode + :PARENT_OF +
    :REPRESENTS + :MUTATED_ON in a single online-Cypher pass. The agent
    must specify a server-side path to a .trees file (paths are not
    transmitted as content over MCP).

    Args:
        tree_sequence_path: Server-side filesystem path to the .trees file.
        run_id: Globally unique identifier for this run.
        source: Inferer that produced the ARG ("tsinfer", "tsdate",
            "singer", "relate", "msprime").
        params: Free-form parameter dict serialised onto the :ARGRun node.
        sample_id_map: Optional dict mapping tskit sample index (int) to
            existing :Sample.sampleId. Defaults to identity ("sample_{i}").

    Returns JSON of the resulting ARGRunSummary.
    """
    try:
        import tskit  # type: ignore
    except ImportError:
        return json.dumps({
            "error": "tskit is not installed in the MCP server environment; "
                     "install graphpop-mcp[arg]"
        })
    from graphpop_import.arg_importer import ARGIngester

    treeseq = tskit.load(tree_sequence_path)
    ingester = ARGIngester(_get_driver())
    sample_id_callable = (
        (lambda i: sample_id_map[i]) if sample_id_map else None
    )
    summary = ingester.ingest(
        treeseq,
        run_id=run_id,
        source=source,
        params=params,
        sample_id_map=sample_id_callable,
    )
    return json.dumps({
        "run_id": summary.run_id,
        "source": summary.source,
        "n_samples": summary.n_samples,
        "sequence_length": summary.sequence_length,
        "n_trees": summary.n_trees,
        "n_nodes": summary.n_nodes,
        "n_edges": summary.n_edges,
        "n_mutations": summary.n_mutations,
        "created_at": summary.created_at.isoformat(),
    })


@mcp.tool()
def graphpop_arg_list_runs() -> str:
    """List every :ARGRun currently stored in the database.

    Returns JSON array of run summaries (run_id, source, n_samples,
    sequence_length, n_trees, n_nodes, n_edges, n_mutations, created_at).
    """
    from graphpop_import.arg_importer import ARGIngester

    ingester = ARGIngester(_get_driver())
    runs = ingester.list_runs()
    return json.dumps([
        {
            "run_id": r.run_id,
            "source": r.source,
            "n_samples": r.n_samples,
            "sequence_length": r.sequence_length,
            "n_trees": r.n_trees,
            "n_nodes": r.n_nodes,
            "n_edges": r.n_edges,
            "n_mutations": r.n_mutations,
            "created_at": r.created_at.isoformat() if r.created_at else None,
        }
        for r in runs
    ])


@mcp.tool()
def graphpop_arg_delete_run(run_id: str) -> str:
    """Delete an :ARGRun and every :TreeNode + relationship that references it.

    Other runs in the database are left untouched. Returns JSON
    {"run_id": ..., "n_tree_nodes_deleted": N}.
    """
    from graphpop_import.arg_importer import ARGIngester

    ingester = ARGIngester(_get_driver())
    n = ingester.delete_run(run_id)
    return json.dumps({"run_id": run_id, "n_tree_nodes_deleted": n})


@mcp.tool()
def graphpop_kinship_branch_grm(
    run_id: str,
    start: int | None = None,
    end: int | None = None,
    include_self: bool = True,
    restrict_to_pathway: str | None = None,
    mutation_filter: str | None = None,
    time_window: list[float] | None = None,
) -> str:
    """Branch GRM via ARG traversal (Fan, Mancuso & Chiang 2022).

    Unconditional mode validates against the ``egrm.varGRM`` Python
    reference to relative error < 1e-6. Conditional predicates close
    open ARG-literature gaps:

    * ``restrict_to_pathway`` (G2) — branches counted only if they
      carry mutations in the named :Pathway.
    * ``mutation_filter`` (G2) — branches counted only if they carry
      mutations with the named :HAS_CONSEQUENCE.consequence.
    * ``time_window`` — geometric branch-time truncation to
      [t_lo, t_hi] in generations (lineage-time stratification).

    All three compose multiplicatively. Requires :ARGRun(run_id) to be
    ingested (see graphpop_arg_ingest).

    Args:
        run_id: :ARGRun.runId identifying which inferred ARG to use.
        start: Region start in bp (default 0).
        end: Region end in bp (default end-of-ARG-sequence).
        include_self: Emit (s, s) self-pair rows (default True).
        restrict_to_pathway: :Pathway.pathwayId to restrict to.
        mutation_filter: HAS_CONSEQUENCE.consequence string to restrict to.
        time_window: Two-element [t_lo, t_hi] in generations.

    Returns JSON array with one row per pair: sample_a, sample_b, phi
    (= B_ij), n_snp (= n_trees that contributed), method ("branch_grm").
    """
    opts: dict = {}
    if start is not None:
        opts["start"] = start
    if end is not None:
        opts["end"] = end
    if not include_self:
        opts["include_self"] = False
    if restrict_to_pathway:
        opts["restrict_to_pathway"] = restrict_to_pathway
    if mutation_filter:
        opts["mutation_filter"] = mutation_filter
    if time_window is not None:
        if len(time_window) != 2:
            return json.dumps({"error": "time_window must be [t_lo, t_hi]"})
        opts["time_window"] = list(time_window)
    cypher = "CALL graphpop.kinship.branch_grm($run_id, $options)"
    results = _run_procedure(
        cypher, {"run_id": run_id, "options": opts}
    )
    return json.dumps(results)


@mcp.tool()
def graphpop_kinship_branch_grm_posterior(
    posterior_run_ids: list[str],
    start: int | None = None,
    end: int | None = None,
    include_self: bool = True,
    restrict_to_pathway: str | None = None,
    mutation_filter: str | None = None,
    time_window: list[float] | None = None,
) -> str:
    """Posterior-aware branch GRM over multiple :ARGRun (closes G1 in full).

    Aggregates element-wise per-run branch_grm via Welford. Returns
    posterior mean and standard error of the mean. Conditional
    predicates (restrict_to_pathway, mutation_filter, time_window)
    compose unchanged — combine with restrict_to_pathway to get
    "posterior-mean kinship through pathway X with credible interval"
    in a single query (the GraphPop v2 paper headline).

    Args:
        posterior_run_ids: List of :ARGRun.runId values from the
            posterior set (e.g. SINGER's M=100 samples).
        start, end: Region bounds (bp); same across all runs.
        include_self: Emit (s, s) self-pair rows (default True).
        restrict_to_pathway: :Pathway.pathwayId (G2).
        mutation_filter: HAS_CONSEQUENCE.consequence string (G2).
        time_window: Two-element [t_lo, t_hi] in generations.

    Returns JSON array with one row per pair:
    sample_a, sample_b, b_ij_mean, b_ij_sd (NaN when n_runs < 2),
    n_runs, method.
    """
    opts: dict = {}
    if start is not None:
        opts["start"] = start
    if end is not None:
        opts["end"] = end
    if not include_self:
        opts["include_self"] = False
    if restrict_to_pathway:
        opts["restrict_to_pathway"] = restrict_to_pathway
    if mutation_filter:
        opts["mutation_filter"] = mutation_filter
    if time_window is not None:
        if len(time_window) != 2:
            return json.dumps({"error": "time_window must be [t_lo, t_hi]"})
        opts["time_window"] = list(time_window)
    cypher = ("CALL graphpop.kinship.branch_grm_posterior("
              "$posterior_run_ids, $options)")
    results = _run_procedure(
        cypher,
        {"posterior_run_ids": list(posterior_run_ids), "options": opts}
    )
    return json.dumps(results)


@mcp.tool()
def graphpop_ancestry_ingest(
    run_id: str,
    painting: list[dict],
    painter: str = "manual",
    replace: bool = False,
) -> str:
    """Ingest per-TreeNode ancestry painting (M4.B).

    Each painting row is a dict with keys ``tskit_node_id`` (int),
    ``population_id`` (str), and optional ``posterior_prob`` (float,
    default 1.0). Multiple rows per node are allowed for probabilistic
    painting (probs should sum to 1).

    Args:
        run_id: :ARGRun.runId to attach the painting to.
        painting: List of {tskit_node_id, population_id, posterior_prob}.
        painter: Identifier of the painting source.
        replace: Delete prior painting for (run_id, painter) before insert.

    Returns JSON with the ingest summary.
    """
    from graphpop_import.ancestry_ingester import (
        AncestryIngester, PaintingRow,
    )

    rows = [
        PaintingRow(
            tskit_node_id=int(p["tskit_node_id"]),
            population_id=str(p["population_id"]),
            posterior_prob=float(p.get("posterior_prob", 1.0)),
        )
        for p in painting
    ]
    ing = AncestryIngester(_get_driver())
    summary = ing.ingest(run_id, rows, painter=painter, replace=replace)
    return json.dumps({
        "run_id": summary.run_id,
        "painter": summary.painter,
        "n_painted_nodes": summary.n_painted_nodes,
        "n_edges": summary.n_edges,
        "n_populations": summary.n_populations,
        "created_at": summary.created_at.isoformat(),
    })


@mcp.tool()
def graphpop_ancestry_ingest_from_samples(
    run_id: str,
    sample_ancestry: dict[str, str],
    painter: str = "majority_vote",
    replace: bool = False,
) -> str:
    """Propagate sample-level ancestry to every TreeNode via majority-vote DFS.

    Args:
        run_id: :ARGRun.runId to paint.
        sample_ancestry: dict mapping :Sample.sampleId -> population_id.
        painter: Identifier of the painting source.
        replace: Delete prior painting for (run_id, painter) before insert.

    Returns JSON with the ingest summary.
    """
    from graphpop_import.ancestry_ingester import AncestryIngester

    ing = AncestryIngester(_get_driver())
    summary = ing.ingest_from_samples(
        run_id, sample_ancestry, painter=painter, replace=replace)
    return json.dumps({
        "run_id": summary.run_id,
        "painter": summary.painter,
        "n_painted_nodes": summary.n_painted_nodes,
        "n_edges": summary.n_edges,
        "n_populations": summary.n_populations,
        "created_at": summary.created_at.isoformat(),
    })


@mcp.tool()
def graphpop_ancestry_list(run_id: str | None = None) -> str:
    """List paintings in the database, optionally filtered by run_id.

    Returns JSON array of {run_id, painter, n_painted_nodes, n_edges,
    n_populations}.
    """
    from graphpop_import.ancestry_ingester import AncestryIngester

    ing = AncestryIngester(_get_driver())
    paintings = ing.list_paintings(run_id=run_id)
    return json.dumps([
        {
            "run_id": p.run_id,
            "painter": p.painter,
            "n_painted_nodes": p.n_painted_nodes,
            "n_edges": p.n_edges,
            "n_populations": p.n_populations,
        }
        for p in paintings
    ])


@mcp.tool()
def graphpop_ancestry_delete(run_id: str, painter: str) -> str:
    """Delete every :HAS_ANCESTRY edge for a (run_id, painter) pair.

    Returns JSON {run_id, painter, n_edges_deleted}.
    """
    from graphpop_import.ancestry_ingester import AncestryIngester

    ing = AncestryIngester(_get_driver())
    n = ing.delete_painting(run_id, painter)
    return json.dumps({"run_id": run_id, "painter": painter,
                        "n_edges_deleted": n})


@mcp.tool()
def graphpop_kinship_branch_grm_by_ancestry(
    run_id: str,
    painter: str | None = None,
    start: int | None = None,
    end: int | None = None,
    include_self: bool = True,
    restrict_to_pathway: str | None = None,
    mutation_filter: str | None = None,
    time_window: list[float] | None = None,
) -> str:
    """Ancestry-decomposed branch GRM (closes G3).

    Returns one row per (sample_a, sample_b, ancestry) triple. Sum
    across ancestries reproduces the unconditional branch_grm matrix
    (modulo unpainted nodes). Conditional predicates compose
    unchanged. Requires :HAS_ANCESTRY edges ingested via the
    graphpop_ancestry_ingest* tools.

    Args:
        run_id: :ARGRun.runId to use.
        painter: Painter id (e.g. "majority_vote", "rfmix") when the
            run has multiple paintings; omit for any.
        start, end: Region bounds (bp).
        include_self: Emit (s, s) rows (default True).
        restrict_to_pathway: G2 predicate.
        mutation_filter: G2 predicate.
        time_window: G1 lineage-time slice [t_lo, t_hi].

    Returns JSON array of rows: sample_a, sample_b, ancestry,
    b_ij_component, n_branches, method.
    """
    opts: dict = {}
    if painter:
        opts["painter"] = painter
    if start is not None:
        opts["start"] = start
    if end is not None:
        opts["end"] = end
    if not include_self:
        opts["include_self"] = False
    if restrict_to_pathway:
        opts["restrict_to_pathway"] = restrict_to_pathway
    if mutation_filter:
        opts["mutation_filter"] = mutation_filter
    if time_window is not None:
        if len(time_window) != 2:
            return json.dumps({"error": "time_window must be [t_lo, t_hi]"})
        opts["time_window"] = list(time_window)
    cypher = "CALL graphpop.kinship.branch_grm_by_ancestry($run_id, $options)"
    results = _run_procedure(
        cypher, {"run_id": run_id, "options": opts}
    )
    return json.dumps(results)


@mcp.tool()
def graphpop_kinship_branch_grm_apply(
    run_id: str,
    vector: list[float],
    start: int | None = None,
    end: int | None = None,
    restrict_to_pathway: str | None = None,
    mutation_filter: str | None = None,
    time_window: list[float] | None = None,
) -> str:
    """Apply branch GRM to a vector via Algorithm V (G·v, biobank-scale).

    Returns G · v without materialising G, in O(T · N log N) per
    multiplication. Conditional predicates from branch_grm compose
    unchanged. Vector length must equal the number of samples in the
    :ARGRun.

    Args:
        run_id: :ARGRun.runId.
        vector: Length-n_samples list of floats (haplotype-ordered).
        start, end: Region bounds (bp).
        restrict_to_pathway: G2 predicate.
        mutation_filter: G2 predicate.
        time_window: [t_lo, t_hi] in generations.

    Returns JSON array: sample_id, col, value, n_branches, method.
    """
    opts: dict = {}
    if start is not None:
        opts["start"] = start
    if end is not None:
        opts["end"] = end
    if restrict_to_pathway:
        opts["restrict_to_pathway"] = restrict_to_pathway
    if mutation_filter:
        opts["mutation_filter"] = mutation_filter
    if time_window is not None:
        opts["time_window"] = list(time_window)
    cypher = ("CALL graphpop.kinship.branch_grm_apply"
              "($run_id, $vector, $options)")
    results = _run_procedure(
        cypher,
        {"run_id": run_id, "vector": list(vector), "options": opts}
    )
    return json.dumps(results)


@mcp.tool()
def graphpop_kinship_ibs(
    chr: str,
    pop: str,
    start: int | None = None,
    end: int | None = None,
    min_snp: int | None = None,
    min_ibs: float | None = None,
    include_self: bool = False,
    samples: list[str] | None = None,
) -> str:
    """Identity-by-state pairwise statistic on packed genotypes (M4.2).

    Matrix-only sibling of kinship.king; no ARG required. IBS is in
    [0, 1]; identical samples = 1, opposite homozygotes average toward 0.

    Args:
        chr: Chromosome ID.
        pop: Population name.
        start, end: Region bounds (bp).
        min_snp: Skip pairs with fewer informative variants (default 1000).
        min_ibs: Skip pairs with IBS below this threshold (default 0.0).
        include_self: Emit (s, s) self-pair rows (default False).
        samples: Restrict to a subset of sample IDs.

    Returns JSON array: sample_a, sample_b, phi (= IBS), ibs0,
    het_het (= ibs2), n_snp (= n variants used), n_aa_min, method.
    """
    opts: dict = {}
    if start is not None:
        opts["start"] = start
    if end is not None:
        opts["end"] = end
    if min_snp is not None:
        opts["min_snp"] = min_snp
    if min_ibs is not None:
        opts["min_ibs"] = min_ibs
    if include_self:
        opts["include_self"] = True
    if samples:
        opts["samples"] = samples
    cypher = "CALL graphpop.kinship.ibs($chr, $pop, $options)"
    results = _run_procedure(
        cypher, {"chr": chr, "pop": pop, "options": opts}
    )
    return json.dumps(results)


@mcp.tool()
def graphpop_ibd_ingest(
    segments: list[dict],
    source: str = "hap_ibd",
    replace: bool = False,
) -> str:
    """Ingest external IBD-caller segments (M4.3).

    Each segment is a dict with keys ``sample_a``, ``sample_b``,
    ``chr``, ``start``, ``end``, optional ``length_cM``.

    Returns JSON with the ingest summary.
    """
    from graphpop_import.ibd_ingester import IBDIngester, IBDSegmentRow
    rows = [
        IBDSegmentRow(
            sample_a=str(s["sample_a"]),
            sample_b=str(s["sample_b"]),
            chr=str(s["chr"]),
            start=int(s["start"]),
            end=int(s["end"]),
            length_cM=(float(s["length_cM"]) if s.get("length_cM") is not None
                       else None),
        )
        for s in segments
    ]
    ing = IBDIngester(_get_driver())
    summary = ing.ingest(rows, source=source, replace=replace)
    return json.dumps({
        "source": summary.source,
        "n_segments": summary.n_segments,
        "n_pairs": summary.n_pairs,
        "created_at": summary.created_at.isoformat(),
    })


@mcp.tool()
def graphpop_ibd_from_arg(
    run_id: str,
    max_tmrca: float | None = None,
    min_length_bp: int = 0,
    chr: str = "chr1",
) -> str:
    """ARG-derived IBD segments (M4.3).

    Walks :PARENT_OF to extract maximal IBD segments per pair under
    an optional TMRCA cap. Equivalent to tskit's ibd_segments.

    Returns JSON array of segments.
    """
    opts: dict = {"chr": chr}
    if max_tmrca is not None:
        opts["max_tmrca"] = max_tmrca
    if min_length_bp:
        opts["min_length_bp"] = min_length_bp
    cypher = "CALL graphpop.ibd.from_arg($run_id, $options)"
    results = _run_procedure(
        cypher, {"run_id": run_id, "options": opts})
    return json.dumps(results)


@mcp.tool()
def graphpop_ibd_kinship(
    source: str,
    min_length_bp: int = 0,
    total_genome_cM: float | None = None,
) -> str:
    """Browning-style kinship from total IBD length (M4.3).

    Returns JSON array of pairs.
    """
    opts: dict = {}
    if min_length_bp:
        opts["min_length_bp"] = min_length_bp
    if total_genome_cM is not None:
        opts["total_genome_cM"] = total_genome_cM
    cypher = "CALL graphpop.ibd.kinship($source, $options)"
    results = _run_procedure(
        cypher, {"source": source, "options": opts})
    return json.dumps(results)


@mcp.tool()
def graphpop_ibd_list() -> str:
    """List every IBD source currently in the database."""
    from graphpop_import.ibd_ingester import IBDIngester
    ing = IBDIngester(_get_driver())
    sources = ing.list_sources()
    return json.dumps([
        {"source": s.source, "n_edges": s.n_edges, "n_pairs": s.n_pairs}
        for s in sources
    ])


@mcp.tool()
def graphpop_ibd_delete(source: str) -> str:
    """Delete every :IBD_SEGMENT edge for the given source."""
    from graphpop_import.ibd_ingester import IBDIngester
    ing = IBDIngester(_get_driver())
    n = ing.delete_source(source)
    return json.dumps({"source": source, "n_edges_deleted": n})


@mcp.tool()
def graphpop_kinship_branch_grm_pca(
    run_id: str,
    k: int,
    n_iter: int = 0,
    seed: int = 42,
    restrict_to_pathway: str | None = None,
    mutation_filter: str | None = None,
    time_window: list[float] | None = None,
) -> str:
    """Top-K PCs of branch GRM via Lanczos (M4.4).

    Tractable at biobank scale (no full matrix materialised).
    Conditional predicates compose. Returns one row per (sample, pc).
    """
    opts: dict = {"seed": seed}
    if n_iter > 0:
        opts["n_iter"] = n_iter
    if restrict_to_pathway:
        opts["restrict_to_pathway"] = restrict_to_pathway
    if mutation_filter:
        opts["mutation_filter"] = mutation_filter
    if time_window is not None:
        opts["time_window"] = list(time_window)
    cypher = ("CALL graphpop.kinship.branch_grm_pca("
              "$run_id, $k, $options)")
    return json.dumps(_run_procedure(
        cypher, {"run_id": run_id, "k": k, "options": opts}))


@mcp.tool()
def graphpop_kinship_branch_grm_he(
    run_id: str,
    phenotype: list[float],
    n_hutchinson: int = 50,
    seed: int = 42,
    restrict_to_pathway: str | None = None,
    mutation_filter: str | None = None,
    time_window: list[float] | None = None,
) -> str:
    """Haseman-Elston heritability via Algorithm V + Hutchinson trace.

    Single-component RHE-mc (Pazokitoroudi et al. 2020).
    Conditional predicates compose. Returns h2, se, num, tr_g_sq.
    """
    opts: dict = {"n_hutchinson": n_hutchinson, "seed": seed}
    if restrict_to_pathway:
        opts["restrict_to_pathway"] = restrict_to_pathway
    if mutation_filter:
        opts["mutation_filter"] = mutation_filter
    if time_window is not None:
        opts["time_window"] = list(time_window)
    cypher = ("CALL graphpop.kinship.branch_grm_he("
              "$run_id, $phenotype, $options)")
    return json.dumps(_run_procedure(
        cypher, {"run_id": run_id, "phenotype": list(phenotype),
                 "options": opts}))


@mcp.tool()
def graphpop_sim_msprime(
    config_path: str,
    output_path: str,
    run_id: str | None = None,
    ingest: bool = False,
) -> str:
    """Run msprime from a YAML config (M12.A).

    Dumps the tree sequence to output_path. With ingest=true and a
    run_id, also persists the result as a :ARGRun via the existing
    ARGIngester. Returns JSON with simulation metadata.
    """
    try:
        from graphpop_sim import MsprimeConfig, MsprimeRunner
    except ImportError:
        return json.dumps({
            "error": "graphpop-sim not installed; "
                     "run `pip install graphpop-sim`",
        })
    cfg = MsprimeConfig.from_yaml(config_path)
    ts = MsprimeRunner(cfg).simulate_and_dump(output_path)
    out = {
        "config_path": config_path,
        "output_path": output_path,
        "n_samples": int(ts.num_samples),
        "n_trees": int(ts.num_trees),
        "n_mutations": int(ts.num_mutations),
        "sequence_length": int(ts.sequence_length),
        "ingested": False,
    }
    if ingest:
        try:
            from graphpop_import.arg_importer import ARGIngester
        except ImportError:
            out["error"] = "graphpop-import not installed; cannot ingest"
            return json.dumps(out)
        if not run_id:
            out["error"] = "run_id is required when ingest=true"
            return json.dumps(out)
        driver = _get_driver()
        ing = ARGIngester(driver)
        summary = ing.ingest(ts, run_id=run_id, source="msprime",
                              params=cfg.to_dict())
        out["ingested"] = True
        out["run_id"] = summary.run_id
        out["n_tree_nodes"] = summary.n_nodes
        out["n_edges"] = summary.n_edges
    return json.dumps(out)


@mcp.tool()
def graphpop_sim_slim(
    template: str,
    params: dict[str, float | int | str],
    output_path: str,
    run_id: str | None = None,
    ingest: bool = False,
) -> str:
    """Run a curated SLiM template (M12.B).

    template ∈ {"neutral", "sweep_genic", "bottleneck"}. Passes
    each (key, value) in params to SLiM as a global variable.
    Requires the SLiM binary on PATH. Returns JSON with run
    metadata.
    """
    try:
        from graphpop_sim import SlimRunner
    except ImportError:
        return json.dumps({
            "error": "graphpop-sim not installed; "
                     "run `pip install graphpop-sim`",
        })
    if not SlimRunner.is_available():
        return json.dumps({
            "error": "SLiM binary not on PATH",
        })
    runner = SlimRunner()
    result = runner.run(template=template, params=dict(params),
                         output_path=output_path)
    out = {
        "template": template,
        "output_path": output_path,
        "return_code": result.return_code,
        "ingested": False,
    }
    if ingest and result.return_code == 0:
        try:
            import tskit
            from graphpop_import.arg_importer import ARGIngester
        except ImportError:
            out["error"] = "graphpop-import + tskit required for ingest"
            return json.dumps(out)
        if not run_id:
            out["error"] = "run_id is required when ingest=true"
            return json.dumps(out)
        ts = tskit.load(output_path)
        driver = _get_driver()
        ing = ARGIngester(driver)
        summary = ing.ingest(ts, run_id=run_id, source="slim",
                              params={"template": template,
                                      "params": dict(params)})
        out["ingested"] = True
        out["run_id"] = summary.run_id
    return json.dumps(out)


@mcp.tool()
def graphpop_sim_abc(
    prior_path: str,
    observed_path: str,
    output_path: str,
    base_config_path: str | None = None,
    n_sim: int = 1000,
    epsilon: float = 0.05,
    summary_keys: str = "pi,theta_w,tajimas_d",
    seed: int = 42,
) -> str:
    """Rejection-ABC for demographic / mutational-rate inference (M12.C).

    Reads priors (YAML) + observed summary stats (TSV) + optional
    base msprime config; simulates n_sim draws; accepts top
    epsilon fraction by Euclidean distance on summary_keys; writes
    posterior TSV to output_path.
    """
    try:
        from graphpop_sim import (
            MsprimeConfig, load_priors, posterior_to_tsv, run_abc,
        )
        from graphpop_sim.cli import _load_observed
    except ImportError:
        return json.dumps({
            "error": "graphpop-sim not installed; "
                     "run `pip install graphpop-sim`",
        })
    priors = load_priors(prior_path)
    obs = _load_observed(observed_path)
    keys = [s.strip() for s in summary_keys.split(",") if s.strip()]
    base = (MsprimeConfig.from_yaml(base_config_path)
            if base_config_path else None)
    result = run_abc(
        priors=priors, observed=obs, n_sim=n_sim, epsilon=epsilon,
        summary_keys=keys, base_config=base, seed=seed,
    )
    n = posterior_to_tsv(result, output_path)
    return json.dumps({
        "n_total": result.n_total,
        "n_accepted": result.n_accepted,
        "epsilon": result.epsilon,
        "output_path": output_path,
        "n_rows_written": n,
    })


@mcp.tool()
def graphpop_recombination_ldhat_mcmc(
    sample_ids: list[str],
    window_size: int = 10_000,
    step: int | None = None,
    max_pair_distance: int = 5_000,
    min_maf: float = 0.05,
    n_iter: int = 5_000,
    burn_in: int = 1_000,
    prop_sd: float = 0.5,
    sigma: float = 0.1,
    seed: int = 42,
) -> str:
    """Per-window LDhat-style MCMC posterior on ρ (M13.A).

    Metropolis-Hastings sampler on log10(ρ) with Hudson 1985
    Gaussian-approx likelihood. Deterministic given seed.
    Returns JSON array with posterior_mean / 2.5 / 97.5 quantiles.
    """
    opts: dict = {
        "window_size": window_size,
        "max_pair_distance": max_pair_distance,
        "min_maf": min_maf,
        "n_iter": n_iter,
        "burn_in": burn_in,
        "prop_sd": prop_sd,
        "sigma": sigma,
        "seed": seed,
    }
    if step is not None:
        opts["step"] = step
    cypher = "CALL graphpop.recombination.ldhat_mcmc($sids, $options)"
    return json.dumps(_run_procedure(
        cypher, {"sids": list(sample_ids), "options": opts}))


@mcp.tool()
def graphpop_recombination_hmm_smooth(
    sample_ids: list[str],
    window_size: int = 10_000,
    step: int | None = None,
    max_pair_distance: int = 5_000,
    n_states: int = 20,
    state_log_lo: float = -10.0,
    state_log_hi: float = -4.0,
    emission_sd: float = 0.5,
    switch_rate: float = 0.1,
) -> str:
    """Pyrho-style HMM smoothing of per-window ρ (M13.B).

    Discrete-state HMM over log10(ρ̂) with Gaussian emissions and a
    banded-walk transition prior. Returns JSON array of per-window
    rows with raw and smoothed ρ + Viterbi state.
    """
    opts: dict = {
        "window_size": window_size,
        "max_pair_distance": max_pair_distance,
        "n_states": n_states,
        "state_log_lo": state_log_lo,
        "state_log_hi": state_log_hi,
        "emission_sd": emission_sd,
        "switch_rate": switch_rate,
    }
    if step is not None:
        opts["step"] = step
    cypher = "CALL graphpop.recombination.hmm_smooth($sids, $options)"
    return json.dumps(_run_procedure(
        cypher, {"sids": list(sample_ids), "options": opts}))


@mcp.tool()
def graphpop_recombination_stratified_ld_decay(
    sample_ids: list[str],
    stratify_by: str = "population",
    window_size: int = 10_000,
    step: int | None = None,
    max_pair_distance: int = 5_000,
    min_maf: float = 0.05,
) -> str:
    """Per-(window, stratum) ρ via :Sample.population or :Sample.sex
    (M13.C).

    Splits the focal sample set by :Sample.<stratify_by> and runs
    the M11 Hudson-Kaplan moment estimator independently per
    stratum. Returns JSON array of per-(window, stratum) rows.
    """
    opts: dict = {
        "stratify_by": stratify_by,
        "window_size": window_size,
        "max_pair_distance": max_pair_distance,
        "min_maf": min_maf,
    }
    if step is not None:
        opts["step"] = step
    cypher = "CALL graphpop.recombination.stratified_ld_decay($sids, $options)"
    return json.dumps(_run_procedure(
        cypher, {"sids": list(sample_ids), "options": opts}))


@mcp.tool()
def graphpop_recombination_hotspots(
    sample_ids: list[str],
    method: str = "ld_decay",
    window_size: int = 10_000,
    step: int | None = None,
    max_pair_distance: int = 5_000,
    min_maf: float = 0.05,
    fdr_q: float = 0.05,
) -> str:
    """Recombination-hotspot detection with BH FDR (M13.D).

    Per-window ρ → z-score against the log10(ρ) genome-wide null →
    one-sided upper-tail normal p-value → Benjamini-Hochberg
    adjustment. Returns JSON array of per-window rows with
    is_hotspot flag at the requested FDR.
    """
    opts: dict = {
        "method": method,
        "window_size": window_size,
        "max_pair_distance": max_pair_distance,
        "min_maf": min_maf,
        "fdr_q": fdr_q,
    }
    if step is not None:
        opts["step"] = step
    cypher = "CALL graphpop.recombination.hotspots($sids, $options)"
    return json.dumps(_run_procedure(
        cypher, {"sids": list(sample_ids), "options": opts}))


@mcp.tool()
def graphpop_recombination_arg_breakpoints(
    run_id: str,
    window_size: int = 10_000,
    step: int | None = None,
) -> str:
    """ARG-derived per-window recombination rate (M11).

    For each sliding window: counts distinct :PARENT_OF breakpoints,
    sums marginal-tree branch lengths, and reports
    rho_per_bp = breakpoints / total_branch_length. Returns JSON
    array of per-window rows.
    """
    opts: dict = {"window_size": window_size}
    if step is not None:
        opts["step"] = step
    cypher = "CALL graphpop.recombination.arg_breakpoints($rid, $options)"
    return json.dumps(_run_procedure(
        cypher, {"rid": run_id, "options": opts}))


@mcp.tool()
def graphpop_recombination_ld_decay(
    sample_ids: list[str],
    window_size: int = 10_000,
    step: int | None = None,
    min_maf: float = 0.05,
    max_pair_distance: int = 5_000,
) -> str:
    """LD-decay Hudson-Kaplan moment estimator for per-window ρ (M11).

    Per window: pulls variants with derived-allele frequency ≥
    min_maf, computes pairwise r² for pairs within max_pair_distance
    bp, fits Hudson 1985 E[r²|n, ρ·d] via bisection to recover
    rho_per_bp. v1 limit: ≤ 63 focal samples. Returns JSON array.
    """
    opts: dict = {
        "window_size": window_size,
        "min_maf": min_maf,
        "max_pair_distance": max_pair_distance,
    }
    if step is not None:
        opts["step"] = step
    cypher = "CALL graphpop.recombination.ld_decay($sids, $options)"
    return json.dumps(_run_procedure(
        cypher, {"sids": list(sample_ids), "options": opts}))


@mcp.tool()
def graphpop_embedding_knn(
    sample_id: str,
    k: int = 10,
) -> str:
    """Top-k nearest neighbours by cosine similarity over
    :Sample.embedding (M10).

    Loads every embedded sample's vector in one Cypher pass and
    computes pairwise cosine against the query. Returns JSON array
    of {query_sample_id, neighbor_sample_id, cosine_similarity, rank}.
    """
    cypher = "CALL graphpop.embedding.knn($sid, $k)"
    return json.dumps(_run_procedure(
        cypher, {"sid": sample_id, "k": k}))


@mcp.tool()
def graphpop_embedding_cluster(
    method: str = "kmeans",
    k: int = 5,
    seed: int = 42,
    max_iter: int = 100,
) -> str:
    """k-means clustering over :Sample.embedding (M10).

    v1: method='kmeans' only (Lloyd's with k-means++ init).
    Returns JSON array of {sample_id, cluster_id,
    distance_to_centroid, n_clusters, method}.
    """
    opts: dict = {"seed": seed, "max_iter": max_iter}
    cypher = "CALL graphpop.embedding.cluster($method, $k, $options)"
    return json.dumps(_run_procedure(
        cypher, {"method": method, "k": k, "options": opts}))


@mcp.tool()
def graphpop_community_louvain(
    source: str,
    edge_weight: str = "unit",
    seed: int = 42,
    persist: bool = True,
) -> str:
    """Pure-Java Louvain modularity community detection (M9).

    Reads :RELATIVE edges with the given source, runs Blondel et al.
    2008 Louvain on the weighted graph, returns one row per sample
    with {sample_id, community_id, modularity, n_communities, source}.

    edge_weight: 'unit' (1.0), 'phi' (RELATIVE.phi), 'degree'
    (1 / (1 + RELATIVE.degree); closer relatives weigh more).
    persist=true also writes :IN_COMMUNITY edges and :Community nodes.
    """
    opts: dict = {
        "edge_weight": edge_weight,
        "seed": seed,
        "persist": persist,
    }
    cypher = "CALL graphpop.community.louvain($src, $options)"
    return json.dumps(_run_procedure(
        cypher, {"src": source, "options": opts}))


@mcp.tool()
def graphpop_pedigree_reconstruct(
    output_path: str,
    family_id: str | None = None,
    source: str = "king",
) -> str:
    """PRIMUS-lite pedigree reconstruction (M9, graphpop-pedigree).

    Reads :IN_FAMILY clusters and :RELATIVE labels from Neo4j,
    reconstructs the most-likely pedigree per family (parent-child
    trios, full-sib cohorts, 3-gen lineages with sex metadata),
    writes a PLINK-compatible PED file. Returns JSON with the
    output path and row count.
    """
    try:
        from graphpop_pedigree import PedExporter, PedigreeReconstructor
    except ImportError:
        return json.dumps({
            "error": "graphpop-pedigree not installed; run "
                     "`pip install graphpop-pedigree`",
        })

    sample_query = (
        "MATCH (s:Sample)-[:IN_FAMILY {source: $src}]->(f:Family) "
        "RETURN f.family_id AS family_id, s.sampleId AS sid, "
        "s.sex AS sex"
    )
    sample_rows = _run_procedure(sample_query, {"src": source})
    if family_id:
        sample_rows = [r for r in sample_rows
                       if r.get("family_id") == family_id]
    if not sample_rows:
        return json.dumps({
            "output_path": output_path,
            "n_rows": 0,
            "warning": "no samples — run graphpop.relate.families first",
        })

    edge_query = (
        "MATCH (a:Sample)-[r:RELATIVE {source: $src}]->(b:Sample) "
        "RETURN a.sampleId AS sa, b.sampleId AS sb, "
        "r.relationship AS rel, r.degree AS deg"
    )
    edge_rows = _run_procedure(edge_query, {"src": source})

    samples = [
        (r["family_id"], r["sid"], r.get("sex"))
        for r in sample_rows
    ]
    edges = [
        (r["sa"], r["sb"], r["rel"], r["deg"])
        for r in edge_rows
    ]
    rec = PedigreeReconstructor(samples, edges, source=source)
    rows = (rec.reconstruct_family(family_id)
            if family_id else rec.reconstruct_all())
    n = PedExporter.write(rows, output_path)
    return json.dumps({
        "output_path": output_path,
        "n_rows": n,
        "source": source,
        "family_id": family_id,
    })


@mcp.tool()
def graphpop_selection_allele_age_scan(
    run_id: str,
    n_freq_bins: int = 20,
    min_freq: float = 0.05,
    max_freq: float = 0.95,
) -> str:
    """Per-variant allele-age z-score conditional on frequency (M8).

    Stratifies all variants by derived-allele frequency (logit-spaced
    bins), computes per-bin (mean, sd) of log-age, emits a z-score
    per variant. Negative z = sweep candidate (too young for its
    frequency). Returns JSON array of per-variant rows.
    """
    opts: dict = {
        "n_freq_bins": n_freq_bins,
        "min_freq": min_freq,
        "max_freq": max_freq,
    }
    cypher = "CALL graphpop.selection.allele_age_scan($rid, $options)"
    return json.dumps(_run_procedure(
        cypher, {"rid": run_id, "options": opts}))


@mcp.tool()
def graphpop_selection_branch_outlier_scan(
    run_id: str,
    sample_ids: list[str],
    window_size: int = 10_000,
    step: int | None = None,
) -> str:
    """Per-window branch-length outlier scan (M8, Speidel et al. 2019).

    Slides a window across the genome; per window computes total
    branch length restricted to the focal sample set; emits z-score
    against the genome-wide null. Strong negative z = sweep candidate.
    Returns JSON array of per-window rows.
    """
    opts: dict = {"window_size": window_size}
    if step is not None:
        opts["step"] = step
    cypher = "CALL graphpop.selection.branch_outlier_scan($rid, $sids, $options)"
    return json.dumps(_run_procedure(
        cypher, {"rid": run_id, "sids": list(sample_ids), "options": opts}))


@mcp.tool()
def graphpop_demography_ne_trajectory(
    run_id: str,
    sample_ids: list[str],
    time_bins: list[float],
    ploidy: int = 2,
    population: str | None = None,
) -> str:
    """Closed-form Ne(t) trajectory (M7).

    Inverts graphpop.arg.coalescence_rate per bin: Ne(bin) =
    1 / (ploidy * rate(bin)). Standard error via delta method on
    Poisson event count. Empty bins return Ne=+Inf with
    flag='no_events'. Returns JSON array of per-bin rows.

    Pass sample_ids=[] together with population='EUR' (etc.) to
    resolve the focal set from :Sample.population.
    """
    opts: dict = {"time_bins": list(time_bins), "ploidy": ploidy}
    if population:
        opts["population"] = population
    cypher = "CALL graphpop.demography.ne_trajectory($rid, $sids, $options)"
    return json.dumps(_run_procedure(
        cypher, {"rid": run_id, "sids": list(sample_ids), "options": opts}))


@mcp.tool()
def graphpop_arg_tmrca(
    run_id: str,
    sample_a: str,
    sample_b: str,
    position: int | None = None,
    window_start: int | None = None,
    window_end: int | None = None,
) -> str:
    """Pairwise TMRCA in an ARG (M6).

    position set → single-position TMRCA (one row per haplotype pair).
    window_start/window_end set → span-weighted mean across overlapping
    marginal trees. No options → genome-wide mean. Returns JSON array.
    """
    opts: dict = {}
    if position is not None:
        opts["position"] = position
    if window_start is not None:
        opts["window_start"] = window_start
    if window_end is not None:
        opts["window_end"] = window_end
    cypher = "CALL graphpop.arg.tmrca($rid, $sa, $sb, $options)"
    return json.dumps(_run_procedure(
        cypher, {"rid": run_id, "sa": sample_a, "sb": sample_b,
                 "options": opts}))


@mcp.tool()
def graphpop_arg_branch_diversity(
    run_id: str,
    sample_ids: list[str],
    mode: str = "pi",
    windows: list[int] | None = None,
) -> str:
    """Branch-mode diversity from an ARG (M6, tskit-equivalent).

    mode='pi' (v1 only) reproduces tskit.TreeSequence.diversity(samples,
    mode='branch'). windows=[bp,bp,...] for per-window output, otherwise
    one row over the full sequence. Returns JSON array.
    """
    opts: dict = {"mode": mode}
    if windows is not None:
        opts["windows"] = list(windows)
    cypher = "CALL graphpop.arg.branch_diversity($rid, $sids, $options)"
    return json.dumps(_run_procedure(
        cypher, {"rid": run_id, "sids": list(sample_ids), "options": opts}))


@mcp.tool()
def graphpop_arg_coalescence_rate(
    run_id: str,
    sample_ids: list[str],
    time_bins: list[float],
) -> str:
    """Per-time-bin coalescence rate from an ARG (M6).

    Speidel/tsdate-style estimator: rate(bin) = events / lineage_pair_time,
    both span-weighted across marginal trees. time_bins must be
    monotonically increasing. Returns JSON array of per-bin rows.
    """
    opts: dict = {"time_bins": list(time_bins)}
    cypher = "CALL graphpop.arg.coalescence_rate($rid, $sids, $options)"
    return json.dumps(_run_procedure(
        cypher, {"rid": run_id, "sids": list(sample_ids), "options": opts}))


@mcp.tool()
def graphpop_arg_allele_age(
    run_id: str,
    variant_id: str,
) -> str:
    """Time bracket for a variant's :MUTATED_ON edge (M6).

    Returns JSON with child_time, parent_time, midpoint_time, and
    n_carriers (descendants of the child TreeNode in the marginal tree
    containing the mutation site).
    """
    cypher = "CALL graphpop.arg.allele_age($rid, $vid)"
    return json.dumps(_run_procedure(
        cypher, {"rid": run_id, "vid": variant_id}))


@mcp.tool()
def graphpop_relate_classify(
    source: str,
    min_phi: float | None = None,
    persist: bool = True,
    identical: float | None = None,
    first_degree: float | None = None,
    second_degree: float | None = None,
    third_degree: float | None = None,
    ibs0_threshold: float | None = None,
) -> str:
    """Per-pair relationship classification (Manichaikul 2010, M5).

    SOURCE selects the pairwise signal: 'king' or any IBD source name
    (e.g. 'hap_ibd', 'arg_derived'). Writes idempotent :RELATIVE edges
    per source. Returns JSON array of {sample_a, sample_b, relationship,
    degree, phi, ibs0_frac, source}.
    """
    opts: dict = {"persist": persist}
    if min_phi is not None:
        opts["min_phi"] = min_phi
    if identical is not None:
        opts["identical"] = identical
    if first_degree is not None:
        opts["first_degree"] = first_degree
    if second_degree is not None:
        opts["second_degree"] = second_degree
    if third_degree is not None:
        opts["third_degree"] = third_degree
    if ibs0_threshold is not None:
        opts["ibs0_threshold"] = ibs0_threshold
    cypher = "CALL graphpop.relate.classify($source, $options)"
    return json.dumps(_run_procedure(
        cypher, {"source": source, "options": opts}))


@mcp.tool()
def graphpop_relate_families(
    source: str,
    max_degree: int | None = None,
    min_total_ibd_bp: int | None = None,
    min_phi: float | None = None,
    persist: bool = True,
) -> str:
    """Family clustering via connected components (M5).

    Edge predicate (mutually exclusive — pick one):
      - max_degree (default 2): link pairs with :RELATIVE.degree <= N
      - min_total_ibd_bp: link pairs with summed :IBD_SEGMENT length >= N
      - min_phi: link pairs with :RELATIVE.phi >= F

    Each connected component is one extended family; isolated samples
    form singleton families. Returns JSON array of {sample_id,
    family_id, family_size, method}.
    """
    n_pred = sum(x is not None for x in (max_degree, min_total_ibd_bp, min_phi))
    if n_pred > 1:
        raise ValueError(
            "max_degree, min_total_ibd_bp, and min_phi are mutually "
            "exclusive; pick one.")
    opts: dict = {"persist": persist}
    if max_degree is not None:
        opts["max_degree"] = max_degree
    if min_total_ibd_bp is not None:
        opts["min_total_ibd_bp"] = min_total_ibd_bp
    if min_phi is not None:
        opts["min_phi"] = min_phi
    cypher = "CALL graphpop.relate.families($source, $options)"
    return json.dumps(_run_procedure(
        cypher, {"source": source, "options": opts}))


def main():
    """Entry point for the graphpop-mcp command."""
    mcp.run()


if __name__ == "__main__":
    main()
