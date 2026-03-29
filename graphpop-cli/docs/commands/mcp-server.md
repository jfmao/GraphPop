# graphpop-mcp — MCP server for AI agent access

Programmatic interface to GraphPop via the Model Context Protocol

## Description

`graphpop-mcp` is a standalone MCP (Model Context Protocol) server that exposes
GraphPop's capabilities as typed tools for AI agents and external applications.
It enables large language model agents to compute population genetics statistics,
query gene annotations, detect convergent selection signals, and explore graph
topology — all without writing Cypher or using the command line.

The MCP server exposes 21 tools that mirror the CLI's functionality:

| Category | Tools | Count |
|----------|-------|-------|
| FAST PATH procedures | diversity, divergence, sfs, joint_sfs, genome_scan, pop_summary | 6 |
| FULL PATH procedures | ihs, xpehh, nsl, roh, garud_h, ld | 6 |
| Multi-stat integration | converge, rank_genes | 2 |
| Annotation lookup | lookup_gene, lookup_pathway, lookup_region | 3 |
| Filtering | filter | 1 |
| Utilities | status, inventory, query | 3 |

## Installation

```bash
cd graphpop-mcp
pip install -e .
```

## Usage

### Start the MCP server

```bash
# Via environment variables (same as CLI)
export GRAPHPOP_URI=bolt://localhost:7687
export GRAPHPOP_USER=neo4j
export GRAPHPOP_PASSWORD=mypassword
export GRAPHPOP_DATABASE=rice3k

graphpop-mcp
```

### Configure for Claude Desktop / Claude Code

Add to your MCP configuration:

```json
{
  "mcpServers": {
    "graphpop": {
      "command": "graphpop-mcp",
      "env": {
        "GRAPHPOP_URI": "bolt://localhost:7687",
        "GRAPHPOP_USER": "neo4j",
        "GRAPHPOP_PASSWORD": "mypassword",
        "GRAPHPOP_DATABASE": "rice3k"
      }
    }
  }
}
```

### Configure for ToolUniverse

The MCP server is registered with ToolUniverse (Zitnik Lab, Harvard) under the
name "GraphPop". AI agent frameworks that query ToolUniverse can discover and
use GraphPop tools automatically.

## Tools Reference

### Population Genetics Procedures

All 12 procedures accept the same parameters as the CLI commands:

- `graphpop_diversity(chr, start, end, pop, consequence?, pathway?, ...)` → JSON
- `graphpop_divergence(chr, start, end, pop1, pop2, pop3?, ...)` → JSON
- `graphpop_sfs(chr, start, end, pop, folded?, ...)` → JSON
- `graphpop_joint_sfs(chr, start, end, pop1, pop2, folded?, ...)` → JSON
- `graphpop_genome_scan(chr, pop, window?, step?, pop2?, consequence?, ...)` → JSON
- `graphpop_pop_summary(chr, pop, consequence?, pathway?, ...)` → JSON
- `graphpop_ld(chr, start, end, pop, max_dist?, r2_threshold?, ...)` → JSON
- `graphpop_ihs(chr, pop, min_af?, ...)` → JSON
- `graphpop_xpehh(chr, pop1, pop2, min_af?, ...)` → JSON
- `graphpop_nsl(chr, pop, min_af?, ...)` → JSON
- `graphpop_roh(chr, pop, method?, min_length?, ...)` → JSON
- `graphpop_garud_h(chr, pop, window?, step?, ...)` → JSON

### Analytical Tools

- `graphpop_converge(pop, stats, thresholds, chr?, pop2?)` — Find convergent multi-statistic signals
- `graphpop_rank_genes(pop, chr?, top?, pop2?)` — Rank genes by composite selection evidence
- `graphpop_filter(statistic, chr, pop, consequence?, pathway?, gene?, min_score?)` — Post-hoc annotation filtering

### Annotation Lookup

- `graphpop_lookup_gene(gene)` — Variant count, consequences, pathways, selection stats for a gene
- `graphpop_lookup_pathway(pathway)` — Member genes and variant counts for a pathway
- `graphpop_lookup_region(chr, start, end)` — Genes and stats in a genomic region

### Utilities

- `graphpop_status()` — Node counts and population list
- `graphpop_inventory()` — Comprehensive database inventory (what's computed, what's loaded)
- `graphpop_query(cypher, params?)` — Run arbitrary Cypher

## Environment Variables

The MCP server shares configuration with the CLI:

| Variable | Default | Description |
|----------|---------|-------------|
| `GRAPHPOP_URI` | `bolt://localhost:7687` | Neo4j bolt URI |
| `GRAPHPOP_USER` | `neo4j` | Neo4j username |
| `GRAPHPOP_PASSWORD` | `graphpop` | Neo4j password |
| `GRAPHPOP_DATABASE` | `neo4j` | Neo4j database name |

For backward compatibility, `GRAPHPOP_NEO4J_URI`, `GRAPHPOP_NEO4J_USER`, and
`GRAPHPOP_NEO4J_PASSWORD` are also accepted (CLI-style names take precedence).

## Examples

### AI agent workflow

An AI agent using the MCP server can autonomously:

1. Check what's available: `graphpop_inventory()`
2. Compute diversity: `graphpop_diversity(chr="chr22", start=1, end=51304566, pop="CEU")`
3. Condition on annotation: `graphpop_diversity(chr="chr22", start=1, end=51304566, pop="CEU", consequence="missense_variant")`
4. Find convergent signals: `graphpop_converge(pop="CEU", stats="ihs,xpehh", thresholds="2.0,2.0")`
5. Look up top genes: `graphpop_rank_genes(pop="CEU", top=10)`
6. Explore a gene: `graphpop_lookup_gene(gene="KCNE1")`
7. Filter by annotation: `graphpop_filter(statistic="ihs", chr="chr22", pop="CEU", consequence="missense_variant")`

### Direct Cypher query

```python
graphpop_query(
    cypher="MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene) WHERE v.ihs_CEU > 2.5 RETURN g.symbol, count(v) ORDER BY count(v) DESC LIMIT 10"
)
```

## See Also

- `graphpop` CLI — command-line interface (60 commands)
- `graphpop-import` — data import pipeline
- GraphPop paper — Supplementary Note 8 (CLI Architecture)
- ToolUniverse registration — `graphpop-mcp/tooluniverse/graphpop_tools.json`
