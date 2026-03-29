# graphpop query

## Title

Run Arbitrary Cypher Queries and Format Output

## Description

Executes a user-supplied Cypher statement against the active Neo4j database and
formats the results as TSV, CSV, or JSON. This is a general-purpose escape hatch
for any query not covered by the built-in GraphPop commands.

## Usage

```
graphpop query CYPHER [OPTIONS]
```

## Arguments

| Argument | Type | Description |
|---|---|---|
| `CYPHER` | TEXT | A Cypher query string. Must be quoted on the command line. |

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `-o`, `--output` | PATH | stdout | Output file path. |
| `--format` | CHOICE | `tsv` | Output format: `tsv`, `csv`, or `json`. |

## Value

Returns the query results formatted as the specified output type. Column names
are derived from the Cypher `RETURN` clause aliases. Floats are formatted to 6
significant figures. Lists are comma-joined. Null values appear as `NA` (TSV/CSV)
or `null` (JSON).

## Details

The Cypher string is passed directly to the Neo4j driver. Any valid Cypher
statement is accepted, including `MATCH`, `CALL`, `WITH`, subqueries, and
aggregations. Read and write queries are both supported (use write queries with
care).

The query runs against the active database as configured in
`~/.graphpop/config.yaml` or overridden with `--database`.

**Tips for shell quoting:** Use double quotes around the Cypher string on the
command line. If the Cypher itself contains quotes, use single quotes inside:

```bash
graphpop query "MATCH (v:Variant) WHERE v.chr = 'chr22' RETURN count(v) AS n"
```

For complex queries, consider writing the Cypher to a file and using shell
substitution:

```bash
graphpop query "$(cat my_query.cypher)" -o results.tsv
```

### Connection Errors

If Neo4j is not running, the command prints a helpful error message with the
current connection URI and suggests checking the configuration.

## Examples

### Count variants per chromosome

```bash
graphpop query "MATCH (v:Variant) RETURN v.chr AS chr, count(v) AS n_variants ORDER BY chr"
```

### List populations and sample counts

```bash
graphpop query \
    "MATCH (p:Population)<-[:IN_POPULATION]-(s:Sample) \
     RETURN p.populationId AS population, count(s) AS n_samples \
     ORDER BY n_samples DESC"
```

### Find variants in a specific gene

```bash
graphpop query \
    "MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene) \
     WHERE g.symbol = 'BRCA1' \
     RETURN v.variantId AS id, v.pos AS pos, v.ref AS ref, v.alt AS alt \
     ORDER BY v.pos" \
    -o brca1_variants.tsv
```

### Export as JSON for programmatic use

```bash
graphpop query \
    "MATCH (v:Variant) WHERE v.ihs_EUR IS NOT NULL AND abs(v.ihs_EUR) > 3 \
     RETURN v.variantId AS id, v.pos AS pos, v.ihs_EUR AS ihs \
     ORDER BY abs(v.ihs_EUR) DESC LIMIT 50" \
    --format json -o top_ihs.json
```

### Run a GraphPop stored procedure directly

```bash
graphpop query \
    "CALL graphpop.diversity('chr22', 1, 50818468, 'EUR') \
     YIELD pi, theta_w, tajima_d"
```

## See Also

- `graphpop filter` -- Structured filtering of persisted statistics.
- `graphpop export-windows` -- Structured export of GenomicWindow nodes.
- `graphpop db info` -- Show database schema (node labels, relationship types).
