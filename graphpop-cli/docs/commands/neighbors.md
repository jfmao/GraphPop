# Explore graph neighborhood of a gene

## Description

`graphpop neighbors` traverses the graph outward from a focal gene, discovering
related genes via shared pathway membership, linkage disequilibrium edges, or
shared Gene Ontology terms. This is the primary tool for graph-native
exploration: rather than relying on a flat gene list, it uses the topology of the
population variation graph to reveal biologically meaningful connections.

The traversal starts at the named Gene node and follows edges of the specified
type(s) up to the requested hop depth. At each hop, the command collects gene
symbols, relationship types, and any persisted statistics. The result is a table
of neighboring genes ranked by graph distance and annotated with the path that
connects them to the focal gene.

Use `neighbors` when you have identified a gene of interest (e.g., from
`graphpop rank-genes` or `graphpop lookup`) and want to understand its
functional and genomic context within the graph.

## Usage

```
graphpop neighbors GENE [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `GENE` | string | yes | -- | Gene symbol to start the traversal from (e.g., `BRCA2`, `Os01g0100100`). |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--hops` | integer | 1 | Maximum traversal depth (1, 2, or 3). Higher values explore further but produce larger result sets. |
| `--via` | string | `IN_PATHWAY` | Edge type(s) to traverse. Comma-separated list from: `IN_PATHWAY`, `LD`, `HAS_GO_TERM`. |
| `--pop` | string | none | Population for which to report persisted statistics (Fst, iHS) on discovered genes. |
| `--min-ld` | float | 0.3 | Minimum r-squared threshold when traversing LD edges. Ignored for other edge types. |
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, `json`. |
| `-o`, `--output` | path | stdout | Write output to a file instead of stdout. |

## Value

Returns one row per discovered neighbor gene:

| Column | Type | Description |
|--------|------|-------------|
| `gene` | string | Neighbor gene symbol. |
| `chr` | string | Chromosome of the neighbor gene. |
| `hops` | integer | Graph distance from the focal gene (1 = directly connected). |
| `via` | string | Edge type(s) used to reach this gene (e.g., `IN_PATHWAY`, `LD`). |
| `shared_entity` | string | Name of the shared pathway, GO term, or LD variant connecting the two genes. |
| `mean_fst` | float | Mean Fst for the neighbor gene (if `--pop` is specified and scores are persisted). |
| `max_abs_ihs` | float | Maximum |iHS| in the neighbor gene (if available). |

## Details

### Traversal semantics

- **IN_PATHWAY**: Gene --> IN_PATHWAY --> Pathway <-- IN_PATHWAY <-- Gene.
  One hop discovers all genes sharing at least one pathway with the focal gene.
  Two hops extend through a second shared pathway (genes that share a pathway
  with a gene that shares a pathway with the focal gene).

- **LD**: Gene --> contains Variant --> LD --> Variant <-- in Gene. One hop
  discovers genes containing variants in LD (r-squared above `--min-ld`) with
  variants in the focal gene. This captures local genomic co-inheritance.

- **HAS_GO_TERM**: Gene --> HAS_GO_TERM --> GOTerm <-- HAS_GO_TERM <-- Gene.
  Discovers genes annotated with the same GO terms.

When multiple `--via` types are specified, the traversal follows all listed edge
types. A gene reachable by multiple paths appears once, with the shortest hop
distance and all connecting paths listed.

### Hop limits

Maximum hops is capped at 3 to prevent combinatorial explosion. At depth 3 via
IN_PATHWAY, the result set may contain hundreds of genes. Use `--pop` to include
selection statistics for prioritization.

### Performance

One-hop IN_PATHWAY queries resolve in under 500 ms for typical pathway sizes
(10-100 genes). LD traversals depend on the density of persisted LD edges and
may take 1-5 seconds for a gene with thousands of variants. Three-hop traversals
may take 5-15 seconds on densely connected subgraphs.

## Examples

```bash
# Find genes sharing a pathway with BRCA2
graphpop neighbors BRCA2 --via IN_PATHWAY --pop EUR

# Explore LD neighborhood of a rice gene, 2 hops deep
graphpop neighbors Os01g0100100 --via LD --hops 2 --min-ld 0.5 --pop GJ-tmp

# Combine pathway and GO term traversal
graphpop neighbors LCT --via IN_PATHWAY,HAS_GO_TERM --hops 2 --pop CEU \
    -o lct_neighbors.tsv

# JSON output for programmatic use
graphpop neighbors EDAR --via IN_PATHWAY --pop EAS --format json

# Deep exploration with 3 hops through pathways
graphpop neighbors SLC24A5 --via IN_PATHWAY --hops 3 --pop EUR \
    -o slc24a5_network.tsv
```

## See Also

- `graphpop lookup gene` -- detailed annotation for a single gene
- `graphpop rank-genes` -- rank genes by composite selection score
- `graphpop converge` -- find regions where multiple statistics converge
- `graphpop ld` -- compute linkage disequilibrium between variant pairs
