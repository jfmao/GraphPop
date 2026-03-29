# graphpop filter

## Title

Query Persisted Statistics with Annotation-Based Filters

## Description

Retrieves already-computed population genetics statistics from graph nodes and
filters them by functional annotation (VEP consequence, pathway, gene). This is
the recommended way to perform **conditioned analysis** for haplotype-based
statistics (iHS, XP-EHH, nSL) and window-level statistics (Fst, pi, Tajima's D,
H12).

This command does **not** re-compute statistics. It queries values that were
previously persisted to graph nodes by running the corresponding procedure with
`--persist`.

## Usage

```
graphpop filter STATISTIC CHR POPULATION [OPTIONS]
```

## Arguments

| Argument | Type | Description |
|---|---|---|
| `STATISTIC` | CHOICE | One of: `ihs`, `xpehh`, `nsl`, `fst`, `pi`, `tajima_d`, `h12`. |
| `CHR` | TEXT | Chromosome identifier (e.g., `chr1`, `chr22`). |
| `POPULATION` | TEXT | Population identifier (e.g., `EUR`, `GJ-tmp`). |

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `-o`, `--output` | PATH | stdout | Output file path. |
| `--format` | CHOICE | `tsv` | Output format: `tsv`, `csv`, or `json`. |
| `--consequence` | TEXT | *(none)* | Filter by VEP consequence type (e.g., `missense_variant`, `synonymous_variant`, `stop_gained`). |
| `--pathway` | TEXT | *(none)* | Filter by pathway name (partial match via `CONTAINS`). |
| `--gene` | TEXT | *(none)* | Filter by gene name or Ensembl gene ID. |
| `--min-score` | FLOAT | *(none)* | Minimum absolute score threshold. |
| `--max-score` | FLOAT | *(none)* | Maximum absolute score threshold. |
| `--pop2` | TEXT | *(none)* | Second population for XP-EHH queries. Required for `xpehh`. |
| `--limit` | INT | `10000` | Maximum number of rows returned. |

## Value

Returns a table (TSV by default) with columns depending on the statistic type:

### Per-variant statistics (ihs, xpehh, nsl)

| Column | Description |
|---|---|
| `variant_id` | Variant identifier (`chr:pos:ref:alt`) |
| `pos` | Genomic position |
| `ihs` / `xpehh` / `nsl` | Standardized score |
| `ihs_unstd` / `xpehh_unstd` / `nsl_unstd` | Unstandardized score |
| `consequence` | VEP consequence (if `--consequence` filter used) |
| `impact` | VEP impact level (if `--consequence` filter used) |
| `gene` | Gene symbol (if `--gene` filter used) |

### Window-level statistics (fst, pi, tajima_d)

| Column | Description |
|---|---|
| `window_id` | Window identifier |
| `start` | Window start position |
| `end` | Window end position |
| `fst` / `pi` / `tajima_d` | Statistic value |
| `n_variants` | Number of variants in window |

### Sweep windows (h12)

| Column | Description |
|---|---|
| `window_id` | Window identifier |
| `chr` | Chromosome |
| `start` | Window start position |
| `end` | Window end position |
| `h12` | Garud's H12 |
| `h2_h1` | H2/H1 ratio |
| `hap_div` | Haplotype diversity |

## Details

### The Compute-Then-Filter Workflow

GraphPop uses two approaches for conditioned analysis, depending on the
statistic type:

**Fast-path statistics** (diversity, divergence, SFS, genome-scan) support
direct conditioning via `--consequence`, `--pathway`, and `--gene` options on
the procedure call itself. These filter variants *before* computing the
statistic:

```bash
# Direct conditioning (built into the procedure call)
graphpop diversity chr1 1 43270923 GJ-tmp --consequence missense_variant
```

**Full-path statistics** (iHS, XP-EHH, nSL) **cannot** be conditioned during
computation because they require genome-wide haplotype structure. The Extended
Haplotype Homozygosity (EHH) decay depends on all variants in the haplotype,
not just those matching a filter. Therefore, full-path statistics must be:

1. **Computed genome-wide** with `--persist` to store scores on Variant nodes.
2. **Filtered post-hoc** using `graphpop filter` to extract variants matching
   annotation criteria.

This two-step workflow is statistically correct because it preserves the
genome-wide standardization of scores (allele-frequency bin normalization for
iHS/nSL, genome-wide normalization for XP-EHH).

### How Annotation Filtering Works

The `--consequence`, `--pathway`, and `--gene` filters work by traversing the
graph from Variant nodes to annotation nodes:

- `--consequence`: Follows `Variant-[:HAS_CONSEQUENCE]->()` edges and matches
  on the `consequence` property.
- `--pathway`: Follows
  `Variant-[:HAS_CONSEQUENCE]->(:Gene)-[:IN_PATHWAY]->(:Pathway)` and matches
  pathway names with `CONTAINS` (partial match).
- `--gene`: Follows `Variant-[:HAS_CONSEQUENCE]->(:Gene)` and matches on
  `geneId` or `symbol`.

Multiple filters can be combined. They are joined with AND logic.

### Score Thresholds

For `--min-score` and `--max-score`, the absolute value of the score is used
for comparison on per-variant statistics (iHS, XP-EHH, nSL). This means
`--min-score 2.0` returns variants with `|score| >= 2.0`, capturing both
positive selection signals (positive scores) and derived allele sweeps
(negative scores).

For window-level statistics (Fst, pi, Tajima's D, H12), the raw value is used
(not absolute).

### Common Consequence Types

| Consequence | Description |
|---|---|
| `missense_variant` | Non-synonymous amino acid change |
| `synonymous_variant` | No amino acid change (neutral proxy) |
| `stop_gained` | Premature stop codon (loss of function) |
| `frameshift_variant` | Insertion/deletion shifting reading frame |
| `splice_donor_variant` | Disrupts splice donor site |
| `intergenic_variant` | Between genes |
| `intron_variant` | Within an intron |
| `5_prime_UTR_variant` | In the 5' untranslated region |
| `3_prime_UTR_variant` | In the 3' untranslated region |

## Examples

### iHS at missense variants (selection on protein-changing mutations)

```bash
# Step 1: Compute iHS genome-wide and persist scores
graphpop ihs chr1 GJ-tmp --persist

# Step 2: Filter to missense variants with strong signals
graphpop filter ihs chr1 GJ-tmp \
    --consequence missense_variant \
    --min-score 2.0 \
    -o ihs_missense_sweeps.tsv
```

### XP-EHH between two rice subpopulations in a specific pathway

```bash
# Step 1: Compute XP-EHH and persist
graphpop xpehh chr1 GJ-tmp XI-1A --persist

# Step 2: Filter to starch biosynthesis pathway
graphpop filter xpehh chr1 GJ-tmp \
    --pop2 XI-1A \
    --pathway "Starch biosynthesis" \
    -o xpehh_starch.tsv
```

### nSL at a specific gene

```bash
# Step 1: Compute nSL and persist
graphpop nsl chr1 GJ-tmp --persist

# Step 2: Filter to the GW5 gene (grain width QTL)
graphpop filter nsl chr1 GJ-tmp \
    --gene GW5 \
    --min-score 2.0 \
    -o nsl_gw5.tsv
```

### High-Fst genomic windows

```bash
# Step 1: Run genome scan to create GenomicWindow nodes
graphpop genome-scan chr1 GJ-tmp 100000 50000 --persist

# Step 2: Filter windows with extreme differentiation
graphpop filter fst chr1 GJ-tmp \
    --min-score 0.5 \
    -o high_fst_windows.tsv
```

### Tajima's D indicating balancing selection

```bash
graphpop filter tajima_d chr22 EUR \
    --min-score 2.0 \
    -o balancing_selection.tsv
```

### Garud's H12 sweep windows at missense variants

```bash
# Step 1: Compute Garud's H
graphpop garud-h chr1 GJ-tmp 100000 50000 --persist

# Step 2: Filter for strong sweeps
graphpop filter h12 chr1 GJ-tmp \
    --min-score 0.3 \
    -o hard_sweeps.tsv
```

### Combined filters: missense variants in a pathway with strong iHS

```bash
graphpop filter ihs chr1 EUR \
    --consequence missense_variant \
    --pathway "Cardiac repolarization" \
    --min-score 2.5 \
    --format json \
    -o cardiac_sweeps.json
```

### piN/piS-style analysis (synonymous vs. missense)

```bash
# Direct conditioning on fast-path statistics:
graphpop diversity chr1 1 43270923 GJ-tmp --consequence missense_variant -o pi_N.tsv
graphpop diversity chr1 1 43270923 GJ-tmp --consequence synonymous_variant -o pi_S.tsv

# Or via filter on persisted genome-scan windows:
graphpop filter pi chr1 GJ-tmp -o pi_all_windows.tsv
```

## See Also

- `graphpop ihs` -- Compute integrated haplotype score.
- `graphpop xpehh` -- Compute cross-population extended haplotype homozygosity.
- `graphpop nsl` -- Compute number of segregating sites by length.
- `graphpop diversity` -- Compute diversity with direct conditioning.
- `graphpop genome-scan` -- Sliding-window scan that creates GenomicWindow nodes.
- `graphpop garud-h` -- Compute Garud's H statistics.
- `graphpop export-windows` -- Alternative way to export GenomicWindow nodes.
