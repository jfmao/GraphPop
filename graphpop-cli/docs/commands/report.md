# Generate automated HTML analysis report

## Description

`graphpop report` produces a self-contained HTML report summarizing the current
state of a GraphPop analysis. The report aggregates persisted statistics from the
graph into a single document with dataset overview, per-population diversity
table, pairwise Fst matrix, top selection signals, and annotation summary. It is
intended as a quick-look deliverable after running a full analysis pipeline.

The command queries the graph directly -- it does not read intermediate TSV files.
All tables and embedded figures are generated at report time from the data
currently stored on Population, GenomicWindow, Variant, Gene, and Pathway nodes.
The HTML output is self-contained (CSS and images are inlined) and can be shared
without a running Neo4j instance.

## Usage

```
graphpop report -o REPORT_FILE [OPTIONS]
```

## Arguments

This command takes no positional arguments.

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o`, `--output` | path | required | Output HTML file path (e.g., `report.html`). |
| `--title` | string | `"GraphPop Analysis Report"` | Report title displayed in the header. |
| `--sections` | string | `all` | Comma-separated list of sections to include: `overview`, `diversity`, `fst`, `selection`, `annotation`, `roh`. |
| `--pop` | string | all | Restrict report to specific populations (comma-separated). |
| `--top-genes` | integer | 20 | Number of top selection candidate genes to display. |
| `--top-windows` | integer | 50 | Number of top outlier windows to display per statistic. |
| `--chr` | string | all | Restrict to specific chromosomes (comma-separated). |

## Value

A self-contained HTML file with the following sections:

### 1. Dataset Overview
- Database name, creation date, GraphPop version.
- Node and relationship counts (mirrors `graphpop inventory`).
- Populations and chromosomes included.

### 2. Diversity Table
- One row per population: pi, theta_w, Tajima's D, H_exp, H_obs, F_IS.
- Color-coded cells for extreme values.
- Sortable columns.

### 3. Fst Matrix
- Pairwise Fst heatmap (Weir-Cockerham) as an embedded SVG.
- Numeric values annotated in cells for matrices with 15 or fewer populations.

### 4. Selection Signals
- Top N genes by composite selection score (same ranking as `graphpop rank-genes`).
- Top N outlier windows for each persisted selection statistic (iHS, XP-EHH, H12).
- Embedded mini-Manhattan plots per chromosome (if genome-scan data exists).

### 5. Annotation Summary
- Consequence type distribution (pie chart).
- Pathway coverage: number of pathways loaded, genes per pathway distribution.
- HIGH-impact variant count per population.

### 6. ROH Summary (if ROH data persisted)
- Per-population FROH distribution (violin plot).
- Top individuals by FROH.

## Details

### Self-contained output

All CSS, JavaScript (for sortable tables), and figures are embedded directly in
the HTML. The report can be opened in any modern browser without an internet
connection. Figures use inline SVG or base64-encoded PNG.

### Section dependencies

Each section is generated only if the required data exists in the graph. If
selection scores have not been persisted, the Selection Signals section is
omitted with a note. Use `--sections` to explicitly request only available
sections.

### Performance

Report generation queries multiple node types and performs aggregations. Typical
runtime is 10-30 seconds for a dataset with ~1M variants, 26 populations, and
full selection scans persisted. The bottleneck is the selection signal ranking
query across all genes.

## Examples

```bash
# Full report with all sections
graphpop report -o analysis_report.html

# Report for a subset of populations
graphpop report -o eur_report.html --pop EUR,CEU,GBR,FIN,TSI,IBS \
    --title "European Population Analysis"

# Only diversity and Fst sections
graphpop report -o quick_report.html --sections overview,diversity,fst

# Rice dataset report with top 50 genes
graphpop report -o rice_report.html --top-genes 50 \
    --title "Rice 3K Analysis Report"

# Single-chromosome report
graphpop report -o chr22_report.html --chr chr22 --top-windows 100
```

## See Also

- `graphpop inventory` -- quick text summary of graph contents
- `graphpop rank-genes` -- gene ranking used in the Selection Signals section
- `graphpop pop-summary` -- per-population diversity statistics
- `graphpop run-all` -- orchestrate full analysis before generating report
- `graphpop plot` -- individual publication-quality figures
