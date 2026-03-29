# GraphPop CLI — Command Tree

```
graphpop
│
├─── SETUP & SERVER ────────────────────────────────────────────────
│    ├── setup              Download, configure, and initialize Neo4j
│    │                        --neo4j-home, --pagecache, --heap, --password
│    │                        --deploy-plugin (GraphPop procedures JAR)
│    ├── start              Start the Neo4j database server
│    ├── stop               Stop the Neo4j database server
│    └── status             Check Neo4j status, version, plugin, config
│
├─── DATABASE MANAGEMENT ───────────────────────────────────────────
│    ├── db list            List all databases with sizes and status
│    ├── db create NAME     Create a new database
│    ├── db switch NAME     Set the active database in config
│    ├── db drop NAME       Drop a database (with confirmation)
│    ├── db info            Show node/edge counts, procedures installed
│    ├── import             Import VCF → Neo4j graph database
│    │                        --vcf, --panel, --database (user-named)
│    │                        --vep, --pathways, --go-terms, --ancestral
│    ├── dump               Dump database to file for sharing
│    │                        --database, -o output.dump, --manifest
│    └── load               Restore database from dump file
│                             --dump-file, --database, --overwrite
│
├─── CONFIGURATION & VALIDATION ────────────────────────────────────
│    ├── config init        Create config file interactively
│    ├── config show        Display current configuration
│    ├── config set K V     Set a configuration value
│    ├── config path        Print config file path
│    ├── validate           Check graph integrity, indexes, procedures
│    │                        --fix (attempt automatic repairs)
│    └── inventory          Show what has been computed on this database
│                             Populations analyzed, chromosomes covered,
│                             which statistics are persisted, annotation
│                             types loaded (VEP, pathways, GO, ancestral)
│
├─── POPULATION GENETICS PROCEDURES (12) ───────────────────────────
│    │
│    ├── FAST PATH (allele count arrays, O(V×K), sample-independent)
│    │   │  All support: --consequence, --pathway, --gene conditioning
│    │   │
│    │   ├── diversity      π, θ_W, Tajima's D, Fay & Wu's H, H_e, H_o, F_IS
│    │   ├── divergence     Hudson Fst, W&C Fst, D_xy, D_a, PBS
│    │   ├── sfs            Site frequency spectrum (folded/unfolded)
│    │   ├── joint-sfs      2D joint SFS between two populations
│    │   ├── genome-scan    Sliding-window scan → GenomicWindow nodes
│    │   │                    --persist writes permanent window records
│    │   └── pop-summary    Whole-chromosome summary → Population nodes
│    │
│    └── FULL PATH (bit-packed haplotypes, individual genotype data)
│        │  Use --persist to write scores to graph; then use 'filter'
│        │
│        ├── ihs            Integrated haplotype score (AF-bin standardized)
│        ├── xpehh          Cross-population EHH (genome-wide standardized)
│        ├── nsl            Number of segregating sites by length
│        ├── roh            Runs of homozygosity (HMM or sliding-window)
│        │                    --method hmm|window  --min-length
│        ├── garud-h        Garud's H1, H12, H2/H1 (sweep detection)
│        │                    window_size, step_size
│        └── ld             Pairwise r², D' (writes LD edges)
│                             --threshold  --persist
│
├─── CONDITIONED ANALYSIS & FILTERING ──────────────────────────────
│    │
│    │  FAST PATH: direct conditioning via --consequence/--pathway/--gene
│    │    graphpop diversity chr1 1 43270923 POP --consequence missense_variant
│    │    graphpop genome-scan chr1 POP 100000 50000 --pathway "Starch biosyn."
│    │
│    │  FULL PATH: compute-then-filter (statistically correct)
│    │    Step 1: graphpop ihs chr1 POP --persist
│    │    Step 2: graphpop filter ihs chr1 POP --consequence missense_variant
│    │
│    └── filter             Query persisted stats with annotation filters
│                             ihs|xpehh|nsl|fst|pi|tajima_d|h12
│                             --consequence, --pathway, --gene
│                             --min-score, --max-score, --pop2
│
├─── ANNOTATION LOOKUP & EXPLORATION ───────────────────────────────
│    ├── lookup gene GENE   Show all variants, pathways, GO terms for a gene
│    │                        Returns: variant count, consequence types,
│    │                        pathway membership, selection statistics
│    ├── lookup pathway PW  Show all genes and variants in a pathway
│    │                        Returns: gene list, mean Fst, variant count
│    ├── lookup variant ID  Show full annotation for a specific variant
│    │                        Returns: allele freqs, Fst, iHS, consequence,
│    │                        gene, pathway, GO terms
│    ├── lookup region      Show genes, variants, and statistics in a region
│    │     CHR START END      Returns: gene list + summary stats per gene
│    └── neighbors GENE     Show genes connected via shared pathways or LD
│         [--hops N]          Graph neighborhood exploration (1-3 hops)
│
├─── MULTI-STATISTIC INTEGRATION ───────────────────────────────────
│    ├── converge           Find loci where multiple statistics agree
│    │     --stats ihs,xpehh,h12,fst  (which statistics to intersect)
│    │     --thresholds 2.0,2.0,0.3,0.5  (per-statistic thresholds)
│    │     --chr CHR --pop POP
│    │     Returns: regions with convergent signals + gene annotations
│    │
│    ├── rank-genes         Rank genes by composite selection evidence
│    │     --pop POP          Combines: max |iHS|, max |XP-EHH|, max H12,
│    │     --top N            mean Fst, consequence enrichment score
│    │     Returns: ranked gene table with per-statistic scores
│    │
│    └── compare            Compare statistics between two populations
│         POP1 POP2 CHR       Returns: per-window delta-pi, delta-Fst,
│         --stat pi|fst|...   PBS, differential iHS, with significance
│
├─── ORCHESTRATION & AGGREGATION ───────────────────────────────────
│    ├── run-all            Full-genome analysis (all pops × all chrs)
│    │                        --phase 1|2|all  --resume  --persist
│    │                        --populations, --chromosomes, --xpehh-pairs
│    │                        -d output_dir/  --json-output results.json
│    ├── aggregate          Generate summary tables from run-all results
│    │                        -d results_dir/  -o tables/
│    │                        → population_summary, fst_matrix, roh_summary,
│    │                          selection_peaks, sweep_windows
│    ├── export-windows     Batch export GenomicWindow nodes with filters
│    │                        --min-fst, --max-tajima-d, --min-pi, --limit
│    ├── batch              Run any command across multiple pops/chrs
│    │                        graphpop batch diversity --pops EUR,AFR,EAS
│    │                        --chrs chr1,chr2,chr3 -d output/
│    │                        Parallelizes with --workers N
│    └── report             Generate automated analysis summary report
│                             --database NAME  -o report.html
│                             Includes: dataset overview, diversity table,
│                             Fst matrix, top selection signals, QC metrics,
│                             key figures (auto-generated)
│
├─── DATA EXTRACTION & EXPORT ──────────────────────────────────────
│    ├── extract variants   Export variant data with filters
│    │     --chr, --start, --end, --pop, --min-af, --max-af
│    │     --consequence, --pathway, --gene
│    │     --fields pos,ref,alt,af,fst,ihs  (select output columns)
│    │     Returns: TSV of variant-level data
│    │
│    ├── extract genotypes  Export genotype matrix for a region
│    │     --chr, --start, --end, --pop, --samples
│    │     --format dosage|gt|haplotype
│    │     Returns: sample × variant matrix
│    │
│    ├── extract samples    Export per-sample summary statistics
│    │     --pop POP          Returns: sampleId, het_rate, n_roh, froh,
│    │                        rare_burden, population assignment
│    │
│    └── export-bed         Export regions as BED file (for bedtools)
│         --stat fst|ihs|h12  --threshold FLOAT
│         --chr CHR --pop POP
│         Returns: BED format (chr, start, end, name, score)
│
├─── VISUALIZATION (graphpop plot) ─────────────────────────────────
│    │
│    ├── Standard plot types
│    │   ├── diversity-bar  Per-population diversity ranking
│    │   ├── fst-heatmap    Pairwise Fst matrix with color scale
│    │   ├── manhattan      Genome-wide statistic scan
│    │   │                    --stat ihs|xpehh|nsl|fst|pi|tajima_d|h12
│    │   │                    --threshold, --abs-value|--raw-value
│    │   ├── pinpis         πN/πS ratios across populations
│    │   ├── sfs-plot       Site frequency spectrum histogram
│    │   │                    --log-scale
│    │   └── roh-landscape  Per-population FROH violin plots
│    │
│    ├── Genomic context plots
│    │   ├── gene-zoom      Multi-track view at a focal gene/region
│    │   │     GENE|CHR:START-END  --pop POP
│    │   │     Shows: Fst track, iHS track, H12 track, gene model,
│    │   │     consequence annotations — stacked alignment
│    │   │
│    │   └── chromosome     Multi-statistic chromosome ideogram
│    │         --chr CHR --pop POP
│    │         --stats fst,ihs,pi  (which tracks to show)
│    │         Alternating chromosome colors, gene annotation band
│    │
│    └── Population structure plots
│        ├── pop-tree       UPGMA/NJ tree from Fst matrix
│        │     INPUT_DIR -o tree.png
│        │     --method upgma|nj  --bootstrap N
│        │
│        ├── pca-scatter    PCA from allele frequency matrix
│        │     INPUT_FILE -o pca.png
│        │     --color-by population|group  --pc 1,2
│        │
│        └── heatmap        General-purpose heatmap from any matrix
│              INPUT_FILE -o heatmap.png
│              --cmap viridis|RdBu|YlOrRd  --annotate
│
├─── UTILITIES ─────────────────────────────────────────────────────
│    └── query              Run arbitrary Cypher → TSV/CSV/JSON output
│                             graphpop query "MATCH (v:Variant) WHERE ..."
│
├─── MCP SERVER (graphpop-mcp) ──────────────────────────────────────
│    │  Programmatic access for AI agents via Model Context Protocol
│    │  21 tools — same procedures + analytical + lookup + utilities
│    │
│    ├── Procedures (12)     Same as CLI: diversity, divergence, sfs, etc.
│    ├── Analytical (2)      converge, rank_genes
│    ├── Lookup (3)          lookup_gene, lookup_pathway, lookup_region
│    ├── Filtering (1)       filter (post-hoc annotation filtering)
│    └── Utilities (3)       status, inventory, query
│
│    Configuration: same env vars as CLI (GRAPHPOP_URI, GRAPHPOP_DATABASE)
│    Start: graphpop-mcp
│    Registration: ToolUniverse (Zitnik Lab, Harvard)
│
└─── GLOBAL OPTIONS ────────────────────────────────────────────────
     ├── --uri              Neo4j bolt URI (or GRAPHPOP_URI env var)
     ├── --user             Neo4j username (or GRAPHPOP_USER)
     ├── --password         Neo4j password (or GRAPHPOP_PASSWORD)
     ├── --database         Neo4j database (or GRAPHPOP_DATABASE)
     ├── --config           Config file path (~/.graphpop/config.yaml)
     └── --version          Show version

     OUTPUT OPTIONS (per command):
     ├── -o, --output FILE  Output file (default: stdout, pipeable)
     └── --format tsv|csv|json  Output format (default: tsv)
```

## New Commands (not yet implemented)

Commands marked with ★ are newly proposed additions to the CLI.

### ★ Annotation Lookup (`graphpop lookup`)
Enables users to explore the graph interactively without writing Cypher.
- `graphpop lookup gene KCNE1` → variants, pathways, GO terms, selection stats
- `graphpop lookup pathway "Cardiac repolarization"` → gene list, Fst
- `graphpop lookup variant chr22:16050075:A:G` → full annotation record
- `graphpop lookup region chr6 9000000 9600000` → genes + stats in region

### ★ Graph Neighborhood (`graphpop neighbors`)
Explores the graph structure around a gene or variant.
- `graphpop neighbors KCNE1 --hops 2` → genes sharing pathways with KCNE1
- Useful for hypothesis generation: "what other genes are functionally related?"

### ★ Multi-Statistic Convergence (`graphpop converge`)
Finds loci where multiple independently computed statistics agree.
- `graphpop converge --stats ihs,xpehh,h12 --thresholds 2,2,0.3 --pop EUR`
- Returns regions with gene annotations — the key GraphPop novelty.

### ★ Gene Ranking (`graphpop rank-genes`)
Ranks genes by composite selection evidence.
- `graphpop rank-genes --pop GJ-tmp --top 50 -o top_genes.tsv`
- Combines max |iHS|, max |XP-EHH|, max H12, mean Fst per gene.

### ★ Population Comparison (`graphpop compare`)
Statistical comparison between two populations.
- `graphpop compare EUR AFR chr22 --stat pi` → per-window delta with significance

### ★ Batch Execution (`graphpop batch`)
Run any command across multiple populations/chromosomes with parallelism.
- `graphpop batch diversity --pops EUR,AFR,EAS --chrs chr1,chr2 --workers 4`
- More flexible than `run-all` for custom analysis combinations.

### ★ Automated Report (`graphpop report`)
One-command summary of an entire database.
- `graphpop report --database rice3k -o report.html`
- Generates: dataset overview, diversity table, Fst heatmap, top selection peaks, QC summary, auto-generated figures.

### ★ Data Extraction (`graphpop extract`)
Flexible variant/genotype/sample extraction with filters.
- `graphpop extract variants --chr chr22 --pop EUR --consequence missense_variant`
- `graphpop extract genotypes --chr chr22 --region 16000000-17000000 --pop EUR`
- `graphpop extract samples --pop EUR` → per-sample statistics

### ★ BED Export (`graphpop export-bed`)
Export high-signal regions in BED format for bedtools integration.
- `graphpop export-bed --stat fst --threshold 0.5 --pop EUR`
- Enables interoperability with the bedtools/UCSC ecosystem.

### ★ Database Inventory (`graphpop inventory`)
Show what analyses have been run on this database.
- Which populations? Which chromosomes? What's persisted (iHS, XP-EHH, windows)?
- What annotations are loaded (VEP, Reactome, GO, ancestral)?
- More detailed than `db info` — answers "what can I query?"

### ★ Additional Plot Types
- `graphpop plot gene-zoom` — multi-track regional plot at a gene (Fst + iHS + H12 + gene model)
- `graphpop plot chromosome` — chromosome-wide ideogram with statistics tracks
- `graphpop plot pop-tree` — phylogenetic tree from Fst matrix
- `graphpop plot pca-scatter` — PCA scatter from allele frequencies
- `graphpop plot heatmap` — general-purpose heatmap from any matrix TSV

---

## End-to-End Workflow

```bash
# ─── First-time setup ───
graphpop setup --password mypass
graphpop start

# ─── Import data ───
graphpop import --vcf rice.vcf.gz --panel panel.txt \
    --database rice3k --vep annotations.vcf \
    --pathways plant_reactome.tsv
graphpop inventory                            # What was imported?

# ─── Explore the graph ───
graphpop lookup gene GW5                      # What do we know about GW5?
graphpop lookup pathway "Starch biosynthesis"  # What genes are in this pathway?
graphpop lookup region Chr5 5200000 5400000   # What's in the GW5 region?

# ─── Quick analysis ───
graphpop diversity Chr1 1 43270923 GJ-tmp     # Diversity → stdout
graphpop diversity Chr1 1 43270923 GJ-tmp \
    --consequence missense_variant -o piN.tsv # Conditioned diversity

# ─── Full-genome run ───
graphpop run-all -d results/                  # All pops × all chrs
graphpop aggregate -d results/ -o tables/     # Summary tables
graphpop report -o rice3k_report.html         # Automated summary

# ─── Visualize ───
graphpop plot diversity-bar results/diversity/ -o fig_diversity.png
graphpop plot fst-heatmap results/divergence/ -o fig_fst.png
graphpop plot pop-tree results/divergence/ -o fig_tree.png
graphpop plot manhattan results/ihs/GJ-tmp.tsv --stat ihs -o fig_ihs.png
graphpop plot gene-zoom GW5 --pop GJ-tmp -o fig_gw5.png

# ─── Selection scan ───
graphpop ihs Chr1 GJ-tmp --persist
graphpop filter ihs Chr1 GJ-tmp \
    --consequence missense_variant \
    --min-score 2.0 -o ihs_sweeps.tsv

# ─── Multi-statistic convergence ───
graphpop converge --stats ihs,xpehh,h12 \
    --thresholds 2.0,2.0,0.3 \
    --pop GJ-tmp -o convergent_signals.tsv
graphpop rank-genes --pop GJ-tmp --top 50 -o top_genes.tsv

# ─── Compare populations ───
graphpop compare GJ-tmp GJ-trp Chr1 --stat pi -o delta_pi.tsv

# ─── Extract and export ───
graphpop extract variants --chr Chr5 --start 5200000 --end 5400000 \
    --pop GJ-tmp --fields pos,ref,alt,af,fst,ihs -o gw5_variants.tsv
graphpop export-bed --stat fst --threshold 0.5 --pop GJ-tmp -o high_fst.bed
graphpop export-windows --min-fst 0.5 -o outlier_windows.tsv

# ─── Batch operations ───
graphpop batch diversity --pops GJ-tmp,GJ-trp,XI-1A \
    --chrs Chr1,Chr2,Chr3 --workers 3 -d batch_results/

# ─── Share database ───
graphpop dump --database rice3k -o rice3k_v1.dump

# ─── Database management ───
graphpop db list
graphpop db info
graphpop validate
```

## Statistics

| Category | Implemented | Proposed | Total |
|----------|------------|----------|-------|
| Setup & Server | 4 | — | 4 |
| Database Management | 8 | — | 8 |
| Configuration | 5 | 1 (inventory) | 6 |
| Procedures | 12 | — | 12 |
| Conditioned Analysis | 1 (filter) | — | 1 |
| Annotation Lookup | — | 5 (lookup ×4, neighbors) | 5 |
| Multi-Stat Integration | — | 3 (converge, rank-genes, compare) | 3 |
| Orchestration | 3 | 2 (batch, report) | 5 |
| Data Extraction | 1 (export-windows) | 4 (extract ×3, export-bed) | 5 |
| Visualization | 6 | 5 (gene-zoom, chromosome, pop-tree, pca-scatter, heatmap) | 11 |
| Utilities | 1 (query) | — | 1 |
| **Total** | **41** | **20** | **61** |

## Implementation Priority

### Phase 1 — High value, enables core analysis workflows
1. `lookup` (gene, pathway, variant, region) — graph exploration without Cypher
2. `converge` — the signature GraphPop analysis (multi-statistic convergence)
3. `inventory` — users need to know what's been computed
4. `plot gene-zoom` — the most requested genomics visualization
5. `rank-genes` — the deliverable users want (ranked gene list)

### Phase 2 — Expands utility and interoperability
6. `extract variants` / `extract samples` — flexible data export
7. `export-bed` — bedtools interoperability
8. `batch` — flexible parallelized execution
9. `compare` — population comparison
10. `plot pop-tree` — population phylogeny

### Phase 3 — Polish and automation
11. `report` — automated analysis summary
12. `plot chromosome` — ideogram visualization
13. `plot pca-scatter` — PCA from frequencies
14. `plot heatmap` — general-purpose heatmap
15. `extract genotypes` — genotype matrix export
16. `neighbors` — graph neighborhood exploration
