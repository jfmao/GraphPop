# GraphPop Human Interpretation — Restart Guide

## Checkpoint Status (as of 2026-03-18)

All progress is saved in `human_interpretation_results.json`. Each subcommand
skips already-completed calls on restart.

| Subcommand | Total Calls | Completed | Status |
|------------|-------------|-----------|--------|
| pinsps | 1,144 | 1,144 | DONE |
| annotate | — | — | DONE |
| fay_wu | 572 | 572 | DONE |
| usfs | 572 | 572 | DONE |
| divergence | 7,150 | 7,150 | DONE |
| tree | — | — | DONE |
| pbs | 132 | 44 | **88 remaining** |
| hwscan | ~132 | 0 | not started |
| roh_hmm | ~572 | 0 | not started |
| daf_enrichment | — | 0 | not started |
| gscan | 528 | 0 | not started |
| report | — | 0 | not started |

## Prerequisites

1. Neo4j running with the human 1000G database loaded
2. Gene annotations loaded (91,973 Gene nodes, 891,689 HAS_CONSEQUENCE edges)
3. `human_full_results.json` in the working directory (~1.7 GB)
4. `human_interpretation_results.json` in the working directory (checkpoint file)
5. Python with `neo4j` driver installed

## Quick Restart

```bash
cd /path/to/GraphPop

# If you moved the project, update these paths in human_interpret.py if needed:
#   RESULTS_FILE = "human_full_results.json"
#   OUTPUT_FILE = "human_interpretation_results.json"
# (both are relative to CWD)

# Verify Neo4j is accessible
cypher-shell -u neo4j -p graphpop "RETURN 1 AS ok;"

# Verify Gene nodes exist
cypher-shell -u neo4j -p graphpop "MATCH (g:Gene) RETURN count(g);"
# Should return 91,973

# Resume remaining subcommands (in order):
python scripts/human_interpret.py pbs          # ~88 calls remaining, ~5 hours
python scripts/human_interpret.py hwscan       # ~132 calls, ~9 hours
python scripts/human_interpret.py roh_hmm      # ~572 calls, hours
python scripts/human_interpret.py daf_enrichment  # fast
python scripts/human_interpret.py gscan        # 528 calls, ~35 hours (slowest)
python scripts/human_interpret.py report       # fast, generates markdown
```

Or use the batch script:
```bash
bash scripts/run_remaining_interpret.sh
```

## Estimated Times (per subcommand)

- **genome_scan-based** (pbs, hwscan, gscan): ~3-4 min/call (WRITE procedure)
- **diversity-based** (fay_wu, usfs, pinsps): ~3 min/pop (~1-2 hours total)
- **divergence**: ~9 sec/call with 6 threads (~18 hours for all 7,150)
- **roh_hmm**: depends on sample count

## If Changing Neo4j Credentials or URI

Edit the constants at the top of `scripts/human_interpret.py`:
```python
NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"
```

## Key Files to Preserve

- `human_interpretation_results.json` — checkpoint (resume depends on this)
- `human_full_results.json` — source data from full analysis
- `data/raw/1000g/full_genome_annotations/` — SnpEff annotation CSVs
  - `gene_nodes.csv` (91,973 genes)
  - `has_consequence_edges.csv` (90M edges, 9.2 GB)

## If Gene Annotations Need Reloading

```bash
python scripts/load_annotations_driver.py data/raw/1000g/full_genome_annotations/
# Loads ~892K functional edges (HIGH/MODERATE/LOW impact)
# Use --all-impacts for all 90M edges (slow, not needed for paper analyses)
```
