# GraphPop CLI

Command-line interface for **GraphPop** — a graph database-native population genomics engine that reduces summary statistic complexity from O(V×N) to O(V×K), independent of sample count.

## Quick Start

```bash
pip install graphpop-cli
graphpop setup --password mypass    # Downloads Neo4j + procedures plugin
graphpop start                      # Start the database
```

**Prerequisites:** Python 3.10+, Java 21+ (for Neo4j runtime).

## Features

- **60 commands** across 11 functional domains
- **12 population genetics procedures**: diversity, Fst, SFS, iHS, XP-EHH, nSL, ROH, Garud's H, LD, genome scan, pop summary, joint SFS
- **Annotation conditioning**: `--consequence`, `--pathway`, `--gene` flags on any procedure
- **Persistent analytical records**: `--persist` writes results to graph nodes
- **Publication-ready plots**: 11 visualization types following Nature Methods guidelines

## Usage

```bash
# Population diversity
graphpop diversity chr1 1 50000000 EUR -o diversity.tsv

# Annotation-conditioned analysis
graphpop diversity chr1 1 43270923 GJ-tmp --consequence missense_variant

# Selection scan
graphpop ihs chr22 EUR --persist -o ihs.tsv

# Multi-statistic convergence
graphpop converge --stats ihs,xpehh,h12 --thresholds 2,2,0.3 --pop EUR
```

## Documentation

- [Full documentation](https://github.com/jfmao/GraphPop)
- [Rice 3K vignette](https://github.com/jfmao/GraphPop/blob/main/graphpop-cli/vignettes/rice-3k-analysis.md)
- [Human 1000G vignette](https://github.com/jfmao/GraphPop/blob/main/graphpop-cli/vignettes/human-1000g-analysis.md)

## License

MIT
