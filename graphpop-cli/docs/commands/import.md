# graphpop import

## Title

Import VCF Data into a Neo4j Graph Database

## Description

Orchestrates the full import pipeline: parses a VCF file with a population panel
to generate CSV files, runs `neo4j-admin database import` for bulk loading, and
optionally loads functional annotations (VEP consequences, pathways, GO terms,
ancestral alleles). The resulting database is a fully populated GraphPop
variation graph ready for analysis.

## Usage

```
graphpop import [OPTIONS]
```

## Arguments

This command takes no positional arguments. All inputs are specified as options.

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `--vcf` | PATH | *(required)* | Input VCF file. Bgzipped (`.vcf.gz`) recommended for large files. |
| `--panel` | PATH | *(required)* | Population panel file. Tab-separated with columns `sample_id` and `population`. |
| `--database` | TEXT | *(required)* | Name for the Neo4j database. User-specified; becomes the active database in config. |
| `--vep` | PATH | *(none)* | VEP or SnpEff annotation file with variant consequence predictions. |
| `--pathways` | PATH | *(none)* | Reactome or Plant Reactome pathway annotation file (TSV). |
| `--go-terms` | PATH | *(none)* | GO term annotation file in UniProt GOA format. |
| `--ancestral` | PATH | *(none)* | Ancestral allele FASTA from Ensembl EPO alignments. Required for unfolded SFS and Fay & Wu's H. |
| `--csv-dir` | PATH | `/tmp/graphpop_csv_<database>` | Directory for intermediate CSV files. Useful for debugging or reusing CSVs. |
| `--neo4j-home` | PATH | *(from config)* | Neo4j installation directory. |
| `--threads` | INT | `4` | Number of threads for CSV generation. |
| `--skip-csv` | FLAG | `false` | Skip Step 1 (CSV generation). Use existing CSVs in `--csv-dir`. |
| `--skip-import` | FLAG | `false` | Skip Step 2 (bulk import). Generate CSVs only. |
| `--skip-annotations` | FLAG | `false` | Skip Step 3 (annotation loading). |

## Value

No return value. Side effects:

- CSV files generated in `--csv-dir`.
- Neo4j database created and populated via `neo4j-admin database import`.
- Functional annotations loaded via Cypher transactions (if annotation files
  provided).
- `~/.graphpop/config.yaml` updated to use the new database.

## Details

### The Three-Step Pipeline

**Step 1: CSV Generation** (`--skip-csv` to bypass)

Parses the VCF and population panel to produce Neo4j bulk-import CSV files:

| CSV File | Content |
|---|---|
| `variants_*.csv` | Variant nodes with `chr`, `pos`, `ref`, `alt`, `pop_ids[]`, `ac[]`, `an[]`, `af[]`, `gt_packed` |
| `samples.csv` | Sample nodes with `sampleId` |
| `populations.csv` | Population nodes with `populationId`, `n_samples` |
| `chromosomes.csv` | Chromosome nodes with `chromosomeId`, `length` |
| `next_*.csv` | NEXT relationships linking variants in positional order |
| `on_chromosome_*.csv` | ON_CHROMOSOME relationships from variants to chromosomes |
| `in_population.csv` | IN_POPULATION relationships from samples to populations |

The command first tries to import the `graphpop_import` Python package directly.
If not installed, it falls back to running the parser as a subprocess.

**Step 2: Bulk Import** (`--skip-import` to bypass)

Runs `neo4j-admin database import full` with the generated CSV header and data
files. This is the fastest way to load large graphs into Neo4j.

**Important:** Neo4j must be stopped before bulk import. If the database already
exists, you will be prompted to confirm overwriting.

**Step 3: Annotation Loading** (`--skip-annotations` to bypass)

Loads functional annotations into the existing graph via Cypher transactions.
This step requires Neo4j to be running. Annotations create:

- `Gene` nodes and `HAS_CONSEQUENCE` edges (from VEP/SnpEff)
- `Pathway` nodes and `IN_PATHWAY` edges (from Reactome)
- `GOTerm` nodes and `HAS_GO_TERM` edges (from GOA)
- `ancestral_allele` and `is_polarized` properties on Variant nodes (from
  ancestral FASTA)

### Panel File Format

The panel file is a tab-separated text file with at least two columns:

```
sample	pop
HG00096	GBR
HG00097	GBR
NA18486	YRI
NA18487	YRI
```

The first row can be a header (detected automatically) or data. Column names
are flexible; the first column is treated as sample ID and the second as
population ID.

### Database Naming Conventions

Choose descriptive names that encode the organism and scope:

- `human_chr22` -- Human chromosome 22
- `rice3k` -- Rice 3K full genome
- `sim_bottleneck_100k` -- Simulated bottleneck scenario

## Examples

### Minimal import (VCF + panel only)

```bash
graphpop stop
graphpop import \
    --vcf data/raw/1000G/chr22.vcf.gz \
    --panel data/raw/1000G/integrated_call_samples_v3.20130502.ALL.panel \
    --database human_chr22
graphpop start
graphpop db info
```

### Full human 1000 Genomes import with annotations

```bash
graphpop stop

graphpop import \
    --vcf data/raw/1000G/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
    --panel data/raw/1000G/integrated_call_samples_v3.20130502.ALL.panel \
    --database human_chr22 \
    --vep data/raw/1000G/chr22_vep.vcf \
    --pathways data/raw/reactome/UniProt2Reactome_All_Levels.txt \
    --go-terms data/raw/goa/goa_human.gaf \
    --ancestral data/raw/ancestral/homo_sapiens_ancestor_22.fa \
    --threads 8

graphpop start
graphpop db info
graphpop validate
```

### Rice 3K import with Plant Reactome

```bash
graphpop stop

graphpop import \
    --vcf data/raw/rice3k/rice3k_chr1.vcf.gz \
    --panel data/raw/rice3k/rice3k_panel.txt \
    --database rice3k \
    --vep data/raw/rice3k/rice3k_vep.vcf \
    --pathways data/raw/plant_reactome/plant_reactome_oryza.tsv \
    --ancestral data/raw/ancestral/oryza_rufipogon_ancestor.fa \
    --csv-dir data/processed/rice_csv \
    --threads 8

graphpop start
graphpop validate
graphpop diversity chr1 1 43270923 GJ-tmp
```

### Two-stage import (CSV on one machine, import on another)

```bash
# Machine A: generate CSVs only
graphpop import \
    --vcf large_dataset.vcf.gz \
    --panel panel.txt \
    --database mydb \
    --csv-dir /shared/csv_output \
    --skip-import --skip-annotations \
    --threads 16

# Machine B: import from existing CSVs
graphpop import \
    --vcf large_dataset.vcf.gz \
    --panel panel.txt \
    --database mydb \
    --csv-dir /shared/csv_output \
    --skip-csv
```

## See Also

- `graphpop setup` -- Install and configure Neo4j before importing.
- `graphpop start` / `graphpop stop` -- Server lifecycle (must stop before bulk import).
- `graphpop db info` -- Verify node/edge counts after import.
- `graphpop validate` -- Check database integrity and completeness.
- `graphpop dump` -- Back up the imported database.
- `graphpop run-all` -- Run full-genome analysis after import.
