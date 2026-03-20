# GraphPop Cluster Deployment Guide

This guide covers deploying and using GraphPop on HPC clusters (SLURM, PBS/Torque) where users typically do not have root access or persistent services.

GraphPop shares its cluster deployment strategy with [GraphMana](https://github.com/your-org/graphmana), the companion data management platform. Both tools use the same Neo4j lifecycle management, filesystem checks, and job script patterns.

## Overview

GraphPop supports four deployment models on clusters, from simplest to most scalable:

| Model | Neo4j Lifecycle | Best For |
|-------|----------------|----------|
| Dedicated node | Persistent service | Labs with sysadmin support |
| Interactive job | Manual start/stop | Exploratory analysis |
| Batch job | Auto start/stop | Automated pipelines |
| Two-step split | CSV without Neo4j | Large imports, max parallelism |

## Model 1: Dedicated Node

A lab secures a persistent VM or node running Neo4j as a service. Users connect via Bolt from any cluster node.

```bash
# On the dedicated node (one-time setup)
python -c "
from graphpop_import.cluster import setup_neo4j
setup_neo4j('/opt/neo4j', memory_auto=True,
            plugin_jar='graphpop-procedures/target/graphpop-procedures-1.0.jar')
"

# Start Neo4j
python -c "
from graphpop_import.cluster import start_neo4j
start_neo4j('/opt/neo4j/neo4j-community-2026.01.4', wait=True)
"

# From any cluster node, run analyses against the remote Neo4j
python scripts/human_full_analysis.py --phase all
# (set NEO4J_URI=bolt://dedicated-node:7687 in the script)
```

## Model 2: Interactive Job

Request an interactive node, start Neo4j in user space, run commands, then stop.

```bash
# Request resources
srun --nodes=1 --cpus-per-task=16 --mem=128G --time=8:00:00 --pty bash

# Set up environment
module load java/21
source ~/graphpop-env/bin/activate

# Start Neo4j (or source the interactive helper)
source scripts/cluster/slurm_interactive.sh $HOME/neo4j /scratch/$USER/graphpop_db

# Run GraphPop analyses
python scripts/human_full_analysis.py --phase all
python scripts/human_interpret.py pinsps
python scripts/human_interpret.py gscan
python scripts/human_interpret.py report

# Stop Neo4j when done
python -c "from graphpop_import.cluster import stop_neo4j; stop_neo4j('$HOME/neo4j/neo4j-community-2026.01.4')"
```

## Model 3: Batch Job with Auto Start/Stop

Submit a SLURM job that automatically manages the Neo4j lifecycle:

```bash
sbatch scripts/cluster/slurm_ingest_single.sh data.vcf.gz populations.tsv
```

This starts Neo4j before the operation and stops it when done, even if an error occurs.

## Model 4: Two-Step Split (Recommended)

Separate CPU-intensive CSV generation from database loading. CSV generation needs no Neo4j and is embarrassingly parallel — ideal for cluster compute nodes.

### Step 1: Generate CSVs (any compute node)

```bash
sbatch scripts/cluster/slurm_prepare_csv.sh data.vcf.gz populations.tsv
```

This runs the GraphPop import pipeline to produce neo4j-admin import CSVs: variants (with packed genotypes), samples, populations, chromosomes, windows, and edges.

### Step 2: Load into Neo4j (needs local SSD)

```bash
sbatch scripts/cluster/slurm_load_csv.sh /scratch/$USER/graphpop_csv $HOME/neo4j
```

Uses `neo4j-admin import` for fast bulk loading. The filesystem check ensures the data directory is on local storage.

### Step 3: Run analyses

```bash
sbatch scripts/cluster/slurm_analysis.sh scripts/human_full_analysis.py --phase all
```

## Installation on a Cluster

### 1. Install Neo4j in User Space

```bash
# Load Java 21 (cluster-specific — check `module avail`)
module load java/21

# Download and configure Neo4j
python -c "
from graphpop_import.cluster import setup_neo4j
result = setup_neo4j(
    '$HOME/neo4j',
    data_dir='/scratch/$USER/graphpop_db',
    plugin_jar='graphpop-procedures/target/graphpop-procedures-1.0.jar',
    memory_auto=True,
)
print(f\"Installed: {result['neo4j_home']}\")
print(f\"Data dir:  {result['data_dir']}\")
print(f\"Java:      {result['java_version']}\")
"
```

This downloads Neo4j Community, extracts it, deploys the GraphPop procedures JAR, configures memory based on available RAM, and sets the data directory. No root access needed.

### 2. Install Python Dependencies

```bash
# With pip (user install)
pip install --user -e graphpop-import/
pip install neo4j cyvcf2 numpy

# Or with a virtual environment
python -m venv ~/graphpop-env
source ~/graphpop-env/bin/activate
pip install -e graphpop-import/
pip install neo4j cyvcf2 numpy
```

### 3. Build GraphPop Procedures JAR

```bash
module load java/21
cd graphpop-procedures
JAVA_HOME=/path/to/java-21 ./mvnw package -DskipTests
# JAR at: target/graphpop-procedures-1.0.jar
```

### 4. Verify Java Version

Neo4j 2026.x requires Java 21+. Check with:

```bash
java -version
```

If the default Java is too old, load the correct module:

```bash
module load java/21    # SLURM clusters
module load jdk/21     # Some PBS systems
```

## Filesystem Guidance

**Critical**: Neo4j performs extremely poorly on network filesystems due to random I/O patterns.

| Component | Filesystem | Notes |
|-----------|-----------|-------|
| Neo4j data directory | **Local SSD/scratch only** | `/scratch`, `/tmp/local`, `/local` |
| CSV output directory | Shared OK | Sequential writes |
| VCF input files | Shared OK | Sequential reads |
| Analysis output (JSON) | Shared OK | Sequential writes |

Use the filesystem check to verify:

```bash
python -c "
from graphpop_import.cluster import check_neo4j_data_dir
result = check_neo4j_data_dir('/scratch/\$USER/graphpop_db')
print(f\"Status: {'WARNING' if result['is_network'] else 'OK'} ({result['fs_type']})\")
if result['warning']:
    print(result['warning'])
"
```

## SLURM Job Scripts

Example scripts are provided in `scripts/cluster/`:

| Script | Purpose | CPUs | RAM | Time |
|--------|---------|------|-----|------|
| `slurm_prepare_csv.sh` | CSV generation (Step 1) | 16 | 64 GB | 4h |
| `slurm_load_csv.sh` | Neo4j load (Step 2) | 4 | 128 GB | 2h |
| `slurm_ingest_single.sh` | Combined import | 16 | 128 GB | 8h |
| `slurm_analysis.sh` | Run analyses | 8 | 64 GB | 12h |
| `slurm_interactive.sh` | Interactive session | — | — | — |

### PBS/Torque

| Script | Purpose |
|--------|---------|
| `pbs_prepare_csv.sh` | CSV generation (Step 1) |
| `pbs_analysis.sh` | Run analyses |

## GraphPop-Specific Considerations

### Two-Phase Workflow: Import + Analysis

Unlike simple data management tools, GraphPop's workflow has two distinct phases:

1. **Import phase**: VCF → CSV → Neo4j database (+ VEP annotations, pathways)
2. **Analysis phase**: Run 12 stored procedures via Cypher (diversity, Fst, SFS, iHS, XP-EHH, nSL, ROH, etc.)

Both phases require Neo4j, but CSV generation (Step 1 of import) does not.

### Procedures JAR Deployment

GraphPop's statistical procedures are Java stored procedures compiled into a JAR. The JAR must be in Neo4j's `plugins/` directory before starting Neo4j:

```bash
cp graphpop-procedures/target/graphpop-procedures-1.0.jar \
   $NEO4J_HOME/plugins/
```

The `setup_neo4j()` function accepts a `plugin_jar` parameter to handle this automatically.

### Memory Profile

GraphPop procedures are designed for low memory consumption:

- **FAST PATH** (diversity, Fst, SFS, etc.): ~150 MB constant, regardless of sample count
- **FULL PATH** (iHS, XP-EHH, nSL, LD): chunked processing, ~100-200 MB per chunk
- **Neo4j page cache**: proportional to database size (50% of RAM recommended)
- **JVM heap**: 4-8 GB typically sufficient (25% of RAM)

For a full human genome (3,202 samples, 85M variants, ~195 GB database):

| Component | Recommended |
|-----------|-------------|
| JVM heap | 4-8 GB |
| Page cache | 16-32 GB |
| Total RAM | 64-128 GB |
| Storage | 250+ GB local SSD |

### Analysis Resume Support

GraphPop analysis scripts (`human_full_analysis.py`, `rice_full_analysis.py`) support automatic resume:

- Results are saved incrementally to JSON
- Already-completed (population, chromosome, procedure) combos are skipped
- If a job times out, resubmit the same script — it picks up where it left off

```bash
# First run (times out after 4 hours)
sbatch --time=4:00:00 slurm_analysis.sh scripts/human_full_analysis.py --phase 1

# Resume (skips completed work)
sbatch --time=4:00:00 slurm_analysis.sh scripts/human_full_analysis.py --phase 1
```

## Resource Allocation

| Operation | CPUs | RAM | Time (chr22, 3K) | Time (WGS, 3K) |
|-----------|------|-----|-------------------|-----------------|
| prepare-csv | 16 | 64 GB | Minutes | 1-2 hours |
| load-csv (neo4j-admin) | 4 | 128 GB | Minutes | 10-30 min |
| Phase 1: per-pop stats | 8 | 64 GB | Minutes | 2-4 hours |
| Phase 2: XP-EHH pairs | 8 | 64 GB | Minutes | 4-8 hours |
| Phase 3: ancestral stats | 8 | 64 GB | Minutes | 2-4 hours |
| Phase 4: pop-specific | 8 | 64 GB | Minutes | 1-2 hours |
| Interpretation (piN/piS) | 8 | 64 GB | Minutes | 4-8 hours |

### Storage Requirements

| Samples | WGS Variants | Estimated DB Size |
|---------|-------------|-------------------|
| 100 | 85M | ~20 GB |
| 1,000 | 85M | ~50 GB |
| 3,200 | 85M | ~130-200 GB |
| 10,000 | 85M | ~400-550 GB |
| 50,000 | 85M | ~2-3 TB |

Request sufficient scratch space before starting a large import.

## Troubleshooting

**Neo4j extremely slow**: Data directory is on a network filesystem (NFS/Lustre/GPFS). Move to local SSD/scratch. Run the filesystem check to verify.

**Neo4j won't start**: Java 21 not found. Check `module load java/21` or equivalent. Also check that Bolt port 7687 is not blocked or in use.

**Procedures not found**: `graphpop.diversity` etc. not recognized. Ensure the GraphPop JAR is in `$NEO4J_HOME/plugins/` and Neo4j was restarted after deployment.

**Out of memory during analysis**: GraphPop procedures are designed for ~150 MB, but Neo4j page cache needs RAM for the database. Increase total RAM allocation or tune page cache size down.

**Out of disk**: Full human WGS (3,200 samples) produces ~195 GB database. Ensure `/scratch` allocation is sufficient.

**prepare-csv succeeds but load-csv fails**: CSV files may have been corrupted on shared storage. Verify file sizes match expectations. Re-run prepare-csv if needed.

**Permission denied on Neo4j binary**: Ensure Neo4j was installed with `setup_neo4j()` (sets correct permissions).

**Port already in use**: Another user may be running Neo4j on the same node. Edit `$NEO4J_HOME/conf/neo4j.conf` to change `server.bolt.listen_address` and `server.http.listen_address` to use different ports.

**Analysis script times out**: Use the resume support — resubmit the same script and it picks up where it left off.
