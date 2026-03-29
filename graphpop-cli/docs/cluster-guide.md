# Cluster Computing Guide

How to deploy and use GraphPop on HPC clusters (SLURM, PBS/Torque)

## Architecture: How GraphPop Works on a Cluster

GraphPop uses Neo4j, a **single-server** database. Unlike Spark-based tools
(Hail) that distribute computation across nodes, Neo4j runs on **one node**.
Analysis clients (`graphpop` CLI) connect to it via the bolt:// protocol.

```
┌─────────────────────────────────────────────────────────┐
│                    HPC Cluster                          │
│                                                         │
│   ┌──────────────────────┐    ┌──────────────────────┐  │
│   │  Database Node        │    │  Compute Nodes       │  │
│   │  (interactive/         │    │  (batch jobs)        │  │
│   │   dedicated)          │    │                      │  │
│   │                       │    │  graphpop diversity   │  │
│   │  Neo4j Server         │◄───│  graphpop ihs        │  │
│   │  pagecache: 20 GB     │    │  graphpop run-all    │  │
│   │  heap: 4 GB           │    │                      │  │
│   │                       │    │  Connects via         │  │
│   │  Local SSD/NVMe       │    │  bolt://<db-node>:7687│  │
│   │  (Neo4j data dir)     │    │                      │  │
│   └──────────────────────┘    └──────────────────────┘  │
│                                                         │
│   ┌──────────────────────┐                              │
│   │  Any Node             │                              │
│   │  (no Neo4j needed)   │                              │
│   │                       │                              │
│   │  graphpop-import      │  CSV generation is           │
│   │  (VCF → CSV)          │  embarrassingly parallel     │
│   └──────────────────────┘                              │
└─────────────────────────────────────────────────────────┘
```

**Key principle:** Neo4j runs on ONE node with fast local storage. Analysis
jobs run on ANY node and connect via network. CSV generation needs no database.

## Resource Requirements

| Component | CPUs | Memory | Storage | Notes |
|-----------|------|--------|---------|-------|
| Neo4j server | 4-8 | 32-64 GB | 50-200 GB SSD | pagecache 20G + heap 4G |
| CSV generation | 4-16 | 32-64 GB | 50-200 GB scratch | I/O bound, parallelizable by chromosome |
| neo4j-admin import | 4 | 64-128 GB | Same as Neo4j | Runs once, offline |
| Analysis (per job) | 4-8 | 8-32 GB | Minimal | Client only, ~160 MB for GraphPop |
| Haplotype export | 4 | 32-64 GB | 2-10 GB per dataset | One-time |

## Step-by-Step Deployment

### 1. Install GraphPop on the cluster

```bash
# Load required modules (cluster-specific)
module load java/21      # or java/openjdk-21
module load miniconda3   # or anaconda3

# Create conda environment
conda create -n graphpop python=3.12 -y
conda activate graphpop

# Clone and install
git clone https://github.com/jfmao/GraphPop.git
cd GraphPop

# Build Java procedures
cd graphpop-procedures && mvn package -DskipTests && cd ..

# Install Python components
pip install -e graphpop-cli/
pip install -e graphpop-import/

# Verify
graphpop --version
```

### 2. Set up Neo4j on a database node

Neo4j should run on a node with **fast local storage** (SSD or NVMe), not on
the shared network filesystem (NFS/Lustre). Most clusters provide local scratch
at `/tmp`, `/scratch`, or `/local`.

**Option A: Interactive session (recommended for first-time setup)**

```bash
# Request an interactive session with local SSD
srun --nodes=1 --cpus-per-task=8 --mem=64G \
     --time=12:00:00 --tmp=200G --pty bash

# Set up Neo4j
graphpop setup --neo4j-home $HOME/neo4j --password mypassword \
    --pagecache 20g --heap 4g \
    --deploy-plugin graphpop-procedures/target/graphpop-procedures-0.1.0-SNAPSHOT.jar

# Move Neo4j data to local SSD for performance
mkdir -p /local/$USER/neo4j_data
ln -sf /local/$USER/neo4j_data $HOME/neo4j/data

# Start Neo4j
graphpop start

# Note the hostname for client connections
echo "Neo4j running on: $(hostname)"
# e.g., Neo4j running on: node042
```

**Option B: Dedicated database node (if your cluster supports it)**

Some clusters have persistent interactive nodes or database servers. Ask your
sysadmin about allocating a dedicated node for Neo4j.

### 3. Import data

CSV generation runs on **any node** (no Neo4j needed):

```bash
# Submit CSV generation as a SLURM job
sbatch scripts/cluster/slurm_prepare_csv.sh /data/rice.vcf.gz /data/panel.txt

# Or use the graphpop CLI
graphpop import --vcf /data/rice.vcf.gz --panel /data/panel.txt \
    --database rice3k --skip-import --csv-dir /scratch/$USER/csv_out
```

For multi-chromosome VCFs, parallelize by chromosome:

```bash
#!/bin/bash
#SBATCH --job-name=graphpop-csv
#SBATCH --array=1-22
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00

CHR=$SLURM_ARRAY_TASK_ID

graphpop import \
    --vcf /data/ALL.chr${CHR}.vcf.gz \
    --panel /data/panel.txt \
    --database human1000g \
    --csv-dir /scratch/$USER/csv_chr${CHR} \
    --skip-import
```

Then bulk-import all CSVs on the database node:

```bash
# On the database node (where Neo4j runs)
graphpop stop  # Neo4j must be stopped for bulk import

graphpop import \
    --database human1000g \
    --csv-dir /scratch/$USER/csv_merged \
    --skip-csv  # CSVs already generated

graphpop start
graphpop validate
```

### 4. Run analyses

Analysis jobs connect to Neo4j via bolt:// from any compute node:

```bash
#!/bin/bash
#SBATCH --job-name=graphpop-analysis
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00

module load java/21
conda activate graphpop

# Point to the database node
export GRAPHPOP_URI=bolt://node042:7687
export GRAPHPOP_PASSWORD=mypassword
export GRAPHPOP_DATABASE=rice3k

# Run analysis
graphpop diversity Chr1 1 43270923 GJ-tmp -o results/diversity_GJtmp_Chr1.tsv
graphpop ihs Chr1 GJ-tmp --persist -o results/ihs_GJtmp_Chr1.tsv
```

For full-genome analysis, use array jobs:

```bash
#!/bin/bash
#SBATCH --job-name=graphpop-fullgenome
#SBATCH --array=1-12        # 12 chromosomes for rice
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00

module load java/21
conda activate graphpop

export GRAPHPOP_URI=bolt://node042:7687
export GRAPHPOP_PASSWORD=mypassword
export GRAPHPOP_DATABASE=rice3k

CHR="Chr${SLURM_ARRAY_TASK_ID}"

# Per-population analysis for this chromosome
for POP in GJ-tmp GJ-trp XI-1A XI-1B XI-2 XI-3 XI-adm cA-Aus cB-Bas GJ-adm GJ-sbtrp admix; do
    graphpop diversity $CHR 1 50000000 $POP \
        -o results/diversity/${POP}_${CHR}.tsv
    graphpop ihs $CHR $POP --persist \
        -o results/ihs/${POP}_${CHR}.tsv
    graphpop garud-h $CHR $POP 100000 50000 \
        -o results/garud/${POP}_${CHR}.tsv
done
```

Or use `graphpop run-all` from a single job:

```bash
#!/bin/bash
#SBATCH --job-name=graphpop-runall
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

export GRAPHPOP_URI=bolt://node042:7687
export GRAPHPOP_PASSWORD=mypassword

graphpop run-all --phase all -d /scratch/$USER/results/ --resume
```

### 5. Post-analysis (any node, no Neo4j needed)

Aggregation, plotting, and reporting don't need Neo4j:

```bash
# These work from the results directory alone
graphpop aggregate -d /scratch/$USER/results/ -o tables/
graphpop plot diversity-bar results/diversity/ -o figures/fig_diversity.png
graphpop plot fst-heatmap results/divergence/ -o figures/fig_fst.png
graphpop report -o analysis_report.html  # This one needs Neo4j for live queries
```

### 6. Share results

```bash
# Dump the database (on the database node, Neo4j must be stopped)
graphpop stop
graphpop dump --database rice3k -o /data/shared/rice3k_v1.dump
graphpop start
```

## SLURM Job Templates

Ready-to-use templates are in `scripts/cluster/`:

| Script | Purpose | Neo4j needed? |
|--------|---------|---------------|
| `slurm_prepare_csv.sh` | VCF → CSV generation | No |
| `slurm_load_csv.sh` | neo4j-admin bulk import | Neo4j stopped |
| `slurm_analysis.sh` | Run analysis scripts | Neo4j running |
| `slurm_ingest_single.sh` | Combined import (CSV + load) | Starts/stops Neo4j |
| `slurm_interactive.sh` | Interactive session setup | Starts Neo4j |

## PBS/Torque Equivalents

| Script | Purpose |
|--------|---------|
| `pbs_prepare_csv.sh` | VCF → CSV generation |
| `pbs_analysis.sh` | Run analysis scripts |

Adapt SLURM scripts to PBS by changing headers:

```bash
# SLURM                          # PBS
#SBATCH --job-name=NAME          #PBS -N NAME
#SBATCH --cpus-per-task=8        #PBS -l ncpus=8
#SBATCH --mem=64G                #PBS -l mem=64GB
#SBATCH --time=12:00:00          #PBS -l walltime=12:00:00
#SBATCH --array=1-22             #PBS -J 1-22
$SLURM_ARRAY_TASK_ID             $PBS_ARRAY_INDEX
```

## Storage Strategy

```
$HOME/                          # Persistent, small
├── neo4j/                      # Neo4j installation (~500 MB)
├── GraphPop/                   # Code repository (~50 MB)
└── .graphpop/config.yaml       # Connection config

/scratch/$USER/                 # Fast scratch, large
├── neo4j_data/                 # Neo4j database (50-200 GB) ← LOCAL SSD
├── csv_out/                    # Import CSVs (temporary, 50-200 GB)
├── results/                    # Analysis output (1-10 GB)
└── hap_cache/                  # Haplotype arrays (2-10 GB)

/data/shared/                   # Shared storage
├── vcf/                        # Input VCF files
├── annotations/                # VEP, Reactome, ancestral alleles
└── dumps/                      # Database snapshots for sharing
```

**Critical:** Put Neo4j data on local SSD (`/local`, `/tmp`, or `$TMPDIR`),
not on NFS/Lustre. Neo4j performance degrades dramatically on network
filesystems due to random I/O patterns.

## Troubleshooting

### Neo4j won't start
- Check Java version: `java -version` (need 21+)
- Check port availability: `ss -tlnp | grep 7687`
- Check data directory permissions: Neo4j needs write access
- Check logs: `$HOME/neo4j/logs/neo4j.log`

### Connection refused from compute nodes
- Verify Neo4j is listening on all interfaces, not just localhost:
  Edit `$HOME/neo4j/conf/neo4j.conf`:
  ```
  server.default_listen_address=0.0.0.0
  ```
- Check cluster firewall rules between nodes
- Verify the hostname: `graphpop --uri bolt://ACTUAL_HOSTNAME:7687 status`

### Slow performance on NFS
- Move Neo4j data to local SSD:
  ```bash
  rsync -a $HOME/neo4j/data/ /local/$USER/neo4j_data/
  ln -sf /local/$USER/neo4j_data $HOME/neo4j/data
  ```

### Out of memory
- Reduce Neo4j pagecache: edit `neo4j.conf`, set `server.memory.pagecache.size=8g`
- For analysis jobs, GraphPop client needs only ~160 MB; reduce `--mem` allocation
- For iHS/XP-EHH chunking, the 5 MB window + 2 MB margin keeps memory under 200 MB

### Job array vs run-all
- **Job arrays** (one job per chromosome): better for cluster utilization, handles failures per-chromosome, easier to resubmit failed tasks
- **run-all** (single long job): simpler, has built-in resume, better for small datasets
- Recommendation: use job arrays for full-genome analysis, run-all for targeted investigations

## Example: Complete Rice 3K Workflow on SLURM

```bash
# 1. Set up (once, interactive)
srun --nodes=1 --cpus-per-task=8 --mem=64G --time=2:00:00 --pty bash
graphpop setup --password mypass --deploy-plugin graphpop-procedures/target/*.jar
graphpop start
echo "DB_NODE=$(hostname)" > ~/.graphpop/cluster_info
exit

# 2. Import (batch job)
sbatch scripts/cluster/slurm_prepare_csv.sh /data/rice.vcf.gz /data/panel.txt
# Wait for completion, then on the DB node:
sbatch scripts/cluster/slurm_load_csv.sh /scratch/$USER/csv_out $HOME/neo4j

# 3. Full-genome analysis (array job, 12 chromosomes)
export GRAPHPOP_URI=bolt://$(cat ~/.graphpop/cluster_info | grep DB_NODE | cut -d= -f2):7687
sbatch --array=1-12 scripts/cluster/slurm_fullgenome_array.sh

# 4. Aggregate and visualize (any node, after all array jobs complete)
graphpop aggregate -d /scratch/$USER/results/ -o tables/
graphpop plot diversity-bar results/diversity/ -o figures/fig_diversity.png
graphpop plot fst-heatmap results/divergence/ -o figures/fig_fst.png

# 5. Dump and share
graphpop dump --database rice3k -o /data/shared/rice3k_v1.dump
```
