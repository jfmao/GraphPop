# slurm_interactive.sh

## Interactive Session Setup with Neo4j (SLURM)

## Description

Sets up an interactive GraphPop session with Neo4j running on a SLURM
compute node. Unlike the other cluster scripts, this file is **sourced**
(not submitted with sbatch). It starts Neo4j, validates the filesystem, and
provides usage instructions. Use this for exploratory analysis, debugging,
and development on the cluster.

Scheduler: **SLURM** (interactive `srun` session)

## Usage

```bash
# Step 1: Allocate an interactive session
srun --nodes=1 --cpus-per-task=16 --mem=128G --time=8:00:00 --pty bash

# Step 2: Source the script (do NOT sbatch it)
source slurm_interactive.sh [neo4j_home] [data_dir]
```

## Arguments

| Position | Name | Required | Default | Description |
|---|---|---|---|---|
| 1 | `neo4j_home` | No | `$HOME/neo4j` | Neo4j installation directory |
| 2 | `data_dir` | No | `/scratch/$USER/graphpop_db` | Neo4j data directory |

## Details

The script performs the following steps when sourced:

1. Prints the current Neo4j home, data directory, and compute node hostname.
2. Validates the data directory filesystem using
   `graphpop_import.cluster.check_neo4j_data_dir`. Prints a warning if a
   network filesystem is detected (NFS, Lustre) but does not abort.
3. Starts Neo4j via `graphpop_import.cluster.start_neo4j` with a 120-second
   timeout. Reports the bolt port on success.
4. Prints example commands for running analyses and stopping Neo4j.

### Why Source Instead of Submit

Sourcing keeps the Neo4j process in the current interactive shell session.
If submitted via sbatch, the job would start Neo4j and then exit, leaving
the database running in the background with no interactive access. By
sourcing, the user retains a live shell on the same node where Neo4j is
running.

### Stopping Neo4j

When the interactive session is complete, stop Neo4j before exiting:

```bash
python -c "from graphpop_import.cluster import stop_neo4j; stop_neo4j('$HOME/neo4j')"
```

Or simply exit the `srun` session -- the process will be terminated when
SLURM reclaims the allocation.

### Typical Interactive Workflow

```bash
# Source the setup
source slurm_interactive.sh

# Run analyses interactively
python scripts/human_full_analysis.py --phase diversity
python scripts/human_interpret.py pinsps

# Explore with ad-hoc queries
graphpop query "MATCH (v:Variant) WHERE v.chr = 'chr22' RETURN count(v)"

# Stop Neo4j when done
python -c "from graphpop_import.cluster import stop_neo4j; stop_neo4j('$HOME/neo4j')"
```

## Examples

```bash
# Default locations
srun --cpus-per-task=16 --mem=128G --time=8:00:00 --pty bash
source slurm_interactive.sh

# Custom Neo4j installation
srun --cpus-per-task=16 --mem=128G --time=8:00:00 --pty bash
source slurm_interactive.sh /opt/neo4j /local/$USER/neo4j_db

# Short session for quick checks
srun --cpus-per-task=4 --mem=32G --time=1:00:00 --pty bash
source slurm_interactive.sh

# Long session for development
srun --cpus-per-task=16 --mem=128G --time=24:00:00 --pty bash
source slurm_interactive.sh ~/neo4j /scratch/$USER/graphpop_db
```

## See Also

- [slurm-setup-neo4j](slurm-setup-neo4j.md) -- initial Neo4j installation
- [slurm-analysis](slurm-analysis.md) -- batch analysis submission
- [cluster-index](cluster-index.md) -- overview of all cluster scripts
- [cluster-guide.md](../../cluster-guide.md) -- full deployment guide
