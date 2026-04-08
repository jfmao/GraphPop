# slurm_setup_neo4j.sh

## One-Time Neo4j Setup on a SLURM Cluster Node

## Description

Downloads, configures, and starts Neo4j on a SLURM-managed cluster node.
Deploys the GraphPop stored-procedures plugin, sets memory parameters,
configures remote access, and optionally symlinks the data directory to
local SSD or scratch storage for performance. This script should be run
once per cluster before any ingestion or analysis jobs.

Scheduler: **SLURM**

## Usage

```bash
sbatch slurm_setup_neo4j.sh

# Or interactively:
srun --cpus-per-task=4 --mem=64G --time=2:00:00 --pty bash
bash slurm_setup_neo4j.sh
```

## SBATCH Directives

| Directive | Value | Description |
|---|---|---|
| `--job-name` | `graphpop-setup` | Job name in SLURM queue |
| `--cpus-per-task` | `4` | CPU cores allocated |
| `--mem` | `64G` | Memory allocation |
| `--time` | `2:00:00` | Wall-time limit |
| `--output` | `graphpop_setup_%j.log` | Log file (`%j` = job ID) |

## Environment Variables

| Variable | Default | Description |
|---|---|---|
| `NEO4J_HOME` | `$HOME/neo4j` | Neo4j installation directory |
| `PAGECACHE` | `20g` | Neo4j page cache size |
| `HEAP` | `4g` | Neo4j JVM heap size |
| `GRAPHPOP_PASSWORD` | `graphpop` | Initial database password |
| `PLUGIN_JAR` | `graphpop-procedures/target/graphpop-procedures-0.1.0-SNAPSHOT.jar` | Path to compiled procedures JAR |

## Details

The script performs the following steps:

1. Loads Java 21 and activates the `graphpop` conda environment.
2. Runs `graphpop setup` to download Neo4j and configure memory settings.
3. Deploys the GraphPop procedures JAR to the plugins directory (if found).
4. Detects fast local storage (`/local/$USER` or `/scratch/$USER`) and
   symlinks the Neo4j data directory there for improved I/O performance.
5. Enables remote bolt connections by setting
   `server.default_listen_address=0.0.0.0` in `neo4j.conf`.
6. Starts Neo4j and verifies it is running.
7. Prints the `GRAPHPOP_URI` connection string for use by other nodes.

If the plugin JAR is not found, setup proceeds without it and prints a
warning with build instructions.

## Prerequisites

- Java 21+ available via `module load java/21` (required by Neo4j runtime)
- `graphpop` CLI installed (`pip install graphpop-cli`)
- The procedures plugin JAR is automatically downloaded during `graphpop setup` (no Maven build required)

## Examples

```bash
# Basic setup with defaults
sbatch slurm_setup_neo4j.sh

# Custom memory and password
PAGECACHE=30g HEAP=8g GRAPHPOP_PASSWORD=secret sbatch slurm_setup_neo4j.sh

# Custom Neo4j location on fast storage
NEO4J_HOME=/scratch/$USER/neo4j sbatch slurm_setup_neo4j.sh

# Interactive setup for debugging
srun --cpus-per-task=4 --mem=64G --time=2:00:00 --pty bash
GRAPHPOP_PASSWORD=mypass bash slurm_setup_neo4j.sh
```

## Output

After completion, the log reports the hostname and connection URI:

```
Database node: compute-042
Connect from other nodes with:
  export GRAPHPOP_URI=bolt://compute-042:7687
  export GRAPHPOP_PASSWORD=graphpop
```

Save the URI for use in subsequent analysis jobs.

## See Also

- [slurm-load-csv](slurm-load-csv.md) -- bulk import after setup
- [slurm-interactive](slurm-interactive.md) -- interactive session with Neo4j
- [cluster-index](cluster-index.md) -- overview of all cluster scripts
- [cluster-guide.md](../../cluster-guide.md) -- full deployment guide
