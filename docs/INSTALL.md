# GraphPop Installation Guide

GraphPop runs on Linux and macOS without admin privileges. Choose the
installation method that suits your environment.

## Quick Install (recommended)

```bash
pip install graphpop-cli                  # Install the CLI
graphpop setup --password mypassword      # Downloads Neo4j + procedures plugin
graphpop start                            # Start the database
graphpop status                           # Verify
```

After installation:

```bash
graphpop --version
graphpop doctor                           # Full health check
```

## Method 1: pip install (requires Python 3.10+ and Java 21+)

```bash
# Create a conda environment with Java (needed by Neo4j runtime)
conda create -n graphpop -c conda-forge python=3.12 openjdk=21 -y
conda activate graphpop

# Install GraphPop CLI
pip install graphpop-cli

# Set up Neo4j (downloads and configures, no admin needed)
graphpop setup --password mypassword
graphpop start
```

## Method 2: Conda install

```bash
git clone https://github.com/jfmao/GraphPop.git
cd GraphPop
conda env create -f environment.yml       # Creates 'graphpop' environment with Java 21 + Python
conda activate graphpop
pip install -e graphpop-cli -e graphpop-import -e graphpop-mcp
graphpop setup --password mypassword
graphpop start
```

## Method 3: Docker

```bash
docker pull ghcr.io/jfmao/graphpop:latest
docker run -p 7474:7474 -p 7687:7687 -v graphpop_data:/data ghcr.io/jfmao/graphpop:latest
```

Or build from source:

```bash
git clone https://github.com/jfmao/GraphPop.git
cd GraphPop
docker compose up --build
```

## Method 4: From source (development)

```bash
git clone https://github.com/jfmao/GraphPop.git
cd GraphPop

# Create conda environment
conda create -n graphpop -c conda-forge \
    python=3.12 numpy click openjdk=21 maven -y
conda activate graphpop

# Build Java procedures
cd graphpop-procedures && ./mvnw package -DskipTests && cd ..

# Install Python packages in development mode
cd graphpop-cli && pip install -e ".[dev]" && cd ..
cd graphpop-import && pip install -e . && cd ..

# Run tests
cd graphpop-cli && pytest -v && cd ..

# Set up Neo4j with locally built plugin
graphpop setup --password mypassword \
    --deploy-plugin graphpop-procedures/target/graphpop-procedures-0.1.0.jar
graphpop start
```

## Method 5: HPC cluster (no Docker)

```bash
# On a login node:
module load java/21  # or use conda: conda install -c conda-forge openjdk=21

graphpop setup --password mypassword
graphpop start
```

See [Cluster Deployment Guide](cluster.md) and the
[SLURM/PBS job scripts](../graphpop-cli/docs/commands/cluster/) for
interactive, batch, and array-job recipes.

## What gets installed

| Component | Size | Source |
|-----------|------|--------|
| GraphPop CLI (Python) | ~2 MB | pip / conda |
| GraphPop procedures (Java JAR) | ~31 KB | GitHub Releases or local build |
| Neo4j Community 5.x | ~130 MB | Downloaded by `graphpop setup` |
| JDK 21 (if via conda) | ~200 MB | conda-forge |
| GraphPop Import Pipeline (Python) | ~1 MB | pip (optional) |
| GraphPop MCP Server (Python) | ~1 MB | pip (optional) |

Total: ~350 MB for a complete installation (CLI + Neo4j + Java).

## No admin privileges needed

GraphPop installs entirely in user space:
- conda/miniforge installs to `~/miniforge3/`
- Neo4j installs to `~/neo4j/` (or any user-writable directory)
- Java (if via conda) installs inside the conda environment
- The GraphPop procedures JAR is automatically deployed to the Neo4j plugins directory
- Configuration is stored in `~/.graphpop/config.yaml`

No `sudo`, no system packages, no Docker daemon required.

## Using an existing Neo4j installation

If Neo4j 2025.x is already installed on your machine, point GraphPop at it:

```bash
graphpop setup --skip-download --neo4j-home /path/to/neo4j --password mypassword
```

This skips the Neo4j download and only deploys the GraphPop procedures JAR,
sets the password, and writes `~/.graphpop/config.yaml`. Subsequent commands
(`start`, `stop`, `status`, etc.) read the config file automatically — you no
longer need `--neo4j-home` on every invocation.

If the existing Neo4j is **running**, use `--adopt` to deploy the plugin and
restart it:

```bash
graphpop setup --adopt --neo4j-home /path/to/neo4j --password mypassword
```

The `--adopt` path stops and restarts the instance (an interactive confirmation
is required unless `--yes` is passed). It is intended for user-owned,
single-user installs only — not shared system Neo4j at `/opt/neo4j` or
container-managed instances. For those, install a separate GraphPop-managed
Neo4j on different ports (see below).

### Skipping plugin deployment

If you only want to configure Neo4j without deploying the GraphPop procedures
plugin (e.g., you already have it deployed):

```bash
graphpop setup --skip-download --skip-plugin --neo4j-home /path/to/neo4j --password mypassword
```

## Offline / air-gapped install

If the target machine cannot reach the internet, download the Neo4j tarball
on a machine that can and transfer it:

1. Download from `https://dist.neo4j.org/neo4j-community-2025.12.1-unix.tar.gz`.
2. Optionally download the GraphPop procedures JAR from the
   [GitHub Releases](https://github.com/jfmao/GraphPop/releases) page.
3. Transfer to the target machine.
4. Run:

```bash
graphpop setup \
    --neo4j-tarball /path/to/neo4j-community-2025.12.1-unix.tar.gz \
    --password mypassword
```

To also deploy a local copy of the procedures JAR (instead of downloading
from GitHub Releases):

```bash
graphpop setup \
    --neo4j-tarball /path/to/neo4j-community-2025.12.1-unix.tar.gz \
    --deploy-plugin /path/to/graphpop-procedures.jar \
    --password mypassword
```

## Port conflicts and multiple Neo4j instances

If another Neo4j is already listening on the default port 7687, `graphpop setup`
will detect the conflict and print instructions. To run GraphPop's Neo4j
alongside another instance, specify alternative ports:

```bash
graphpop setup --bolt-port 7688 --http-port 7475 --password mypassword
```

The chosen ports are stored in `~/.graphpop/config.yaml` and used by all
subsequent commands automatically.

To check the full installation health at any time:

```bash
graphpop doctor
```

This verifies Java, Neo4j home, running process, port reachability, plugin
deployment, config file, and password connectivity.

## Verifying the installation

```bash
graphpop --version                # Shows CLI version
graphpop status                   # Shows Neo4j status, version, config, plugin
graphpop doctor                   # Full health check
graphpop --help                   # Full command list
```

## First import

```bash
# Start Neo4j
graphpop start

# Import a VCF
graphpop import \
    --vcf your_data.vcf.gz \
    --panel your_panel.txt \
    --database mydb

# Check status
graphpop status

# Run an analysis
graphpop diversity chr1 1 50000000 ALL
```

## Install the MCP server (optional, for AI agent access)

```bash
pip install -e graphpop-mcp

# Configure for Claude Desktop or other MCP clients
export GRAPHPOP_URI=bolt://localhost:7687
export GRAPHPOP_PASSWORD=mypassword
graphpop-mcp
```

## Troubleshooting

**Java not found**: Install via conda (no admin needed):
`conda install -c conda-forge openjdk=21`. On HPC clusters: `module load java/21`.

**Neo4j extremely slow**: Ensure the data directory is on local SSD, not NFS.

**Port 7687 in use**: See "Port conflicts" section above, or run
`graphpop doctor` for a full diagnostic.

**Plugin not loaded after setup**: Restart Neo4j with `graphpop stop && graphpop start`.
Verify with `graphpop doctor`.

**Permission denied on Neo4j binary**: Run `chmod +x ~/neo4j/bin/neo4j` and
`chmod +x ~/neo4j/bin/neo4j-admin`.

**graphpop setup fails mid-download**: The partial download is cached at
`/tmp/neo4j-community-2025.12.1-unix.tar.gz`. Delete it and retry, or use
`--neo4j-tarball` with a manually downloaded copy.
