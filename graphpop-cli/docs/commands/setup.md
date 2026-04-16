# graphpop setup

## Title

Set Up Neo4j for GraphPop

## Description

Downloads Neo4j Community Edition, automatically downloads and deploys the
pre-compiled GraphPop stored procedures plugin, configures memory settings,
sets the initial database password, and creates the GraphPop configuration
file at `~/.graphpop/config.yaml`.

No Java or Maven installation is required — the procedures plugin JAR is
downloaded as a pre-compiled binary from GitHub Releases. If installed via
conda, the bundled JAR is used automatically.

This is the first command you should run on a new machine. It is safe to re-run;
if Neo4j is already installed it will prompt before overwriting.

## Usage

```
graphpop setup [OPTIONS]
```

## Arguments

This command takes no positional arguments.

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `--neo4j-home` | PATH | `~/neo4j` | Neo4j installation directory. |
| `--pagecache` | TEXT | `16g` | Neo4j page cache size (e.g. `16g`, `20g`). Controls how much of the graph is kept in memory. |
| `--heap` | TEXT | `4g` | JVM heap size for Neo4j (initial and max). |
| `--password` | TEXT | *(prompted)* | Password for the `neo4j` user. Prompted interactively with confirmation if not supplied on the command line. |
| `--skip-download` | FLAG | `false` | Skip downloading Neo4j. Use this when Neo4j is already installed at `--neo4j-home`. |
| `--deploy-plugin` | PATH | *(none)* | Path to a local `graphpop-procedures.jar` file. Overrides the automatic download from GitHub Releases. |
| `--skip-plugin` | FLAG | `false` | Skip deploying the GraphPop procedures plugin entirely. |
| `--bolt-port` | INT | `7687` | Bolt protocol port. Use a non-default port to run alongside another Neo4j instance. |
| `--http-port` | INT | `7474` | HTTP port for the Neo4j Browser interface. |
| `--neo4j-tarball` | PATH | *(none)* | Path to a local Neo4j tarball for offline/air-gapped installation. Skips the download step. |
| `--adopt` | FLAG | `false` | Adopt a running Neo4j instance: stops it, deploys the plugin, configures, and restarts. Requires confirmation (or `--yes`). |
| `--yes` | FLAG | `false` | Skip interactive confirmations. Intended for use with `--adopt` in scripts. |

## Value

No return value. Side effects:

- Neo4j Community Edition 2025.12.1 downloaded and extracted to `--neo4j-home`.
- `neo4j.conf` updated with the specified memory settings and GraphPop-specific
  configuration (`dbms.security.procedures.unrestricted=graphpop.*`).
- Initial password set via `neo4j-admin dbms set-initial-password`.
- GraphPop procedures JAR deployed (auto-downloaded from GitHub Releases,
  or from conda bundle, or from `--deploy-plugin` path).
- Configuration file written to `~/.graphpop/config.yaml` with URI, user,
  password, database, and neo4j_home.

## Details

The setup process runs six steps in order:

1. **Download Neo4j** -- Fetches the `neo4j-community-2025.12.1-unix.tar.gz`
   tarball from `dist.neo4j.org` into `/tmp` (cached across runs) and extracts
   it to `--neo4j-home`. Skipped if `--skip-download` is set or the user
   declines to overwrite an existing installation.

2. **Verify installation** -- Confirms that `<neo4j-home>/bin/neo4j` exists.

3. **Configure Neo4j** -- Writes the following settings to `neo4j.conf`:
   - `server.memory.pagecache.size` (from `--pagecache`)
   - `server.memory.heap.initial_size` and `server.memory.heap.max_size` (from `--heap`)
   - `server.directories.import=import`
   - `db.tx_log.rotation.retention_policy=2 days 2G`
   - `dbms.security.procedures.unrestricted=graphpop.*`

4. **Set password** -- Runs `neo4j-admin dbms set-initial-password`. If the
   password has already been set, a warning is printed.

5. **Deploy plugin** -- Copies the JAR to `<neo4j-home>/plugins/graphpop-procedures.jar`.

6. **Create config** -- Writes `~/.graphpop/config.yaml` with bolt URI
   (`bolt://localhost:<bolt-port>`), user (`neo4j`), the supplied password,
   database (`neo4j`), and the neo4j_home path. Custom ports are stored in the
   config and used by all subsequent commands automatically.

### Adopt mode (`--adopt`)

When `--adopt` is specified, setup targets a **running** Neo4j instance:

1. Stops Neo4j (`neo4j stop`)
2. Deploys the GraphPop procedures plugin
3. Applies configuration (memory, ports)
4. Restarts Neo4j (`neo4j start`)
5. Writes config file

An interactive confirmation is required before stopping Neo4j. Pass `--yes`
to skip it (for scripting). This mode is intended for user-owned, single-user
installs — not shared system Neo4j or container-managed instances.

### Offline install (`--neo4j-tarball`)

For air-gapped or HPC environments, download the Neo4j tarball on a connected
machine and transfer it. The `--neo4j-tarball` flag extracts the local tarball
instead of downloading. Combine with `--deploy-plugin` to also supply a local
copy of the procedures JAR.

### Port configuration (`--bolt-port`, `--http-port`)

Non-default ports are written to both `neo4j.conf` and `~/.graphpop/config.yaml`.
All subsequent CLI commands (`start`, `stop`, `status`, `doctor`, analysis
commands) read the config file automatically and use the correct ports.

Use custom ports when another Neo4j instance is already running on the default
ports. The setup command checks for port conflicts before proceeding.

### Memory Sizing Guidelines

| System RAM | Recommended `--pagecache` | Recommended `--heap` |
|---|---|---|
| 16 GB | `8g` | `2g` |
| 32 GB | `16g` | `4g` |
| 64 GB | `20g` | `4g` |
| 128 GB | `40g` | `8g` |

The page cache should be large enough to hold the graph store files. For a
1000 Genomes chr22 dataset (~1.1M variants, 2504 samples), the store is about
4 GB. For Rice 3K chr1 (~2M variants, 3000 samples), plan for 6-8 GB.

## Examples

### Basic setup with defaults

```bash
graphpop setup --password mySecurePass123
```

Sets up Neo4j at `~/neo4j` with 16 GB page cache and 4 GB heap.

### Custom installation directory and memory

```bash
graphpop setup \
    --neo4j-home /opt/neo4j \
    --pagecache 20g \
    --heap 8g \
    --password mySecurePass123
```

### Use an existing Neo4j installation

```bash
# Skip download, just deploy the plugin and write config:
graphpop setup --skip-download --password mySecurePass123

# Deploy a locally built JAR (for developers):
cd graphpop-procedures && ./mvnw package -DskipTests
graphpop setup \
    --skip-download \
    --deploy-plugin target/graphpop-procedures-0.1.0-SNAPSHOT.jar \
    --password mySecurePass123
```

### Adopt a running Neo4j instance

```bash
graphpop setup --adopt --neo4j-home /path/to/neo4j --password mySecurePass123

# Non-interactive (for scripts):
graphpop setup --adopt --yes --neo4j-home /path/to/neo4j --password mySecurePass123
```

### Run alongside another Neo4j instance

```bash
graphpop setup --bolt-port 7688 --http-port 7475 --password mySecurePass123
```

### Offline / air-gapped install

```bash
graphpop setup \
    --neo4j-tarball /path/to/neo4j-community-2025.12.1-unix.tar.gz \
    --deploy-plugin /path/to/graphpop-procedures.jar \
    --password mySecurePass123
```

### Non-interactive password (CI/scripting)

```bash
graphpop setup --password "$NEO4J_PASSWORD" --skip-download
```

## See Also

- [graphpop doctor](doctor.md) -- Verify the installation health after setup.
- [graphpop start](start.md) -- Start the Neo4j server after setup.
- [graphpop status](status.md) -- Verify that Neo4j is running and the plugin is deployed.
- [graphpop config init](config-init.md) -- Create or re-create the config file interactively.
- [graphpop import](import.md) -- Import VCF data into the graph database.
