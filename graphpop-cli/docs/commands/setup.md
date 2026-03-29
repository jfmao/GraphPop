# graphpop setup

## Title

Set Up Neo4j for GraphPop

## Description

Downloads Neo4j Community Edition, configures memory settings, sets the initial
database password, optionally deploys the GraphPop stored procedures plugin, and
creates the GraphPop configuration file at `~/.graphpop/config.yaml`.

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
| `--deploy-plugin` | PATH | *(none)* | Path to a `graphpop-procedures.jar` file to copy into the Neo4j `plugins/` directory. |

## Value

No return value. Side effects:

- Neo4j Community Edition 5.26.0 downloaded and extracted to `--neo4j-home`.
- `neo4j.conf` updated with the specified memory settings and GraphPop-specific
  configuration (`dbms.security.procedures.unrestricted=graphpop.*`).
- Initial password set via `neo4j-admin dbms set-initial-password`.
- GraphPop procedures JAR deployed (if `--deploy-plugin` given).
- Configuration file written to `~/.graphpop/config.yaml` with URI, user,
  password, database, and neo4j_home.

## Details

The setup process runs six steps in order:

1. **Download Neo4j** -- Fetches the `neo4j-community-5.26.0-unix.tar.gz`
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
   (`bolt://localhost:7687`), user (`neo4j`), the supplied password, database
   (`neo4j`), and the neo4j_home path.

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

### Re-deploy plugin without re-downloading Neo4j

```bash
cd graphpop-procedures && mvn package -DskipTests
graphpop setup \
    --skip-download \
    --deploy-plugin target/graphpop-procedures-1.0-SNAPSHOT.jar \
    --password mySecurePass123
```

### Non-interactive password (CI/scripting)

```bash
graphpop setup --password "$NEO4J_PASSWORD" --skip-download
```

## See Also

- `graphpop start` -- Start the Neo4j server after setup.
- `graphpop status` -- Verify that Neo4j is running and the plugin is deployed.
- `graphpop config init` -- Create or re-create the config file interactively.
- `graphpop import` -- Import VCF data into the graph database.
