# graphpop dump

## Title

Dump a Neo4j Database to a File for Sharing

## Description

Creates a `neo4j-admin dump` file from a named database, suitable for sharing,
backup, or transferring to another machine. Optionally generates a JSON manifest
file containing database metadata (node/edge counts, populations, chromosomes).

## Usage

```
graphpop dump [OPTIONS]
```

## Arguments

This command takes no positional arguments.

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `--database` | TEXT | *(required)* | Name of the database to dump. |
| `-o`, `--output` | PATH | `graphpop_<database>_<YYYYMMDD>.dump` | Output dump file path. |
| `--neo4j-home` | PATH | *(from config)* | Neo4j installation directory. |
| `--manifest` / `--no-manifest` | FLAG | `--manifest` | Generate a JSON manifest file alongside the dump. |

## Value

No return value. Creates two files:

- `<output>.dump` -- The Neo4j database dump (binary format).
- `<output>.manifest.json` -- Metadata including node counts, edge counts,
  populations, chromosomes, and dump size (only if `--manifest` is enabled).

### Manifest Format

```json
{
  "database": "rice3k",
  "date": "2026-03-27T14:30:00",
  "dump_file": "graphpop_rice3k_20260327.dump",
  "dump_size_bytes": 2147483648,
  "node_counts": {
    "Variant": 2103456,
    "Sample": 3000,
    "Population": 12,
    "Chromosome": 12,
    "Gene": 28236
  },
  "edge_counts": {
    "CARRIES": 18456789,
    "NEXT": 2103455
  },
  "populations": {
    "GJ-tmp": 840,
    "GJ-trp": 498,
    "XI-1A": 352
  },
  "chromosomes": {
    "chr1": 43270923,
    "chr2": 35937250
  }
}
```

If Neo4j is not running at dump time, the manifest will contain only basic
file metadata and a note that detailed counts are unavailable.

## Details

Runs `neo4j-admin database dump --to-path=<dir> <database>`.

**Prerequisites:**

- **Neo4j must be stopped** before dumping. The `neo4j-admin dump` command
  requires exclusive access to the database files. Run `graphpop stop` first.
- The `neo4j-admin` binary must exist at `<neo4j-home>/bin/neo4j-admin`.

The manifest generation attempts to query the database for metadata. If Neo4j
is not running (the normal case during dump), the manifest will contain only
the database name, date, and file size. To get full metadata in the manifest,
start Neo4j temporarily, run the dump with manifest, then stop again -- or
generate the manifest separately before stopping.

Dump files are portable across machines running the same major version of Neo4j.

## Examples

### Basic dump with default naming

```bash
graphpop stop
graphpop dump --database rice3k
# Creates: graphpop_rice3k_20260327.dump
#          graphpop_rice3k_20260327.manifest.json
graphpop start
```

### Dump with custom output path

```bash
graphpop stop
graphpop dump --database human_chr22 -o backups/human_chr22_v1.dump
graphpop start
```

### Dump without manifest

```bash
graphpop stop
graphpop dump --database rice3k --no-manifest
graphpop start
```

### Full backup workflow

```bash
graphpop stop
graphpop dump --database rice3k -o /backup/rice3k_$(date +%Y%m%d).dump
graphpop dump --database human_chr22 -o /backup/human_chr22_$(date +%Y%m%d).dump
graphpop start
```

## See Also

- `graphpop load` -- Restore a database from a dump file.
- `graphpop stop` -- Stop Neo4j before dumping (required).
- `graphpop start` -- Restart Neo4j after dumping.
- `graphpop db list` -- List databases to choose which to dump.
