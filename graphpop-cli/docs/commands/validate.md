# graphpop validate

## Title

Validate Graph Database Integrity and Completeness

## Description

Runs a series of integrity checks on the current GraphPop database to verify
that required node labels, indexes, stored procedures, and variant properties
are present. Reports issues with `[OK]`, `[FAIL]`, and `[--]` indicators.
Optionally attempts automatic fixes for missing indexes.

## Usage

```
graphpop validate [OPTIONS]
```

## Arguments

This command takes no positional arguments.

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `--fix` | FLAG | `false` | Attempt to automatically fix issues (currently: create missing indexes). |

## Value

Prints a structured validation report to stdout. Exits with code 0 regardless
of findings (check the output for `[FAIL]` markers). Example:

```
Validating database: rice3k

Node labels:
  [OK]   Variant: 2,103,456
  [OK]   Sample: 3,000
  [OK]   Population: 12
  [OK]   Chromosome: 12
  [OK]   Gene: 28,236
  [OK]   Pathway: 312
  [--]   GOTerm: not present (optional)
  [OK]   GenomicWindow: 4,800

Indexes:
  [OK]   Variant(chr, pos)
  [OK]   Population(populationId)
  [OK]   Sample(sampleId)

GraphPop procedures:
  [OK]   graphpop.diversity
  [OK]   graphpop.divergence
  ...

Variant node properties:
  [OK]   chr
  [OK]   pos
  [OK]   ref
  [OK]   alt
  [OK]   pop_ids
  [OK]   ac
  [OK]   an
  [OK]   af
  [OK]   gt_packed
  [--]   phase_packed: not present (optional)
  [OK]   ancestral_allele
  [OK]   is_polarized

========================================
VALIDATION: All checks passed
```

## Details

### Checks Performed

**1. Node Labels**

Verifies that the following required labels exist and have at least one node:

| Label | Status | Description |
|---|---|---|
| `Variant` | Required | Variant nodes with allele count arrays |
| `Sample` | Required | Sample nodes linked to populations |
| `Population` | Required | Population nodes with sample counts |
| `Chromosome` | Required | Chromosome nodes with length |
| `Gene` | Optional | Gene nodes from VEP/SnpEff annotations |
| `Pathway` | Optional | Pathway nodes from Reactome |
| `GOTerm` | Optional | GO term nodes from GOA |
| `GenomicWindow` | Optional | Window nodes from genome-scan |

Missing required labels are reported as `[FAIL]`. Missing optional labels are
reported as `[--]` (informational, not an error).

**2. Indexes**

Checks for the following indexes required for query performance:

| Index | Properties | Purpose |
|---|---|---|
| `Variant(chr, pos)` | Composite | Fast variant lookup by genomic position |
| `Population(populationId)` | Single | Population lookup in procedure calls |
| `Sample(sampleId)` | Single | Sample lookup for genotype queries |

Missing indexes are reported as `[FAIL]`. With `--fix`, the command attempts to
create them using `CREATE INDEX ... IF NOT EXISTS`.

**3. GraphPop Procedures**

Checks that all 12 expected stored procedures are installed:

`graphpop.diversity`, `graphpop.divergence`, `graphpop.sfs`,
`graphpop.joint_sfs`, `graphpop.genome_scan`, `graphpop.pop_summary`,
`graphpop.ld`, `graphpop.ihs`, `graphpop.xpehh`, `graphpop.nsl`,
`graphpop.roh`, `graphpop.garud_h`.

Missing procedures indicate that the `graphpop-procedures.jar` plugin is not
deployed or is an older version.

**4. Variant Node Properties**

Samples one Variant node and checks for required and optional properties:

| Property | Status | Description |
|---|---|---|
| `chr` | Required | Chromosome identifier |
| `pos` | Required | Genomic position |
| `ref` | Required | Reference allele |
| `alt` | Required | Alternate allele |
| `pop_ids` | Required | Array of population IDs |
| `ac` | Required | Allele count array |
| `an` | Required | Allele number array |
| `af` | Required | Allele frequency array |
| `gt_packed` | Optional | Bit-packed genotype array (full-path stats) |
| `phase_packed` | Optional | Bit-packed phase array |
| `ancestral_allele` | Optional | Ancestral allele (for unfolded SFS) |
| `is_polarized` | Optional | Whether ancestral state is known |

### The --fix Flag

Currently, `--fix` only creates missing indexes. It does not create missing
nodes, install procedures, or add missing properties. For those issues:

- Missing nodes: Re-run `graphpop import`.
- Missing procedures: Re-deploy with `graphpop setup --deploy-plugin`.
- Missing properties: Re-run the annotation loading step of import.

## Examples

### Basic validation

```bash
graphpop validate
```

### Validate and fix missing indexes

```bash
graphpop validate --fix
```

### Validate a specific database

```bash
graphpop --database human_chr22 validate
```

### Post-import validation

```bash
graphpop import --vcf data.vcf.gz --panel panel.txt --database mydb
graphpop start
graphpop validate
# Review output, then:
graphpop validate --fix  # if indexes are missing
```

## See Also

- `graphpop db info` -- Show node/edge counts (less detailed than validate).
- `graphpop import` -- Re-import if required nodes are missing.
- `graphpop setup` -- Re-deploy the plugin if procedures are missing.
- `graphpop status` -- Check server and plugin status.
