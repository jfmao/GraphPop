#!/usr/bin/env python3
"""Build ancestral allele FASTA files for rice from Ensembl Plants 8-rice EPO.

Downloads EMF files from Ensembl Plants FTP, parses the Ortheus-reconstructed
ancestral sequences, and writes per-chromosome FASTA files compatible with
annotate_ancestral.py.

Usage:
    python scripts/build_rice_ancestral.py download   # Download EMF files
    python scripts/build_rice_ancestral.py build       # Parse EMF → FASTA
    python scripts/build_rice_ancestral.py annotate    # Annotate Neo4j Variant nodes
    python scripts/build_rice_ancestral.py all         # Run all steps
"""

import argparse
import gzip
import logging
import os
import re
import subprocess
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-5s %(message)s",
)
log = logging.getLogger(__name__)

# ── Constants ─────────────────────────────────────────────────────────

EMF_BASE_URL = (
    "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-62/"
    "emf/ensembl-compara/multiple_alignments/8_rice.epo"
)

EMF_DIR = Path("data/raw/rice_epo/emf")
FASTA_DIR = Path("data/raw/rice_epo")

# Number of EMF segments per chromosome (from Ensembl Plants release 62)
CHR_SEGMENTS = {
    1: 9, 2: 8, 3: 7, 4: 8, 5: 6, 6: 7,
    7: 7, 8: 7, 9: 5, 10: 5, 11: 7, 12: 7,
}

# Chromosome lengths (IRGSP-1.0)
CHR_LENGTHS = {
    1: 43270923, 2: 35937250, 3: 36413819, 4: 35502694,
    5: 29958434, 6: 31248787, 7: 29697621, 8: 28443022,
    9: 23012720, 10: 23207287, 11: 29021106, 12: 27531856,
}

COMP = str.maketrans("ACGTacgtNn.-", "TGCAtgcaNn.-")

SATIVA_SPECIES = "oryza_sativa"


# ── Download ──────────────────────────────────────────────────────────


def cmd_download(args):
    """Download 8-rice EPO EMF files from Ensembl Plants FTP."""
    EMF_DIR.mkdir(parents=True, exist_ok=True)

    files_to_download = []
    for chr_num, n_segs in CHR_SEGMENTS.items():
        for seg in range(1, n_segs + 1):
            fname = f"8_rice.epo.{chr_num}_{seg}.emf.gz"
            local_path = EMF_DIR / fname
            if local_path.exists() and local_path.stat().st_size > 0:
                continue
            files_to_download.append((fname, local_path))

    if not files_to_download:
        log.info("All EMF files already downloaded (%d files)", sum(CHR_SEGMENTS.values()))
        return

    log.info("Downloading %d EMF files to %s", len(files_to_download), EMF_DIR)
    for i, (fname, local_path) in enumerate(files_to_download, 1):
        url = f"{EMF_BASE_URL}/{fname}"
        log.info("  [%d/%d] %s", i, len(files_to_download), fname)
        try:
            subprocess.run(
                ["wget", "-q", url, "-O", str(local_path)],
                check=True, timeout=120,
            )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            log.error("Failed to download %s: %s", fname, e)
            if local_path.exists():
                local_path.unlink()

    n_ok = sum(1 for c, n in CHR_SEGMENTS.items()
               for s in range(1, n + 1)
               if (EMF_DIR / f"8_rice.epo.{c}_{s}.emf.gz").exists())
    log.info("Download complete: %d/%d files", n_ok, sum(CHR_SEGMENTS.values()))


# ── EMF Parser ────────────────────────────────────────────────────────


def _find_root_ancestor_id(tree_str):
    """Extract the root ancestor ID from an EMF TREE Newick string.

    The root label is the text after the final ')' before the ';'.
    Example: ...))Aseq_Ancestor_9910_30249_1_55285[+]:0.006;
    → 'Ancestor_9910_30249'
    """
    tree_str = tree_str.strip().rstrip(";")
    # Root label is after the last ')'
    idx = tree_str.rfind(")")
    if idx < 0:
        return None
    root_label = tree_str[idx + 1:]
    # Extract ancestor ID: Aseq_Ancestor_XXXX_YYYYY_start_end[strand]:brlen
    m = re.match(r"Aseq_(Ancestor_\d+_\d+)_", root_label)
    if m:
        return m.group(1)
    return None


def _parse_emf_blocks(emf_gz_path):
    """Parse an EMF file and yield alignment blocks.

    Each block yields:
        (seq_headers, tree_str, data_columns)

    seq_headers: list of dicts with keys:
        species, region, start, end, strand, chr_length, ancestor_id (if ancestral)
    data_columns: list of strings, each string = one alignment column
    """
    with gzip.open(emf_gz_path, "rt") as f:
        seq_headers = []
        tree_str = None
        data_lines = []
        in_data = False

        for line in f:
            line = line.rstrip("\n")

            if line.startswith("SEQ "):
                in_data = False
                parts = line.split()
                # SEQ species region start end strand (chr_length=N)
                species = parts[1]
                region = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                strand = int(parts[5])
                chr_len_match = re.search(r"chr_length=(\d+)", line)
                chr_length = int(chr_len_match.group(1)) if chr_len_match else 0

                header = {
                    "species": species,
                    "region": region,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "chr_length": chr_length,
                }
                # Tag ancestral sequences with their ID
                if species == "ancestral_sequences":
                    header["ancestor_id"] = region  # e.g. "Ancestor_9910_30249"
                seq_headers.append(header)

            elif line.startswith("TREE "):
                tree_str = line[5:]

            elif line == "DATA":
                in_data = True
                data_lines = []

            elif line == "//":
                if seq_headers and data_lines:
                    yield seq_headers, tree_str, data_lines
                seq_headers = []
                tree_str = None
                data_lines = []
                in_data = False

            elif in_data and not line.startswith("#") and not line.startswith("ID "):
                data_lines.append(line)


def _process_block(seq_headers, tree_str, data_lines, ancestral_dict):
    """Process one alignment block and update ancestral_dict.

    ancestral_dict: {1-based genomic position: ancestral_char}
    """
    # Find oryza_sativa index
    sativa_idx = None
    sativa_header = None
    for i, h in enumerate(seq_headers):
        if h["species"] == SATIVA_SPECIES:
            sativa_idx = i
            sativa_header = h
            break

    if sativa_idx is None:
        return 0  # No sativa in this block

    # Find root ancestral index
    root_ancestor_id = _find_root_ancestor_id(tree_str) if tree_str else None
    root_idx = None
    if root_ancestor_id:
        for i, h in enumerate(seq_headers):
            if h.get("ancestor_id") == root_ancestor_id:
                root_idx = i
                break

    if root_idx is None:
        # Fallback: use the immediate parent of sativa (index sativa_idx + 1,
        # which in EMF is typically the first ancestral node after sativa)
        for i, h in enumerate(seq_headers):
            if h["species"] == "ancestral_sequences" and i == sativa_idx + 1:
                root_idx = i
                break

    if root_idx is None:
        return 0

    # Walk through alignment columns
    strand = sativa_header["strand"]
    if strand == 1 or strand == -1:
        pos = sativa_header["start"] if strand == 1 else sativa_header["end"]
    else:
        return 0

    n_added = 0
    for col_str in data_lines:
        if len(col_str) <= max(sativa_idx, root_idx):
            continue

        sativa_char = col_str[sativa_idx]
        anc_char = col_str[root_idx]

        # Skip gaps in sativa (no reference position)
        if sativa_char == "-":
            continue

        # Get the ancestral base
        if anc_char not in ("-", ".", "N", "n"):
            # For negative strand, complement back to forward strand
            if strand == -1:
                anc_forward = anc_char.translate(COMP)
            else:
                anc_forward = anc_char

            # Only store if position not already filled (first block wins)
            if pos not in ancestral_dict:
                ancestral_dict[pos] = anc_forward
                n_added += 1

        # Advance position
        if strand == 1:
            pos += 1
        else:
            pos -= 1

    return n_added


def build_chromosome_fasta(chr_num):
    """Parse all EMF segments for one chromosome and build ancestral FASTA.

    Returns the ancestral sequence string (1-based positions mapped to 0-based
    indices, unknown positions filled with 'N').
    """
    chr_len = CHR_LENGTHS[chr_num]
    ancestral_dict = {}  # {1-based pos: char}

    n_segs = CHR_SEGMENTS[chr_num]
    for seg in range(1, n_segs + 1):
        emf_path = EMF_DIR / f"8_rice.epo.{chr_num}_{seg}.emf.gz"
        if not emf_path.exists():
            log.warning("Missing EMF file: %s", emf_path)
            continue

        n_blocks = 0
        n_bases = 0
        for seq_headers, tree_str, data_lines in _parse_emf_blocks(emf_path):
            added = _process_block(seq_headers, tree_str, data_lines, ancestral_dict)
            n_bases += added
            n_blocks += 1

        log.info("  Segment %d_%d: %d blocks, %d ancestral bases",
                 chr_num, seg, n_blocks, n_bases)

    # Build FASTA string (0-based index = 1-based position - 1)
    seq_chars = ["N"] * chr_len
    for pos, char in ancestral_dict.items():
        if 1 <= pos <= chr_len:
            seq_chars[pos - 1] = char

    coverage = sum(1 for c in seq_chars if c != "N")
    log.info("  Chr%d: %d/%d positions covered (%.1f%%)",
             chr_num, coverage, chr_len, 100 * coverage / chr_len)

    return "".join(seq_chars)


def write_fasta(chr_num, seq, out_dir):
    """Write ancestral sequence as FASTA file."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"oryza_sativa_ancestor_Chr{chr_num}.fa"

    with open(out_path, "w") as f:
        f.write(f">Chr{chr_num} Ancestral sequence from 8-rice EPO (Ensembl Plants r62)\n")
        # Write sequence in 80-char lines
        for i in range(0, len(seq), 80):
            f.write(seq[i:i + 80] + "\n")

    log.info("  Written: %s (%d bp)", out_path, len(seq))
    return out_path


# ── Build ─────────────────────────────────────────────────────────────


def cmd_build(args):
    """Parse EMF files and build per-chromosome ancestral FASTA files."""
    chromosomes = list(range(1, 13))
    if args.chr:
        chromosomes = [int(c.replace("Chr", "").replace("chr", "")) for c in args.chr]

    log.info("Building ancestral FASTA for %d chromosomes", len(chromosomes))

    for chr_num in chromosomes:
        log.info("Processing Chr%d...", chr_num)
        seq = build_chromosome_fasta(chr_num)
        write_fasta(chr_num, seq, FASTA_DIR)

    log.info("All FASTA files written to %s", FASTA_DIR)


# ── Annotate ──────────────────────────────────────────────────────────


def cmd_annotate(args):
    """Annotate Neo4j Variant nodes using the built ancestral FASTA files."""
    chromosomes = list(range(1, 13))
    if args.chr:
        chromosomes = [int(c.replace("Chr", "").replace("chr", "")) for c in args.chr]

    for chr_num in chromosomes:
        fasta_path = FASTA_DIR / f"oryza_sativa_ancestor_Chr{chr_num}.fa"
        if not fasta_path.exists():
            log.error("FASTA not found: %s (run 'build' first)", fasta_path)
            continue

        chr_name = f"Chr{chr_num}"
        log.info("Annotating %s from %s...", chr_name, fasta_path)

        cmd = [
            "python", "scripts/annotate_ancestral.py",
            "--chr", chr_name,
            "--fasta", str(fasta_path),
            "--uri", args.uri,
            "--user", args.user,
            "--password", args.password,
            "--batch-size", str(args.batch_size),
        ]
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            log.error("Annotation failed for %s: %s", chr_name, e)

    log.info("Annotation complete for all chromosomes")


# ── All ───────────────────────────────────────────────────────────────


def cmd_all(args):
    """Run download + build + annotate."""
    cmd_download(args)
    cmd_build(args)
    cmd_annotate(args)


# ── Main ──────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Build rice ancestral allele FASTA from Ensembl 8-rice EPO")
    sub = parser.add_subparsers(dest="command")

    sub.add_parser("download", help="Download EMF files from Ensembl Plants FTP")

    p_build = sub.add_parser("build", help="Parse EMF → FASTA")
    p_build.add_argument("--chr", nargs="*", help="Specific chromosomes (e.g. Chr1 Chr2)")

    p_ann = sub.add_parser("annotate", help="Annotate Neo4j Variant nodes")
    p_ann.add_argument("--chr", nargs="*", help="Specific chromosomes")
    p_ann.add_argument("--uri", default="bolt://localhost:7687")
    p_ann.add_argument("--user", default="neo4j")
    p_ann.add_argument("--password", default="graphpop")
    p_ann.add_argument("--batch-size", type=int, default=10000)

    p_all = sub.add_parser("all", help="Download + build + annotate")
    p_all.add_argument("--chr", nargs="*", help="Specific chromosomes")
    p_all.add_argument("--uri", default="bolt://localhost:7687")
    p_all.add_argument("--user", default="neo4j")
    p_all.add_argument("--password", default="graphpop")
    p_all.add_argument("--batch-size", type=int, default=10000)

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        sys.exit(1)

    {"download": cmd_download, "build": cmd_build,
     "annotate": cmd_annotate, "all": cmd_all}[args.command](args)


if __name__ == "__main__":
    main()
