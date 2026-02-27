"""Parse VCF files and extract variant/genotype data for Neo4j import."""


class VCFParser:
    """Read a VCF file and yield per-variant records with population allele counts.

    TODO: Implement using cyvcf2.
    """
