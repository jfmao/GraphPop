# GraphPop

GraphPop is a graph-native population genomics platform built on Neo4j. It models
variants, samples, populations, and genomic annotations as a property graph and
provides SIMD-accelerated stored procedures for computing population genetics
statistics (allele frequencies, F-statistics, LD, selection scans) directly inside
the database. See [`docs/GraphPop_Compiled_Holistic_Design.md`](docs/GraphPop_Compiled_Holistic_Design.md)
for the full design document.
