package org.graphpop.procedures;

import org.graphpop.VectorOps;
import org.neo4j.graphdb.RelationshipType;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Computes pairwise linkage disequilibrium (LD) within a genomic region
 * and writes LD relationship edges between Variant nodes.
 *
 * <p>Architecture: load → parallel compute → sequential write.
 * <ol>
 *   <li>Load dense haplotype matrix via {@link HaplotypeMatrix} (single-threaded)</li>
 *   <li>Compute pairwise r² and D' in parallel (ForkJoinPool, no DB access)</li>
 *   <li>Write LD edges sequentially</li>
 * </ol>
 * </p>
 *
 * <p>Optimizations (M1.4):
 * <ul>
 *   <li>D' computed directly from byte[] haplotypes — zero array allocation per pair</li>
 *   <li>Thread-local List collection via flatMap — no CAS contention</li>
 * </ul>
 * </p>
 */
public class LDProcedure {

    private static final RelationshipType LD_REL = RelationshipType.withName("LD");

    @Context
    public Transaction tx;

    public static class LDResult {
        public String variant1;
        public String variant2;
        public double r2;
        public double dprime;
        public long distance;

        public LDResult() {}
    }

    @Procedure(name = "graphpop.ld", mode = Mode.WRITE)
    @SuppressWarnings("unchecked")
    public Stream<LDResult> ld(
            @Name("chr") String chr,
            @Name("start") long start,
            @Name("end") long end,
            @Name("pop") String pop,
            @Name(value = "max_dist", defaultValue = "500000") long maxDist,
            @Name(value = "r2_threshold", defaultValue = "0.2") double r2Threshold,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        VariantFilter filter = VariantFilter.fromOptions(options);
        String variantType = filter.variantType;
        boolean writeEdges = true;
        if (options != null && options.containsKey("write_edges")) {
            writeEdges = Boolean.TRUE.equals(options.get("write_edges"));
        }

        List<String> sampleList = options != null ? (List<String>) options.get("samples") : null;
        boolean useSubset = sampleList != null && !sampleList.isEmpty();

        // Phase 1: Load dense haplotype matrix (single-threaded)
        HaplotypeMatrix matrix;
        if (useSubset) {
            Map<String, Integer> sampleIndex = GenotypeLoader.buildSampleIndex(tx, sampleList);
            matrix = HaplotypeMatrix.load(tx, chr, sampleIndex, start, end, variantType);
        } else {
            matrix = HaplotypeMatrix.load(tx, chr, pop, start, end, variantType);
        }
        if (matrix == null || matrix.nVariants < 2) return Stream.empty();

        // Identify polymorphic variants (skip monomorphic, apply filters)
        int[] polyIdx = identifyPolymorphic(matrix, filter);
        if (polyIdx.length < 2) return Stream.empty();

        // Pre-compute dosage vectors for all polymorphic variants
        double[][] dosages = new double[polyIdx.length][];
        for (int i = 0; i < polyIdx.length; i++) {
            dosages[i] = matrix.dosage(polyIdx[i]);
        }

        int nHaps = matrix.nHaplotypes;

        // Phase 2: Parallel pairwise r² + D' computation (no DB access)
        List<LDHit> hits = IntStream.range(0, polyIdx.length).parallel()
                .mapToObj(i -> {
                    int vi = polyIdx[i];
                    long posI = matrix.positions[vi];
                    byte[] hapI = matrix.haplotypes[vi];
                    List<LDHit> local = new ArrayList<>();

                    for (int j = i + 1; j < polyIdx.length; j++) {
                        int vj = polyIdx[j];
                        long posJ = matrix.positions[vj];
                        long dist = posJ - posI;
                        if (dist > maxDist) break;

                        double r2 = VectorOps.pearsonR2(dosages[i], dosages[j]);
                        if (r2 >= r2Threshold) {
                            double dprime = VectorOps.dPrimePacked(hapI, matrix.haplotypes[vj], nHaps);
                            local.add(new LDHit(i, j, vi, vj, r2, dprime, dist));
                        }
                    }
                    return local;
                })
                .flatMap(List::stream)
                .collect(Collectors.toList());

        // Phase 3: Sequential write of LD edges (skipped if write_edges=false)
        List<LDResult> results = new ArrayList<>();
        for (LDHit hit : hits) {
            if (writeEdges) {
                var rel = matrix.nodes[hit.vi].createRelationshipTo(matrix.nodes[hit.vj], LD_REL);
                rel.setProperty("r2", hit.r2);
                rel.setProperty("dprime", hit.dprime);
                rel.setProperty("population", pop);
            }

            LDResult lr = new LDResult();
            lr.variant1 = matrix.variantIds[hit.vi];
            lr.variant2 = matrix.variantIds[hit.vj];
            lr.r2 = hit.r2;
            lr.dprime = hit.dprime;
            lr.distance = hit.dist;
            results.add(lr);
        }

        return results.stream();
    }

    private static int[] identifyPolymorphic(HaplotypeMatrix matrix, VariantFilter filter) {
        List<Integer> poly = new ArrayList<>();
        for (int i = 0; i < matrix.nVariants; i++) {
            int ac = matrix.acs[i];
            int an = matrix.ans[i];
            if (ac <= 0 || ac >= an || an < 2) continue;
            double af = (double) ac / an;

            if (filter.isActive()) {
                if (af < filter.minAf || af > filter.maxAf) continue;
                // variant_type and HWE filters would require Node access;
                // for FULL PATH, we skip those (nodes may be null in test matrices)
            } else {
                // No filter active — backward compatible
            }

            poly.add(i);
        }
        return poly.stream().mapToInt(Integer::intValue).toArray();
    }

    private static class LDHit {
        final int idxI, idxJ;
        final int vi, vj;
        final double r2, dprime;
        final long dist;

        LDHit(int idxI, int idxJ, int vi, int vj, double r2, double dprime, long dist) {
            this.idxI = idxI;
            this.idxJ = idxJ;
            this.vi = vi;
            this.vj = vj;
            this.r2 = r2;
            this.dprime = dprime;
            this.dist = dist;
        }
    }
}
