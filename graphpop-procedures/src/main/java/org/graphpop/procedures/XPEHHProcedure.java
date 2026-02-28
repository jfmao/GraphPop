package org.graphpop.procedures;

import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.*;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Computes Cross-Population Extended Haplotype Homozygosity (XP-EHH).
 *
 * <p>Architecture: load → parallel compute → standardize → sequential write.
 * <ol>
 *   <li>Load paired haplotype matrices via {@link HaplotypeMatrix#loadPair}
 *       (single-threaded, one pass over variants)</li>
 *   <li>Compute iHH per population in parallel per variant (ForkJoinPool)</li>
 *   <li>Standardize genome-wide</li>
 *   <li>Write XP-EHH properties on Variant nodes (single-threaded)</li>
 * </ol>
 * </p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.xpehh('chr22', 'EUR', 'AFR', {min_af: 0.05, start: 16000000, end: 17000000})
 * YIELD variantId, pos, af_pop1, af_pop2, xpehh_unstd, xpehh
 * </pre>
 * </p>
 */
public class XPEHHProcedure {

    private static final double DEFAULT_MIN_AF = 0.05;

    @Context
    public Transaction tx;

    public static class XPEHHResult {
        public String variantId;
        public long pos;
        public double af_pop1;
        public double af_pop2;
        public double xpehh_unstd;
        public double xpehh;

        public XPEHHResult() {}
    }

    @Procedure(name = "graphpop.xpehh", mode = Mode.WRITE)
    public Stream<XPEHHResult> xpehh(
            @Name("chr") String chr,
            @Name("pop1") String pop1,
            @Name("pop2") String pop2,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        double minAf = DEFAULT_MIN_AF;
        Long regionStart = null;
        Long regionEnd = null;
        if (options != null) {
            if (options.containsKey("min_af"))
                minAf = ((Number) options.get("min_af")).doubleValue();
            if (options.containsKey("start"))
                regionStart = ((Number) options.get("start")).longValue();
            if (options.containsKey("end"))
                regionEnd = ((Number) options.get("end")).longValue();
        }

        // Phase 1: Load paired haplotype matrices (single-threaded, one variant pass)
        HaplotypeMatrix[] pair = HaplotypeMatrix.loadPair(tx, chr, pop1, pop2,
                regionStart, regionEnd);
        if (pair == null) return Stream.empty();
        HaplotypeMatrix m1 = pair[0];
        HaplotypeMatrix m2 = pair[1];

        // Identify focal variants (polymorphic in at least one population)
        final double minAfFinal = minAf;
        int[] focalIndices = identifyFocal(m1, m2, minAfFinal);
        if (focalIndices.length == 0) return Stream.empty();

        // Phase 2: Parallel iHH computation per variant (no DB access)
        double[] xpehhUnstd = new double[focalIndices.length];
        boolean[] valid = new boolean[focalIndices.length];

        IntStream.range(0, focalIndices.length).parallel().forEach(i -> {
            int vidx = focalIndices[i];

            // For XP-EHH: compute EHH for ALL haplotypes in each population
            int[] allHaps1 = allHaplotypeIndices(m1.nHaplotypes);
            int[] allHaps2 = allHaplotypeIndices(m2.nHaplotypes);

            double iHH1 = EHHComputer.computeIHH(m1, vidx, allHaps1);
            double iHH2 = EHHComputer.computeIHH(m2, vidx, allHaps2);

            if (iHH1 > 0 && iHH2 > 0) {
                xpehhUnstd[i] = Math.log(iHH1 / iHH2);
                valid[i] = true;
            }
        });

        // Phase 3: Genome-wide standardization
        double sum = 0.0;
        int count = 0;
        for (int i = 0; i < focalIndices.length; i++) {
            if (valid[i]) { sum += xpehhUnstd[i]; count++; }
        }
        if (count == 0) return Stream.empty();

        double mean = sum / count;
        double varSum = 0.0;
        for (int i = 0; i < focalIndices.length; i++) {
            if (valid[i]) {
                double d = xpehhUnstd[i] - mean;
                varSum += d * d;
            }
        }
        double std = count > 1 ? Math.sqrt(varSum / (count - 1)) : 1.0;

        // Phase 4: Sequential write to graph
        String propName = "xpehh_" + pop1 + "_" + pop2;
        String propNameUnstd = "xpehh_unstd_" + pop1 + "_" + pop2;

        List<XPEHHResult> results = new ArrayList<>();
        for (int i = 0; i < focalIndices.length; i++) {
            if (!valid[i]) continue;
            int vidx = focalIndices[i];
            double xpehhStd = std > 0 ? (xpehhUnstd[i] - mean) / std : 0.0;

            m1.nodes[vidx].setProperty(propName, xpehhStd);
            m1.nodes[vidx].setProperty(propNameUnstd, xpehhUnstd[i]);

            XPEHHResult xr = new XPEHHResult();
            xr.variantId = m1.variantIds[vidx];
            xr.pos = m1.positions[vidx];
            xr.af_pop1 = m1.afs[vidx];
            xr.af_pop2 = m2.afs[vidx];
            xr.xpehh_unstd = xpehhUnstd[i];
            xr.xpehh = xpehhStd;
            results.add(xr);
        }

        results.sort(Comparator.comparingLong(r -> r.pos));
        return results.stream();
    }

    private static int[] identifyFocal(HaplotypeMatrix m1, HaplotypeMatrix m2, double minAf) {
        List<Integer> focal = new ArrayList<>();
        for (int i = 0; i < m1.nVariants; i++) {
            int an1 = m1.ans[i], an2 = m2.ans[i];
            int ac1 = m1.acs[i], ac2 = m2.acs[i];
            double af1 = m1.afs[i], af2 = m2.afs[i];

            if (an1 < 2 || an2 < 2) continue;
            if (ac1 == 0 && ac2 == 0) continue;
            if (ac1 == an1 && ac2 == an2) continue;

            boolean poly1 = ac1 > 0 && ac1 < an1 && af1 >= minAf && af1 <= (1.0 - minAf);
            boolean poly2 = ac2 > 0 && ac2 < an2 && af2 >= minAf && af2 <= (1.0 - minAf);
            if (!poly1 && !poly2) continue;

            focal.add(i);
        }
        return focal.stream().mapToInt(Integer::intValue).toArray();
    }

    private static int[] allHaplotypeIndices(int nHaplotypes) {
        int[] indices = new int[nHaplotypes];
        for (int i = 0; i < nHaplotypes; i++) indices[i] = i;
        return indices;
    }
}
