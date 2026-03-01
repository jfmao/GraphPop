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
 * Computes Runs of Homozygosity (ROH) for each sample in a population.
 *
 * <p>Architecture: FULL PATH, 4-phase:
 * <ol>
 *   <li>Load HaplotypeMatrix (single-threaded)</li>
 *   <li>Per-sample parallel scanning — sliding window heterozygosity check</li>
 *   <li>Merge windows into ROH segments, filter by min_length/min_snps</li>
 *   <li>Return per-sample summary as result stream</li>
 * </ol>
 * </p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.roh('chr22', 'EUR', {min_length: 1000000, min_snps: 50, max_het: 1, window_snps: 50})
 * YIELD sampleId, n_roh, total_length, froh, mean_length, max_length
 * </pre>
 * </p>
 */
public class ROHProcedure {

    private static final int DEFAULT_WINDOW_SNPS = 50;
    private static final int DEFAULT_MAX_HET = 1;
    private static final long DEFAULT_MIN_LENGTH = 1_000_000;
    private static final int DEFAULT_MIN_SNPS = 50;

    @Context
    public Transaction tx;

    public static class ROHResult {
        public String sampleId;
        public long n_roh;
        public long total_length;
        public double froh;
        public double mean_length;
        public long max_length;

        public ROHResult() {}
    }

    @Procedure(name = "graphpop.roh", mode = Mode.READ)
    @SuppressWarnings("unchecked")
    public Stream<ROHResult> roh(
            @Name("chr") String chr,
            @Name("pop") String pop,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        int windowSnps = getInt(options, "window_snps", DEFAULT_WINDOW_SNPS);
        int maxHet = getInt(options, "max_het", DEFAULT_MAX_HET);
        long minLength = getLong(options, "min_length", DEFAULT_MIN_LENGTH);
        int minSnps = getInt(options, "min_snps", DEFAULT_MIN_SNPS);

        List<String> sampleList = options != null ? (List<String>) options.get("samples") : null;
        boolean useSubset = sampleList != null && !sampleList.isEmpty();

        // Phase 1: Load HaplotypeMatrix
        HaplotypeMatrix matrix;
        if (useSubset) {
            Map<String, Integer> sampleIndex = GenotypeLoader.buildSampleIndex(tx, sampleList);
            matrix = HaplotypeMatrix.load(tx, chr, sampleIndex, null, null);
        } else {
            matrix = HaplotypeMatrix.load(tx, chr, pop, null, null);
        }
        if (matrix == null || matrix.nVariants < windowSnps) return Stream.empty();

        int nSamples = matrix.nSamples;
        int nVariants = matrix.nVariants;

        // Chromosome span for FROH computation
        long chrStart = matrix.positions[0];
        long chrEnd = matrix.positions[nVariants - 1];
        long chrLength = chrEnd - chrStart;
        if (chrLength <= 0) chrLength = 1;

        // Sample IDs in order
        String[] sampleIds = new String[nSamples];
        for (Map.Entry<String, Integer> entry : matrix.sampleIndex.entrySet()) {
            sampleIds[entry.getValue()] = entry.getKey();
        }

        // Phase 2 + 3: Per-sample parallel ROH detection
        final long chrLenFinal = chrLength;
        ROHResult[] results = new ROHResult[nSamples];

        IntStream.range(0, nSamples).parallel().forEach(s -> {
            // Build heterozygosity flag array for this sample
            boolean[] isHet = new boolean[nVariants];
            for (int v = 0; v < nVariants; v++) {
                int h0 = matrix.haplotypes[v][2 * s];
                int h1 = matrix.haplotypes[v][2 * s + 1];
                isHet[v] = (h0 != h1);
            }

            // Sliding window scan: mark each SNP as "in homozygous window"
            boolean[] homWindow = new boolean[nVariants];
            int hetCount = 0;

            // Initialize first window
            for (int v = 0; v < windowSnps && v < nVariants; v++) {
                if (isHet[v]) hetCount++;
            }
            if (hetCount <= maxHet) {
                for (int v = 0; v < windowSnps && v < nVariants; v++) {
                    homWindow[v] = true;
                }
            }

            // Slide the window
            for (int v = 1; v <= nVariants - windowSnps; v++) {
                // Remove leftmost SNP of previous window
                if (isHet[v - 1]) hetCount--;
                // Add rightmost SNP of new window
                if (isHet[v + windowSnps - 1]) hetCount++;

                if (hetCount <= maxHet) {
                    for (int w = v; w < v + windowSnps; w++) {
                        homWindow[w] = true;
                    }
                }
            }

            // Phase 3: Merge contiguous homWindow=true segments into ROH segments
            List<long[]> segments = new ArrayList<>();
            int segStart = -1;
            for (int v = 0; v < nVariants; v++) {
                if (homWindow[v]) {
                    if (segStart < 0) segStart = v;
                } else {
                    if (segStart >= 0) {
                        segments.add(new long[]{segStart, v - 1});
                        segStart = -1;
                    }
                }
            }
            if (segStart >= 0) {
                segments.add(new long[]{segStart, nVariants - 1});
            }

            // Filter segments by min_length and min_snps
            long totalLength = 0;
            long maxLen = 0;
            int nRoh = 0;

            for (long[] seg : segments) {
                int startIdx = (int) seg[0];
                int endIdx = (int) seg[1];
                int nSnpsInSeg = endIdx - startIdx + 1;
                long segLength = matrix.positions[endIdx] - matrix.positions[startIdx];

                if (segLength >= minLength && nSnpsInSeg >= minSnps) {
                    nRoh++;
                    totalLength += segLength;
                    if (segLength > maxLen) maxLen = segLength;
                }
            }

            ROHResult r = new ROHResult();
            r.sampleId = sampleIds[s];
            r.n_roh = nRoh;
            r.total_length = totalLength;
            r.froh = (double) totalLength / chrLenFinal;
            r.mean_length = nRoh > 0 ? (double) totalLength / nRoh : 0.0;
            r.max_length = maxLen;
            results[s] = r;
        });

        return Arrays.stream(results)
                .sorted(Comparator.comparing(r -> r.sampleId));
    }

    private static int getInt(Map<String, Object> options, String key, int defaultValue) {
        if (options == null) return defaultValue;
        Object val = options.get(key);
        if (val == null) return defaultValue;
        return ((Number) val).intValue();
    }

    private static long getLong(Map<String, Object> options, String key, long defaultValue) {
        if (options == null) return defaultValue;
        Object val = options.get(key);
        if (val == null) return defaultValue;
        return ((Number) val).longValue();
    }
}
