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
 * <p>Supports two methods via the {@code method} option:</p>
 * <ul>
 *   <li>{@code "hmm"} (default) — 2-state Viterbi HMM (Narasimhan et al. 2016,
 *       bcftools roh). Uses allele-frequency-weighted emissions and
 *       recombination-scaled transitions.</li>
 *   <li>{@code "sliding_window"} — PLINK-style sliding window with
 *       configurable window_snps, max_het, density, and gap filters.</li>
 * </ul>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.roh('chr22', 'EUR', {min_length: 1000000})
 * YIELD sampleId, n_roh, total_length, froh, mean_length, max_length
 *
 * CALL graphpop.roh('chr22', 'EUR', {method: 'hmm', az_length: 2000000})
 * CALL graphpop.roh('chr22', 'EUR', {method: 'sliding_window', window_snps: 50})
 * </pre>
 * </p>
 */
public class ROHProcedure {

    private static final int DEFAULT_WINDOW_SNPS = 50;
    private static final int DEFAULT_MAX_HET = 1;
    private static final long DEFAULT_MIN_LENGTH = 1_000_000;
    private static final int DEFAULT_MIN_SNPS = 100;
    /** Max gap (bp) between consecutive SNPs within a ROH. PLINK default = 1000 kb. */
    private static final long DEFAULT_MAX_GAP = 1_000_000;
    /** Min SNP density: at least 1 SNP per this many bp. PLINK default = 50 kb. */
    private static final long DEFAULT_DENSITY_KB = 50;

    /** GRCh38/hg38 chromosome lengths (bp). Used for accurate FROH computation. */
    private static final java.util.Map<String, Long> CHR_LENGTHS_HG38 =
            java.util.Map.ofEntries(
                java.util.Map.entry("chr1",  248956422L),
                java.util.Map.entry("chr2",  242193529L),
                java.util.Map.entry("chr3",  198295559L),
                java.util.Map.entry("chr4",  190214555L),
                java.util.Map.entry("chr5",  181538259L),
                java.util.Map.entry("chr6",  170805979L),
                java.util.Map.entry("chr7",  159345973L),
                java.util.Map.entry("chr8",  145138636L),
                java.util.Map.entry("chr9",  138394717L),
                java.util.Map.entry("chr10", 133797422L),
                java.util.Map.entry("chr11", 135086622L),
                java.util.Map.entry("chr12", 133275309L),
                java.util.Map.entry("chr13", 114364328L),
                java.util.Map.entry("chr14", 107043718L),
                java.util.Map.entry("chr15", 101991189L),
                java.util.Map.entry("chr16",  90338345L),
                java.util.Map.entry("chr17",  83257441L),
                java.util.Map.entry("chr18",  80373285L),
                java.util.Map.entry("chr19",  58617616L),
                java.util.Map.entry("chr20",  64444167L),
                java.util.Map.entry("chr21",  46709983L),
                java.util.Map.entry("chr22",  50818468L),
                java.util.Map.entry("chrX",  156040895L),
                java.util.Map.entry("chrY",   57227415L)
            );

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
        // Default to SNPs only — PLINK 1.9 --homozyg uses SNPs only when reading
        // VCF. Indels inflate variant density, making windows physically shorter and
        // detecting spurious ROH. Users can override with variant_type: 'ALL'.
        if (options == null) {
            options = Map.of("variant_type", "SNP");
        } else if (!options.containsKey("variant_type")) {
            options = new HashMap<>(options);
            options.put("variant_type", "SNP");
        } else if ("ALL".equals(options.get("variant_type"))) {
            options = new HashMap<>(options);
            options.remove("variant_type");  // no type filter
        }

        long minLength = getLong(options, "min_length", DEFAULT_MIN_LENGTH);
        int minSnps = getInt(options, "min_snps", DEFAULT_MIN_SNPS);
        String method = getString(options, "method", "hmm");

        List<String> sampleList = options != null ? (List<String>) options.get("samples") : null;
        boolean useSubset = sampleList != null && !sampleList.isEmpty();

        // Determine variant type filter for matrix loading
        String variantType = null;
        Object vtOpt = options.get("variant_type");
        if (vtOpt instanceof String) variantType = (String) vtOpt;

        // Minimum variants required depends on method
        int minVariants = "sliding_window".equals(method)
                ? getInt(options, "window_snps", DEFAULT_WINDOW_SNPS)
                : 2;

        // Phase 1: Load HaplotypeMatrix (SNP-only by default)
        HaplotypeMatrix matrix;
        if (useSubset) {
            Map<String, Integer> sampleIndex = GenotypeLoader.buildSampleIndex(tx, sampleList);
            matrix = HaplotypeMatrix.load(tx, chr, sampleIndex, null, null, variantType);
        } else {
            matrix = HaplotypeMatrix.load(tx, chr, pop, null, null, variantType);
        }
        if (matrix == null || matrix.nVariants < minVariants) return Stream.empty();

        int nSamples = matrix.nSamples;
        int nVariants = matrix.nVariants;

        // Chromosome length for FROH computation.
        // Use known hg38 lengths; fall back to variant span if chromosome unknown.
        long chrLength = CHR_LENGTHS_HG38.getOrDefault(chr,
                matrix.positions[nVariants - 1] - matrix.positions[0]);
        if (chrLength <= 0) chrLength = 1;

        // Sample IDs in order
        String[] sampleIds = new String[nSamples];
        for (Map.Entry<String, Integer> entry : matrix.sampleIndex.entrySet()) {
            sampleIds[entry.getValue()] = entry.getKey();
        }

        final long chrLenFinal = chrLength;
        ROHResult[] results = new ROHResult[nSamples];

        if ("hmm".equals(method)) {
            // ---- HMM-based ROH (Narasimhan et al. 2016) ----
            ROHHmm hmm = ROHHmm.fromOptions(options);

            IntStream.range(0, nSamples).parallel().forEach(s -> {
                // Extract diploid genotype for this sample (bit-packed access)
                int[] genotypes = new int[nVariants];
                for (int v = 0; v < nVariants; v++) {
                    genotypes[v] = matrix.hap(v, 2 * s) + matrix.hap(v, 2 * s + 1);
                }

                // Run Viterbi
                List<long[]> segments = hmm.viterbi(
                        matrix.positions, matrix.afs, genotypes, minLength, minSnps);

                // Compute summary stats
                long totalLen = 0, maxLen = 0;
                for (long[] seg : segments) {
                    long len = seg[1] - seg[0];
                    totalLen += len;
                    if (len > maxLen) maxLen = len;
                }

                ROHResult r = new ROHResult();
                r.sampleId = sampleIds[s];
                r.n_roh = segments.size();
                r.total_length = totalLen;
                r.froh = (double) totalLen / chrLenFinal;
                r.mean_length = segments.isEmpty() ? 0 : (double) totalLen / segments.size();
                r.max_length = maxLen;
                results[s] = r;
            });
        } else {
            // ---- Sliding window ROH (PLINK-style) ----
            int windowSnps = getInt(options, "window_snps", DEFAULT_WINDOW_SNPS);
            int maxHet = getInt(options, "max_het", DEFAULT_MAX_HET);
            long maxGap = getLong(options, "max_gap", DEFAULT_MAX_GAP);
            long densityKb = getLong(options, "density_kb", DEFAULT_DENSITY_KB);
            double windowThreshold = getDouble(options, "window_threshold", 0.05);

            final long maxGapFinal = maxGap;
            final long densityKbFinal = densityKb;
            final double windowThresholdFinal = windowThreshold;

            IntStream.range(0, nSamples).parallel().forEach(s -> {
                // Build heterozygosity flag array for this sample (bit-packed)
                boolean[] isHet = new boolean[nVariants];
                for (int v = 0; v < nVariants; v++) {
                    int h0 = matrix.hap(v, 2 * s);
                    int h1 = matrix.hap(v, 2 * s + 1);
                    isHet[v] = (h0 != h1);
                }

                // Sliding window scan: count how many overlapping windows each SNP
                // participates in, and how many of those are homozygous "hits".
                int[] hitCount = new int[nVariants];
                int[] totalCount = new int[nVariants];
                int hetCount = 0;

                // Initialize first window
                for (int v = 0; v < windowSnps && v < nVariants; v++) {
                    if (isHet[v]) hetCount++;
                }
                boolean isHit = hetCount <= maxHet;
                for (int v = 0; v < windowSnps && v < nVariants; v++) {
                    totalCount[v]++;
                    if (isHit) hitCount[v]++;
                }

                // Slide the window
                for (int v = 1; v <= nVariants - windowSnps; v++) {
                    if (isHet[v - 1]) hetCount--;
                    if (isHet[v + windowSnps - 1]) hetCount++;

                    isHit = hetCount <= maxHet;
                    for (int w = v; w < v + windowSnps; w++) {
                        totalCount[w]++;
                        if (isHit) hitCount[w]++;
                    }
                }

                // A SNP is in a ROH if its window hit rate >= threshold
                boolean[] homWindow = new boolean[nVariants];
                for (int v = 0; v < nVariants; v++) {
                    if (totalCount[v] > 0) {
                        double rate = (double) hitCount[v] / totalCount[v];
                        homWindow[v] = rate >= windowThresholdFinal;
                    }
                }

                // Merge contiguous homWindow=true segments into ROH segments
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

                // Split segments at internal gaps > maxGap
                List<long[]> splitSegments = new ArrayList<>();
                for (long[] seg : segments) {
                    int sStart = (int) seg[0];
                    int sEnd = (int) seg[1];
                    int subStart = sStart;
                    for (int v = sStart + 1; v <= sEnd; v++) {
                        long gap = matrix.positions[v] - matrix.positions[v - 1];
                        if (gap > maxGapFinal) {
                            if (v - 1 > subStart) {
                                splitSegments.add(new long[]{subStart, v - 1});
                            }
                            subStart = v;
                        }
                    }
                    splitSegments.add(new long[]{subStart, sEnd});
                }

                // Filter segments by min_length, min_snps, and density
                long totalLength = 0;
                long maxLen = 0;
                int nRoh = 0;

                for (long[] seg : splitSegments) {
                    int startIdx = (int) seg[0];
                    int endIdx = (int) seg[1];
                    int nSnpsInSeg = endIdx - startIdx + 1;
                    long segLength = matrix.positions[endIdx] - matrix.positions[startIdx];

                    if (segLength < minLength || nSnpsInSeg < minSnps) continue;

                    // Density check
                    if (densityKbFinal > 0 && segLength > 0) {
                        double kbPerSnp = (segLength / 1000.0) / nSnpsInSeg;
                        if (kbPerSnp > densityKbFinal) continue;
                    }

                    nRoh++;
                    totalLength += segLength;
                    if (segLength > maxLen) maxLen = segLength;
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
        }

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

    private static double getDouble(Map<String, Object> options, String key, double defaultValue) {
        if (options == null) return defaultValue;
        Object val = options.get(key);
        if (val == null) return defaultValue;
        return ((Number) val).doubleValue();
    }

    private static String getString(Map<String, Object> options, String key, String defaultValue) {
        if (options == null) return defaultValue;
        Object val = options.get(key);
        if (val == null) return defaultValue;
        return val.toString();
    }
}
