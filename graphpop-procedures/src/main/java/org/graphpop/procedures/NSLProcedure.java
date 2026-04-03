package org.graphpop.procedures;

import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Description;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.*;
import java.util.stream.Stream;

/**
 * Computes the nSL (number of Segregating sites by Length) statistic
 * for detecting recent positive selection.
 *
 * <p>Implements the original algorithm from Ferrer-Admetlla et al. (2014):
 * for each focal variant, nSL = log(SL1 / SL0) where SL1 is the mean
 * pairwise Shared Suffix Length among haplotypes carrying allele 1, and
 * SL0 is the mean among allele 0 carriers. Shared suffix length for a pair
 * is the number of contiguous sites (from the focal outward) where the two
 * haplotypes carry identical alleles.</p>
 *
 * <p>Architecture: load &rarr; forward+backward pairwise SSL scan &rarr;
 * AF-bin standardize &rarr; sequential write.</p>
 *
 * <p>Reference: Ferrer-Admetlla A, Liang M, Korneliussen T, Nielsen R.
 * On detecting incomplete soft or hard selective sweeps using haplotype structure.
 * Mol Biol Evol. 2014;31(5):1275-91.</p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.nsl('chr22', 'EUR', {min_af: 0.05, start: 16000000, end: 17000000})
 * YIELD variantId, pos, af, nsl_unstd, nsl
 * </pre>
 * </p>
 */
public class NSLProcedure {

    private static final double DEFAULT_MIN_AF = 0.05;
    private static final int DEFAULT_AF_BINS = 20;

    @Context
    public Transaction tx;

    public static class NSLResult {
        public String variantId;
        public long pos;
        public double af;
        public double nsl_unstd;
        public double nsl;

        public NSLResult() {}
    }

    @Procedure(name = "graphpop.nsl", mode = Mode.WRITE)
    @Description("Compute the number of segregating sites by length (nSL) statistic for detecting selection.")
    @SuppressWarnings("unchecked")
    public Stream<NSLResult> nsl(
            @Name("chr") String chr,
            @Name("pop") String pop,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        // Default to SNPs only — nSL was designed for SNP data and indels in
        // the scan path alter SSL accumulation. Users can override with
        // variant_type: 'ALL' to include all biallelic variants.
        if (options == null) {
            options = Map.of("variant_type", "SNP");
        } else if (!options.containsKey("variant_type")) {
            options = new HashMap<>(options);
            options.put("variant_type", "SNP");
        } else if ("ALL".equals(options.get("variant_type"))) {
            options = new HashMap<>(options);
            options.remove("variant_type");  // no type filter
        }
        VariantFilter filter = VariantFilter.fromOptions(options);
        double minAf = filter.minAf > 0 ? filter.minAf : DEFAULT_MIN_AF;
        String variantType = filter.variantType;
        Long regionStart = null;
        Long regionEnd = null;
        if (options != null) {
            if (options.containsKey("start"))
                regionStart = ((Number) options.get("start")).longValue();
            if (options.containsKey("end"))
                regionEnd = ((Number) options.get("end")).longValue();
        }

        List<String> sampleList = options != null ? (List<String>) options.get("samples") : null;
        boolean useSubset = sampleList != null && !sampleList.isEmpty();
        int afBins = options != null && options.containsKey("n_af_bins")
                ? ((Number) options.get("n_af_bins")).intValue() : DEFAULT_AF_BINS;

        // Phase 1: Load dense haplotype matrix (single-threaded)
        HaplotypeMatrix matrix;
        if (useSubset) {
            Map<String, Integer> sampleIndex = GenotypeLoader.buildSampleIndex(tx, sampleList);
            matrix = HaplotypeMatrix.load(tx, chr, sampleIndex, regionStart, regionEnd, variantType);
        } else {
            matrix = HaplotypeMatrix.load(tx, chr, pop, regionStart, regionEnd, variantType);
        }
        if (matrix == null || matrix.nVariants == 0) return Stream.empty();

        // Identify focal variants: common, polymorphic, passing filters
        int[] focalIndices = identifyFocal(matrix, minAf, filter);
        if (focalIndices.length == 0) return Stream.empty();

        int nFocal = focalIndices.length;
        int nHaps = matrix.nHaplotypes;

        // Phase 2: Pairwise SSL scan (Ferrer-Admetlla 2014)
        // Forward and backward scans run in parallel (independent state)
        double[] sl0Fwd = new double[nFocal];
        double[] sl1Fwd = new double[nFocal];
        double[] sl0Rev = new double[nFocal];
        double[] sl1Rev = new double[nFocal];

        Thread fwdThread = new Thread(() ->
                scanDirection(matrix, focalIndices, nHaps, sl0Fwd, sl1Fwd, +1), "nSL-fwd");
        Thread revThread = new Thread(() ->
                scanDirection(matrix, focalIndices, nHaps, sl0Rev, sl1Rev, -1), "nSL-rev");

        fwdThread.start();
        revThread.start();
        try {
            fwdThread.join();
            revThread.join();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            return Stream.empty();
        }

        // Combine forward + backward: nSL = log(SL1 / SL0)
        double[] nslUnstd = new double[nFocal];
        boolean[] valid = new boolean[nFocal];

        for (int i = 0; i < nFocal; i++) {
            double sl0 = sl0Fwd[i] + sl0Rev[i];
            double sl1 = sl1Fwd[i] + sl1Rev[i];
            if (sl0 > 0 && sl1 > 0) {
                nslUnstd[i] = Math.log(sl1 / sl0);
                valid[i] = true;
            }
        }

        // Phase 3: Standardize within allele-frequency bins
        Map<Integer, List<Integer>> bins = new HashMap<>();
        for (int i = 0; i < nFocal; i++) {
            if (!valid[i]) continue;
            int bin = Math.min((int) (matrix.afs[focalIndices[i]] * afBins), afBins - 1);
            bins.computeIfAbsent(bin, k -> new ArrayList<>()).add(i);
        }

        double[] nslStd = new double[nFocal];
        for (List<Integer> binList : bins.values()) {
            double sum = 0.0;
            for (int idx : binList) sum += nslUnstd[idx];
            double mean = sum / binList.size();

            double varSum = 0.0;
            for (int idx : binList) {
                double d = nslUnstd[idx] - mean;
                varSum += d * d;
            }
            double std = binList.size() > 1 ? Math.sqrt(varSum / (binList.size() - 1)) : 1.0;

            for (int idx : binList) {
                nslStd[idx] = std > 0 ? (nslUnstd[idx] - mean) / std : 0.0;
            }
        }

        // Phase 4: Sequential write to graph
        List<NSLResult> results = new ArrayList<>();
        for (int i = 0; i < nFocal; i++) {
            if (!valid[i]) continue;
            int vidx = focalIndices[i];

            matrix.nodes[vidx].setProperty("nsl_" + pop, nslStd[i]);
            matrix.nodes[vidx].setProperty("nsl_unstd_" + pop, nslUnstd[i]);

            NSLResult nr = new NSLResult();
            nr.variantId = matrix.variantIds[vidx];
            nr.pos = matrix.positions[vidx];
            nr.af = matrix.afs[vidx];
            nr.nsl_unstd = nslUnstd[i];
            nr.nsl = nslStd[i];
            results.add(nr);
        }

        results.sort(Comparator.comparingLong(r -> r.pos));
        return results.stream();
    }

    /**
     * Scan focal variants in one direction computing mean pairwise SSL per allele class.
     *
     * <p>Maintains an int[] ssl of size nHaps*(nHaps-1)/2 tracking the shared suffix
     * length for every haplotype pair. At each focal variant, for each pair (j,k):
     * if both carry the same allele at the current site, ssl[pair] += 1 and the
     * accumulated value is added to the appropriate allele-class sum. If they differ,
     * ssl[pair] is reset to 0.</p>
     *
     * <p>This faithfully reproduces the algorithm from Ferrer-Admetlla et al. (2014)
     * and matches scikit-allel's nsl01_scan implementation.</p>
     *
     * @param matrix       dense haplotype matrix
     * @param focalIndices indices of focal variants in position order
     * @param nHaps        number of haplotypes
     * @param sl0Out       output: mean SSL for allele-0 pairs at each focal variant
     * @param sl1Out       output: mean SSL for allele-1 pairs at each focal variant
     * @param direction    +1 for forward (ascending index), -1 for backward
     */
    static void scanDirection(HaplotypeMatrix matrix, int[] focalIndices,
                              int nHaps, double[] sl0Out, double[] sl1Out,
                              int direction) {
        int nFocal = focalIndices.length;
        long nPairs = (long) nHaps * (nHaps - 1) / 2;
        int[] ssl = new int[(int) nPairs];

        // Build a lookup: matrix variant index -> focal array index (or -1)
        // This lets us quickly check if a variant is focal during the scan
        int maxIdx = focalIndices[nFocal - 1];
        int minIdx = focalIndices[0];

        // Iterate focal variants in scan order
        int startF, endF, stepF;
        if (direction > 0) {
            startF = 0;
            endF = nFocal;
            stepF = 1;
        } else {
            startF = nFocal - 1;
            endF = -1;
            stepF = -1;
        }

        for (int fi = startF; fi != endF; fi += stepF) {
            int vidx = focalIndices[fi];
            byte[] haps = matrix.haplotypes[vidx];

            // Update SSL for all pairs and accumulate per-allele sums
            long ssl00Sum = 0, ssl11Sum = 0;
            long n00 = 0, n11 = 0;
            int u = 0;

            for (int j = 0; j < nHaps; j++) {
                int a1 = (haps[j >> 3] >> (j & 7)) & 1;
                for (int k = j + 1; k < nHaps; k++) {
                    int a2 = (haps[k >> 3] >> (k & 7)) & 1;
                    if (a1 == a2) {
                        int l = ssl[u] + 1;
                        ssl[u] = l;
                        if (a1 == 0) {
                            ssl00Sum += l;
                            n00++;
                        } else {
                            ssl11Sum += l;
                            n11++;
                        }
                    } else {
                        ssl[u] = 0;
                    }
                    u++;
                }
            }

            // Mean SSL per allele class
            sl0Out[fi] = n00 > 0 ? (double) ssl00Sum / n00 : 0.0;
            sl1Out[fi] = n11 > 0 ? (double) ssl11Sum / n11 : 0.0;
        }
    }

    private static int[] identifyFocal(HaplotypeMatrix matrix, double minAf, VariantFilter filter) {
        List<Integer> focal = new ArrayList<>();
        String variantType = filter.variantType;
        for (int i = 0; i < matrix.nVariants; i++) {
            int ac = matrix.acs[i];
            int an = matrix.ans[i];
            double af = matrix.afs[i];
            if (an < 2 || ac == 0 || ac == an) continue;
            if (af < minAf || af > (1.0 - minAf)) continue;

            if (filter.isActive()) {
                if (af < filter.minAf || af > filter.maxAf) continue;
            }

            // Check variant_type from node property (crucial for nSL scan path)
            if (variantType != null && matrix.nodes != null) {
                Object vt = matrix.nodes[i].getProperty("variant_type", null);
                if (vt == null || !variantType.equals(vt)) continue;
            }

            focal.add(i);
        }
        return focal.stream().mapToInt(Integer::intValue).toArray();
    }
}
