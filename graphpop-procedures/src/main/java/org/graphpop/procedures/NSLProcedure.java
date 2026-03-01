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
 * Computes the nSL (number of Segregating sites by Length) statistic
 * for detecting recent positive selection.
 *
 * <p>nSL is similar to iHS but counts SNP positions traversed instead of
 * integrating over physical distance. This makes it robust to variation in
 * recombination rate and does not require a genetic map.</p>
 *
 * <p>Architecture: load → parallel SL compute → AF-bin standardize → sequential write.
 * Follows the same 4-phase pattern as {@link IHSProcedure}.</p>
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
    private static final int AF_BINS = 20;

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
    @SuppressWarnings("unchecked")
    public Stream<NSLResult> nsl(
            @Name("chr") String chr,
            @Name("pop") String pop,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        VariantFilter filter = VariantFilter.fromOptions(options);
        double minAf = filter.minAf > 0 ? filter.minAf : DEFAULT_MIN_AF;
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

        // Phase 1: Load dense haplotype matrix (single-threaded)
        HaplotypeMatrix matrix;
        if (useSubset) {
            Map<String, Integer> sampleIndex = GenotypeLoader.buildSampleIndex(tx, sampleList);
            matrix = HaplotypeMatrix.load(tx, chr, sampleIndex, regionStart, regionEnd);
        } else {
            matrix = HaplotypeMatrix.load(tx, chr, pop, regionStart, regionEnd);
        }
        if (matrix == null || matrix.nVariants == 0) return Stream.empty();

        // Identify focal variants: common, polymorphic
        int[] focalIndices = identifyFocal(matrix, minAf, filter);
        if (focalIndices.length == 0) return Stream.empty();

        // Phase 2: Parallel SL computation (no DB access)
        double[] nslUnstd = new double[focalIndices.length];
        boolean[] valid = new boolean[focalIndices.length];

        IntStream.range(0, focalIndices.length).parallel().forEach(i -> {
            int vidx = focalIndices[i];

            byte[] focalHaps = matrix.haplotypes[vidx];
            int nHaps = matrix.nHaplotypes;

            int nDerived = 0;
            for (int h = 0; h < nHaps; h++) {
                if (focalHaps[h] == 1) nDerived++;
            }
            int nAncestral = nHaps - nDerived;

            if (nDerived < 2 || nAncestral < 2) return;

            int[] derivedIdx = new int[nDerived];
            int[] ancestralIdx = new int[nAncestral];
            int di = 0, ai = 0;
            for (int h = 0; h < nHaps; h++) {
                if (focalHaps[h] == 1) {
                    derivedIdx[di++] = h;
                } else {
                    ancestralIdx[ai++] = h;
                }
            }

            int slDerived = EHHComputer.computeSL(matrix, vidx, derivedIdx);
            int slAncestral = EHHComputer.computeSL(matrix, vidx, ancestralIdx);

            if (slDerived > 0 && slAncestral > 0) {
                nslUnstd[i] = Math.log((double) slAncestral / slDerived);
                valid[i] = true;
            }
        });

        // Phase 3: Standardize within allele-frequency bins
        Map<Integer, List<Integer>> bins = new HashMap<>();
        for (int i = 0; i < focalIndices.length; i++) {
            if (!valid[i]) continue;
            int bin = Math.min((int) (matrix.afs[focalIndices[i]] * AF_BINS), AF_BINS - 1);
            bins.computeIfAbsent(bin, k -> new ArrayList<>()).add(i);
        }

        double[] nslStd = new double[focalIndices.length];
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
        for (int i = 0; i < focalIndices.length; i++) {
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

    private static int[] identifyFocal(HaplotypeMatrix matrix, double minAf, VariantFilter filter) {
        List<Integer> focal = new ArrayList<>();
        for (int i = 0; i < matrix.nVariants; i++) {
            int ac = matrix.acs[i];
            int an = matrix.ans[i];
            double af = matrix.afs[i];
            if (an < 2 || ac == 0 || ac == an) continue;
            if (af < minAf || af > (1.0 - minAf)) continue;

            if (filter.isActive()) {
                if (af < filter.minAf || af > filter.maxAf) continue;
            }

            focal.add(i);
        }
        return focal.stream().mapToInt(Integer::intValue).toArray();
    }
}
