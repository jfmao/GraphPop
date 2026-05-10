package org.graphpop.procedures.recombination;

import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.Transaction;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Shared helper for the M13 procedures that consume per-window
 * pairwise LD data over {@code :CARRIES} edges. Extracted from
 * {@link LdDecayProcedure} so multiple consumers (MCMC, HMM,
 * stratified, hotspots) reuse the same Cypher load.
 *
 * <p>v1 constraint: ≤ 63 focal samples per call (64-bit bitmask).</p>
 */
public final class LdPairLoader {

    private LdPairLoader() {}

    /** A pair of in-window variants with their physical distance and r². */
    public static final class Pair {
        public final long distance;
        public final double r2;

        Pair(long distance, double r2) {
            this.distance = distance;
            this.r2 = r2;
        }
    }

    /**
     * Load every in-window variant pair (r², d) over the focal
     * sample set.
     *
     * @param tx               read transaction
     * @param sampleIds        focal sample IDs (≤ 63)
     * @param windowStart      window start (bp, inclusive)
     * @param windowEnd        window end   (bp, exclusive)
     * @param minMaf           minimum derived-allele frequency
     * @param maxPairDistance  maximum pair distance (bp)
     */
    public static List<Pair> loadPairs(Transaction tx,
                                        List<String> sampleIds,
                                        long windowStart, long windowEnd,
                                        double minMaf, long maxPairDistance) {
        Set<String> focalSet = new HashSet<>(sampleIds);
        if (focalSet.size() > 63) {
            throw new IllegalArgumentException(
                    "v1 requires ≤ 63 focal samples; got "
                  + focalSet.size());
        }
        if (focalSet.isEmpty()) return List.of();

        List<String> ordered = new ArrayList<>(focalSet);
        ordered.sort(null);
        Map<String, Integer> bit = new HashMap<>();
        for (int i = 0; i < ordered.size(); i++) bit.put(ordered.get(i), i);

        Result r = tx.execute(
                "MATCH (s:Sample)-[c:CARRIES]->(v:Variant) "
              + "WHERE v.position >= $lo AND v.position < $hi "
              + "  AND s.sampleId IN $sids "
              + "RETURN v.variantId AS vid, v.position AS pos, "
              + "s.sampleId AS sid",
                Map.of("lo", windowStart, "hi", windowEnd, "sids", ordered));
        Map<String, long[]> byVariant = new HashMap<>();
        try {
            while (r.hasNext()) {
                Map<String, Object> row = r.next();
                String vid = (String) row.get("vid");
                long pos = ((Number) row.get("pos")).longValue();
                String sid = (String) row.get("sid");
                Integer b = bit.get(sid);
                if (b == null) continue;
                long[] entry = byVariant.computeIfAbsent(vid,
                        k -> new long[]{pos, 0L});
                entry[1] |= 1L << b;
            }
        } finally {
            r.close();
        }

        int nSamples = ordered.size();
        double minCarriers = minMaf * nSamples;
        double maxCarriers = (1.0 - minMaf) * nSamples;
        List<long[]> variants = new ArrayList<>(byVariant.size());
        for (long[] v : byVariant.values()) {
            int k = Long.bitCount(v[1]);
            if (k <= minCarriers || k >= maxCarriers) continue;
            variants.add(v);
        }
        variants.sort((a, b2) -> Long.compare(a[0], b2[0]));

        List<Pair> pairs = new ArrayList<>();
        for (int i = 0; i < variants.size(); i++) {
            long posI = variants.get(i)[0];
            long carriersI = variants.get(i)[1];
            int kI = Long.bitCount(carriersI);
            for (int j = i + 1; j < variants.size(); j++) {
                long posJ = variants.get(j)[0];
                long d = posJ - posI;
                if (d > maxPairDistance) break;
                long carriersJ = variants.get(j)[1];
                int kJ = Long.bitCount(carriersJ);
                int kAB = Long.bitCount(carriersI & carriersJ);
                double r2 = rSquared(kI, kJ, kAB, nSamples);
                if (Double.isFinite(r2)) {
                    pairs.add(new Pair(d, r2));
                }
            }
        }
        return pairs;
    }

    /** Hill 1968 r² formula (extracted from {@link LdDecayProcedure}). */
    public static double rSquared(int kA, int kB, int kAB, int n) {
        if (n <= 1) return Double.NaN;
        double pA = (double) kA / n;
        double pB = (double) kB / n;
        double pAB = (double) kAB / n;
        double d = pAB - pA * pB;
        double denom = pA * (1.0 - pA) * pB * (1.0 - pB);
        if (denom <= 0.0) return Double.NaN;
        return d * d / denom;
    }
}
