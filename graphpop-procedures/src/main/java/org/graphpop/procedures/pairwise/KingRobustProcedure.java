package org.graphpop.procedures.pairwise;

import org.graphpop.procedures.GenotypeLoader;
import org.graphpop.procedures.PackedGenotypeReader;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Description;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

/**
 * KING-robust pairwise kinship (Manichaikul et al. 2010).
 *
 * <p>Validation baseline for the GraphPop kinship suite — agrees with PLINK2
 * {@code --make-king} to {@code r ≥ 0.999} on standard cohorts. This is the
 * matrix-side procedure; for graph-native branch GRMs (Tang &amp; Chiang 2025),
 * see {@code graphpop.kinship.branch_grm} (depends on the Phase 4 ARG layer).</p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.kinship.king('chr22', 'EUR', {min_snp: 1000, min_phi: -0.5})
 *   YIELD sample_a, sample_b, phi, ibs0, het_het, n_snp, n_aa_min, method
 * </pre>
 *
 * <p>Options (defaults shown):
 * <ul>
 *   <li>{@code start} (1) / {@code end} (Long.MAX_VALUE) — region in bp.</li>
 *   <li>{@code min_snp} (1000) — skip pairs with fewer informative variants.</li>
 *   <li>{@code min_phi} (-0.5) — skip pairs with kinship below this threshold.</li>
 *   <li>{@code samples} — explicit list of sample IDs (overrides {@code pop}).</li>
 *   <li>{@code include_self} (false) — emit {@code (s, s)} self-pairs.</li>
 * </ul>
 */
public class KingRobustProcedure {

    private static final int DEFAULT_MIN_SNP = 1000;
    private static final double DEFAULT_MIN_PHI = -0.5;

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.kinship.king", mode = Mode.READ)
    @Description("KING-robust pairwise kinship (Manichaikul 2010) on packed genotypes. "
            + "Validation baseline for the GraphPop kinship suite.")
    @SuppressWarnings("unchecked")
    public Stream<KinshipResult> king(
            @Name("chr") String chr,
            @Name("pop") String pop,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();

        long start = getLong(options, "start", 1L);
        long end = getLong(options, "end", Long.MAX_VALUE);
        int minSnp = (int) getLong(options, "min_snp", DEFAULT_MIN_SNP);
        double minPhi = getDouble(options, "min_phi", DEFAULT_MIN_PHI);
        boolean includeSelf = getBoolean(options, "include_self", false);
        List<String> sampleList = (List<String>) options.get("samples");

        // 1. Build sample index for the analysis subset.
        Map<String, Integer> sampleIndex = (sampleList != null && !sampleList.isEmpty())
                ? GenotypeLoader.buildSampleIndex(tx, sampleList)
                : GenotypeLoader.buildSampleIndex(tx, pop);

        if (sampleIndex.isEmpty()) return Stream.empty();

        int nSamples = sampleIndex.size();
        String[] sampleIds = new String[nSamples];
        for (Map.Entry<String, Integer> e : sampleIndex.entrySet()) {
            sampleIds[e.getValue()] = e.getKey();
        }
        int[] packedIndices = GenotypeLoader.buildPackedIndices(tx, sampleIndex);

        // 2. Stream variants on the chromosome and accumulate per-pair stats.
        int nPairs = nSamples * (nSamples + 1) / 2;
        KingRobustComputer.PairAccumulator[] acc =
                new KingRobustComputer.PairAccumulator[nPairs];
        for (int p = 0; p < nPairs; p++) {
            acc[p] = new KingRobustComputer.PairAccumulator();
        }
        int[] gtBuf = new int[nSamples];

        Result result = tx.execute(
                "MATCH (v:Variant) "
              + "WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end "
              + "RETURN v ORDER BY v.pos",
                Map.of("chr", chr, "start", start, "end", end));
        try {
            while (result.hasNext()) {
                Map<String, Object> row = result.next();
                Node variant = (Node) row.get("v");
                Object gtPackedObj = variant.getProperty("gt_packed", null);
                if (!(gtPackedObj instanceof byte[])) continue;
                byte[] gtPacked = (byte[]) gtPackedObj;

                for (int s = 0; s < nSamples; s++) {
                    int pi = packedIndices[s];
                    gtBuf[s] = (pi >= 0)
                            ? PackedGenotypeReader.genotype(gtPacked, pi)
                            : PackedGenotypeReader.GT_MISSING;
                }
                int p = 0;
                for (int b = 0; b < nSamples; b++) {
                    for (int a = 0; a <= b; a++) {
                        acc[p++].accumulate(gtBuf[a], gtBuf[b]);
                    }
                }
            }
        } finally {
            result.close();
        }

        // 3. Build result rows; apply min_snp and min_phi filters.
        List<KinshipResult> rows = new ArrayList<>();
        int p = 0;
        for (int b = 0; b < nSamples; b++) {
            for (int a = 0; a <= b; a++, p++) {
                if (a == b && !includeSelf) continue;
                KingRobustComputer.PairAccumulator s = acc[p];
                if (s.nUsed < minSnp) continue;
                double phi = s.phi();
                // NaN means undefined (one sample has zero het sites).
                if (Double.isNaN(phi)) continue;
                if (phi < minPhi) continue;
                int nAaMin = Math.min(s.nAaI, s.nAaJ);
                rows.add(new KinshipResult(
                        sampleIds[a], sampleIds[b],
                        phi, s.ibs0, s.hetHet, s.nUsed, nAaMin, "king-robust"));
            }
        }
        return rows.stream();
    }

    private static long getLong(Map<String, Object> opts, String key, long def) {
        Object v = opts.get(key);
        if (v instanceof Number) return ((Number) v).longValue();
        return def;
    }

    private static double getDouble(Map<String, Object> opts, String key, double def) {
        Object v = opts.get(key);
        if (v instanceof Number) return ((Number) v).doubleValue();
        return def;
    }

    private static boolean getBoolean(Map<String, Object> opts, String key, boolean def) {
        Object v = opts.get(key);
        if (v instanceof Boolean) return (Boolean) v;
        return def;
    }
}
