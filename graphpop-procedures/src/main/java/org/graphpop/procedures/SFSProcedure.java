package org.graphpop.procedures;

import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

/**
 * Computes the Site Frequency Spectrum (SFS) for a population in a genomic region.
 *
 * <p>FAST PATH: uses only allele-count arrays on Variant nodes.</p>
 *
 * <p>When {@code folded == false} and ancestral allele annotations are present,
 * produces a polarized unfolded SFS using derived allele counts.</p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.sfs('chr22', 16000000, 17000000, 'AFR', false, {min_af: 0.01})
 * YIELD sfs, n_variants, max_ac, n_polarized
 * </pre>
 * </p>
 */
public class SFSProcedure {

    @Context
    public Transaction tx;

    public static class SFSResult {
        public List<Long> sfs;
        public long n_variants;
        public long max_ac;
        public long n_polarized;

        public SFSResult() {}
    }

    @Procedure(name = "graphpop.sfs", mode = Mode.READ)
    @SuppressWarnings("unchecked")
    public Stream<SFSResult> sfs(
            @Name("chr") String chr,
            @Name("start") long start,
            @Name("end") long end,
            @Name("pop") String pop,
            @Name(value = "folded", defaultValue = "false") boolean folded,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        VariantFilter filter = VariantFilter.fromOptions(options);
        List<String> sampleList = options != null ? (List<String>) options.get("samples") : null;
        boolean useSubset = sampleList != null && !sampleList.isEmpty();

        Map<String, Integer> subsetIndex = null;
        if (useSubset) {
            subsetIndex = GenotypeLoader.buildSampleIndex(tx, sampleList);
        }

        int popIndex = -1;
        int maxAn = 0;
        long nPolarized = 0;

        List<int[]> acAnPairs = new ArrayList<>();  // [ac, an] or [derivedCount, an]

        var result = tx.execute(
                VariantQuery.build(options),
                Map.of("chr", chr, "start", start, "end", end)
        );

        try {
            while (result.hasNext()) {
                Map<String, Object> row = result.next();
                Node variant = (Node) row.get("v");

                int ac, an, hetCount, homAltCount;
                double af;

                if (useSubset) {
                    SampleSubsetComputer.SubsetStats ss =
                            SampleSubsetComputer.compute(variant, subsetIndex);
                    ac = ss.ac; an = ss.an; af = ss.af;
                    hetCount = ss.hetCount; homAltCount = ss.homAltCount;
                } else {
                    if (popIndex < 0) {
                        String[] popIds = ArrayUtil.toStringArray(variant.getProperty("pop_ids"));
                        for (int i = 0; i < popIds.length; i++) {
                            if (popIds[i].equals(pop)) { popIndex = i; break; }
                        }
                        if (popIndex < 0) throw new RuntimeException("Population not found: " + pop);
                    }
                    int[] acArr = ArrayUtil.toIntArray(variant.getProperty("ac"));
                    int[] anArr = ArrayUtil.toIntArray(variant.getProperty("an"));
                    double[] afArr = ArrayUtil.toDoubleArray(variant.getProperty("af"));
                    int[] hetArr = ArrayUtil.toIntArray(variant.getProperty("het_count"));
                    int[] homAltArr = ArrayUtil.toIntArray(variant.getProperty("hom_alt_count"));

                    ac = acArr[popIndex]; an = anArr[popIndex]; af = afArr[popIndex];
                    hetCount = hetArr[popIndex]; homAltCount = homAltArr[popIndex];
                }

                if (an < 2) continue;

                if (filter.isActive() && !filter.passes(ac, an, af, hetCount, homAltCount, variant)) {
                    continue;
                }

                if (an > maxAn) maxAn = an;

                // For unfolded SFS, use derived allele count if ancestral info is available
                int countForSfs = ac;
                if (!folded) {
                    Object ancestralObj = variant.getProperty("ancestral_allele", null);
                    if (ancestralObj != null) {
                        String ancestral = (String) ancestralObj;
                        if ("REF".equals(ancestral)) {
                            countForSfs = ac;  // ALT is derived
                        } else if ("ALT".equals(ancestral)) {
                            countForSfs = an - ac;  // REF is derived
                        }
                        nPolarized++;
                    }
                    // If no ancestral info, use raw AC (unpolarized)
                }

                acAnPairs.add(new int[]{countForSfs, an});
            }
        } finally {
            result.close();
        }

        // Build histogram
        long[] histogram;
        if (folded) {
            int foldedSize = maxAn / 2 + 1;
            histogram = new long[foldedSize];
            for (int[] pair : acAnPairs) {
                int acVal = pair[0];
                int anVal = pair[1];
                int foldedAc = Math.min(acVal, anVal - acVal);
                if (foldedAc >= 0 && foldedAc < histogram.length) {
                    histogram[foldedAc]++;
                }
            }
        } else {
            histogram = new long[maxAn + 1];
            for (int[] pair : acAnPairs) {
                int acVal = pair[0];
                if (acVal >= 0 && acVal < histogram.length) {
                    histogram[acVal]++;
                }
            }
        }

        List<Long> sfsList = new ArrayList<>(histogram.length);
        for (long count : histogram) sfsList.add(count);

        SFSResult out = new SFSResult();
        out.sfs = sfsList;
        out.n_variants = acAnPairs.size();
        out.max_ac = maxAn;
        out.n_polarized = nPolarized;

        return Stream.of(out);
    }
}
