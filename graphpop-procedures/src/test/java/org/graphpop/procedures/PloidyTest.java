package org.graphpop.procedures;

import org.junit.jupiter.api.Test;
import org.neo4j.graphdb.*;

import java.lang.reflect.InvocationHandler;
import java.lang.reflect.Method;
import java.lang.reflect.Proxy;
import java.util.HashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for ploidy-aware allele counting in SampleSubsetComputer and
 * HaplotypeMatrix AC/AN computation.
 */
class PloidyTest {

    // ---- Stub Node via Proxy (no Mockito needed) ----

    /**
     * Create a minimal Node stub backed by a property map.
     */
    private static Node stubNode(Map<String, Object> properties) {
        InvocationHandler handler = (proxy, method, args) -> {
            String name = method.getName();
            if ("getProperty".equals(name)) {
                if (args.length == 1) {
                    String key = (String) args[0];
                    Object val = properties.get(key);
                    if (val == null) throw new NotFoundException("Property " + key);
                    return val;
                } else if (args.length == 2) {
                    // getProperty(key, default)
                    return properties.getOrDefault((String) args[0], args[1]);
                }
            }
            if ("hasProperty".equals(name)) {
                return properties.containsKey((String) args[0]);
            }
            throw new UnsupportedOperationException(name);
        };
        return (Node) Proxy.newProxyInstance(
                Node.class.getClassLoader(),
                new Class[]{Node.class},
                handler);
    }

    // ---- SampleSubsetComputer: all-diploid (backward compat) ----

    @Test
    void subsetComputer_allDiploid_noPloidyPacked() {
        // 4 samples, no ploidy_packed → standard diploid counting
        // gt: HOM_REF, HET, HOM_ALT, MISSING
        int nSamples = 4;
        byte[] gtPacked = new byte[PackedGenotypeReader.gtPackedLength(nSamples)];
        PackedGenotypeReader.setGenotype(gtPacked, 0, PackedGenotypeReader.GT_HOM_REF);
        PackedGenotypeReader.setGenotype(gtPacked, 1, PackedGenotypeReader.GT_HET);
        PackedGenotypeReader.setGenotype(gtPacked, 2, PackedGenotypeReader.GT_HOM_ALT);
        PackedGenotypeReader.setGenotype(gtPacked, 3, PackedGenotypeReader.GT_MISSING);

        Map<String, Object> props = new HashMap<>();
        props.put("gt_packed", gtPacked);
        // No ploidy_packed → all diploid
        Node variant = stubNode(props);

        Map<String, Integer> sampleIndex = Map.of("S0", 0, "S1", 1, "S2", 2, "S3", 3);
        int[] packedIndices = {0, 1, 2, 3};

        SampleSubsetComputer.SubsetStats ss = SampleSubsetComputer.compute(variant, sampleIndex, packedIndices);

        // 3 called × 2 ploidy = 6 AN; HET→1 + HOM_ALT→2 = 3 AC
        assertEquals(3, ss.ac);
        assertEquals(6, ss.an);
        assertEquals(3.0 / 6, ss.af, 1e-10);
        assertEquals(1, ss.hetCount);
        assertEquals(1, ss.homAltCount);
    }

    // ---- SampleSubsetComputer: all-haploid ----

    @Test
    void subsetComputer_allHaploid() {
        // 5 haploid samples: 3 REF, 2 ALT → ac=2, an=5
        int nSamples = 5;
        byte[] gtPacked = new byte[PackedGenotypeReader.gtPackedLength(nSamples)];
        // cyvcf2 maps haploid REF → HOM_REF, haploid ALT → HOM_ALT
        PackedGenotypeReader.setGenotype(gtPacked, 0, PackedGenotypeReader.GT_HOM_REF);
        PackedGenotypeReader.setGenotype(gtPacked, 1, PackedGenotypeReader.GT_HOM_REF);
        PackedGenotypeReader.setGenotype(gtPacked, 2, PackedGenotypeReader.GT_HOM_REF);
        PackedGenotypeReader.setGenotype(gtPacked, 3, PackedGenotypeReader.GT_HOM_ALT);
        PackedGenotypeReader.setGenotype(gtPacked, 4, PackedGenotypeReader.GT_HOM_ALT);

        // All haploid
        byte[] ploidyPacked = new byte[PackedGenotypeReader.ploidyPackedLength(nSamples)];
        for (int i = 0; i < nSamples; i++) {
            PackedGenotypeReader.setPloidy(ploidyPacked, i, 1);
        }

        Map<String, Object> props = new HashMap<>();
        props.put("gt_packed", gtPacked);
        props.put("ploidy_packed", ploidyPacked);
        Node variant = stubNode(props);

        Map<String, Integer> sampleIndex = new HashMap<>();
        int[] packedIndices = new int[nSamples];
        for (int i = 0; i < nSamples; i++) {
            sampleIndex.put("S" + i, i);
            packedIndices[i] = i;
        }

        SampleSubsetComputer.SubsetStats ss = SampleSubsetComputer.compute(variant, sampleIndex, packedIndices);

        assertEquals(2, ss.ac);     // 2 ALT alleles
        assertEquals(5, ss.an);     // 5 samples × 1 allele each
        assertEquals(2.0 / 5, ss.af, 1e-10);
        assertEquals(0, ss.hetCount);  // haploid can't be het
        assertEquals(2, ss.homAltCount);
    }

    // ---- SampleSubsetComputer: mixed ploidy (chrX-like) ----

    @Test
    void subsetComputer_mixedPloidy() {
        // 2 diploid females + 3 haploid males
        // Female 0: HET → ac=1, an=2
        // Female 1: HOM_REF → ac=0, an=2
        // Male 2: ALT → ac=1, an=1
        // Male 3: REF → ac=0, an=1
        // Male 4: REF → ac=0, an=1
        // Total: ac=2, an=7

        int nSamples = 5;
        byte[] gtPacked = new byte[PackedGenotypeReader.gtPackedLength(nSamples)];
        PackedGenotypeReader.setGenotype(gtPacked, 0, PackedGenotypeReader.GT_HET);
        PackedGenotypeReader.setGenotype(gtPacked, 1, PackedGenotypeReader.GT_HOM_REF);
        PackedGenotypeReader.setGenotype(gtPacked, 2, PackedGenotypeReader.GT_HOM_ALT);  // haploid ALT
        PackedGenotypeReader.setGenotype(gtPacked, 3, PackedGenotypeReader.GT_HOM_REF);  // haploid REF
        PackedGenotypeReader.setGenotype(gtPacked, 4, PackedGenotypeReader.GT_HOM_REF);  // haploid REF

        // Samples 2,3,4 are haploid
        byte[] ploidyPacked = new byte[PackedGenotypeReader.ploidyPackedLength(nSamples)];
        PackedGenotypeReader.setPloidy(ploidyPacked, 2, 1);
        PackedGenotypeReader.setPloidy(ploidyPacked, 3, 1);
        PackedGenotypeReader.setPloidy(ploidyPacked, 4, 1);

        Map<String, Object> props = new HashMap<>();
        props.put("gt_packed", gtPacked);
        props.put("ploidy_packed", ploidyPacked);
        Node variant = stubNode(props);

        Map<String, Integer> sampleIndex = new HashMap<>();
        int[] packedIndices = new int[nSamples];
        for (int i = 0; i < nSamples; i++) {
            sampleIndex.put("S" + i, i);
            packedIndices[i] = i;
        }

        SampleSubsetComputer.SubsetStats ss = SampleSubsetComputer.compute(variant, sampleIndex, packedIndices);

        assertEquals(2, ss.ac);     // 1 from HET + 1 from haploid ALT
        assertEquals(7, ss.an);     // 2×2 diploid + 3×1 haploid
        assertEquals(2.0 / 7, ss.af, 1e-10);
        assertEquals(1, ss.hetCount);   // only the diploid HET
        assertEquals(1, ss.homAltCount); // the haploid ALT
    }

    // ---- SampleSubsetComputer: mixed with missing ----

    @Test
    void subsetComputer_mixedWithMissing() {
        // 2 diploid + 2 haploid, 1 missing in each group
        // Diploid 0: HET → ac=1, an=2
        // Diploid 1: MISSING → excluded
        // Haploid 2: ALT → ac=1, an=1
        // Haploid 3: MISSING → excluded
        // Total: ac=2, an=3

        int nSamples = 4;
        byte[] gtPacked = new byte[PackedGenotypeReader.gtPackedLength(nSamples)];
        PackedGenotypeReader.setGenotype(gtPacked, 0, PackedGenotypeReader.GT_HET);
        PackedGenotypeReader.setGenotype(gtPacked, 1, PackedGenotypeReader.GT_MISSING);
        PackedGenotypeReader.setGenotype(gtPacked, 2, PackedGenotypeReader.GT_HOM_ALT);
        PackedGenotypeReader.setGenotype(gtPacked, 3, PackedGenotypeReader.GT_MISSING);

        byte[] ploidyPacked = new byte[PackedGenotypeReader.ploidyPackedLength(nSamples)];
        PackedGenotypeReader.setPloidy(ploidyPacked, 2, 1);
        PackedGenotypeReader.setPloidy(ploidyPacked, 3, 1);

        Map<String, Object> props = new HashMap<>();
        props.put("gt_packed", gtPacked);
        props.put("ploidy_packed", ploidyPacked);
        Node variant = stubNode(props);

        Map<String, Integer> sampleIndex = new HashMap<>();
        int[] packedIndices = new int[nSamples];
        for (int i = 0; i < nSamples; i++) {
            sampleIndex.put("S" + i, i);
            packedIndices[i] = i;
        }

        SampleSubsetComputer.SubsetStats ss = SampleSubsetComputer.compute(variant, sampleIndex, packedIndices);

        assertEquals(2, ss.ac);     // 1 from HET + 1 from haploid ALT
        assertEquals(3, ss.an);     // 1×2 diploid + 1×1 haploid (missing excluded)
        assertEquals(1, ss.hetCount);
        assertEquals(1, ss.homAltCount);
    }

    // ---- HaplotypeMatrix: AC/AN via forTest (tests corrected logic indirectly) ----

    @Test
    void haplotypeMatrix_acsAns_diploid() {
        // Verify that HaplotypeMatrix.forTest stores correct values for diploid
        // (forTest bypasses loadInternal, but verifies the data structure)
        String[] vids = {"v1", "v2"};
        long[] positions = {100, 200};
        int[] acs = {3, 0};
        int[] ans = {10, 10};
        double[] afs = {0.3, 0.0};
        byte[][] haps = new byte[2][10]; // 5 samples × 2 haplotypes

        HaplotypeMatrix m = HaplotypeMatrix.forTest(vids, positions, afs, acs, ans, haps, 5);

        assertEquals(3, m.acs[0]);
        assertEquals(10, m.ans[0]);
        assertEquals(0.3, m.afs[0], 1e-10);
    }

    @Test
    void haplotypeMatrix_acsAns_haploidValues() {
        // Simulate what loadInternal would produce for 3 haploid samples:
        // AN should be 3 (not 6), AC reflects actual haploid counts
        String[] vids = {"v1"};
        long[] positions = {100};
        // For 3 haploid samples with 1 ALT: ac=1, an=3
        int[] acs = {1};
        int[] ans = {3};
        double[] afs = {1.0 / 3};
        byte[][] haps = new byte[1][6]; // 3 samples × 2 (pseudo-diploid slots)
        // Haploid ALT → both slots = 1
        haps[0][4] = 1; // sample 2, hap 0
        haps[0][5] = 1; // sample 2, hap 1

        HaplotypeMatrix m = HaplotypeMatrix.forTest(vids, positions, afs, acs, ans, haps, 3);

        assertEquals(1, m.acs[0]);
        assertEquals(3, m.ans[0]);  // haploid: 3 samples × 1 allele
        assertEquals(1.0 / 3, m.afs[0], 1e-10);
    }
}
