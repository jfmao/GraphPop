package org.graphpop.procedures;

import org.graphpop.VectorOps;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for LD-related VectorOps methods: pearsonR2, dPrime, and dPrimePacked.
 *
 * <p>These test the static helper methods directly without Neo4j, since
 * LDProcedure's DB interaction (HaplotypeMatrix loading) is hard to mock
 * in a unit test.</p>
 */
class LDHelperTest {

    private static final double EPS = 1e-8;

    // -----------------------------------------------------------------------
    // pearsonR2 tests
    // -----------------------------------------------------------------------

    @Test
    void pearsonR2_perfectCorrelation() {
        // Identical dosage vectors -> r^2 = 1.0
        double[] geno = {0, 1, 2, 0, 1, 2, 0, 1};
        assertEquals(1.0, VectorOps.pearsonR2(geno, geno), EPS);
    }

    @Test
    void pearsonR2_perfectNegativeCorrelation() {
        // Perfectly anti-correlated dosages -> r^2 = 1.0 (r=-1, squared)
        double[] geno1 = {0, 0, 2, 2, 0, 0, 2, 2};
        double[] geno2 = {2, 2, 0, 0, 2, 2, 0, 0};
        assertEquals(1.0, VectorOps.pearsonR2(geno1, geno2), EPS);
    }

    @Test
    void pearsonR2_uncorrelated() {
        // Orthogonal pattern -> r^2 near 0
        double[] geno1 = {0, 2, 0, 2, 0, 2, 0, 2};
        double[] geno2 = {0, 0, 2, 2, 0, 0, 2, 2};
        assertEquals(0.0, VectorOps.pearsonR2(geno1, geno2), EPS);
    }

    @Test
    void pearsonR2_monomorphicReturnsZero() {
        // One variant is monomorphic (variance = 0) -> r^2 = 0
        double[] geno1 = {1, 1, 1, 1, 1, 1, 1, 1};
        double[] geno2 = {0, 1, 2, 0, 1, 2, 0, 1};
        assertEquals(0.0, VectorOps.pearsonR2(geno1, geno2), EPS);
    }

    @Test
    void pearsonR2_knownValue() {
        // Hand-computed example:
        // geno1 = [0, 1, 0, 1, 2], mean1 = 4/5 = 0.8
        // geno2 = [0, 0, 1, 1, 2], mean2 = 4/5 = 0.8
        // cov = sum((x-0.8)*(y-0.8)) = (-0.8)(-0.8) + (0.2)(-0.8) + (-0.8)(0.2) + (0.2)(0.2) + (1.2)(1.2)
        //     = 0.64 - 0.16 - 0.16 + 0.04 + 1.44 = 1.80
        // var1 = sum((x-0.8)^2) = 0.64 + 0.04 + 0.64 + 0.04 + 1.44 = 2.80
        // var2 = sum((y-0.8)^2) = 0.64 + 0.64 + 0.04 + 0.04 + 1.44 = 2.80
        // r = 1.80 / sqrt(2.80*2.80) = 1.80/2.80 = 0.642857...
        // r^2 = 0.413265...
        double[] geno1 = {0, 1, 0, 1, 2};
        double[] geno2 = {0, 0, 1, 1, 2};
        double expected = (1.80 / 2.80) * (1.80 / 2.80);
        assertEquals(expected, VectorOps.pearsonR2(geno1, geno2), EPS);
    }

    @Test
    void pearsonR2_emptyReturnsZero() {
        assertEquals(0.0, VectorOps.pearsonR2(new double[0], new double[0]), EPS);
    }

    // -----------------------------------------------------------------------
    // dPrime (int[] version) tests
    // -----------------------------------------------------------------------

    @Test
    void dPrime_perfectLD() {
        // All haplotypes: both loci have same allele -> D' = 1
        int[] hap1 = {0, 0, 0, 1, 1, 1, 1, 1};
        int[] hap2 = {0, 0, 0, 1, 1, 1, 1, 1};
        assertEquals(1.0, VectorOps.dPrime(hap1, hap2, 8), EPS);
    }

    @Test
    void dPrime_noLD() {
        // Perfect independence: 2x2 table is balanced
        // hap1: 4 ref, 4 alt; hap2: 4 ref, 4 alt
        // AB=2, Ab=2, aB=2, ab=2 -> D = 2/8 - (4/8)(4/8) = 0.25 - 0.25 = 0
        int[] hap1 = {0, 0, 1, 1, 0, 0, 1, 1};
        int[] hap2 = {0, 1, 0, 1, 0, 1, 0, 1};
        assertEquals(0.0, VectorOps.dPrime(hap1, hap2, 8), EPS);
    }

    @Test
    void dPrime_monomorphicReturnsZero() {
        int[] hap1 = {0, 0, 0, 0};
        int[] hap2 = {0, 1, 0, 1};
        assertEquals(0.0, VectorOps.dPrime(hap1, hap2, 4), EPS);
    }

    @Test
    void dPrime_knownValue() {
        // 8 haplotypes:
        // hap1: [1,1,1,0,0,0,0,0] -> pA = 3/8
        // hap2: [1,1,0,0,0,0,0,0] -> pB = 2/8 = 1/4
        // n11=2 (both alt at positions 0,1), n10=1, n01=0
        // pAB = 2/8 = 0.25
        // D = pAB - pA*pB = 0.25 - (3/8)(1/4) = 0.25 - 0.09375 = 0.15625
        // D >= 0: Dmax = min(pA*(1-pB), (1-pA)*pB) = min(3/8 * 3/4, 5/8 * 1/4)
        //       = min(0.28125, 0.15625) = 0.15625
        // D' = |0.15625| / 0.15625 = 1.0
        int[] hap1 = {1, 1, 1, 0, 0, 0, 0, 0};
        int[] hap2 = {1, 1, 0, 0, 0, 0, 0, 0};
        assertEquals(1.0, VectorOps.dPrime(hap1, hap2, 8), EPS);
    }

    @Test
    void dPrime_partialLD() {
        // 8 haplotypes with partial LD
        // hap1: [1,1,1,1,0,0,0,0] -> pA = 0.5
        // hap2: [1,1,0,0,1,0,0,0] -> pB = 3/8
        // n11=2, n10=2, n01=1, n00=3
        // pAB = 2/8 = 0.25
        // D = 0.25 - 0.5*0.375 = 0.25 - 0.1875 = 0.0625
        // Dmax = min(0.5*0.625, 0.5*0.375) = min(0.3125, 0.1875) = 0.1875
        // D' = 0.0625 / 0.1875 = 1/3
        int[] hap1 = {1, 1, 1, 1, 0, 0, 0, 0};
        int[] hap2 = {1, 1, 0, 0, 1, 0, 0, 0};
        assertEquals(1.0 / 3.0, VectorOps.dPrime(hap1, hap2, 8), EPS);
    }

    // -----------------------------------------------------------------------
    // dPrime (byte[] version) tests
    // -----------------------------------------------------------------------

    @Test
    void dPrime_byteVersion_matchesIntVersion() {
        int[] intHap1 = {1, 0, 1, 0, 1, 1, 0, 0, 1, 0};
        int[] intHap2 = {1, 1, 0, 0, 1, 0, 1, 0, 1, 0};
        byte[] byteHap1 = {1, 0, 1, 0, 1, 1, 0, 0, 1, 0};
        byte[] byteHap2 = {1, 1, 0, 0, 1, 0, 1, 0, 1, 0};

        double intResult = VectorOps.dPrime(intHap1, intHap2, 10);
        double byteResult = VectorOps.dPrime(byteHap1, byteHap2, 10);
        assertEquals(intResult, byteResult, EPS,
            "byte[] and int[] dPrime should produce identical results");
    }

    // -----------------------------------------------------------------------
    // dPrimePacked (bit-packed) tests
    // -----------------------------------------------------------------------

    @Test
    void dPrimePacked_perfectLD() {
        // 8 haplotypes packed into 1 byte each: bits 0-7
        // hap1: 0b00011111 = haps 0-4 are alt (5 alt, 3 ref)
        // hap2: 0b00011111 = same pattern -> D' = 1
        byte[] packed1 = {(byte) 0b00011111};
        byte[] packed2 = {(byte) 0b00011111};
        assertEquals(1.0, VectorOps.dPrimePacked(packed1, packed2, 8), EPS);
    }

    @Test
    void dPrimePacked_noLD() {
        // 8 haplotypes: balanced independence
        // hap1: 0b00001111 -> haps 0-3 alt, 4-7 ref (pA=0.5)
        // hap2: 0b01010101 -> haps 0,2,4,6 alt (pB=0.5)
        // n11 = bitcount(0b00001111 & 0b01010101) = bitcount(0b00000101) = 2
        // pAB = 2/8 = 0.25, pA*pB = 0.25 -> D=0 -> D'=0
        byte[] packed1 = {(byte) 0b00001111};
        byte[] packed2 = {(byte) 0b01010101};
        assertEquals(0.0, VectorOps.dPrimePacked(packed1, packed2, 8), EPS);
    }

    @Test
    void dPrimePacked_matchesUnpackedDPrime() {
        // Cross-validate packed vs unpacked on the same data
        // 10 haplotypes: 2 bytes (8 + 2 bits)
        int[] intHap1 = {1, 0, 1, 0, 1, 1, 0, 0, 1, 0};
        int[] intHap2 = {1, 1, 0, 0, 1, 0, 1, 0, 1, 0};

        // Pack into bytes: byte[0] = bits 0-7, byte[1] = bits 8-9
        // hap1: bits 0-7 = 10101100 = 0b00110101, bits 8-9 = 10 = 0b01
        byte[] packed1 = new byte[2];
        byte[] packed2 = new byte[2];
        for (int i = 0; i < 10; i++) {
            if (intHap1[i] == 1) packed1[i >> 3] |= (byte) (1 << (i & 7));
            if (intHap2[i] == 1) packed2[i >> 3] |= (byte) (1 << (i & 7));
        }

        double unpackedResult = VectorOps.dPrime(intHap1, intHap2, 10);
        double packedResult = VectorOps.dPrimePacked(packed1, packed2, 10);
        assertEquals(unpackedResult, packedResult, EPS,
            "Packed and unpacked dPrime should produce identical results");
    }

    @Test
    void dPrimePacked_multipleBytes() {
        // 16 haplotypes: 2 full bytes
        // hap1: all alt in first byte, all ref in second
        // hap2: alternating
        byte[] packed1 = {(byte) 0xFF, (byte) 0x00}; // 8 alt + 8 ref -> pA = 0.5
        byte[] packed2 = {(byte) 0x55, (byte) 0x55}; // alternating 01010101 -> pB = 0.5
        // n11 = bitcount(0xFF & 0x55) + bitcount(0x00 & 0x55) = bitcount(0x55) + 0 = 4
        // pAB = 4/16 = 0.25, pA*pB = 0.25 -> D = 0 -> D' = 0
        assertEquals(0.0, VectorOps.dPrimePacked(packed1, packed2, 16), EPS);
    }

    @Test
    void dPrimePacked_monomorphicReturnsZero() {
        byte[] packed1 = {(byte) 0x00}; // all ref
        byte[] packed2 = {(byte) 0x0F}; // 4 alt, 4 ref
        assertEquals(0.0, VectorOps.dPrimePacked(packed1, packed2, 8), EPS);
    }

    @Test
    void dPrimePacked_emptyReturnsZero() {
        assertEquals(0.0, VectorOps.dPrimePacked(new byte[0], new byte[0], 0), EPS);
    }

    // -----------------------------------------------------------------------
    // r^2 and D' relationship sanity check
    // -----------------------------------------------------------------------

    @Test
    void r2AndDPrime_perfectLDImpliesBothOne() {
        // When two biallelic loci are in perfect LD, both r^2 and |D'| should be 1
        // (assuming equal allele frequencies)
        // 10 haplotypes: 5 carry alt at both, 5 carry ref at both
        double[] dosage1 = {0, 0, 0, 0, 0, 2, 2, 2, 2, 2};
        double[] dosage2 = {0, 0, 0, 0, 0, 2, 2, 2, 2, 2};
        int[] hap1 = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
        int[] hap2 = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};

        assertEquals(1.0, VectorOps.pearsonR2(dosage1, dosage2), EPS);
        assertEquals(1.0, VectorOps.dPrime(hap1, hap2, 10), EPS);
    }
}
