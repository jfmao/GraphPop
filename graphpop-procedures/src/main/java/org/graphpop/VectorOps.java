package org.graphpop;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorSpecies;

/**
 * SIMD-accelerated numeric operations for population genomics.
 *
 * <p>Uses the Java Vector API ({@code jdk.incubator.vector}) to perform
 * vectorised arithmetic on allele-frequency arrays and genotype vectors.</p>
 */
public final class VectorOps {

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;

    private VectorOps() {
        // utility class
    }

    /**
     * Compute the dot product of two double arrays using SIMD lanes.
     *
     * @param a first array
     * @param b second array (must be same length as {@code a})
     * @return dot product
     */
    public static double dotProduct(double[] a, double[] b) {
        // TODO: implement SIMD dot product
        throw new UnsupportedOperationException("not yet implemented");
    }
}
