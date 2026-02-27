package org.graphpop;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for {@link VectorOps}.
 */
class VectorOpsTest {

    @Test
    void dotProductPlaceholder() {
        // TODO: replace with real test once dotProduct is implemented
        assertThrows(UnsupportedOperationException.class,
                () -> VectorOps.dotProduct(new double[]{1.0}, new double[]{2.0}));
    }
}
