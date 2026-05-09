package org.graphpop.procedures.pairwise;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class WelfordMatrixTest {

    private static final double EPS = 1e-12;

    @Test
    void singleUpdate_meanEqualsX_seIsNaN() {
        WelfordMatrix w = new WelfordMatrix(2);
        w.update(new double[][]{{1.0, 2.0}, {3.0, 4.0}});

        assertEquals(1, w.nSeen());
        double[][] m = w.mean();
        assertArrayEquals(new double[]{1.0, 2.0}, m[0], EPS);
        assertArrayEquals(new double[]{3.0, 4.0}, m[1], EPS);

        double[][] se = w.stderr();
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                assertTrue(Double.isNaN(se[i][j]),
                    "se[" + i + "][" + j + "] should be NaN");
    }

    @Test
    void twoUpdates_meanIsAverage_seFormula() {
        WelfordMatrix w = new WelfordMatrix(1);
        w.update(new double[][]{{2.0}});
        w.update(new double[][]{{4.0}});

        // mean = 3
        assertEquals(3.0, w.mean()[0][0], EPS);
        // sd (ddof=1) = sqrt(((2-3)^2 + (4-3)^2) / 1) = sqrt(2)
        assertEquals(Math.sqrt(2.0), w.sampleStd()[0][0], EPS);
        // SE = sd / sqrt(2)
        assertEquals(Math.sqrt(2.0) / Math.sqrt(2.0), w.stderr()[0][0], EPS);
        // Equivalently: SE = sqrt(M2 / (n*(n-1))) = sqrt(2 / 2) = 1
        assertEquals(1.0, w.stderr()[0][0], EPS);
    }

    @Test
    void identicalUpdates_seIsZero() {
        WelfordMatrix w = new WelfordMatrix(2);
        double[][] x = {{1.5, -0.5}, {-0.5, 1.5}};
        for (int i = 0; i < 5; i++) w.update(x);
        double[][] m = w.mean();
        double[][] se = w.stderr();
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                assertEquals(x[i][j], m[i][j], EPS);
                assertEquals(0.0, se[i][j], 1e-15);
            }
        }
    }

    @Test
    void tenUpdates_matchesNumpyDdofOne() {
        // numpy.array([1, 2, ..., 10]).std(ddof=1) = sqrt(55/3) for the
        // sequence [1..10] with n=10, mean=5.5, var = 110/9? Let me
        // compute directly: var = sum((x-5.5)^2) / 9.
        double[] xs = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        double mean = 5.5;
        double sumSq = 0;
        for (double x : xs) sumSq += (x - mean) * (x - mean);
        double expectedSd = Math.sqrt(sumSq / 9.0);
        double expectedSe = expectedSd / Math.sqrt(10.0);

        WelfordMatrix w = new WelfordMatrix(1);
        for (double x : xs) w.update(new double[][]{{x}});

        assertEquals(mean, w.mean()[0][0], 1e-12);
        assertEquals(expectedSd, w.sampleStd()[0][0], 1e-12);
        assertEquals(expectedSe, w.stderr()[0][0], 1e-12);
    }

    @Test
    void rejectsMismatchedRowCount() {
        WelfordMatrix w = new WelfordMatrix(3);
        assertThrows(IllegalArgumentException.class,
            () -> w.update(new double[][]{{1.0}, {2.0}}));
    }

    @Test
    void rejectsMismatchedColumnCount() {
        WelfordMatrix w = new WelfordMatrix(2);
        assertThrows(IllegalArgumentException.class,
            () -> w.update(new double[][]{{1.0, 2.0}, {3.0}}));
    }

    @Test
    void preservesSymmetry() {
        WelfordMatrix w = new WelfordMatrix(3);
        double[][] x1 = {{1, 2, 3}, {2, 4, 5}, {3, 5, 6}};
        double[][] x2 = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
        w.update(x1);
        w.update(x2);
        double[][] m = w.mean();
        double[][] se = w.stderr();
        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 3; j++) {
                assertEquals(m[i][j], m[j][i], EPS);
                assertEquals(se[i][j], se[j][i], EPS);
            }
        }
    }
}
