package org.graphpop.procedures.pairwise;

import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Pure-Java unit tests for {@link BranchIbdComputer}. Exercises hand-
 * built ARGs to verify segment boundaries, TMRCA caps, and min-length
 * filters before the headline tskit-reference integration test.
 */
class BranchIbdComputerTest {

    /** Star tree of 3 samples sharing one root (MRCA tmrca=1.0). */
    private static ARG starThree(long L) {
        int nNodes = 4;
        double[] time = {0, 0, 0, 1.0};
        boolean[] isSample = {true, true, true, false};
        long[] flags = {1, 1, 1, 0};
        int[] tskitId = {0, 1, 2, 3};
        int[] edgeParent = {3, 3, 3};
        int[] edgeChild = {0, 1, 2};
        long[] edgeStart = {0, 0, 0};
        long[] edgeEnd = {L, L, L};
        int[] sampleNodes = {0, 1, 2};
        return new ARG("test", nNodes, 3, time, isSample, flags, tskitId,
                edgeParent, edgeChild, edgeStart, edgeEnd, sampleNodes, L);
    }

    /**
     * Two samples that switch MRCA at position 50:
     * <pre>
     *  [0, 50): both share node 4 (tmrca=1.0)
     *  [50, 100): both share node 5 (tmrca=2.0)
     * </pre>
     * 4 → 0; 4 → 1; 5 → 0; 5 → 1.
     */
    private static ARG twoSamplesTwoIntervals() {
        int nNodes = 4;
        double[] time = {0, 0, 1.0, 2.0};
        boolean[] isSample = {true, true, false, false};
        long[] flags = {1, 1, 0, 0};
        int[] tskitId = {0, 1, 2, 3};
        int[] edgeParent = {2, 2, 3, 3};
        int[] edgeChild = {0, 1, 0, 1};
        long[] edgeStart = {0, 0, 50, 50};
        long[] edgeEnd = {50, 50, 100, 100};
        int[] sampleNodes = {0, 1};
        return new ARG("test", nNodes, 4, time, isSample, flags, tskitId,
                edgeParent, edgeChild, edgeStart, edgeEnd, sampleNodes, 100L);
    }

    @Test
    void starTree_threeSamples_emitsThreePairsOneSegmentEach() {
        ARG arg = starThree(100L);
        List<IbdSegment> segs = BranchIbdComputer.compute(
                arg, 0L, 100L,
                BranchIbdComputer.NO_MAX_TMRCA,
                BranchIbdComputer.NO_MIN_LENGTH);
        assertEquals(3, segs.size(), "3 unordered pairs => 3 segments");
        for (IbdSegment s : segs) {
            assertEquals(0L, s.start);
            assertEquals(100L, s.end);
            assertEquals(3, s.mrcaNodeId);     // tskit id = 3 (the root)
            assertEquals(1.0, s.tmrca, 1e-12);
            assertEquals(100L, s.lengthBp());
        }
    }

    @Test
    void twoSamples_twoIntervals_emitsTwoSegmentsAtBreakpoint() {
        ARG arg = twoSamplesTwoIntervals();
        List<IbdSegment> segs = BranchIbdComputer.compute(
                arg, 0L, 100L,
                BranchIbdComputer.NO_MAX_TMRCA,
                BranchIbdComputer.NO_MIN_LENGTH);
        assertEquals(2, segs.size());
        // Sort by start.
        segs.sort((s1, s2) -> Long.compare(s1.start, s2.start));
        assertEquals(0L, segs.get(0).start);
        assertEquals(50L, segs.get(0).end);
        assertEquals(2, segs.get(0).mrcaNodeId);
        assertEquals(1.0, segs.get(0).tmrca, 1e-12);

        assertEquals(50L, segs.get(1).start);
        assertEquals(100L, segs.get(1).end);
        assertEquals(3, segs.get(1).mrcaNodeId);
        assertEquals(2.0, segs.get(1).tmrca, 1e-12);
    }

    @Test
    void tmrcaCap_dropsSegmentsAboveThreshold() {
        ARG arg = twoSamplesTwoIntervals();
        // max_tmrca = 1.5: drops the second interval (tmrca=2.0).
        List<IbdSegment> segs = BranchIbdComputer.compute(
                arg, 0L, 100L, /*maxTmrca=*/ 1.5,
                BranchIbdComputer.NO_MIN_LENGTH);
        assertEquals(1, segs.size());
        assertEquals(0L, segs.get(0).start);
        assertEquals(50L, segs.get(0).end);
        assertEquals(1.0, segs.get(0).tmrca, 1e-12);
    }

    @Test
    void minLengthFilter_dropsShortSegments() {
        ARG arg = twoSamplesTwoIntervals();
        // min_length_bp = 60 -> both 50-bp segments dropped.
        List<IbdSegment> segs = BranchIbdComputer.compute(
                arg, 0L, 100L,
                BranchIbdComputer.NO_MAX_TMRCA,
                /*minLengthBp=*/ 60L);
        assertEquals(0, segs.size());
    }

    @Test
    void minLengthFilter_keepsLongSegments() {
        ARG arg = starThree(100L);
        List<IbdSegment> segs = BranchIbdComputer.compute(
                arg, 0L, 100L,
                BranchIbdComputer.NO_MAX_TMRCA,
                /*minLengthBp=*/ 50L);
        // 100-bp segments survive.
        assertEquals(3, segs.size());
    }

    @Test
    void regionRestriction_clipsOutput() {
        ARG arg = twoSamplesTwoIntervals();
        // Region [40, 80) intersects both intervals partially.
        List<IbdSegment> segs = BranchIbdComputer.compute(
                arg, 40L, 80L,
                BranchIbdComputer.NO_MAX_TMRCA,
                BranchIbdComputer.NO_MIN_LENGTH);
        assertEquals(2, segs.size());
        segs.sort((s1, s2) -> Long.compare(s1.start, s2.start));
        assertEquals(40L, segs.get(0).start);
        assertEquals(50L, segs.get(0).end);
        assertEquals(50L, segs.get(1).start);
        assertEquals(80L, segs.get(1).end);
    }

    @Test
    void emptyArg_returnsEmpty() {
        ARG arg = new ARG("empty", 0, 0,
                new double[0], new boolean[0], new long[0], new int[0],
                new int[0], new int[0], new long[0], new long[0],
                new int[0], 100L);
        assertEquals(0, BranchIbdComputer.compute(
                arg, 0L, 100L,
                BranchIbdComputer.NO_MAX_TMRCA,
                BranchIbdComputer.NO_MIN_LENGTH).size());
    }

    @Test
    void canonicalOrdering_aLessThanB() {
        ARG arg = starThree(100L);
        List<IbdSegment> segs = BranchIbdComputer.compute(
                arg, 0L, 100L,
                BranchIbdComputer.NO_MAX_TMRCA,
                BranchIbdComputer.NO_MIN_LENGTH);
        for (IbdSegment s : segs) {
            assertTrue(s.sampleA < s.sampleB,
                "segment " + s + " not in canonical order");
        }
    }
}
