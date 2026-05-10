package org.graphpop.procedures.pairwise;

/**
 * Identity-by-descent segment between two samples.
 *
 * <p>Emitted by {@link BranchIbdComputer} for each maximal genomic
 * interval where two haplotypes share the same most-recent common
 * ancestor (under an optional TMRCA cap). Mirrors the per-segment
 * tuples returned by {@code tskit.TreeSequence.ibd_segments(...,
 * store_pairs=True, store_segments=True)}.</p>
 */
public final class IbdSegment {

    public final int sampleA;     // packed sample index (0..n-1) -- canonical (min)
    public final int sampleB;     // packed sample index (0..n-1) -- canonical (max)
    public final long start;       // bp inclusive
    public final long end;         // bp exclusive
    public final int mrcaNodeId;  // tskit node id of the segment's MRCA
    public final double tmrca;     // generations

    public IbdSegment(int sampleA, int sampleB, long start, long end,
                      int mrcaNodeId, double tmrca) {
        // Canonical ordering: lower index first.
        if (sampleA <= sampleB) {
            this.sampleA = sampleA;
            this.sampleB = sampleB;
        } else {
            this.sampleA = sampleB;
            this.sampleB = sampleA;
        }
        this.start = start;
        this.end = end;
        this.mrcaNodeId = mrcaNodeId;
        this.tmrca = tmrca;
    }

    public long lengthBp() { return end - start; }

    @Override
    public String toString() {
        return "IbdSegment{" + sampleA + "," + sampleB
            + " [" + start + "," + end + ") mrca=" + mrcaNodeId
            + " tmrca=" + tmrca + "}";
    }
}
