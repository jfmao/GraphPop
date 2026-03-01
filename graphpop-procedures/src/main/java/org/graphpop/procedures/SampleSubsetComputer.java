package org.graphpop.procedures;

import org.neo4j.graphdb.Direction;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.RelationshipType;

import java.util.Map;

/**
 * Recomputes population statistics on the fly for a custom sample subset.
 *
 * <p>When a user provides a {@code samples} option, FAST PATH procedures
 * cannot use the pre-computed per-population arrays (ac, an, af, het_count, etc.)
 * because those are for the full population. This class traverses CARRIES edges
 * for the specified samples and recomputes AC/AN/AF/het per variant.</p>
 *
 * <p>This is inherently slower than FAST PATH (O(V × S_subset) rather than O(V × K)),
 * but gives full flexibility for user-defined sample subsets.</p>
 */
final class SampleSubsetComputer {

    private static final RelationshipType CARRIES_REL = RelationshipType.withName("CARRIES");

    private SampleSubsetComputer() {}

    /**
     * Recompute allele statistics for a single variant from CARRIES edges.
     *
     * @param variant     the Variant node
     * @param sampleIndex mapping from sampleId to array position (only samples in this map are counted)
     * @return computed statistics for the variant
     */
    static SubsetStats compute(Node variant, Map<String, Integer> sampleIndex) {
        int nSamples = sampleIndex.size();
        int ac = 0;
        int hetCount = 0;
        int homAltCount = 0;

        for (var rel : variant.getRelationships(Direction.INCOMING, CARRIES_REL)) {
            Node sampleNode = rel.getStartNode();
            String sid = (String) sampleNode.getProperty("sampleId");
            if (!sampleIndex.containsKey(sid)) continue;

            int gt = ((Number) rel.getProperty("gt")).intValue();
            if (gt == 1) {
                ac += 1;
                hetCount++;
            } else if (gt == 2) {
                ac += 2;
                homAltCount++;
            }
        }

        int an = nSamples * 2;
        double af = an > 0 ? (double) ac / an : 0.0;

        return new SubsetStats(ac, an, af, hetCount, homAltCount);
    }

    /**
     * Container for recomputed per-variant statistics from a sample subset.
     */
    static class SubsetStats {
        final int ac;
        final int an;
        final double af;
        final int hetCount;
        final int homAltCount;

        SubsetStats(int ac, int an, double af, int hetCount, int homAltCount) {
            this.ac = ac;
            this.an = an;
            this.af = af;
            this.hetCount = hetCount;
            this.homAltCount = homAltCount;
        }
    }
}
