package org.graphpop.procedures.pairwise;

/**
 * Result row for {@code graphpop.ibd.from_arg}: one row per emitted
 * {@code :IBD_SEGMENT}. The procedure also writes the same edges
 * into Neo4j for downstream querying.
 */
public class IbdSegmentResult {

    public String sample_a;
    public String sample_b;
    public String chr;
    public long start;
    public long end;
    public long length_bp;
    public long mrca_node_id;
    public double tmrca;
    public String source;
    public String runId;

    public IbdSegmentResult() {}

    public IbdSegmentResult(String a, String b, String chr,
                             long start, long end, long lengthBp,
                             long mrcaNodeId, double tmrca,
                             String source, String runId) {
        this.sample_a = a;
        this.sample_b = b;
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.length_bp = lengthBp;
        this.mrca_node_id = mrcaNodeId;
        this.tmrca = tmrca;
        this.source = source;
        this.runId = runId;
    }
}
