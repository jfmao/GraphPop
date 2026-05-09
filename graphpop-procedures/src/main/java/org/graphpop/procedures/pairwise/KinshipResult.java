package org.graphpop.procedures.pairwise;

/**
 * Result row emitted by pairwise kinship procedures.
 *
 * <p>Fields are public to satisfy Neo4j's stored-procedure POJO contract.</p>
 */
public class KinshipResult {

    public String sample_a;
    public String sample_b;
    public double phi;
    public long ibs0;
    public long het_het;
    public long n_snp;
    public long n_aa_min;
    public String method;

    public KinshipResult() {}

    public KinshipResult(String a, String b, double phi, long ibs0, long hetHet,
                         long nSnp, long nAaMin, String method) {
        this.sample_a = a;
        this.sample_b = b;
        this.phi = phi;
        this.ibs0 = ibs0;
        this.het_het = hetHet;
        this.n_snp = nSnp;
        this.n_aa_min = nAaMin;
        this.method = method;
    }
}
