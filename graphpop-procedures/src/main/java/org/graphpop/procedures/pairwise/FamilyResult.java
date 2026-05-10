package org.graphpop.procedures.pairwise;

/**
 * Per-sample family assignment row for {@code graphpop.relate.families}.
 */
public class FamilyResult {

    public String sample_id;
    public String family_id;
    public long family_size;
    public String method;

    public FamilyResult() {}

    public FamilyResult(String sampleId, String familyId,
                        long familySize, String method) {
        this.sample_id = sampleId;
        this.family_id = familyId;
        this.family_size = familySize;
        this.method = method;
    }
}
