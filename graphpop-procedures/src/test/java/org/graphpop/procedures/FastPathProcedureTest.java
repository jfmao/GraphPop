package org.graphpop.procedures;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.neo4j.driver.Driver;
import org.neo4j.driver.GraphDatabase;
import org.neo4j.driver.Record;
import org.neo4j.driver.Result;
import org.neo4j.driver.Session;
import org.neo4j.harness.Neo4j;
import org.neo4j.harness.Neo4jBuilders;

import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Integration tests for the fast-path population genomics procedures
 * (Diversity, Divergence, SFS, JointSFS) using an embedded Neo4j instance.
 *
 * <p>Test data: 5 Variant nodes on chr1 (positions 100-500) with two populations
 * (POP1: 5 diploid samples / an=10, POP2: 4 diploid samples / an=8).
 * All expected values are hand-computed and validated against formulas.</p>
 */
class FastPathProcedureTest {

    private static final double EPS = 1e-4;

    private static Neo4j embeddedNeo4j;
    private static Driver driver;

    @BeforeAll
    static void setUp() {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(DiversityProcedure.class)
                .withProcedure(DivergenceProcedure.class)
                .withProcedure(SFSProcedure.class)
                .withProcedure(JointSFSProcedure.class)
                .build();

        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        try (Session session = driver.session()) {
            // Create Population nodes (required by PopulationContext.resolve)
            // POP1: 5 diploid samples, n=10 sequences
            //   a_n = H(9) = 2.828968254..., a_n2 = H2(9) = 1.539767731...
            session.run(
                "CREATE (:Population {populationId: 'POP1', n_samples: 5, " +
                "a_n: 2.8289682539682538, a_n2: 1.5397677311665408})"
            );
            // POP2: 4 diploid samples, n=8 sequences
            //   a_n = H(7) = 2.592857142..., a_n2 = H2(7) = 1.511797052...
            session.run(
                "CREATE (:Population {populationId: 'POP2', n_samples: 4, " +
                "a_n: 2.5928571428571425, a_n2: 1.5117970521541951})"
            );

            // Create 5 Variant nodes with known allele count arrays
            session.run(
                "CREATE (:Variant {variantId: 'chr1:100:A:T', chr: 'chr1', pos: 100, " +
                "ref: 'A', alt: 'T', " +
                "pop_ids: ['POP1','POP2'], ac: [2,1], an: [10,8], af: [0.2, 0.125], " +
                "het_count: [2,1], hom_alt_count: [0,0], ancestral_allele: 'REF'})"
            );
            session.run(
                "CREATE (:Variant {variantId: 'chr1:200:C:G', chr: 'chr1', pos: 200, " +
                "ref: 'C', alt: 'G', " +
                "pop_ids: ['POP1','POP2'], ac: [4,3], an: [10,8], af: [0.4, 0.375], " +
                "het_count: [2,1], hom_alt_count: [1,1], ancestral_allele: 'REF'})"
            );
            session.run(
                "CREATE (:Variant {variantId: 'chr1:300:G:A', chr: 'chr1', pos: 300, " +
                "ref: 'G', alt: 'A', " +
                "pop_ids: ['POP1','POP2'], ac: [6,5], an: [10,8], af: [0.6, 0.625], " +
                "het_count: [2,1], hom_alt_count: [2,2], ancestral_allele: 'ALT'})"
            );
            session.run(
                "CREATE (:Variant {variantId: 'chr1:400:T:C', chr: 'chr1', pos: 400, " +
                "ref: 'T', alt: 'C', " +
                "pop_ids: ['POP1','POP2'], ac: [8,2], an: [10,8], af: [0.8, 0.25], " +
                "het_count: [2,2], hom_alt_count: [3,0], ancestral_allele: 'REF'})"
            );
            session.run(
                "CREATE (:Variant {variantId: 'chr1:500:A:G', chr: 'chr1', pos: 500, " +
                "ref: 'A', alt: 'G', " +
                "pop_ids: ['POP1','POP2'], ac: [0,0], an: [10,8], af: [0.0, 0.0], " +
                "het_count: [0,0], hom_alt_count: [0,0], ancestral_allele: 'REF'})"
            );
        }
    }

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    // -----------------------------------------------------------------------
    // Diversity tests — POP1
    // -----------------------------------------------------------------------

    @Test
    void diversityPop1_pi() {
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.diversity('chr1', 1, 600, 'POP1') " +
                "YIELD pi, theta_w, tajima_d, het_exp, het_obs, fis, n_variants, n_segregating"
            ).single();

            // pi = sum(2*af*(1-af)*an/(an-1)) / L = 1.77778 / 5 = 0.35556
            assertEquals(0.3556, rec.get("pi").asDouble(), EPS);
            assertEquals(5L, rec.get("n_variants").asLong());
            assertEquals(4L, rec.get("n_segregating").asLong());
        }
    }

    @Test
    void diversityPop1_thetaW() {
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.diversity('chr1', 1, 600, 'POP1') " +
                "YIELD theta_w"
            ).single();

            // theta_w = S / a_n / L = 4 / 2.82897 / 5 = 0.28279
            assertEquals(0.2828, rec.get("theta_w").asDouble(), EPS);
        }
    }

    @Test
    void diversityPop1_tajimaD() {
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.diversity('chr1', 1, 600, 'POP1') " +
                "YIELD tajima_d"
            ).single();

            // Tajima's D computed from piTotal=1.77778, S=4, n=10
            assertEquals(0.9879, rec.get("tajima_d").asDouble(), EPS);
        }
    }

    @Test
    void diversityPop1_heterozygosity() {
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.diversity('chr1', 1, 600, 'POP1') " +
                "YIELD het_exp, het_obs, fis"
            ).single();

            // het_exp = mean(2*af*(1-af)) = 0.32
            assertEquals(0.32, rec.get("het_exp").asDouble(), EPS);
            // het_obs = mean(het_count / nDiploid) = mean(2/5, 2/5, 2/5, 2/5, 0/5) = 0.32
            assertEquals(0.32, rec.get("het_obs").asDouble(), EPS);
            // fis = 1 - het_obs/het_exp = 0
            assertEquals(0.0, rec.get("fis").asDouble(), EPS);
        }
    }

    @Test
    void diversityPop1_fayWuH() {
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.diversity('chr1', 1, 600, 'POP1') " +
                "YIELD fay_wu_h, fay_wu_h_norm, n_polarized"
            ).single();

            assertEquals(5L, rec.get("n_polarized").asLong());
            // H = (piPolarized - thetaH) / nPolarized = (1.77778 - 2.22222) / 5 = -0.08889
            assertEquals(-0.0889, rec.get("fay_wu_h").asDouble(), EPS);
            // Normalized H
            assertEquals(-0.9364, rec.get("fay_wu_h_norm").asDouble(), EPS);
        }
    }

    @Test
    void diversityPop1_subregion() {
        // Only variants at positions 100-300 (3 variants, 3 segregating)
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.diversity('chr1', 100, 300, 'POP1') " +
                "YIELD n_variants, n_segregating"
            ).single();

            assertEquals(3L, rec.get("n_variants").asLong());
            assertEquals(3L, rec.get("n_segregating").asLong());
        }
    }

    @Test
    void diversityPop1_emptyRegion() {
        // No variants in region 700-800
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.diversity('chr1', 700, 800, 'POP1') " +
                "YIELD n_variants, pi"
            ).single();

            assertEquals(0L, rec.get("n_variants").asLong());
            assertEquals(0.0, rec.get("pi").asDouble(), EPS);
        }
    }

    // -----------------------------------------------------------------------
    // Divergence tests — POP1 vs POP2
    // -----------------------------------------------------------------------

    @Test
    void divergence_hudsonFst() {
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.divergence('chr1', 1, 600, 'POP1', 'POP2') " +
                "YIELD fst_hudson, fst_wc, dxy, da, n_variants"
            ).single();

            assertEquals(5L, rec.get("n_variants").asLong());
            assertEquals(0.0593, rec.get("fst_hudson").asDouble(), EPS);
        }
    }

    @Test
    void divergence_wcFst() {
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.divergence('chr1', 1, 600, 'POP1', 'POP2') " +
                "YIELD fst_wc"
            ).single();

            assertEquals(0.0328, rec.get("fst_wc").asDouble(), EPS);
        }
    }

    @Test
    void divergence_dxy() {
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.divergence('chr1', 1, 600, 'POP1', 'POP2') " +
                "YIELD dxy, da"
            ).single();

            // dxy = mean(p1*(1-p2) + p2*(1-p1)) = 0.375
            assertEquals(0.3750, rec.get("dxy").asDouble(), EPS);
            // da = dxy - (piWithin1 + piWithin2) / 2 = 0.0222
            assertEquals(0.0222, rec.get("da").asDouble(), EPS);
        }
    }

    @Test
    void divergence_fstPositive() {
        // Fst should be non-negative with these different populations
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.divergence('chr1', 1, 600, 'POP1', 'POP2') " +
                "YIELD fst_hudson, fst_wc"
            ).single();

            assertTrue(rec.get("fst_hudson").asDouble() >= 0.0);
            assertTrue(rec.get("fst_wc").asDouble() >= 0.0);
        }
    }

    @Test
    void divergence_emptyRegion() {
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.divergence('chr1', 700, 800, 'POP1', 'POP2') " +
                "YIELD n_variants, fst_hudson"
            ).single();

            assertEquals(0L, rec.get("n_variants").asLong());
            assertEquals(0.0, rec.get("fst_hudson").asDouble(), EPS);
        }
    }

    // -----------------------------------------------------------------------
    // SFS tests — folded and unfolded
    // -----------------------------------------------------------------------

    @Test
    void sfsPop1_folded() {
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.sfs('chr1', 1, 600, 'POP1', true) " +
                "YIELD sfs, n_variants, max_ac"
            ).single();

            assertEquals(5L, rec.get("n_variants").asLong());
            assertEquals(10L, rec.get("max_ac").asLong());

            // Folded SFS: ac -> min(ac, an-ac)
            // [2,4,4,2,0] -> folded: [2,4,4,2,0] -> bins 0..5: [1,0,2,0,2,0]
            List<Object> sfs = rec.get("sfs").asList();
            assertEquals(6, sfs.size()); // bins 0..5
            assertEquals(1L, sfs.get(0)); // bin 0: monomorphic variant (ac=0)
            assertEquals(0L, sfs.get(1)); // bin 1
            assertEquals(2L, sfs.get(2)); // bin 2: two variants with folded ac=2
            assertEquals(0L, sfs.get(3)); // bin 3
            assertEquals(2L, sfs.get(4)); // bin 4: two variants with folded ac=4
            assertEquals(0L, sfs.get(5)); // bin 5
        }
    }

    @Test
    void sfsPop1_unfolded() {
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.sfs('chr1', 1, 600, 'POP1', false) " +
                "YIELD sfs, n_variants, n_polarized"
            ).single();

            assertEquals(5L, rec.get("n_variants").asLong());
            assertEquals(5L, rec.get("n_polarized").asLong());

            // Unfolded SFS: derived counts
            // v1: REF ancestral -> derived=ac=2
            // v2: REF ancestral -> derived=ac=4
            // v3: ALT ancestral -> derived=an-ac=10-6=4
            // v4: REF ancestral -> derived=ac=8
            // v5: REF ancestral -> derived=ac=0
            // histogram[0]=1, histogram[2]=1, histogram[4]=2, histogram[8]=1
            List<Object> sfs = rec.get("sfs").asList();
            assertEquals(11, sfs.size()); // bins 0..10
            assertEquals(1L, sfs.get(0));  // bin 0
            assertEquals(0L, sfs.get(1));  // bin 1
            assertEquals(1L, sfs.get(2));  // bin 2
            assertEquals(0L, sfs.get(3));  // bin 3
            assertEquals(2L, sfs.get(4));  // bin 4
            assertEquals(0L, sfs.get(5));  // bin 5
            assertEquals(0L, sfs.get(6));  // bin 6
            assertEquals(0L, sfs.get(7));  // bin 7
            assertEquals(1L, sfs.get(8));  // bin 8
            assertEquals(0L, sfs.get(9));  // bin 9
            assertEquals(0L, sfs.get(10)); // bin 10
        }
    }

    @Test
    void sfsPop2_folded() {
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.sfs('chr1', 1, 600, 'POP2', true) " +
                "YIELD sfs, n_variants, max_ac"
            ).single();

            assertEquals(5L, rec.get("n_variants").asLong());
            assertEquals(8L, rec.get("max_ac").asLong());

            // POP2 ac: [1,3,5,2,0] -> folded min(ac, 8-ac): [1,3,3,2,0]
            // bins 0..4: [1,1,1,2,0]
            List<Object> sfs = rec.get("sfs").asList();
            assertEquals(5, sfs.size()); // bins 0..4
            assertEquals(1L, sfs.get(0)); // bin 0: ac=0
            assertEquals(1L, sfs.get(1)); // bin 1: ac=1
            assertEquals(1L, sfs.get(2)); // bin 2: ac=2
            assertEquals(2L, sfs.get(3)); // bin 3: two variants with folded ac=3
            assertEquals(0L, sfs.get(4)); // bin 4
        }
    }

    @Test
    void sfs_emptyRegion() {
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.sfs('chr1', 700, 800, 'POP1', true) " +
                "YIELD sfs, n_variants"
            ).single();

            assertEquals(0L, rec.get("n_variants").asLong());
        }
    }

    // -----------------------------------------------------------------------
    // Joint SFS tests
    // -----------------------------------------------------------------------

    @Test
    void jointSfs_folded() {
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.joint_sfs('chr1', 1, 600, 'POP1', 'POP2', true) " +
                "YIELD joint_sfs, n_variants, max_ac1, max_ac2, dim1, dim2"
            ).single();

            assertEquals(5L, rec.get("n_variants").asLong());
            assertEquals(10L, rec.get("max_ac1").asLong());
            assertEquals(8L, rec.get("max_ac2").asLong());
            assertEquals(6L, rec.get("dim1").asLong()); // 10/2+1
            assertEquals(5L, rec.get("dim2").asLong()); // 8/2+1

            // Folded pairs: (pop1, pop2) = (2,1), (4,3), (4,3), (2,2), (0,0)
            // 2D histogram (6 x 5):
            //   [0,0]=1  (monomorphic variant)
            //   [2,1]=1  (v1)
            //   [2,2]=1  (v4)  -- wait, v4 folded pop1=min(8,2)=2, pop2=min(2,6)=2
            //   [4,3]=2  (v2 and v3)
            List<Object> jsfs = rec.get("joint_sfs").asList();
            int dim2 = 5;

            // [0][0] = 1 (monomorphic)
            assertEquals(1L, jsfs.get(0 * dim2 + 0));
            // [2][1] = 1 (v1: folded=(2,1))
            assertEquals(1L, jsfs.get(2 * dim2 + 1));
            // [2][2] = 1 (v4: folded=(2,2))
            assertEquals(1L, jsfs.get(2 * dim2 + 2));
            // [4][3] = 2 (v2 and v3: folded=(4,3))
            assertEquals(2L, jsfs.get(4 * dim2 + 3));

            // Verify total counts sum to 5
            long total = 0;
            for (Object v : jsfs) total += (Long) v;
            assertEquals(5L, total);
        }
    }

    @Test
    void jointSfs_unfolded() {
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.joint_sfs('chr1', 1, 600, 'POP1', 'POP2', false) " +
                "YIELD joint_sfs, n_variants, dim1, dim2"
            ).single();

            assertEquals(5L, rec.get("n_variants").asLong());
            assertEquals(11L, rec.get("dim1").asLong()); // 10+1
            assertEquals(9L, rec.get("dim2").asLong());  // 8+1

            // Unfolded: raw ac pairs (pop1, pop2) = (2,1), (4,3), (6,5), (8,2), (0,0)
            List<Object> jsfs = rec.get("joint_sfs").asList();
            int dim2 = 9;

            assertEquals(1L, jsfs.get(0 * dim2 + 0)); // (0,0)
            assertEquals(1L, jsfs.get(2 * dim2 + 1)); // (2,1)
            assertEquals(1L, jsfs.get(4 * dim2 + 3)); // (4,3)
            assertEquals(1L, jsfs.get(6 * dim2 + 5)); // (6,5)
            assertEquals(1L, jsfs.get(8 * dim2 + 2)); // (8,2)

            long total = 0;
            for (Object v : jsfs) total += (Long) v;
            assertEquals(5L, total);
        }
    }

    @Test
    void jointSfs_emptyRegion() {
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.joint_sfs('chr1', 700, 800, 'POP1', 'POP2', false) " +
                "YIELD n_variants"
            ).single();

            assertEquals(0L, rec.get("n_variants").asLong());
        }
    }

    // -----------------------------------------------------------------------
    // Cross-procedure consistency checks
    // -----------------------------------------------------------------------

    @Test
    void diversityAndSfs_segregatingSitesAgree() {
        // The number of segregating sites from diversity should match
        // the sum of non-zero bins in the folded SFS (excluding bin 0)
        try (Session session = driver.session()) {
            Record divRec = session.run(
                "CALL graphpop.diversity('chr1', 1, 600, 'POP1') " +
                "YIELD n_segregating"
            ).single();

            Record sfsRec = session.run(
                "CALL graphpop.sfs('chr1', 1, 600, 'POP1', true) " +
                "YIELD sfs"
            ).single();

            long nSeg = divRec.get("n_segregating").asLong();
            List<Object> sfs = sfsRec.get("sfs").asList();

            // Sum all bins except bin 0 (monomorphic)
            long sfsSegregating = 0;
            for (int i = 1; i < sfs.size(); i++) {
                sfsSegregating += (Long) sfs.get(i);
            }

            assertEquals(nSeg, sfsSegregating,
                "Segregating site count from diversity should match non-zero SFS bins");
        }
    }

    @Test
    void divergence_dxyGreaterThanOrEqualToDa() {
        // By definition da = dxy - mean(piWithin), so dxy >= da
        try (Session session = driver.session()) {
            Record rec = session.run(
                "CALL graphpop.divergence('chr1', 1, 600, 'POP1', 'POP2') " +
                "YIELD dxy, da"
            ).single();

            assertTrue(rec.get("dxy").asDouble() >= rec.get("da").asDouble(),
                "dxy should be >= da (since da = dxy - mean piWithin)");
        }
    }
}
