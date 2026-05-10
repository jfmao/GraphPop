package org.graphpop.procedures.recombination;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.neo4j.driver.Driver;
import org.neo4j.driver.GraphDatabase;
import org.neo4j.driver.Record;
import org.neo4j.driver.Session;
import org.neo4j.harness.Neo4j;
import org.neo4j.harness.Neo4jBuilders;

import java.util.List;
import java.util.Map;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Integration + unit tests for {@link LdhatMcmcProcedure}.
 */
class LdhatMcmcProcedureTest {

    private static Neo4j embeddedNeo4j;
    private static Driver driver;

    @BeforeAll
    static void setUp() {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(LdhatMcmcProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());
        try (Session session = driver.session()) {
            for (int i = 0; i < 10; i++) {
                session.run(
                    "CREATE (:Sample {sampleId: $sid})",
                    Map.of("sid", "S" + i));
            }
            // Two variants close together with high r² → posterior should
            // favour low ρ.
            createVariant(session, "v1", 100L);
            createVariant(session, "v2", 200L);
            for (int i = 0; i < 5; i++) {
                createCarries(session, "S" + i, "v1");
                createCarries(session, "S" + i, "v2");
            }
            // Two variants far apart with low r² → posterior should
            // favour higher ρ.
            createVariant(session, "v3", 1_100L);
            createVariant(session, "v4", 5_900L);
            for (String s : new String[]{"S0", "S1", "S2", "S3", "S4"}) {
                createCarries(session, s, "v3");
            }
            for (String s : new String[]{"S2", "S3", "S5", "S6", "S7"}) {
                createCarries(session, s, "v4");
            }
        }
    }

    private static void createVariant(Session session, String vid, long pos) {
        session.run("CREATE (:Variant {variantId: $vid, position: $pos})",
            Map.of("vid", vid, "pos", pos));
    }

    private static void createCarries(Session session, String sid, String vid) {
        session.run(
            "MATCH (s:Sample {sampleId: $sid}), (v:Variant {variantId: $vid}) "
          + "CREATE (s)-[:CARRIES {gt: 1}]->(v)",
            Map.of("sid", sid, "vid", vid));
    }

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    private static List<String> sids10() {
        return List.of("S0", "S1", "S2", "S3", "S4",
                       "S5", "S6", "S7", "S8", "S9");
    }

    @Test
    void mcmc_per_window_returns_posterior_summary() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.recombination.ldhat_mcmc($sids, "
              + "{window_size: 10000, step: 10000, "
              + " max_pair_distance: 10000, "
              + " n_iter: 500, burn_in: 100, seed: 42}) "
              + "YIELD start, end, n_variant_pairs, rho_posterior_mean, "
              + "rho_lower_2_5, rho_upper_97_5, n_iter, n_accepted "
              + "RETURN start, end, n_variant_pairs, rho_posterior_mean, "
              + "rho_lower_2_5, rho_upper_97_5, n_iter, n_accepted "
              + "ORDER BY start",
                Map.of("sids", sids10())).list();
            assertTrue(rs.size() >= 1);
            for (Record r : rs) {
                long nPairs = r.get("n_variant_pairs").asLong();
                if (nPairs == 0) {
                    assertTrue(Double.isNaN(
                            r.get("rho_posterior_mean").asDouble()));
                    continue;
                }
                double mean = r.get("rho_posterior_mean").asDouble();
                double lo = r.get("rho_lower_2_5").asDouble();
                double hi = r.get("rho_upper_97_5").asDouble();
                assertTrue(Double.isFinite(mean) && mean > 0.0);
                assertTrue(lo <= mean,
                        "lower " + lo + " > mean " + mean);
                assertTrue(mean <= hi,
                        "mean " + mean + " > upper " + hi);
                assertEquals(500L, r.get("n_iter").asLong());
            }
        }
    }

    @Test
    void mcmc_deterministic_under_same_seed() {
        try (Session session = driver.session()) {
            List<Record> r1 = session.run(
                "CALL graphpop.recombination.ldhat_mcmc($sids, "
              + "{window_size: 10000, step: 10000, "
              + " max_pair_distance: 10000, "
              + " n_iter: 200, burn_in: 50, seed: 1}) "
              + "YIELD start, rho_posterior_mean "
              + "RETURN start, rho_posterior_mean ORDER BY start",
                Map.of("sids", sids10())).list();
            List<Record> r2 = session.run(
                "CALL graphpop.recombination.ldhat_mcmc($sids, "
              + "{window_size: 10000, step: 10000, "
              + " max_pair_distance: 10000, "
              + " n_iter: 200, burn_in: 50, seed: 1}) "
              + "YIELD start, rho_posterior_mean "
              + "RETURN start, rho_posterior_mean ORDER BY start",
                Map.of("sids", sids10())).list();
            assertEquals(r1.size(), r2.size());
            for (int i = 0; i < r1.size(); i++) {
                assertEquals(
                        r1.get(i).get("rho_posterior_mean").asDouble(),
                        r2.get(i).get("rho_posterior_mean").asDouble(),
                        1e-12);
            }
        }
    }

    @Test
    void mcmc_invalid_options_throw() {
        try (Session session = driver.session()) {
            assertThrows(Exception.class, () -> session.run(
                "CALL graphpop.recombination.ldhat_mcmc($sids, "
              + "{window_size: 0}) YIELD start RETURN start",
                Map.of("sids", sids10())).list());
            assertThrows(Exception.class, () -> session.run(
                "CALL graphpop.recombination.ldhat_mcmc($sids, "
              + "{prior_log_lo: 0.0, prior_log_hi: -2.0}) "
              + "YIELD start RETURN start",
                Map.of("sids", sids10())).list());
        }
    }

    @Test
    void mcmc_empty_sample_list_returns_empty() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.recombination.ldhat_mcmc([], {}) "
              + "YIELD start RETURN start").list();
            assertEquals(0, rs.size());
        }
    }

    @Test
    void run_mcmc_for_window_unit_test_recovers_known_rho() {
        // Build a deterministic synthetic dataset: ρ = 1e-3 per bp,
        // pairs at distances {100, 200, ... , 1000} bp with r² =
        // E[r²|n, ρ·d] + tiny noise. The posterior should concentrate
        // near ρ = 1e-3.
        double trueRho = 1e-3;
        Random data_rng = new Random(0);
        List<LdPairLoader.Pair> pairs = new java.util.ArrayList<>();
        for (long d = 100; d <= 1000; d += 100) {
            double mu = HudsonRecombination.expectedR2(trueRho * d);
            for (int rep = 0; rep < 5; rep++) {
                double r2 = mu + data_rng.nextGaussian() * 0.005;
                pairs.add(new LdPairLoader.Pair(d, Math.max(0.0, r2)));
            }
        }
        LdhatMcmcProcedure.McmcSummary s =
                LdhatMcmcProcedure.runMcmcForWindow(
                        pairs, 0.05, 0.3, -12.0, -2.0,
                        2000, 500, new Random(42));
        // Expect posterior mean within an order of magnitude of truth.
        assertTrue(s.posteriorMean > 1e-5 && s.posteriorMean < 1e-1,
                "posterior mean " + s.posteriorMean + " out of range");
        // CI should bracket the posterior mean.
        assertTrue(s.lower2_5 <= s.posteriorMean);
        assertTrue(s.posteriorMean <= s.upper97_5);
    }
}
