package org.graphpop.procedures.recombination;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.graphpop.procedures.arg.ArgFixtureLoader;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.neo4j.driver.Driver;
import org.neo4j.driver.GraphDatabase;
import org.neo4j.driver.Record;
import org.neo4j.driver.Session;
import org.neo4j.harness.Neo4j;
import org.neo4j.harness.Neo4jBuilders;

import java.io.InputStream;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Integration test for {@link ArgBreakpointsProcedure}. Validates
 * per-window breakpoint counts + Hudson-scaled ρ on the 20-sample
 * fixture (6 marginal trees, ~5 internal breakpoints across
 * 50 kb).
 */
class ArgBreakpointsProcedureTest {

    private static final String RUN_ID = "egrm_fixture_20samples";

    private static Neo4j embeddedNeo4j;
    private static Driver driver;

    @BeforeAll
    static void setUp() throws Exception {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(ArgBreakpointsProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        ObjectMapper mapper = new ObjectMapper();
        JsonNode arg;
        try (InputStream is = ArgBreakpointsProcedureTest.class
                .getResourceAsStream("/egrm_fixture_20samples_arg.json")) {
            arg = mapper.readTree(is);
        }
        ArgFixtureLoader.loadInto(driver, arg, RUN_ID);
    }

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    @Test
    void arg_breakpoints_emits_one_row_per_window() {
        try (Session session = driver.session()) {
            // 50_000 bp sequence, non-overlapping 10 kb → 5 windows.
            List<Record> rs = session.run(
                "CALL graphpop.recombination.arg_breakpoints($rid, "
              + "{window_size: 10000, step: 10000}) "
              + "YIELD start, end, n_breakpoints, n_marginal_trees, "
              + "total_branch_length, rho_per_bp "
              + "RETURN start, end, n_breakpoints, n_marginal_trees, "
              + "total_branch_length, rho_per_bp ORDER BY start",
                Map.of("rid", RUN_ID)).list();
            assertEquals(5, rs.size());
            for (int i = 0; i < rs.size(); i++) {
                Record r = rs.get(i);
                assertEquals(i * 10_000L, r.get("start").asLong());
                assertEquals((i + 1) * 10_000L, r.get("end").asLong());
                // total branch length must be positive on neutral data.
                assertTrue(r.get("total_branch_length").asDouble() > 0.0,
                        "window " + i + " has no branches");
                // ρ ≥ 0 by construction; finite.
                double rho = r.get("rho_per_bp").asDouble();
                assertTrue(rho >= 0.0 && Double.isFinite(rho),
                        "ρ_per_bp = " + rho + " in window " + i);
                // n_marginal_trees consistent with n_breakpoints + 1.
                assertEquals(r.get("n_breakpoints").asLong() == 0 ? 1
                                : r.get("n_breakpoints").asLong() + 1,
                        r.get("n_marginal_trees").asLong());
            }
        }
    }

    @Test
    void arg_breakpoints_overlapping_windows_yields_more_rows() {
        try (Session session = driver.session()) {
            // 50_000 bp, 10 kb window, 5 kb step → 10 windows.
            List<Record> rs = session.run(
                "CALL graphpop.recombination.arg_breakpoints($rid, "
              + "{window_size: 10000, step: 5000}) "
              + "YIELD start RETURN start ORDER BY start",
                Map.of("rid", RUN_ID)).list();
            assertEquals(10, rs.size());
        }
    }

    @Test
    void arg_breakpoints_total_count_matches_fixture_n_trees_minus_one() {
        try (Session session = driver.session()) {
            // 6 marginal trees on the 20-sample fixture → 5 breakpoints
            // across the genome. With non-overlapping 10 kb windows
            // covering the full 50 kb, the sum of n_breakpoints across
            // all windows should equal the number of breakpoints with
            // strict-interior rule (boundaries between windows are
            // missed, so the count can be ≤ 5).
            Long sum = session.run(
                "CALL graphpop.recombination.arg_breakpoints($rid, "
              + "{window_size: 50000, step: 50000}) "
              + "YIELD n_breakpoints "
              + "RETURN sum(n_breakpoints) AS s",
                Map.of("rid", RUN_ID)).single().get("s").asLong();
            // With one window covering [0, 50000), all 5 internal
            // breakpoints are captured.
            assertEquals(5L, sum, "expected 5 breakpoints for 6 trees");
        }
    }

    @Test
    void arg_breakpoints_invalid_options_throw() {
        try (Session session = driver.session()) {
            assertThrows(Exception.class, () -> session.run(
                "CALL graphpop.recombination.arg_breakpoints($rid, "
              + "{window_size: 0}) YIELD start RETURN start",
                Map.of("rid", RUN_ID)).list());
        }
    }

    @Test
    void arg_breakpoints_unknown_run_returns_empty() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.recombination.arg_breakpoints("
              + "'does_not_exist', {}) "
              + "YIELD start RETURN start").list();
            assertEquals(0, rs.size());
        }
    }
}
