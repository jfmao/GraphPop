package org.graphpop.procedures.pairwise;

import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

class FamilyClustererTest {

    private static Map<String, FamilyClusterer.Assignment> byId(
            List<FamilyClusterer.Assignment> ass) {
        Map<String, FamilyClusterer.Assignment> m = new HashMap<>();
        for (FamilyClusterer.Assignment a : ass) m.put(a.sampleId, a);
        return m;
    }

    @Test
    void empty_graph_each_sample_is_own_family() {
        List<String> samples = List.of("S0", "S1", "S2", "S3");
        List<String[]> edges = List.of();
        List<FamilyClusterer.Assignment> ass =
                FamilyClusterer.cluster(samples, edges);
        for (FamilyClusterer.Assignment a : ass) {
            assertEquals(a.sampleId, a.familyId);
            assertEquals(1, a.familySize);
        }
    }

    @Test
    void disjoint_trios_and_singletons() {
        // 3 trios (S0-S1-S2, S3-S4-S5, S6-S7-S8) + 3 singletons (S9, S10, S11).
        List<String> samples = List.of(
                "S0", "S1", "S2", "S3", "S4", "S5",
                "S6", "S7", "S8", "S9", "S10", "S11");
        List<String[]> edges = List.of(
                new String[]{"S0", "S1"}, new String[]{"S1", "S2"},
                new String[]{"S3", "S4"}, new String[]{"S4", "S5"},
                new String[]{"S6", "S7"}, new String[]{"S7", "S8"});
        Map<String, FamilyClusterer.Assignment> ass =
                byId(FamilyClusterer.cluster(samples, edges));

        // Each trio has size 3.
        for (String s : List.of("S0", "S1", "S2")) {
            assertEquals(3, ass.get(s).familySize);
            assertEquals("S0", ass.get(s).familyId);
        }
        for (String s : List.of("S3", "S4", "S5")) {
            assertEquals(3, ass.get(s).familySize);
            assertEquals("S3", ass.get(s).familyId);
        }
        for (String s : List.of("S6", "S7", "S8")) {
            assertEquals(3, ass.get(s).familySize);
            assertEquals("S6", ass.get(s).familyId);
        }
        // Singletons.
        for (String s : List.of("S9", "S10", "S11")) {
            assertEquals(1, ass.get(s).familySize);
            assertEquals(s, ass.get(s).familyId);
        }
    }

    @Test
    void one_big_family() {
        List<String> samples = List.of("A", "B", "C", "D");
        // Chain of edges connecting all four.
        List<String[]> edges = List.of(
                new String[]{"A", "B"}, new String[]{"B", "C"},
                new String[]{"C", "D"});
        Map<String, FamilyClusterer.Assignment> ass =
                byId(FamilyClusterer.cluster(samples, edges));
        for (FamilyClusterer.Assignment a : ass.values()) {
            assertEquals(4, a.familySize);
            assertEquals("A", a.familyId, "smallest sample-id wins");
        }
    }

    @Test
    void edge_to_unknown_sample_is_skipped() {
        List<String> samples = List.of("A", "B");
        List<String[]> edges = List.of(
                new String[]{"A", "B"},
                new String[]{"A", "Z"});  // Z not in samples
        Map<String, FamilyClusterer.Assignment> ass =
                byId(FamilyClusterer.cluster(samples, edges));
        assertEquals(2, ass.get("A").familySize);
        assertEquals(2, ass.get("B").familySize);
    }

    @Test
    void family_id_is_lex_smallest_sample() {
        List<String> samples = List.of("ZZ", "AA", "MM");
        List<String[]> edges = List.of(
                new String[]{"ZZ", "AA"},
                new String[]{"AA", "MM"});
        Map<String, FamilyClusterer.Assignment> ass =
                byId(FamilyClusterer.cluster(samples, edges));
        for (FamilyClusterer.Assignment a : ass.values()) {
            assertEquals("AA", a.familyId);
            assertEquals(3, a.familySize);
        }
    }

    @Test
    void path_compression_correctness_on_long_chain() {
        // Chain of 100 samples connected in a line.
        int n = 100;
        List<String> samples = new java.util.ArrayList<>(n);
        for (int i = 0; i < n; i++) samples.add(String.format("S%03d", i));
        List<String[]> edges = new java.util.ArrayList<>(n - 1);
        for (int i = 0; i < n - 1; i++) {
            edges.add(new String[]{samples.get(i), samples.get(i + 1)});
        }
        Map<String, FamilyClusterer.Assignment> ass =
                byId(FamilyClusterer.cluster(samples, edges));
        for (FamilyClusterer.Assignment a : ass.values()) {
            assertEquals(n, a.familySize);
            assertEquals("S000", a.familyId);
        }
    }
}
