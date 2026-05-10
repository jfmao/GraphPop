package org.graphpop.procedures.community;

import org.junit.jupiter.api.Test;

import java.util.HashSet;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class LouvainComputerTest {

    /**
     * Three disjoint K4 cliques bridged by 3 weak edges. Louvain must
     * recover the three cliques as the 3 communities; modularity is
     * high (~0.55).
     */
    @Test
    void three_disjoint_cliques_resolved() {
        int n = 12;  // nodes 0-3, 4-7, 8-11
        // Within-clique edges (weight 10 each)
        // Cross-clique bridges (weight 1)
        int[] u = new int[6 * 3 + 3];
        int[] v = new int[6 * 3 + 3];
        double[] w = new double[6 * 3 + 3];
        int idx = 0;
        for (int base : new int[]{0, 4, 8}) {
            for (int i = 0; i < 4; i++) {
                for (int j = i + 1; j < 4; j++) {
                    u[idx] = base + i;
                    v[idx] = base + j;
                    w[idx] = 10.0;
                    idx++;
                }
            }
        }
        // Bridges
        u[idx] = 0;  v[idx] = 4;  w[idx] = 1.0;  idx++;
        u[idx] = 4;  v[idx] = 8;  w[idx] = 1.0;  idx++;
        u[idx] = 8;  v[idx] = 0;  w[idx] = 1.0;  idx++;

        LouvainComputer.Result r = LouvainComputer.run(
                n, u, v, w, 20, 50, 1e-9, 42L);
        assertEquals(3, r.nCommunities, "must recover 3 communities");
        // Same community within each clique.
        for (int base : new int[]{0, 4, 8}) {
            int c0 = r.community[base];
            for (int k = 1; k < 4; k++) {
                assertEquals(c0, r.community[base + k],
                        "clique " + base + " split: node " + (base + k));
            }
        }
        // Cross-clique nodes belong to different communities.
        Set<Integer> seen = new HashSet<>();
        seen.add(r.community[0]);
        seen.add(r.community[4]);
        seen.add(r.community[8]);
        assertEquals(3, seen.size());
        // Modularity at least 0.55 on this fixture (analytical lower bound).
        assertTrue(r.modularity >= 0.55,
                "modularity " + r.modularity + " below 0.55 threshold");
    }

    /**
     * Single edge (2 nodes connected) — Louvain places them in one
     * community. Modularity is exactly 0 by Newman's formula: the
     * sum of weights inside the community (= 2) equals
     * {@code Σ_tot² / (2m) = 4 / 2 = 2}, so Q = 0.
     */
    @Test
    void single_edge_one_community() {
        LouvainComputer.Result r = LouvainComputer.run(
                2, new int[]{0}, new int[]{1}, new double[]{1.0},
                10, 20, 1e-9, 7L);
        assertEquals(1, r.nCommunities);
        assertEquals(r.community[0], r.community[1]);
        assertEquals(0.0, r.modularity, 1e-12);
    }

    /**
     * Empty graph (no edges). Each node forms its own (trivial)
     * community; modularity = 0.
     */
    @Test
    void empty_graph_each_node_singleton() {
        LouvainComputer.Result r = LouvainComputer.run(
                5, new int[0], new int[0], new double[0],
                10, 20, 1e-9, 1L);
        assertEquals(5, r.nCommunities);
        Set<Integer> seen = new HashSet<>();
        for (int c : r.community) seen.add(c);
        assertEquals(5, seen.size());
        assertEquals(0.0, r.modularity);
    }

    /**
     * Single isolated node (n=1). Sanity.
     */
    @Test
    void single_node_no_edges() {
        LouvainComputer.Result r = LouvainComputer.run(
                1, new int[0], new int[0], new double[0],
                10, 20, 1e-9, 0L);
        assertEquals(1, r.nCommunities);
        assertEquals(0, r.community[0]);
    }

    /**
     * Determinism: same seed → identical communities.
     */
    @Test
    void deterministic_under_same_seed() {
        int n = 8;
        int[] u = {0, 1, 2, 3, 0, 4, 5, 6};
        int[] v = {1, 2, 3, 0, 4, 5, 6, 7};
        double[] w = {1, 1, 1, 1, 0.1, 1, 1, 1};
        LouvainComputer.Result r1 = LouvainComputer.run(
                n, u, v, w, 10, 20, 1e-9, 42L);
        LouvainComputer.Result r2 = LouvainComputer.run(
                n, u, v, w, 10, 20, 1e-9, 42L);
        assertArrayEquals(r1.community, r2.community);
        assertEquals(r1.modularity, r2.modularity, 1e-12);
    }

    /**
     * Weights affect the partition: pump up an inter-community edge
     * and the two communities merge.
     */
    @Test
    void heavy_inter_edge_merges_communities() {
        // Two K3 cliques bridged with a tiny edge → 2 communities.
        int n = 6;
        int[] uLight = {0, 0, 1, 3, 3, 4, 0};
        int[] vLight = {1, 2, 2, 4, 5, 5, 3};
        double[] wLight = {10, 10, 10, 10, 10, 10, 0.01};
        LouvainComputer.Result light = LouvainComputer.run(
                n, uLight, vLight, wLight, 10, 20, 1e-9, 9L);
        assertEquals(2, light.nCommunities);

        // Same topology but the bridge weight dominates → single community.
        double[] wHeavy = {1, 1, 1, 1, 1, 1, 100};
        LouvainComputer.Result heavy = LouvainComputer.run(
                n, uLight, vLight, wHeavy, 10, 20, 1e-9, 9L);
        assertEquals(1, heavy.nCommunities);
    }
}
