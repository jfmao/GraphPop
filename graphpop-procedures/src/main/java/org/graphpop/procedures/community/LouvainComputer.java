package org.graphpop.procedures.community;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Pure-Java Louvain modularity optimisation
 * (Blondel, Guillaume, Lambiotte &amp; Lefebvre, 2008).
 *
 * <p>Greedy, multi-level, deterministic given a seed. Operates on a
 * packed undirected weighted graph; returns per-node community ids
 * (re-labelled to {@code [0, n_communities)}) plus the final
 * modularity.</p>
 *
 * <p>Skipped: Leiden's refinement phase (Traag et al. 2019) — gives
 * stricter community guarantees but adds ~30 % code; defer until
 * downstream needs it.</p>
 */
public final class LouvainComputer {

    private LouvainComputer() {}

    /** Result bundle. */
    public static final class Result {
        public final int[] community;
        public final double modularity;
        public final int nCommunities;

        Result(int[] community, double modularity, int nCommunities) {
            this.community = community;
            this.modularity = modularity;
            this.nCommunities = nCommunities;
        }
    }

    /**
     * Run Louvain on the given undirected weighted graph.
     *
     * @param n         number of nodes (packed indices {@code [0, n)}).
     * @param edges     undirected edges as parallel arrays {@code (u[i], v[i], w[i])};
     *                  each edge stored once, self-loops permitted, weights {@code > 0}.
     * @param u         edge endpoint a (length = m)
     * @param v         edge endpoint b (length = m)
     * @param w         edge weight    (length = m)
     * @param maxOuter  cap on outer (aggregation) iterations
     * @param maxInner  cap on inner (local-move) iterations per level
     * @param tol       minimum modularity gain to accept a move
     * @param seed      RNG seed (deterministic node-iteration order)
     */
    public static Result run(int n, int[] u, int[] v, double[] w,
                              int maxOuter, int maxInner, double tol,
                              long seed) {
        if (n <= 0) {
            return new Result(new int[0], 0.0, 0);
        }
        // Phase 0: build the level-0 packed graph.
        Level lvl = Level.fromEdges(n, u, v, w);

        // Outer loop: iterate until no further aggregation improves modularity.
        int[] currentToOriginal = identityArray(n);
        Random rng = new Random(seed);
        for (int outer = 0; outer < maxOuter; outer++) {
            // Phase 1: local moves on current level.
            boolean moved = optimiseLocally(lvl, maxInner, tol, rng);
            if (!moved && outer > 0) break;

            // Compress communities to dense ids [0, k).
            int[] dense = denseRelabel(lvl.community);
            int k = uniqueCount(dense);
            if (k == lvl.n) break;  // every node alone — no further compression

            // Map level-0 → current → next-level community id.
            for (int i = 0; i < currentToOriginal.length; i++) {
                int curr = currentToOriginal[i];  // index into lvl
                currentToOriginal[i] = dense[lvl.community[curr]];
            }
            lvl = aggregate(lvl, dense, k);
        }

        // Final relabel + modularity.
        int[] finalDense = denseRelabel(currentToOriginal);
        int finalK = uniqueCount(finalDense);
        double q = modularityFromOriginal(n, u, v, w, finalDense);
        return new Result(finalDense, q, finalK);
    }

    private static boolean optimiseLocally(Level lvl, int maxInner,
                                            double tol, Random rng) {
        boolean anyMove = false;
        int[] order = new int[lvl.n];
        for (int i = 0; i < lvl.n; i++) order[i] = i;
        for (int iter = 0; iter < maxInner; iter++) {
            shuffle(order, rng);
            boolean moved = false;
            for (int idx : order) {
                int bestC = chooseCommunity(lvl, idx, tol);
                if (bestC != lvl.community[idx]) {
                    move(lvl, idx, bestC);
                    moved = true;
                    anyMove = true;
                }
            }
            if (!moved) break;
        }
        return anyMove;
    }

    /** Pick the community that maximises modularity gain (or the current one). */
    private static int chooseCommunity(Level lvl, int idx, double tol) {
        int currentC = lvl.community[idx];
        double m2 = lvl.totalWeight * 2.0;  // 2m
        double k_i = lvl.degree[idx];

        // Aggregate per-community connections from this node.
        Map<Integer, Double> linkToCommunity = new HashMap<>();
        double selfLoop = 0.0;
        int[] adj = lvl.adj[idx];
        double[] adjW = lvl.adjW[idx];
        for (int e = 0; e < adj.length; e++) {
            int nb = adj[e];
            double we = adjW[e];
            if (nb == idx) {
                selfLoop += we;
                continue;
            }
            int cN = lvl.community[nb];
            linkToCommunity.merge(cN, we, Double::sum);
        }
        // Subtract this node from its own community first.
        double sigmaTotCurrent = lvl.communityTot[currentC] - k_i;

        double bestGain = 0.0;
        int bestC = currentC;
        for (Map.Entry<Integer, Double> e : linkToCommunity.entrySet()) {
            int c = e.getKey();
            double k_i_in = e.getValue();
            double sigmaTot = lvl.communityTot[c];
            if (c == currentC) sigmaTot = sigmaTotCurrent;
            // Blondel et al. 2008 simplified gain: ΔQ = k_i_in / m
            //                                     - Σ_tot · k_i / (2m²)
            double gain = k_i_in / lvl.totalWeight
                          - (sigmaTot * k_i) / (m2 * lvl.totalWeight);
            if (gain > bestGain + tol) {
                bestGain = gain;
                bestC = c;
            }
        }
        // Self-loop is handled by leaving the node in its current community
        // when no neighbour beats it; bestC defaults to currentC.
        return bestC;
    }

    private static void move(Level lvl, int idx, int newC) {
        int oldC = lvl.community[idx];
        double k = lvl.degree[idx];
        lvl.communityTot[oldC] -= k;
        lvl.communityTot[newC] += k;
        lvl.community[idx] = newC;
    }

    private static Level aggregate(Level lvl, int[] dense, int k) {
        // Group edges by (denseCommunity[u], denseCommunity[v]).
        Map<Long, Double> edgeMap = new HashMap<>();
        for (int idx = 0; idx < lvl.n; idx++) {
            int cu = dense[lvl.community[idx]];
            int[] adj = lvl.adj[idx];
            double[] adjW = lvl.adjW[idx];
            for (int e = 0; e < adj.length; e++) {
                int nb = adj[e];
                if (nb < idx && nb != idx) continue;  // each undirected edge once
                int cv = dense[lvl.community[nb]];
                int a = Math.min(cu, cv);
                int b = Math.max(cu, cv);
                long key = ((long) a << 32) | (b & 0xFFFFFFFFL);
                edgeMap.merge(key, adjW[e], Double::sum);
            }
        }
        int[] uA = new int[edgeMap.size()];
        int[] vA = new int[edgeMap.size()];
        double[] wA = new double[edgeMap.size()];
        int i = 0;
        for (Map.Entry<Long, Double> e : edgeMap.entrySet()) {
            long key = e.getKey();
            uA[i] = (int) (key >>> 32);
            vA[i] = (int) (key & 0xFFFFFFFFL);
            wA[i] = e.getValue();
            i++;
        }
        return Level.fromEdges(k, uA, vA, wA);
    }

    private static int[] denseRelabel(int[] arr) {
        Map<Integer, Integer> map = new HashMap<>();
        int[] out = new int[arr.length];
        int next = 0;
        for (int i = 0; i < arr.length; i++) {
            Integer d = map.get(arr[i]);
            if (d == null) {
                d = next++;
                map.put(arr[i], d);
            }
            out[i] = d;
        }
        return out;
    }

    private static int uniqueCount(int[] dense) {
        int max = -1;
        for (int x : dense) if (x > max) max = x;
        return max + 1;
    }

    private static int[] identityArray(int n) {
        int[] a = new int[n];
        for (int i = 0; i < n; i++) a[i] = i;
        return a;
    }

    private static void shuffle(int[] arr, Random rng) {
        for (int i = arr.length - 1; i > 0; i--) {
            int j = rng.nextInt(i + 1);
            int tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
        }
    }

    /**
     * Newman 2006 modularity computed against the level-0 graph:
     *
     * <pre>
     * Q = (1/2m) Σ_c [Σ_in(c) − Σ_tot(c)² / (2m)]
     * </pre>
     *
     * <p>where {@code Σ_in(c)} = sum of weights of edges with both
     * endpoints in {@code c} (off-diagonal entries counted twice,
     * self-loops counted once — i.e. literal Σ_{i,j∈c} A_ij), and
     * {@code Σ_tot(c)} = Σ_{i∈c} k_i.</p>
     */
    private static double modularityFromOriginal(int n, int[] u, int[] v,
                                                  double[] w, int[] comm) {
        double[] degree = new double[n];
        double m = 0.0;
        for (int e = 0; e < u.length; e++) {
            degree[u[e]] += w[e];
            if (u[e] != v[e]) degree[v[e]] += w[e];
            m += w[e];
        }
        if (m <= 0) return 0.0;

        int nComm = uniqueCount(comm);
        double[] internal = new double[nComm];
        double[] total = new double[nComm];
        for (int i = 0; i < n; i++) total[comm[i]] += degree[i];
        for (int e = 0; e < u.length; e++) {
            if (comm[u[e]] != comm[v[e]]) continue;
            int c = comm[u[e]];
            // A_ij + A_ji for i≠j; A_ii for i=j.
            internal[c] += (u[e] == v[e]) ? w[e] : 2.0 * w[e];
        }

        double q = 0.0;
        for (int c = 0; c < nComm; c++) {
            q += internal[c] - (total[c] * total[c]) / (2.0 * m);
        }
        return q / (2.0 * m);
    }

    /**
     * Per-level packed state: adjacency lists, per-node degree,
     * per-community totals, current community assignment.
     */
    private static final class Level {
        final int n;
        final int[][] adj;
        final double[][] adjW;
        final double[] degree;        // Σ_j A_ij  (includes 2× self-loops)
        final int[] community;
        final double[] communityTot;  // Σ_tot per community
        final double totalWeight;     // m = (1/2) Σ_ij A_ij

        Level(int n, int[][] adj, double[][] adjW,
              double[] degree, int[] community,
              double[] communityTot, double totalWeight) {
            this.n = n;
            this.adj = adj;
            this.adjW = adjW;
            this.degree = degree;
            this.community = community;
            this.communityTot = communityTot;
            this.totalWeight = totalWeight;
        }

        static Level fromEdges(int n, int[] uA, int[] vA, double[] wA) {
            int m = uA.length;
            int[] deg = new int[n];
            for (int e = 0; e < m; e++) {
                deg[uA[e]]++;
                if (uA[e] != vA[e]) deg[vA[e]]++;
            }
            int[][] adj = new int[n][];
            double[][] adjW = new double[n][];
            int[] cursor = new int[n];
            for (int i = 0; i < n; i++) {
                adj[i] = new int[deg[i]];
                adjW[i] = new double[deg[i]];
            }
            for (int e = 0; e < m; e++) {
                int a = uA[e], b = vA[e];
                double we = wA[e];
                adj[a][cursor[a]] = b;
                adjW[a][cursor[a]] = we;
                cursor[a]++;
                if (a != b) {
                    adj[b][cursor[b]] = a;
                    adjW[b][cursor[b]] = we;
                    cursor[b]++;
                }
            }
            double[] degree = new double[n];
            double total = 0.0;
            for (int e = 0; e < m; e++) {
                degree[uA[e]] += wA[e];
                if (uA[e] != vA[e]) degree[vA[e]] += wA[e];
                total += wA[e];
            }
            int[] comm = new int[n];
            double[] commTot = new double[n];
            for (int i = 0; i < n; i++) {
                comm[i] = i;
                commTot[i] = degree[i];
            }
            return new Level(n, adj, adjW, degree, comm, commTot, total);
        }
    }
}
