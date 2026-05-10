package org.graphpop.procedures.pairwise;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Disjoint-set / union-find for connected-components family
 * clustering. Each connected component on the relatedness graph is
 * one extended family.
 */
public final class FamilyClusterer {

    private FamilyClusterer() {}

    /**
     * Per-sample assignment.
     */
    public static final class Assignment {
        public final String sampleId;
        public final String familyId;     // representative sample of the component
        public final int familySize;

        Assignment(String sampleId, String familyId, int familySize) {
            this.sampleId = sampleId;
            this.familyId = familyId;
            this.familySize = familySize;
        }
    }

    /**
     * Compute connected components.
     *
     * @param samples all sample IDs in the cohort (length n).
     *                Singletons (no edges) get their own family of size 1.
     * @param edges   undirected (sa, sb) pairs joining two samples.
     */
    public static List<Assignment> cluster(List<String> samples,
                                            List<String[]> edges) {
        int n = samples.size();
        Map<String, Integer> idx = new HashMap<>(n * 2);
        for (int i = 0; i < n; i++) idx.put(samples.get(i), i);

        int[] parent = new int[n];
        int[] rank = new int[n];
        for (int i = 0; i < n; i++) parent[i] = i;

        for (String[] e : edges) {
            Integer a = idx.get(e[0]);
            Integer b = idx.get(e[1]);
            if (a == null || b == null) continue;
            union(parent, rank, a, b);
        }

        // Find canonical sample per component (smallest sampleId).
        Map<Integer, String> rootToFamilyId = new HashMap<>();
        Map<Integer, Integer> rootSize = new HashMap<>();
        for (int i = 0; i < n; i++) {
            int r = find(parent, i);
            String sid = samples.get(i);
            String existing = rootToFamilyId.get(r);
            if (existing == null || sid.compareTo(existing) < 0) {
                rootToFamilyId.put(r, sid);
            }
            rootSize.merge(r, 1, Integer::sum);
        }

        List<Assignment> out = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            int r = find(parent, i);
            out.add(new Assignment(
                    samples.get(i),
                    rootToFamilyId.get(r),
                    rootSize.get(r)));
        }
        return out;
    }

    private static int find(int[] parent, int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];  // path compression
            x = parent[x];
        }
        return x;
    }

    private static void union(int[] parent, int[] rank, int a, int b) {
        int ra = find(parent, a), rb = find(parent, b);
        if (ra == rb) return;
        if (rank[ra] < rank[rb]) { int t = ra; ra = rb; rb = t; }
        parent[rb] = ra;
        if (rank[ra] == rank[rb]) rank[ra]++;
    }
}
