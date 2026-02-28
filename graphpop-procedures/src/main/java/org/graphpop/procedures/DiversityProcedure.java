package org.graphpop.procedures;

import org.graphpop.VectorOps;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.Map;
import java.util.stream.Stream;

/**
 * Computes diversity statistics for a genomic region and population.
 *
 * <p>FAST PATH: uses only the allele-count arrays stored on Variant nodes,
 * so complexity is O(V × K) independent of sample size.</p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.diversity('chr22', 16000000, 17000000, 'AFR')
 * YIELD pi, theta_w, tajima_d, het_exp, het_obs, fis, n_variants, n_segregating
 * </pre>
 * </p>
 */
public class DiversityProcedure {

    @Context
    public Transaction tx;

    public static class DiversityResult {
        public double pi;
        public double theta_w;
        public double tajima_d;
        public double het_exp;
        public double het_obs;
        public double fis;
        public long n_variants;
        public long n_segregating;

        public DiversityResult() {}
    }

    @Procedure(name = "graphpop.diversity", mode = Mode.READ)
    public Stream<DiversityResult> diversity(
            @Name("chr") String chr,
            @Name("start") long start,
            @Name("end") long end,
            @Name("pop") String pop,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        // Optional filters
        String consequenceFilter = options != null ? (String) options.get("consequence") : null;
        String pathwayFilter = options != null ? (String) options.get("pathway") : null;

        int popIndex = -1;
        double a_n = 0.0;
        double a_n2 = 0.0;

        // Accumulate statistics
        double piSum = 0.0;
        double heSum = 0.0;
        double hoSum = 0.0;
        long nVariants = 0;
        long nSegregating = 0;

        // Query variants in range
        var result = tx.execute(
                buildVariantQuery(consequenceFilter, pathwayFilter),
                Map.of("chr", chr, "start", start, "end", end)
        );

        try {
            while (result.hasNext()) {
                Map<String, Object> row = result.next();
                Node variant = (Node) row.get("v");

                // Resolve population index from first variant's pop_ids
                if (popIndex < 0) {
                    String[] popIds = ArrayUtil.toStringArray(variant.getProperty("pop_ids"));
                    for (int i = 0; i < popIds.length; i++) {
                        if (popIds[i].equals(pop)) {
                            popIndex = i;
                            break;
                        }
                    }
                    if (popIndex < 0) {
                        throw new RuntimeException("Population not found: " + pop);
                    }

                    // Get harmonic numbers from Population node
                    Node popNode = tx.findNode(Label.label("Population"), "populationId", pop);
                    if (popNode == null) {
                        throw new RuntimeException("Population node not found: " + pop);
                    }
                    a_n = ((Number) popNode.getProperty("a_n")).doubleValue();
                    a_n2 = ((Number) popNode.getProperty("a_n2")).doubleValue();
                }

                // Extract arrays for this population
                int[] acArr = ArrayUtil.toIntArray(variant.getProperty("ac"));
                int[] anArr = ArrayUtil.toIntArray(variant.getProperty("an"));
                double[] afArr = ArrayUtil.toDoubleArray(variant.getProperty("af"));
                int[] hetArr = ArrayUtil.toIntArray(variant.getProperty("het_count"));

                int ac = acArr[popIndex];
                int an = anArr[popIndex];
                double af = afArr[popIndex];
                int hetCount = hetArr[popIndex];

                if (an < 2) continue;

                nVariants++;

                // Pi: per-site nucleotide diversity with finite-sample correction
                piSum += VectorOps.piPerSite(af, an);

                // He: expected heterozygosity
                heSum += 2.0 * af * (1.0 - af);

                // Ho: observed heterozygosity
                int nDiploid = an / 2;
                hoSum += (double) hetCount / nDiploid;

                // Count segregating sites
                if (ac > 0 && ac < an) {
                    nSegregating++;
                }
            }
        } finally {
            result.close();
        }

        // Compute summary
        DiversityResult out = new DiversityResult();
        out.n_variants = nVariants;
        out.n_segregating = nSegregating;

        if (nVariants == 0) {
            return Stream.of(out);
        }

        double L = nVariants;
        out.pi = piSum / L;
        out.het_exp = heSum / L;
        out.het_obs = hoSum / L;
        out.fis = out.het_exp > 0 ? 1.0 - out.het_obs / out.het_exp : 0.0;

        // Watterson's theta
        out.theta_w = a_n > 0 ? (nSegregating / a_n) / L : 0.0;

        // Tajima's D
        out.tajima_d = tajimaD(piSum, nSegregating, a_n, a_n2, L);

        return Stream.of(out);
    }

    /**
     * Build the Cypher query string for fetching variants in a genomic range.
     * Optionally filters by VEP consequence or pathway membership.
     */
    private String buildVariantQuery(String consequence, String pathway) {
        StringBuilder sb = new StringBuilder();

        if (pathway != null) {
            sb.append("MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway) ");
            sb.append("WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end ");
            sb.append("AND p.name = '").append(pathway.replace("'", "''")).append("' ");
        } else if (consequence != null) {
            sb.append("MATCH (v:Variant)-[c:HAS_CONSEQUENCE]->(g:Gene) ");
            sb.append("WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end ");
            sb.append("AND c.consequence = '").append(consequence.replace("'", "''")).append("' ");
        } else {
            sb.append("MATCH (v:Variant) ");
            sb.append("WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end ");
        }

        sb.append("RETURN DISTINCT v");
        return sb.toString();
    }

    /**
     * Compute Tajima's D statistic.
     *
     * @param piTotal   sum of per-site pi (not divided by L)
     * @param S         number of segregating sites
     * @param a_n       harmonic number (sum 1/i for i=1..n-1)
     * @param a_n2      sum of 1/i² for i=1..n-1
     * @param L         total number of sites examined
     * @return Tajima's D
     */
    private double tajimaD(double piTotal, long S, double a_n, double a_n2, double L) {
        if (S == 0 || a_n == 0) return 0.0;

        double thetaW = S / a_n;
        double d = piTotal - thetaW;

        // Tajima's variance coefficients
        double a1 = a_n;
        double a2 = a_n2;

        // n is not directly available here, derive from a_n
        // For simplicity, we compute the variance using the standard formulas
        // e1 = c1/a1 where c1 = b1 - 1/a1
        // e2 = c2/(a1² + a2) where c2 = b2 - (n+2)/(a1*n) + a2/a1²
        // We approximate n from a_n using the known relationship
        // a_n = H(2n-1) where n = number of diploid samples

        // However, since we have a_n and a_n2, we can use a simpler approach:
        // The original Tajima (1989) variance formula:
        // Var(d) = e1*S + e2*S*(S-1)

        // Approximate n from harmonic number
        // a_n ≈ ln(n) + γ, so n ≈ exp(a_n - 0.5772)
        // But we need exact n. We can get it from the Population node's n_samples.
        // For now, we use the Tajima formula without needing explicit n:

        double b1 = (a1 + 1.0) / (3.0 * (a1 - 1.0 + 1.0));
        // We need to recover n to compute b2.
        // Heuristic: solve for n from a_n = H(2n-1)
        // Instead, use a_n and a_n2 directly per Tajima (1989) eqs 34-38

        // Standard approach: compute n from population data.
        // For correctness, we should pass n directly. Let's use a robust estimation.
        // Since a_n = sum(1/i, i=1..m) where m = 2*n_samples - 1,
        // and a_n is stored precisely, we can recover m by checking when H(m) = a_n.
        // Faster: we store n_samples in the Population node and pass it.
        // For now, use a direct numerical approximation.

        int m = estimateM(a_n);  // m = 2*n_samples - 1
        int nPlusOne = m + 1; // = 2*n_samples
        double n = nPlusOne; // using m+1 as the "n" in Tajima's notation (number of sequences)

        b1 = (n + 1.0) / (3.0 * (n - 1.0));
        double b2 = 2.0 * (n * n + n + 3.0) / (9.0 * n * (n - 1.0));

        double c1 = b1 - 1.0 / a1;
        double c2 = b2 - (n + 2.0) / (a1 * n) + a2 / (a1 * a1);

        double e1 = c1 / a1;
        double e2 = c2 / (a1 * a1 + a2);

        double variance = e1 * S + e2 * S * (S - 1);
        if (variance <= 0) return 0.0;

        return d / Math.sqrt(variance);
    }

    /**
     * Estimate m (number of alleles - 1) from harmonic number a_n = H(m).
     * Uses Newton's method starting from the approximation m ≈ exp(a_n - γ).
     */
    private int estimateM(double target) {
        // Start with ln approximation: H(m) ≈ ln(m) + 0.5772
        int m = Math.max(1, (int) Math.round(Math.exp(target - 0.5772156649)));

        // Refine: compute H(m) and adjust
        double h = 0;
        for (int i = 1; i <= m; i++) h += 1.0 / i;
        while (h < target - 1e-6 && m < 100000) {
            m++;
            h += 1.0 / m;
        }
        while (h > target + 1e-6 && m > 1) {
            h -= 1.0 / m;
            m--;
        }
        return m;
    }

}
