package org.graphpop.procedures;

import org.graphpop.VectorOps;
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
 * so complexity is O(V x K) independent of sample size.</p>
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
        String consequenceFilter = options != null ? (String) options.get("consequence") : null;
        String pathwayFilter = options != null ? (String) options.get("pathway") : null;

        PopulationContext ctx = null;

        double piSum = 0.0;
        double heSum = 0.0;
        double hoSum = 0.0;
        long nVariants = 0;
        long nSegregating = 0;

        var result = tx.execute(
                buildVariantQuery(consequenceFilter, pathwayFilter),
                Map.of("chr", chr, "start", start, "end", end)
        );

        try {
            while (result.hasNext()) {
                Map<String, Object> row = result.next();
                Node variant = (Node) row.get("v");

                if (ctx == null) {
                    ctx = PopulationContext.resolve(tx, variant, pop);
                }

                int[] acArr = ArrayUtil.toIntArray(variant.getProperty("ac"));
                int[] anArr = ArrayUtil.toIntArray(variant.getProperty("an"));
                double[] afArr = ArrayUtil.toDoubleArray(variant.getProperty("af"));
                int[] hetArr = ArrayUtil.toIntArray(variant.getProperty("het_count"));

                int ac = acArr[ctx.index];
                int an = anArr[ctx.index];
                double af = afArr[ctx.index];
                int hetCount = hetArr[ctx.index];

                if (an < 2) continue;

                nVariants++;
                piSum += VectorOps.piPerSite(af, an);
                heSum += 2.0 * af * (1.0 - af);

                int nDiploid = an / 2;
                hoSum += (double) hetCount / nDiploid;

                if (ac > 0 && ac < an) {
                    nSegregating++;
                }
            }
        } finally {
            result.close();
        }

        DiversityResult out = new DiversityResult();
        out.n_variants = nVariants;
        out.n_segregating = nSegregating;

        if (nVariants == 0 || ctx == null) {
            return Stream.of(out);
        }

        double L = nVariants;
        out.pi = piSum / L;
        out.het_exp = heSum / L;
        out.het_obs = hoSum / L;
        out.fis = out.het_exp > 0 ? 1.0 - out.het_obs / out.het_exp : 0.0;
        out.theta_w = ctx.a_n > 0 ? (nSegregating / ctx.a_n) / L : 0.0;
        out.tajima_d = TajimaD.compute(piSum, nSegregating, ctx);

        return Stream.of(out);
    }

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
}
