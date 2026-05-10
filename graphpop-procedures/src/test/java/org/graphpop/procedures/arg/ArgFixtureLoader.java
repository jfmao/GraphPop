package org.graphpop.procedures.arg;

import com.fasterxml.jackson.databind.JsonNode;
import org.neo4j.driver.Driver;
import org.neo4j.driver.Session;

import java.util.HashMap;
import java.util.Map;

/**
 * Shared helper for the M6 integration tests: loads the
 * {@code egrm_fixture_20samples_arg.json} ARG topology into a fresh
 * embedded Neo4j harness as {@code :ARGRun}, {@code :TreeNode},
 * {@code :Sample}, {@code :PARENT_OF}, {@code :REPRESENTS}, and
 * {@code :MUTATED_ON} nodes/edges so the procedures can run.
 */
public final class ArgFixtureLoader {

    private ArgFixtureLoader() {}

    public static void loadInto(Driver driver, JsonNode arg, String runId) {
        try (Session session = driver.session()) {
            session.run(
                "CREATE (:ARGRun {runId: $rid, source: 'msprime', "
              + "n_samples: $ns, sequence_length: $sl})",
                Map.of("rid", runId,
                       "ns", arg.get("n_samples").asInt(),
                       "sl", arg.get("sequence_length").asLong()));

            for (JsonNode s : arg.get("samples")) {
                int sn = s.asInt();
                session.run(
                    "CREATE (:Sample {sampleId: $sid, packed_index: $pi})",
                    Map.of("sid", "hap_" + sn, "pi", sn));
            }
            for (JsonNode n : arg.get("nodes")) {
                Map<String, Object> p = new HashMap<>();
                p.put("tid", runId + ":" + n.get("id").asInt());
                p.put("rid", runId);
                p.put("nid", n.get("id").asInt());
                p.put("t", n.get("time").asDouble());
                p.put("is", n.get("is_sample").asBoolean());
                p.put("f", n.get("flags").asLong());
                session.run(
                    "CREATE (:TreeNode {treeNodeId: $tid, runId: $rid, "
                  + "nodeId: $nid, time: $t, is_sample: $is, flags: $f})", p);
            }
            for (JsonNode e : arg.get("edges")) {
                Map<String, Object> p = new HashMap<>();
                p.put("pid", runId + ":" + e.get("parent").asInt());
                p.put("cid", runId + ":" + e.get("child").asInt());
                p.put("rid", runId);
                p.put("s", e.get("start").asLong());
                p.put("en", e.get("end").asLong());
                session.run(
                    "MATCH (p:TreeNode {treeNodeId: $pid}), "
                  + "(c:TreeNode {treeNodeId: $cid}) "
                  + "CREATE (p)-[:PARENT_OF {runId: $rid, "
                  + "start: $s, end: $en}]->(c)", p);
            }
            for (JsonNode s : arg.get("samples")) {
                int sn = s.asInt();
                session.run(
                    "MATCH (n:TreeNode {treeNodeId: $tid}), "
                  + "(s:Sample {sampleId: $sid}) "
                  + "CREATE (n)-[:REPRESENTS {runId: $rid, "
                  + "haplotype: 0}]->(s)",
                    Map.of("tid", runId + ":" + sn,
                           "sid", "hap_" + sn,
                           "rid", runId));
            }
            // :Variant + :MUTATED_ON for allele-age tests.
            JsonNode muts = arg.get("mutations");
            if (muts != null) {
                int idx = 0;
                for (JsonNode m : muts) {
                    String varId = runId + ":mut:" + idx;
                    Map<String, Object> p = new HashMap<>();
                    p.put("vid", varId);
                    p.put("pos", m.get("position").asLong());
                    p.put("ds", m.get("derived_state").asText());
                    session.run(
                        "CREATE (:Variant {variantId: $vid, "
                      + "position: $pos, derived_state: $ds})", p);

                    Map<String, Object> e = new HashMap<>();
                    e.put("vid", varId);
                    e.put("tid", runId + ":" + m.get("child_node_id").asInt());
                    e.put("rid", runId);
                    e.put("ds", m.get("derived_state").asText());
                    e.put("pid", m.get("parent_node_id").asInt());
                    session.run(
                        "MATCH (v:Variant {variantId: $vid}), "
                      + "(t:TreeNode {treeNodeId: $tid}) "
                      + "CREATE (v)-[:MUTATED_ON {runId: $rid, "
                      + "parent_node_id: $pid, derived_state: $ds}]->(t)", e);
                    idx++;
                }
            }
        }
    }
}
