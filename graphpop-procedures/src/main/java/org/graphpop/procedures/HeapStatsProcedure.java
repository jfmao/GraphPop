package org.graphpop.procedures;

import org.neo4j.procedure.Description;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Procedure;

import java.util.stream.Stream;

/**
 * JVM heap statistics for benchmarking memory pressure.
 *
 * <p>Two procedures:
 * <ul>
 *   <li>{@code graphpop.heap_stats()} — point-in-time snapshot (no GC)</li>
 *   <li>{@code graphpop.gc()} — forces GC then returns heap stats
 *       (gives a clean baseline for measuring memory delta)</li>
 * </ul>
 * </p>
 */
public class HeapStatsProcedure {

    public static class HeapResult {
        public long total_mb;
        public long used_mb;
        public long free_mb;
        public long max_mb;

        public HeapResult() {}
    }

    @Procedure(name = "graphpop.heap_stats", mode = Mode.READ)
    @Description("Report JVM heap memory usage for the Neo4j database.")
    public Stream<HeapResult> heapStats() {
        return Stream.of(snapshot());
    }

    @Procedure(name = "graphpop.gc", mode = Mode.READ)
    @Description("Trigger JVM garbage collection.")
    public Stream<HeapResult> gc() {
        System.gc();
        System.gc();   // two rounds for best effort
        return Stream.of(snapshot());
    }

    private static HeapResult snapshot() {
        Runtime rt = Runtime.getRuntime();
        HeapResult r = new HeapResult();
        r.total_mb = rt.totalMemory() / (1024 * 1024);
        r.used_mb = (rt.totalMemory() - rt.freeMemory()) / (1024 * 1024);
        r.free_mb = rt.freeMemory() / (1024 * 1024);
        r.max_mb = rt.maxMemory() / (1024 * 1024);
        return r;
    }
}
