#!/bin/bash
# Run remaining human_interpret.py subcommands sequentially.
# PBS is already running separately — this handles everything after it.
# Each subcommand has resume support, so this script is safe to re-run.

set -euo pipefail
cd /mnt/data/GraphPop

log() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*"; }

# Wait for PBS to finish (check if process exists)
PBS_PID=$(ps aux | grep "human_interpret.py pbs" | grep -v grep | awk '{print $2}' | head -1)
if [ -n "$PBS_PID" ]; then
    log "Waiting for PBS (PID $PBS_PID) to finish..."
    while kill -0 "$PBS_PID" 2>/dev/null; do
        sleep 60
    done
    log "PBS finished."
fi

for cmd in hwscan roh_hmm daf_enrichment gscan report; do
    log "=== Starting: $cmd ==="
    python scripts/human_interpret.py "$cmd" 2>&1 | tee -a "human_interpret_${cmd}.log"
    log "=== Finished: $cmd ==="
done

log "=== All interpretation subcommands complete ==="
