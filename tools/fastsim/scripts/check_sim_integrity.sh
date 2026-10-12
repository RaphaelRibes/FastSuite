#!/bin/bash
# Check sim output integrity + gzip stream completion.
FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
cd "$FASTSUITE_ROOT"
echo "--- files ---"
ls -la local_work/sim/
echo "--- gzip integrity (gzip -t, may take a moment) ---"
for f in local_work/sim/sim_R1.fastq.gz local_work/sim/sim_R2.fastq.gz; do
  if gzip -t "$f" 2>/dev/null; then
    echo "OK: $f"
  else
    echo "TRUNCATED: $f"
  fi
done
echo "--- gt tail ---"
tail -2 local_work/sim/sim.gt.tsv 2>/dev/null
