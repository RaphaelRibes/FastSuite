#!/bin/bash
# Local parse check: the REAL measured mandril profile against the fixed binary.
set -e
cd /home/raphael/hermes_scratch/novaseq-sim

echo "=== parse + small run with the real profile ==="
./target/release/novaseq-sim \
  --pairs 100000 \
  --profile profiles/mandril_E.json \
  --out /tmp/nvsq_real \
  --seed 7 \
  --no-gzip
echo SIM_DONE
