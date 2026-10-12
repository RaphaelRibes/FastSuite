#!/bin/bash
# Local parse check: the REAL measured mandril profile against the fixed binary.
set -e
FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
cd "$FASTSUITE_ROOT"

echo "=== parse + small run with the real profile ==="
./target/release/fastsim \
  --pairs 100000 \
  --profile profiles/mandril_E.json \
  --out /tmp/nvsq_real \
  --seed 7 \
  --no-gzip
echo SIM_DONE
