#!/bin/bash
# Re-simulate 20M with the offset-centering fix, then validate.
set -e
FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
cd "$FASTSUITE_ROOT"
./target/release/fastsim --pairs 20000000 --profile profiles/mandril_E.json --out local_work/sim/sim --seed 7 --no-gzip
echo "SIM_DONE"
ls -la local_work/sim/
