#!/bin/bash
# Re-simulate 20M with the offset-centering fix, then validate.
set -e
cd /home/raphael/hermes_scratch/novaseq-sim
./target/release/fastsim --pairs 20000000 --profile profiles/mandril_E.json --out local_work/sim/sim --seed 7 --no-gzip
echo "SIM_DONE"
ls -la local_work/sim/
