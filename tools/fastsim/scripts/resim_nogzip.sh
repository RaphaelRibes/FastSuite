#!/bin/bash
# Simulate 20M pairs WITHOUT gzip (fast), then validate (gzip.open handles plain).
set -e
cd /home/raphael/hermes_scratch/novaseq-sim
./target/release/fastsim --pairs 20000000 --profile profiles/mandril_E.json --out local_work/sim/sim --seed 7 --no-gzip
echo "SIM_DONE"
ls -la local_work/sim/
/usr/bin/python3 scripts/validate_sim.py profiles/mandril_E.json \
  local_work/sim/sim_R1.fastq local_work/sim/sim_R2.fastq \
  local_work/sim/sim.gt.tsv local_work/sim/val
echo "VALIDATION_EXIT=$?"
