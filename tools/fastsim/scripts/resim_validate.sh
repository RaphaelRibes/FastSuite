#!/bin/bash
# Re-run the 20M sim (previous run was killed -> truncated gzip), then validate.
set -e
cd /home/raphael/hermes_scratch/novaseq-sim
rm -f local_work/sim/sim_*.fastq.gz local_work/sim/sim.gt.tsv
./target/release/fastsim --pairs 20000000 --profile profiles/mandril_E.json --out local_work/sim/sim --seed 7
echo "SIM_DONE"
ls -la local_work/sim/
/usr/bin/python3 scripts/validate_sim.py profiles/mandril_E.json \
  local_work/sim/sim_R1.fastq.gz local_work/sim/sim_R2.fastq.gz \
  local_work/sim/sim.gt.tsv local_work/sim/val
echo "VALIDATION_EXIT=$?"
