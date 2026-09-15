#!/bin/bash
# Validate the 20M sim against the deconvolved profile.
cd /home/raphael/hermes_scratch/novaseq-sim
/usr/bin/python3 scripts/validate_sim.py profiles/mandril_E.json \
  local_work/sim/sim_R1.fastq local_work/sim/sim_R2.fastq \
  local_work/sim/sim.gt.tsv local_work/sim/val
echo "VALIDATION_EXIT=$?"
ls -la local_work/sim/val* 2>/dev/null
