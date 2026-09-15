#!/bin/bash
# Local closed-loop sanity: simulate 20k pairs from the small test profile,
# then validate the simulated reads against the profile.
set -e
cd /home/raphael/hermes_scratch/novaseq-sim
mkdir -p /tmp/nvsq_local

./target/release/novaseq-sim \
  --pairs 20000 \
  --profile profiles/test_small.json \
  --out /tmp/nvsq_local/sim \
  --seed 42 \
  --no-gzip

echo "=== validation ==="
python3 scripts/validate_sim.py \
  profiles/test_small.json \
  /tmp/nvsq_local/sim_R1.fastq \
  /tmp/nvsq_local/sim_R2.fastq \
  /tmp/nvsq_local/sim.gt.tsv \
  /tmp/nvsq_local/val
