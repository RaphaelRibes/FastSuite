#!/bin/bash
# Check gzip state after the timeout.
FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
cd "$FASTSUITE_ROOT/tools/fastsim/local_work/sim"
ls -la sim_R1.fastq.gz sim_R2.fastq.gz
echo "--- testing ---"
gzip -t sim_R1.fastq.gz 2>&1 && echo "R1 OK" || echo "R1 BAD"
gzip -t sim_R2.fastq.gz 2>&1 && echo "R2 OK" || echo "R2 BAD"
