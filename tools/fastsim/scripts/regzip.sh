#!/bin/bash
# Re-gzip the plain sim FASTQs (overwrite truncated gz), keep originals.
FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
cd "$FASTSUITE_ROOT/tools/fastsim/local_work/sim"
gzip -kf sim_R1.fastq
gzip -kf sim_R2.fastq
echo "GZIP_DONE"
ls -la sim_R1.fastq.gz sim_R2.fastq.gz
