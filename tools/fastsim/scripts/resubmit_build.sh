#!/bin/bash
# Resubmit profile build via sbatch, then poll out/err until done.
set -e
timeout 60 ssh -o BatchMode=yes -o ConnectTimeout=15 Io 'cd /scratch/users/ribesr/novaseq-sim && sbatch build_profile.sbatch 2>&1 | tail -1'
echo SUBMITTED
