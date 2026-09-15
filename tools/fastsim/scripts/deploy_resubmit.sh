#!/bin/bash
# Deploy fixed binary + scripts to the cluster, then resubmit validation.
set -e
cd /home/raphael/hermes_scratch/novaseq-sim
timeout 120 scp -q -o BatchMode=yes target/release/novaseq-sim Io:/scratch/users/ribesr/novaseq-sim/bin/
timeout 60 scp -q -o BatchMode=yes scripts/validate_sim.py Io:/scratch/users/ribesr/novaseq-sim/scripts/
timeout 60 ssh -o BatchMode=yes -o ConnectTimeout=15 Io 'cd /scratch/users/ribesr/novaseq-sim && sbatch validate.sbatch 2>&1 | tail -1'
echo DEPLOY_DONE
