#!/bin/bash
# Rebuild profile with deconvolved pos_q, redeploy binary, resubmit validation.
set -e
cd /home/raphael/hermes_scratch/novaseq-sim

git add -A && git -c user.name="RaphaelRibes" -c user.email="raphael.ribes@gmail.com" commit -q -m "fix: deconvolve pos_q to out-of-run curve, collapse only tail runs (closed-loop catch #3: R1 pos-0 +5.4, 3' ends -4)"

timeout 60 scp -q -o BatchMode=yes target/release/novaseq-sim Io:/scratch/users/ribesr/novaseq-sim/bin/
timeout 60 scp -q -o BatchMode=yes scripts/build_profile.py scripts/validate_sim.py Io:/scratch/users/ribesr/novaseq-sim/scripts/
timeout 60 scp -q -o BatchMode=yes build_profile.sbatch validate.sbatch Io:/scratch/users/ribesr/novaseq-sim/

# rebuild the profile (700 JSONs, ~2 min on compute node)
timeout 60 ssh -o BatchMode=yes -o ConnectTimeout=15 Io 'cd /scratch/users/ribesr/novaseq-sim && sbatch build_profile.sbatch 2>&1 | tail -1'
echo BUILD_SUBMITTED
