#!/bin/bash
# Rebuild profile with deconvolved pos_q, redeploy binary, resubmit validation.
set -e
FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
cd "$FASTSUITE_ROOT"

git add -A && git -c user.name="RaphaelRibes" -c user.email="raphael.ribes@gmail.com" commit -q -m "fix: deconvolve pos_q to out-of-run curve, collapse only tail runs (closed-loop catch #3: R1 pos-0 +5.4, 3' ends -4)"

timeout 60 scp -q -o BatchMode=yes target/release/fastsim Io:/scratch/users/$USER/FastSuite/tools/fastsim/bin/
timeout 60 scp -q -o BatchMode=yes scripts/build_profile.py scripts/validate_sim.py Io:/scratch/users/$USER/FastSuite/tools/fastsim/scripts/
timeout 60 scp -q -o BatchMode=yes build_profile.sbatch validate.sbatch Io:/scratch/users/$USER/FastSuite/tools/fastsim/

# rebuild the profile (700 JSONs, ~2 min on compute node)
timeout 60 ssh -o BatchMode=yes -o ConnectTimeout=15 Io "cd $FASTSUITE_ROOT && sbatch build_profile.sbatch 2>&1 | tail -1"
echo BUILD_SUBMITTED
