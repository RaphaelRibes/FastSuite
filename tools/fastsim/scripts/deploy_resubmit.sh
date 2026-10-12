#!/bin/bash
# Deploy fixed binary + scripts to the cluster, then resubmit validation.
set -e
FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
cd "$FASTSUITE_ROOT"
timeout 120 scp -q -o BatchMode=yes target/release/fastsim Io:/scratch/users/$USER/FastSuite/tools/fastsim/bin/
timeout 60 scp -q -o BatchMode=yes scripts/validate_sim.py Io:/scratch/users/$USER/FastSuite/tools/fastsim/scripts/
timeout 60 ssh -o BatchMode=yes -o ConnectTimeout=15 Io "cd $FASTSUITE_ROOT && sbatch validate.sbatch 2>&1 | tail -1"
echo DEPLOY_DONE
