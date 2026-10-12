#!/bin/bash
# Run the Q-diagnosis on the cluster.
set -e
FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
cd "$FASTSUITE_ROOT"
timeout 30 scp -q -o BatchMode=yes scripts/diag_q.py Io:/scratch/users/$USER/FastSuite/tools/fastsim/scripts/
timeout 90 ssh -o BatchMode=yes -o ConnectTimeout=15 Io "python3 $FASTSUITE_ROOT/tools/fastsim/scripts/diag_q.py"
echo DIAG_DONE
