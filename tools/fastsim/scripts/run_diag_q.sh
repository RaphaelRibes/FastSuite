#!/bin/bash
# Run the Q-diagnosis on the cluster.
set -e
cd /home/raphael/hermes_scratch/novaseq-sim
timeout 30 scp -q -o BatchMode=yes scripts/diag_q.py Io:/scratch/users/ribesr/novaseq-sim/scripts/
timeout 90 ssh -o BatchMode=yes -o ConnectTimeout=15 Io 'python3 /scratch/users/ribesr/novaseq-sim/scripts/diag_q.py'
echo DIAG_DONE
