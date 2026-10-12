#!/bin/bash
# Inspect local_work state and run validation if sim files exist.
set -e
FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
cd "$FASTSUITE_ROOT"
echo "--- local_work tree ---"
ls local_work/
echo "--- venv? ---"
ls local_work/venv/bin/python 2>/dev/null && echo VENV_OK || echo NO_VENV
echo "--- validation ---"
ls -la local_work/sim/val* 2>/dev/null || echo "no val files"
