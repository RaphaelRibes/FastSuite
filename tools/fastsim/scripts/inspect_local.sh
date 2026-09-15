#!/bin/bash
# Inspect local_work state and run validation if sim files exist.
set -e
cd /home/raphael/hermes_scratch/novaseq-sim
echo "--- local_work tree ---"
ls local_work/
echo "--- venv? ---"
ls local_work/venv/bin/python 2>/dev/null && echo VENV_OK || echo NO_VENV
echo "--- validation ---"
ls -la local_work/sim/val* 2>/dev/null || echo "no val files"
