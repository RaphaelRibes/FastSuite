#!/bin/bash
# Install the analysis env (numpy, matplotlib, pandas) via pixi.
FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
cd "$FASTSUITE_ROOT"
/home/raphael/.pixi/bin/pixi install -e analysis > /tmp/pixi_local_install.log 2>&1
echo "PIXI_EXIT=$?"
tail -3 /tmp/pixi_local_install.log
ls .pixi/envs/analysis/bin/python3 2>/dev/null && echo PYTHON_OK
