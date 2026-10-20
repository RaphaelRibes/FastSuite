#!/bin/bash
# Run the FastSuite optimisation bench on the desktop.
#
# Usage: bash run_bench.sh [snakemake targets...]   (default: all)
#
# No pixi wrapper: the bench env's bin goes on PATH and every tool is invoked by
# absolute path from config.yaml, so the only thing that matters here is which
# python/snakemake runs the DAG.
set -euo pipefail

export PATH=/home/raphael/FastSuite/.pixi/envs/default/bin:$PATH
cd /home/raphael/FastSuite/bench/optim

TARGETS=("$@")
if [ ${#TARGETS[@]} -eq 0 ]; then
  TARGETS=(all)
fi

exec snakemake \
  --cores "${BENCH_CORES:-8}" \
  --rerun-triggers=mtime \
  --configfile config.yaml \
  --keep-going \
  --printshellcmds \
  "${TARGETS[@]}"
