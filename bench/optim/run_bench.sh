#!/bin/bash
# Run the FastSuite optimisation bench on the desktop.
#
# Usage: bash run_bench.sh [snakemake targets...]   (default: all)
#
# No pixi wrapper: the bench env's bin goes on PATH and every tool is invoked by
# absolute path from config.yaml, so the only thing that matters here is which
# python/snakemake runs the DAG.
#
# Two defaults are deliberately conservative, because the first campaign of this
# round got both wrong and produced 90 rows that had to be thrown away:
#
#   BENCH_CORES=1        every cell is pinned with taskset, so two cells at once
#                        queue on the same core. A `-j 1` cell then reports
#                        cpu_s/wall_s well under 1.00 and the number is a
#                        queueing artefact, not a measurement. Only raise this if
#                        the cells you are running are pinned to disjoint cores.
#
#   config-papers.yaml   the targeted config that produced the table the dossier
#                        and the six papers quote: one size, one thread count,
#                        both arms. config.yaml is the full matrix (three sizes,
#                        three thread counts) and takes hours serially; it is
#                        still available as BENCH_CONFIG=config.yaml.
set -euo pipefail

export PATH=/home/raphael/FastSuite/.pixi/envs/default/bin:$PATH
cd /home/raphael/FastSuite/bench/optim

CORES="${BENCH_CORES:-1}"
CONFIG="${BENCH_CONFIG:-config-papers.yaml}"

TARGETS=("$@")
if [ ${#TARGETS[@]} -eq 0 ]; then
  TARGETS=(all)
fi

echo "run_bench: config=${CONFIG} cores=${CORES} targets=${TARGETS[*]}"
if [ "$CORES" != "1" ]; then
  echo "run_bench: WARNING cores=${CORES} - only sound if these cells are pinned to disjoint cores"
fi

exec snakemake \
  --cores "$CORES" \
  --rerun-triggers=mtime \
  --configfile "$CONFIG" \
  --keep-going \
  --printshellcmds \
  "${TARGETS[@]}"
