#!/bin/bash
#SBATCH --job-name=fdedup-bench
#SBATCH --output=logs/slurm/fdedup-bench-%j.out
#SBATCH --error=logs/slurm/fdedup-bench-%j.err
#SBATCH --time=06:00:00          # Max 6h for the master job
#SBATCH --mem=4G                 # Master job needs little memory
#SBATCH --cpus-per-task=1        # Master job is just coordinating
#SBATCH --ntasks=1
#SBATCH --partition=defq         # <-- CHANGE to your SLURM partition
#SBATCH --account=my_account     # <-- CHANGE to your account
#
# FastDedup Snakemake SLURM Launcher
# ====================================
# Submits Snakemake with the SLURM executor. Snakemake then submits each
# of the 324 benchmark jobs as individual SLURM jobs (1-6 cpus, 16 GB max).
#
# Budget: 64 cores max / 16 GB per task
#
# Usage:
#   sbatch slurm.sh                           # Full benchmark
#   DRY_RUN=1 sbatch slurm.sh                 # Dry-run (check DAG)
#   RESUME=1 sbatch slurm.sh                  # Resume after interruption
#   JOBS=20 sbatch slurm.sh                   # Limit concurrent submits
#   PARTITION=fast qos=high sbatch slurm.sh   # Tune SLURM params
#
# After completion:
#   pixi run python scripts/visualize.py      # Generate graphs

set -euo pipefail

BENCH_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${BENCH_DIR}"

# Create log dir
mkdir -p logs/slurm

# Env overrides
DRY_RUN="${DRY_RUN:-0}"
RESUME="${RESUME:-0}"
JOBS="${JOBS:-100}"

SNAKEMAKE=".pixi/envs/default/bin/snakemake"

echo "============================================"
echo " FastDedup Benchmark — SLURM"
echo "============================================"
echo " Profile:      slurm-profile/"
echo " Max jobs:     ${JOBS}"
echo " Dry run:      ${DRY_RUN}"
echo " Resume:       ${RESUME}"
echo " Partition:    ${SLURM_JOB_PARTITION:-N/A}"
echo " Account:      ${SLURM_JOB_ACCOUNT:-N/A}"
echo " Start:        $(date)"
echo " Work dir:     ${BENCH_DIR}"
echo "============================================"

# Build Snakemake args
ARGS=(
    --profile slurm-profile/
    --jobs "${JOBS}"
    --keep-going
    --latency-wait 60
    --rerun-incomplete
)

if [ "${DRY_RUN}" = "1" ]; then
    ARGS+=("-n")
fi

# Run Snakemake (this submits each benchmark rule to SLURM)
${SNAKEMAKE} "${ARGS[@]}"

EXIT_CODE=$?

echo ""
echo "============================================"
echo " Benchmark finished at $(date)"
echo " Exit code: ${EXIT_CODE}"
echo "============================================"

if [ "${DRY_RUN}" = "1" ]; then
    echo " This was a dry-run. Remove DRY_RUN=1 to execute."
else
    echo ""
    echo " Results: benchmark.tsv"
    echo " Graphs:  pixi run python scripts/visualize.py"
fi

exit ${EXIT_CODE}
