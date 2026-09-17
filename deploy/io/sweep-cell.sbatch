#!/bin/bash
#SBATCH --job-name=fs-cell
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00          # one cell fits the 1 h cap (600 cfg × 100 k pairs ≈ 10 min)
#SBATCH --output=slurm/logs/cell-%A_%a.out
#SBATCH --error=slurm/logs/cell-%A_%a.err
#
# FastSuite campaign cell: one dataset × the 600-config reference grid, scored
# against the simulator's ground truth (steps 1+2 of the validation).
#
# Invoked by sweep-array.sbatch with SLURM_ARRAY_TASK_ID → one dataset row.
set -euo pipefail
WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$WORKDIR"
mkdir -p slurm/logs results/sweep

RUN=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" deploy/io/corpus-manifest.tsv | cut -f1)
say() { echo "[$(date +%H:%M:%S)] cell #${SLURM_ARRAY_TASK_ID} $RUN: $*"; }
die() { echo "CELL FAIL ($RUN): $*" >&2; exit 1; }
say "workdir $WORKDIR"

# --- env by directory (pixi deadlocks on Io) -----------------------------------
BUILD_ENV="${BUILD_ENV:-$WORKDIR/.pixi/envs/build}"
[[ -x "$BUILD_ENV/bin/cargo" ]] || die "missing $BUILD_ENV/bin/cargo"
export PATH="$BUILD_ENV/bin:$PATH"
say "cargo $(cargo --version)"

# --- build once per cell; identical binary for the whole cell -------------------
FS="$WORKDIR/FastSuite"
for t in fastsim fastprep; do
    cargo build --release --manifest-path "$FS/tools/$t/Cargo.toml" || die "$t build failed"
done

# --- the sweep: simulate from the measured profile, score configs vs gt ----------
export FASTSIM_BIN="$FS/tools/fastsim/target/release/fastsim"
export FASTPREP_BIN="$FS/tools/fastprep/target/release/fastprep"
say "simulating + sweeping $RUN (600 configs)"
python3 bench/scripts/resweep_reference.py \
    --run "$RUN" \
    --out "results/sweep/${RUN}.tsv" \
    --pairs 100000
say "done -> results/sweep/${RUN}.tsv"
