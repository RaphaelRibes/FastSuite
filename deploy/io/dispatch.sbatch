#!/bin/bash
#SBATCH --job-name=fs-sweep
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-long
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm/logs/dispatch-%j.out
#SBATCH --error=slurm/logs/dispatch-%j.err
#
# FastSuite PLAN-ML — phase 3 dispatcher. Holds the Snakemake controller while
# it submits the ~570 shard cells as separate SLURM jobs. It does no work
# itself: 1 CPU, 2 GB.
#
# The account is `cirad-long` and the wall clock is 7 days because this job must
# OUTLIVE the whole matrix (~5.5 h at -j 100, but queue delays and retries mean
# the controller needs days). Under `cirad-normal` it can be killed at the QoS
# max runtime and the sweep would stop mid-wave. The CELLS keep cirad-normal —
# that split is the whole point of the two accounts.
#
# Submit from the project root:
#   mkdir -p slurm/logs && sbatch slurm/dispatch.sbatch
#
# Resume, don't replay: it runs `--rerun-triggers=mtime`, so resubmitting after
# a crash/stop redoes only missing shards. To force a real rerun after a code
# change, name the rules (equals form — the option is variadic and would
# otherwise swallow the targets):
#   sbatch slurm/dispatch.sbatch --forcerun sweep_shard
set -euo pipefail

WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$WORKDIR"
mkdir -p slurm/logs

say() { echo "[$(date +%H:%M:%S)] $*"; }
die() { echo "DISPATCH FAIL: $*" >&2; exit 1; }

say "dispatcher job ${SLURM_JOB_ID:-none} on $(hostname)"
say "workdir    $WORKDIR"

# --- environment: by directory, never pixi -------------------------------------
# `pixi install` (and even `pixi shell-hook`) blocks forever on this cluster:
# the rattler package-cache lock is held by stale processes. A populated env
# directory is fully usable on its own.
BENCH_ENV="${BENCH_ENV:-$WORKDIR/.pixi/envs/bench}"
BUILD_ENV="${BUILD_ENV:-$WORKDIR/.pixi/envs/build}"
[[ -x "$BENCH_ENV/bin/snakemake" ]] || die "missing $BENCH_ENV"
export PATH="$BENCH_ENV/bin:$PATH"
[[ -x "$BUILD_ENV/bin/cargo" ]] && export PATH="$BUILD_ENV/bin:$PATH"
say "snakemake  $(snakemake --version)"
say "python     $(python3 -V)"
command -v cargo >/dev/null || die "cargo missing"
say "cargo      $(cargo --version)"

# --- build once, here, so every cell runs the identical binary ------------------
# A rule that rebuilds mid-matrix could change the binary under the sweep and
# invalidate a mixed run. fastdedup/fastmaps/fasthelper are deliberately NOT
# built: fastdedup needs rust 1.85 and the cluster build feature pins 1.83, and
# PLAN-ML uses only these three.
FS="$WORKDIR/FastSuite"
for t in fastsim fastprep fastcheck; do
    if grep -rEq "/[^ ]*(novaseq-bench|nb-staging|fs-staging)/[^ ]*\.pixi/envs" "$FS/tools/$t/target" 2>/dev/null \
       && ! grep -rq "$BENCH_ENV" "$FS/tools/$t/target" 2>/dev/null; then
        say "stale cmake cache in $t/target — cleaning"
        cargo clean --manifest-path "$FS/tools/$t/Cargo.toml"
    fi
    cargo build --release --manifest-path "$FS/tools/$t/Cargo.toml" || die "$t build failed"
done
export FASTSIM_BIN="$FS/tools/fastsim/target/release/fastsim"
export FASTPREP_BIN="$FS/tools/fastprep/target/release/fastprep"
export FASTCHECK_BIN="$FS/tools/fastcheck/target/release/fastcheck"
say "binaries: $($FASTPREP_BIN --version) $($FASTCHECK_BIN --version)"

# --- run -------------------------------------------------------------------------
cd "$WORKDIR/bench"
say "bench      $PWD"
df -h . | tail -1 | awk '{print "[disk] " $4 " free"}'

# --nolock: the controller is the only writer; a stale lock from a killed
# dispatcher would otherwise block resubmission.
# Concurrency is capped by the profile (jobs: 100, and a summed
# cpus_per_task: 100 ceiling over concurrent jobs, so the workflow cannot
# monopolise the partition).
snakemake --profile slurm --nolock --rerun-triggers=mtime "$@"

say "finished"
