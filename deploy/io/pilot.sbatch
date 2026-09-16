#!/bin/bash
#SBATCH --job-name=fs-pilot
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=48:00:00
#SBATCH --output=slurm/logs/pilot-%j.out
#SBATCH --error=slurm/logs/pilot-%j.err
#
# FastSuite PLAN-ML — phase 0 pilot. One job, minutes, and it settles every
# unknown that cannot be checked from a workstation:
#
#   1. does the cluster's rust (1.83 in the build feature) build the three
#      tools PLAN-ML needs? (fastdedup wants 1.85 and is NOT built here)
#   2. does a compute node reach the outside world (ENA for the corpus,
#      crates.io for cargo)?
#   3. does the env copy actually work (snakemake 9 + polars, no pixi)
#   4. one real dataset end to end: fastsim -> fastprep -> fastcheck -> score
#   5. the shard's real wall time and peak RSS on cluster hardware
#
# Submit from the project root (slurm/logs must exist — sbatch creates the
# files, not the directory):
#   mkdir -p slurm/logs && sbatch slurm/pilot.sbatch
#
# Same form as the interactive path (ssh Io + srun --time=01:00:00) so the
# pilot can also be run by hand for debugging.
#
# 1 CPU because every step is single-threaded (PLAN-ML.md §1: -jN buys wall time
# and costs core-hours). 1 GB because cargo/cmake need the most.
set -euo pipefail

WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$WORKDIR"
mkdir -p slurm/logs results scratch/pilot

say() { echo "[$(date +%H:%M:%S)] $*"; }
die() { echo "PILOT FAIL: $*" >&2; exit 1; }
ok()  { echo "  ok   $*"; }

# No /usr/bin/time dependency (not guaranteed on the node): sample VmHWM from
# /proc, the same method used for the local measurements PLAN-ML quotes, so the
# cluster numbers are directly comparable.
peak_run() {                       # peak_run <label> <cmd...>
    local label="$1"; shift
    local t0 t1 hi=0 v pid rc
    t0=$(date +%s.%N)
    "$@" >/dev/null 2>&1 &
    pid=$!
    while kill -0 "$pid" 2>/dev/null; do
        v=$(awk '/VmHWM/{print $2}' "/proc/$pid/status" 2>/dev/null)
        [[ -n "$v" && "$v" -gt "$hi" ]] && hi=$v
        sleep 0.05
    done
    wait "$pid"; rc=$?
    t1=$(date +%s.%N)
    awk -v l="$label" -v k="$hi" -v a="$t0" -v b="$t1" -v r="$rc" \
        'BEGIN{printf "  %-9s wall %6.2f s   peak RSS %6.1f MB   exit %d\n", l, b-a, k/1024, r}'
    return $rc
}

say "pilot on $(hostname), job ${SLURM_JOB_ID:-none}"
say "workdir $WORKDIR"
df -h . | tail -1 | awk '{print "  disk: " $4 " free"}'

# --- 1. environment -------------------------------------------------------------
# The env is a copied directory: `pixi install` deadlocks cluster-wide on the
# rattler cache lock, so pixi is never invoked here.
BENCH_ENV="${BENCH_ENV:-$WORKDIR/.pixi/envs/bench}"
BUILD_ENV="${BUILD_ENV:-$WORKDIR/.pixi/envs/build}"
[[ -x "$BENCH_ENV/bin/snakemake" ]] || die "missing $BENCH_ENV (copy it from novaseq-bench)"
export PATH="$BENCH_ENV/bin:$PATH"
say "env $BENCH_ENV"
say "snakemake $(snakemake --version)"
python3 -c "import polars, sys; print('  ok   python', sys.version.split()[0], 'polars', polars.__version__)" \
    || die "python/polars missing"
[[ -x "$BUILD_ENV/bin/cargo" ]] && PATH="$BUILD_ENV/bin:$PATH"
command -v cargo >/dev/null || die "cargo not found (build env on PATH?)"
say "cargo $(cargo --version)"
rustc --version
ldd --version | head -1

# --- 2. egress from a COMPUTE node ---------------------------------------------
# crates.io is implied by the existing dispatcher pattern; ENA is not.
ena=$(curl -sI -m 15 https://ftp.sra.ebi.ac.uk/ 2>/dev/null | head -1 || true)
if [[ -n "$ena" ]]; then ok "ENA reachable ($ena)"; else echo "  !!   ENA NOT reachable from a compute node"; fi
crates=$(curl -sI -m 15 https://static.crates.io/ 2>/dev/null | head -1 || true)
if [[ -n "$crates" ]]; then ok "crates.io reachable"; else echo "  !!   crates.io NOT reachable"; fi

# --- 3. build the three tools this pipeline uses --------------------------------
# Never ship a workstation binary to the cluster: the local glibc is newer than
# Io's and fastprep links zlib-ng built from source.
FS="$WORKDIR/FastSuite"
for t in fastsim fastprep fastcheck; do
    # Stale-CMake-cache guard (from the proven novaseq-bench dispatcher): cargo
    # re-passes CMAKE_C_COMPILER on every run but reads CMAKE_AR/CMAKE_RANLIB
    # from the cached CMakeCache.txt, so an env that ever lived elsewhere dies
    # at the archive step with "no such file or directory".
    if grep -rEq "/[^ ]*(novaseq-bench|nb-staging|fs-staging)/[^ ]*\.pixi/envs" "$FS/tools/$t/target" 2>/dev/null \
       && ! grep -rq "$BENCH_ENV" "$FS/tools/$t/target" 2>/dev/null; then
        say "  stale cmake cache in $t/target — cleaning"
        cargo clean --manifest-path "$FS/tools/$t/Cargo.toml"
    fi
    say "building $t"
    cargo build --release --manifest-path "$FS/tools/$t/Cargo.toml" \
        || die "$t failed to build (if the error names a rust version, that is the 1.83 pin)"
    [[ -x "$FS/tools/$t/target/release/$t" ]] || die "$t produced no binary"
done
export FASTSIM_BIN="$FS/tools/fastsim/target/release/fastsim"
export FASTPREP_BIN="$FS/tools/fastprep/target/release/fastprep"
export FASTCHECK_BIN="$FS/tools/fastcheck/target/release/fastcheck"
ok "binaries built: $($FASTPREP_BIN --version), $($FASTCHECK_BIN --version)"

# --- 4. one dataset, end to end (100k pairs = the sweep work unit) --------------
W="$WORKDIR/scratch/pilot"
say "smoke: fastsim 100k pairs -> fastprep -> fastcheck -> tier-1 score"
peak_run fastsim  "$FASTSIM_BIN" --pairs 100000 \
    --profile "$FS/tools/fastsim/profiles/test_small.json" --out "$W/sim" --seed 1 \
    || die "fastsim failed"
peak_run fastprep "$FASTPREP_BIN" -1 "$W/sim_R1.fastq.gz" -2 "$W/sim_R2.fastq.gz" \
    -o "$W/out_R1.fastq.gz" -p "$W/out_R2.fastq.gz" \
    --adapter-preset nextera --metadata "$W/pxmt.bin" \
    --json "$W/fastprep.json" -j 1 \
    || die "fastprep failed"
peak_run fastcheck "$FASTCHECK_BIN" "$W/out_R1.fastq.gz" "$W/out_R2.fastq.gz" \
    --json "$W/fastcheck.json" -q \
    || die "fastcheck failed"

# The only Python of the run: PXMT decisions vs the simulator's ground truth.
peak_run score python3 bench/scripts/score_gt.py --gt "$W/sim.gt.tsv" \
    --pxmt "$W/pxmt.bin" --weights config/score_weights.toml --json "$W/score.json" \
    || die "scoring failed"
cat "$W/score.json"

# --- 5. verdict ------------------------------------------------------------------
echo
say "PILOT PASS — record before scaling: the wall/RSS numbers above are what make"
say "the 256 MB per cell and the ~52 min per shard real instead of projected."
