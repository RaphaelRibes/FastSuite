#!/bin/bash
# stage-to-io.sh — stage FastSuite + the sweep campaign to Io (NO upload: dry-run).
#
# The campaign needs only code + profiles + config grid on the cluster — it
# SIMULATES each dataset from its measured profile, it does not ship reads
# (there is no 1.4 TB transfer; profiles are ~10-20 KB each).
#
# Rule (user): scp for transfer, never rsync. Build on the cluster, never ship
# binaries. `pixi install` deadlocks on Io — envs are used by directory.
#
# Usage: bash deploy/io/stage-to-io.sh          # print the plan (dry run)
#        bash deploy/io/stage-to-io.sh --run    # actually scp (ONLY with user OK)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"   # FastSuite/
DEST="Io:/scratch/users/ribesr/fastsuite-bench/"
RUN="${1:-}"
say(){ echo "[stage] $*"; }

# --- FastSuite code (no target/) ----------------------------------------------
CODE=(
  "FastSuite/tools/fastsim/src" "FastSuite/tools/fastsim/Cargo.toml" "FastSuite/tools/fastsim/Cargo.lock"
  "FastSuite/tools/fastprep/src" "FastSuite/tools/fastprep/Cargo.toml" "FastSuite/tools/fastprep/Cargo.lock"
  "FastSuite/tools/fastcheck/src" "FastSuite/tools/fastcheck/Cargo.toml" "FastSuite/tools/fastcheck/Cargo.lock"
  "FastSuite/deploy/io/corpus-manifest.tsv" "FastSuite/deploy/io/sweep-cell.sbatch"
)
# --- bench scripts + config grid ----------------------------------------------
BENCH=(
  "FastSuite-bench/bench/scripts/resweep_reference.py"
  "FastSuite-bench/bench/scripts/score_gt.py"
  "FastSuite-bench/bench/config/score_weights.toml"
  "FastSuite-bench/results/sweep/complete_ERR10123282.tsv"   # the 600-config grid (shared)
)
# --- measured profiles (the shape to simulate) ---------------------------------
PROFILES=()
for r in $(cut -f1 FastSuite/deploy/io/corpus-manifest.tsv | grep -v '^#' | grep -v '^run'); do
  PROFILES+=( "FastSuite-bench/data/complete/$r/profiles/profile_as-measured.json" )
done

say "plan: $(( ${#CODE[@]} + ${#BENCH[@]} + ${#PROFILES[@]} )) files to $DEST (no binaries, no reads)"
[ "$RUN" = "--run" ] || { say "dry run — pass --run to actually scp (needs user OK)"; exit 0; }

say "creating remote layout"
ssh Io "mkdir -p /scratch/users/ribesr/fastsuite-bench/{FastSuite/tools,bench/scripts,bench/config,bench/results/sweep,deploy/io}"
say "scp FastSuite code"
scp -r "${CODE[@]/#/}" "$DEST"
say "scp bench scripts + config grid"
scp -r "${BENCH[@]/#/}" "$DEST"
say "scp measured profiles"
for p in "${PROFILES[@]}"; do scp "$p" "$DEST"; done
say "done. Next (with user OK): ssh Io -> srun -> sbatch --array=0-17 sweep-cell.sbatch"
