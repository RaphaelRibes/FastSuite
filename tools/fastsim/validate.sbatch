#!/bin/bash
#SBATCH --job-name=fastsim-validate
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:30:00
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-normal
#SBATCH --output=logs/validate_%j.out
#SBATCH --error=logs/validate_%j.err

set -e
# Cluster-side: no $HOME here, so the root is the Io scratch checkout
# (override with FASTSUITE_ROOT if the checkout lives elsewhere).
FASTSUITE_ROOT="${FASTSUITE_ROOT:-/scratch/users/$USER/FastSuite}"
cd "$FASTSUITE_ROOT/tools/fastsim"
mkdir -p logs profiles validation

BIN="$FASTSUITE_ROOT/tools/fastsim/target/release/fastsim"
PROFILE=profiles/mandril_E.json
OUT=validation/sim_mandril

echo "=== 1. simulate 20M pairs from the measured profile ==="
"$BIN" --pairs 20000000 --profile "$PROFILE" --out "$OUT" --seed 7
ls -la "$OUT"_R1.fastq.gz "$OUT"_R2.fastq.gz "$OUT".gt.tsv

echo "=== 2. re-profile the simulated reads (closed loop) ==="
# use the deployed polyg-profiler from the mandril project
PROFILER=/scratch/users/ribesr/mandril-polyg-2026/scanner/target/release/polyg-profiler
if [ -x "$PROFILER" ]; then
  "$PROFILER" --r1 "$OUT"_R1.fastq.gz --r2 "$OUT"_R2.fastq.gz -o validation/sim_profiled.json
  echo "profiled: $(python3 -c "import json; v=json.load(open('validation/sim_profiled.json')); print(v['n_pairs'], 'pairs')")"
else
  echo "polyg-profiler not found at $PROFILER"
fi

echo "=== 3. validate simulated vs profile ==="
python3 scripts/validate_sim.py \
  "$PROFILE" \
  "$OUT"_R1.fastq.gz \
  "$OUT"_R2.fastq.gz \
  "$OUT".gt.tsv \
  validation/closed_loop

echo "=== done ==="
ls -la validation/
