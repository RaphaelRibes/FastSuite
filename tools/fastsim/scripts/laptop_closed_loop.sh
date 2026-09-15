#!/bin/bash
# Run the ENTIRE closed loop on the laptop, working inside
# ~/hermes_scratch/novaseq-sim/local_work/ (not /tmp):
# 1. pull 700 stats JSONs from cluster (one-time, ~100MB)
# 2. build_profile.py -> profiles/mandril_E.json (pure python, no numpy)
# 3. simulate 20M pairs with the local release binary
# 4. validate with validate_sim.py (needs numpy+matplotlib)
set -e
cd /home/raphael/hermes_scratch/novaseq-sim
WORK=local_work

# --- which python has numpy+matplotlib? ---
PY=""
for py in /usr/bin/python3 "$HOME/.pixi/bin/python" python3; do
  if [ -x "$py" ]; then
    if "$py" -c "import numpy, matplotlib" 2>/dev/null; then
      PY="$py"
      echo "numpy python: $PY"
      break
    fi
  fi
done
if [ -z "$PY" ]; then
  echo "no numpy python found; creating venv"
  python3 -m venv "$WORK/venv"
  "$WORK/venv/bin/pip" -q install numpy matplotlib 2>&1 | tail -1
  PY="$WORK/venv/bin/python"
fi
echo "PY=$PY"

# --- 1. pull stats (only if missing) ---
if [ -z "$(ls "$WORK"/stats/*.json 2>/dev/null | head -1)" ]; then
  echo "pulling 700 stats JSONs from cluster..."
  timeout 180 rsync -az -e "ssh -o BatchMode=yes" Io:/scratch/users/ribesr/mandril-polyg-2026/stats/ "$WORK/stats/"
  echo "stats: $(ls "$WORK"/stats/*.json | wc -l) files"
else
  echo "stats already present: $(ls "$WORK"/stats/*.json | wc -l) files"
fi

# --- 2. build profile ---
python3 scripts/build_profile.py "$WORK/stats" profiles/mandril_E.json

# --- 3. simulate ---
./target/release/fastsim --pairs 20000000 --profile profiles/mandril_E.json --out "$WORK/sim/sim" --seed 7
ls -la "$WORK/sim/"

# --- 4. validate ---
"$PY" scripts/validate_sim.py profiles/mandril_E.json \
  "$WORK/sim/sim_R1.fastq.gz" "$WORK/sim/sim_R2.fastq.gz" \
  "$WORK/sim/sim.gt.tsv" "$WORK/sim/val"

echo "LOCAL_CLOSED_LOOP_DONE"
