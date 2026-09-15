#!/bin/bash
# Plot head and tail polyG CDFs from the sim GT vs the measured profile.
# Usage: bash scripts/run_plot_heads_tails.sh [PROFILE_JSON] [SIM_GT] [OUT_PREFIX]
cd /home/raphael/hermes_scratch/novaseq-sim
PROFILE="${1:-profiles/mandril_E.json}"
GT="${2:-local_work/sim/sim.gt.tsv}"
PREFIX="${3:-local_work/sim/ht}"
PY=""
for py in /usr/bin/python3 "$HOME/.pixi/bin/python" python3; do
  if [ -x "$py" ]; then
    PY="$py"
    break
  fi
done
"$PY" scripts/plot_heads_tails.py "$PROFILE" "$GT" "$PREFIX"
echo "PLOT_EXIT=$?"
