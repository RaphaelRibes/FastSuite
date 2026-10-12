#!/bin/bash
# Commit the closed-loop validation state.
set -e
FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
cd "$FASTSUITE_ROOT"
git add -A
git -c user.name="RaphaelRibes" -c user.email="raphael.ribes@gmail.com" commit -q -m "feat: deconvolved out-of-run pos_q + laptop closed-loop validation (20M pairs, all metrics pass)"
git log --oneline | head -6
echo COMMIT_DONE
