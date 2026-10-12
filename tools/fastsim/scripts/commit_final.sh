#!/bin/bash
# Commit the final quality model (cap compensation + head deconvolution).
set -e
FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
cd "$FASTSUITE_ROOT"
git add -A
git -c user.name="RaphaelRibes" -c user.email="raphael.ribes@gmail.com" commit -q -m "fix: cap-compensated noise + first-base-aware head deconvolution (validation: R1 mean|diff| 0.127, R2 0.007)"
git log --oneline | head -4
echo COMMIT_DONE
