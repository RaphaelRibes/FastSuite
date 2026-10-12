#!/bin/bash
# Rebuild profile with first-base-aware head deconvolution + qual_noise,
# rebuild binary, resim 20M, validate.
set -e
FASTSUITE_ROOT="${FASTSUITE_ROOT:-$HOME/hermes_scratch/FastSuite}"
[ -d "$FASTSUITE_ROOT" ] || FASTSUITE_ROOT="/scratch/users/$USER/FastSuite"
cd "$FASTSUITE_ROOT"
python3 scripts/build_profile.py local_work/stats profiles/mandril_E.json
echo "PROFILE_REBUILT"
cargo fmt
cargo test 2>&1 | grep -E 'test result|panicked|FAILED' | head -6
cargo clippy --all-targets 2>&1 | grep -cE '^warning: [a-z]' || true
cargo build --release 2>&1 | tail -1
./target/release/fastsim --pairs 20000000 --profile profiles/mandril_E.json --out local_work/sim/sim --seed 7 --no-gzip
echo "SIM_DONE"
