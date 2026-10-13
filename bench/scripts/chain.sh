#!/usr/bin/env bash
# C1 — the full chain, one dataset per worker (PLAN-ML-LOCAL C1).
#
#   fastcheck --json            the real unit's QC                    (C1)
#   fastcheck --profile         the fitter IS fastcheck --profile     (A1, Rust)
#   fastsim --profile           simulate a matching dataset w/ ground truth
#   fastprep sweep              N configs x 3 variants, scored vs gt.tsv
#   winner                      the best config per dataset            (C1)
#
# A2 is satisfied by construction: the recovery test passed 8/10 on this
# desktop with the Task-10 insert fix (censored_share 52.75% vs true 55.77%),
# so the fitter's insert block flows --insert-hist -> --profile -> fastsim.
#
# Usage:  chain.sh <unit_dir_prefix> <work_dir> <out_row.tsv>
#   e.g. chain.sh /mnt/.../units/human_gut/DRR425474 /mnt/.../work/human_gut/DRR425474 row.tsv
#
# Env knobs: THREADS (default 4), PAIRS (default 100000), VARIANT (default 0),
#            CONFIGS (default 60), SIM_PAIRS (default 100000).
#
# Idempotent: every artefact is written under work_dir and the row is appended
# only on success. Re-running recomputes; a partial dir is simply redone.

set -uo pipefail

UNIT="$1"     # e.g. .../units/human_gut/DRR425474  (files: DRR425474.R1.fastq.gz, .R2, .sha256)
WORK="$2"
ROW="$3"      # append one TSV row per config evaluated

THREADS="${THREADS:-4}"
PAIRS="${PAIRS:-100000}"
VARIANT="${VARIANT:-0}"
CONFIGS="${CONFIGS:-60}"
SIM_PAIRS="${SIM_PAIRS:-100000}"

BIN="${BIN:-$HOME/FastSuite/tools}"
FASTCHECK="$BIN/fastcheck/target/release/fastcheck"
FASTSIM="$BIN/fastsim/target/release/fastsim"
FASTPREP="$BIN/fastprep/target/release/fastprep"

[[ -x "$FASTCHECK" && -x "$FASTSIM" && -x "$FASTPREP" ]] || { echo "chain: missing a binary under $BIN" >&2; exit 1; }
[[ -s "$UNIT.R1.fastq.gz" ]] || { echo "chain: no unit at $UNIT" >&2; exit 1; }

DOMAIN=$(basename "$(dirname "$UNIT")")
ACC=$(basename "$UNIT")
mkdir -p "$WORK"
cd "$WORK" || exit 1

# ---------------------------------------------------------------- 1. real QC
# The report IS the fitter's input; --insert-hist makes the profile's insert
# distribution measured rather than invented (20000 pairs sampled).
if [[ ! -s real.report.json ]]; then
    if [[ -s "$UNIT.R2.fastq.gz" ]]; then
        "$FASTCHECK" --json real.report.json --profile real.profile.json \
            --insert-hist --insert-pairs 20000 \
            -j "$THREADS" -- "$UNIT.R1.fastq.gz" "$UNIT.R2.fastq.gz" \
            > real.check.log 2>&1 || { echo "chain: fastcheck failed on $ACC" >&2; tail -5 real.check.log >&2; exit 1; }
    else
        echo "chain: unit $ACC is single-end; the C chain is paired-only" >&2
        exit 1
    fi
fi

# The profile must carry a measured insert block or the simulation is invented.
python3 - "$WORK/real.profile.json" <<'PY' || exit 1
import json, sys
p = json.load(open(sys.argv[1]))
ins = p.get("insert")
if not ins:
    print("chain: profile has no insert block", file=sys.stderr); sys.exit(1)
cdf = ins.get("cdf") or ins.get("cum") or []
if len(cdf) < 2:
    print("chain: profile insert block is empty", file=sys.stderr); sys.exit(1)
print(f"chain: insert cdf {len(cdf)} bins, pairs_sampled={ins.get('pairs_sampled','?')}")
PY

# ------------------------------------------------------- 2. simulate from it
# A variant changes ONLY the non-identifiable fields (A4): the seed, so the
# three variants are three independent draws from the same fitted profile —
# the sweep then measures a config's mean score across draws, which is the
# honest estimator of "does this config trim this data well".
SEED=$(( 42 + VARIANT * 1000 ))
if [[ ! -s "sim.v$VARIANT.gt.tsv" ]]; then
    "$FASTSIM" --pairs "$SIM_PAIRS" --profile real.profile.json \
        --out "sim.v$VARIANT" --seed "$SEED" --force \
        > "sim.v$VARIANT.log" 2>&1
    rc=$?
    # fastsim prints a citation line on stderr at the end of every run: EXPECTED.
    if [[ $rc -ne 0 || ! -s "sim.v$VARIANT.gt.tsv" ]]; then
        echo "chain: fastsim failed on $ACC variant $VARIANT" >&2
        tail -20 "sim.v$VARIANT.log" >&2; exit 1
    fi
fi

# ------------------------------------------------------------- 3. the sweep
# The grid lives in sweep.py (bench/scripts/sweep.py): 60 configs x 3 variants.
# Each config is run once per variant; the row carries the variant id.
SIM_R1="sim.v$VARIANT_R1.fastq.gz"
SIM_R2="sim.v$VARIANT_R2.fastq.gz"
python3 "$(dirname "$0")/sweep.py" \
    --acc "$ACC" --domain "$DOMAIN" --variant "$VARIANT" \
    --fastprep "$FASTPREP" \
    --r1 "$SIM_R1" --r2 "$SIM_R2" \
    --gt "sim.v$VARIANT.gt.tsv" \
    --work "$WORK" --row "$ROW" --threads "$THREADS" --configs "$CONFIGS" \
    > "sweep.v$VARIANT.log" 2>&1
rc=$?
if [[ $rc -ne 0 ]]; then
    echo "chain: sweep failed on $ACC v$VARIANT" >&2
    tail -30 "sweep.v$VARIANT.log" >&2; exit 1
fi

echo "chain: $ACC v$VARIANT done"
