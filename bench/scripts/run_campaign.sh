#!/usr/bin/env bash
# C1-C3 driver: run the full chain over every unit, then aggregate (PLAN-ML-LOCAL C).
#
# Per dataset, for each of the 3 variants:
#   fastcheck --json + --profile (the fitter)  ->  fastsim  ->  the sweep
# Then:
#   aggregate.py   -> results/pilot.tsv, results/winners.tsv   (C1)
#   closed_loop.py -> results/closed_loop.tsv                  (C2, GATE 2)
#   fastp_compare  -> results/fastp-comparison.tsv              (C3)
#
# Usage:
#   run_campaign.sh --units /mnt/Fast/fastsuite-scratch/units \
#                   --results /mnt/Fast/fastsuite-scratch/results \
#                   [--jobs 6] [--variants 3] [--configs 60]
#
# Resumable: each (dataset, variant) writes a status file; finished ones are
# skipped, so the driver can be re-run freely.

set -uo pipefail

UNITS=""
RESULTS=""
JOBS=6
VARIANTS=3
CONFIGS=60
THREADS=4

while [[ $# -gt 0 ]]; do
    case "$1" in
        --units)    UNITS="$2"; shift 2 ;;
        --results)  RESULTS="$2"; shift 2 ;;
        --jobs)     JOBS="$2"; shift 2 ;;
        --variants) VARIANTS="$2"; shift 2 ;;
        --configs)  CONFIGS="$2"; shift 2 ;;
        --threads)  THREADS="$2"; shift 2 ;;
        *)          echo "run_campaign: unknown arg $1" >&2; exit 1 ;;
    esac
done

[[ -n "$UNITS" && -n "$RESULTS" ]] || { echo "need --units --results" >&2; exit 1; }

SCRIPTS="$(cd "$(dirname "$0")" && pwd)"
BIN="$HOME/FastSuite/tools"
WORK="$RESULTS/../work"
mkdir -p "$RESULTS" "$WORK" "$RESULTS/status"

ROWS="$RESULTS/sweep_rows.tsv"
[[ -s "$ROWS" ]] || # Same columns as sweep.py's ROW_FIELDS, in the same order.
printf 'domain\trun_accession\tvariant\tconfig_id\tmin_run\tmercy\tphred_sat\tgate_shape\tmin_len\tmean_q\tadapter_err\tmin_overlap\tr2_reverse\tinterior\tkeep_precision\tkeep_recall\ttrim_mae\tover_trim_rate\tunder_trim_rate\tfalse_keep_rate\tdrop_accuracy\tkept\tshould_keep\tscore\tseconds\tread_len\tmin_len_effective\n' > "$ROWS"

# ---------------------------------------------------------------- the workers
run_one() {
    local unit="$1" variant="$2"
    local domain acc tag status
    domain=$(basename "$(dirname "$unit")")
    acc=$(basename "$unit")
    tag="${domain}_${acc}.v${variant}"
    status="$RESULTS/status/$tag.ok"

    if [[ -s "$status" ]]; then
        return 0
    fi
    THREADS="$THREADS" CONFIGS="$CONFIGS" VARIANT="$variant" \
        bash "$SCRIPTS/chain.sh" "$unit" "$WORK/$domain/$acc" "$ROWS" \
        >> "$RESULTS/status/$tag.log" 2>&1
    if [[ $? -eq 0 ]]; then
        date -Iseconds > "$status"
    else
        echo "run_campaign: $tag FAILED (see status/$tag.log)" >&2
        touch "$RESULTS/status/$tag.failed"
    fi
}
export -f run_one
export SCRIPTS RESULTS WORK ROWS CONFIGS THREADS

# dispatch: every (unit, variant) is independent.
# A unit lives at $UNITS/<domain>/<acc>.R1.fastq.gz (flat per accession, not a
# per-accession directory), so UNIT is the path with the .R1.fastq.gz suffix
# stripped — that is what chain.sh consumes.
for r1 in "$UNITS"/*/*.R1.fastq.gz; do
    [ -e "$r1" ] || continue
    unit="${r1%.R1.fastq.gz}"
    [ -s "$unit.R2.fastq.gz" ] || continue   # the C chain is paired-only
    for v in $(seq 0 $(( VARIANTS - 1 ))); do
        run_one "$unit" "$v" &
        while [[ "$(jobs -rp | wc -l)" -ge "$JOBS" ]]; do sleep 2; done
    done
done
wait

n_ok=$(find "$RESULTS/status" -name '*.ok' | wc -l)
n_fail=$(find "$RESULTS/status" -name '*.failed' | wc -l)
echo "run_campaign: chains done — $n_ok ok, $n_fail failed"

# ---------------------------------------------------------------- C1 aggregate
python3 "$SCRIPTS/aggregate.py" --rows "$ROWS" --outdir "$RESULTS" || exit 1

# ---------------------------------------------------------------- C2 the loop
python3 "$SCRIPTS/closed_loop.py" --fastcheck "$BIN/fastcheck/target/release/fastcheck" \
    --work "$WORK" --row "$RESULTS/closed_loop.tsv" --threads "$THREADS" || exit 1

# ---------------------------------------------------------------- C3 vs fastp
FASTP="$(command -v fastp || true)"
if [[ -n "$FASTP" ]]; then
    python3 "$SCRIPTS/fastp_compare.py" \
        --fastcheck "$BIN/fastcheck/target/release/fastcheck" \
        --fastprep "$BIN/fastprep/target/release/fastprep" \
        --fastp "$FASTP" --units "$UNITS" \
        --winners "$RESULTS/winners.tsv" \
        --row "$RESULTS/fastp-comparison.tsv" --threads "$THREADS" || exit 1
else
    echo "run_campaign: fastp not on PATH — C3 skipped" >&2
fi

echo "run_campaign: done. Results in $RESULTS"
