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
CONFIGS=55  # the authored grid length (see sweep.py); the plan's "60" was never authored
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
# Derive the header from sweep.py's ROW_FIELDS. It used to be restated here and
# went stale when the head_dropout axis landed (c53dbd6): one field short, so
# every column from head_dropout onward shifted and `score` read `should_keep`.
# sweep.py writes the header only when it creates the file, so a wrong header
# here silently corrupts every row of the run.
if [[ ! -s "$ROWS" ]]; then
    python3 -c "
import sys; sys.path.insert(0, '$SCRIPTS')
import sweep; print('\t'.join(sweep.ROW_FIELDS))
" > "$ROWS" || { echo "run_campaign: cannot read sweep.ROW_FIELDS" >&2; exit 1; }
fi

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
    # each variant gets its OWN work dir: with --jobs > 1 two variants of the
    # same dataset run concurrently and would race on the shared sweep dir
    # (fastprep's .part -> rename commit gets unlinked mid-flight by the other
    # variant's cleanup — measured as "committing c41_R1.fastq.gz.part ... No
    # such file or directory"). The per-variant files (sim.vN, sweep.vN) mean
    # no output was ever actually shared.
    vwork="$WORK/$domain/$acc/v$variant"
    mkdir -p "$vwork"
    THREADS="$THREADS" CONFIGS="$CONFIGS" VARIANT="$variant" \
        bash "$SCRIPTS/chain.sh" "$unit" "$vwork" "$ROWS" \
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
