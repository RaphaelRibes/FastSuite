#!/bin/bash
# Cell B, remote: FRESH end-to-end verification.
#
# Same as run_cellB_remote.sh but writes to a fresh --out dir (the report refuses
# to overwrite an existing TSV), so this proves the whole chain reproduces from
# scratch: simulate -> cohort layout -> winners.tsv -> audit -> two reports.
set -e

ROOT=/mnt/Fast/fastsuite-scratch/bench-work/20261019-pe-overlap
CELLB=$ROOT/cellB
OUT=$ROOT/results-fresh
SCRATCH=$ROOT/trim-cellB-fresh
MIRROR=~/FastSuite

V0=$CELLB/cohort/work/synthetic/cellB-pe-overlap/v0
mkdir -p "$V0" "$CELLB/cohort/results" "$OUT" "$SCRATCH"

echo "=== 1. simulate (seed 4242) ==="
BEFORE=$(sha256sum "$V0/sim.v0.gt.tsv" | cut -d' ' -f1)
"$MIRROR/tools/fastsim/target/release/fastsim" \
    --pairs 400000 \
    --profile "$MIRROR/tools/fastsim/profiles/cellB-pe-overlap.json" \
    --out "$V0/sim.v0" \
    --seed 4242 --force --drop-rate 0 --dup-rate 0
AFTER=$(sha256sum "$V0/sim.v0.gt.tsv" | cut -d' ' -f1)
echo "  GT sha256 before $BEFORE"
echo "  GT sha256 after  $AFTER"
[ "$BEFORE" = "$AFTER" ] && echo "  REPRODUCED: YES" || { echo "  REPRODUCED: NO"; exit 1; }

echo "=== 2. winners.tsv ==="
printf 'domain\trun_accession\tconfig_id\tmin_run\tmercy\tphred_sat\tgate_shape\tmin_len\tmean_q\tadapter_err\tmin_overlap\tr2_reverse\tinterior\tmin_len_effective\nsynthetic\tcellB-pe-overlap\tcellB\t12\t3\t20\t1.0\t0.9\t25.0\t0.2\t14\t0\t0\t0.9L\n' \
    > "$CELLB/cohort/results/winners.tsv"

echo "=== 3. audit ==="
python3 "$MIRROR/bench/work/cellB/audit_cohort.py" \
    "$V0/sim.v0.gt.tsv" "$V0/sim.v0_R1.fastq.gz" "$V0/sim.v0_R2.fastq.gz"

echo "=== 4a. WITH --pe-overlap ==="
FASTSCRATCH="$SCRATCH/with-pe" \
FASTPREP_VARIANTS='[["fp-pe",["--pe-overlap"]]]' \
python3 "$MIRROR/bench/scripts/trim_report.py" --mode decomp \
    --work "$CELLB/cohort/work" --results "$CELLB/cohort/results" \
    --out "$OUT/cellB-with-pe.tsv" --tol 0

echo "=== 4b. WITHOUT ==="
FASTSCRATCH="$SCRATCH/without-pe" \
FASTPREP_VARIANTS='' \
python3 "$MIRROR/bench/scripts/trim_report.py" --mode decomp \
    --work "$CELLB/cohort/work" --results "$CELLB/cohort/results" \
    --out "$OUT/cellB-without-pe.tsv" --tol 0

echo "=== 5. the two headline tables ==="
for f in cellB-with-pe.tsv cellB-with-pe-pe.tsv cellB-without-pe.tsv; do
    echo "--- $f ---"
    cat "$OUT/$f"
done
