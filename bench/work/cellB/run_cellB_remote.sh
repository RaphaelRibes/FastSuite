#!/bin/bash
# Cell B, remote: build the synthetic cohort and run the with/without trim report.
#
# The cohort generator is the repo's tools/fastsim; the report is the repo's
# bench/scripts/trim_report.py. Both are taken from ~/FastSuite (the mirror).
#
# Layout produced (the layout trim_report.py discovers):
#   <root>/cohort/work/synthetic/cellB-pe-overlap/v0/sim.v0.{gt.tsv,_R1.fastq.gz,_R2.fastq.gz}
#   <root>/cohort/results/winners.tsv
#   <out>/cellB-{with-pe,without-pe}[-permate|-pe|-permate-pe].tsv
set -e

ROOT=/mnt/Fast/fastsuite-scratch/bench-work/20261019-pe-overlap
CELLB=$ROOT/cellB
OUT=$ROOT/results
SCRATCH=$ROOT/trim-cellB
MIRROR=~/FastSuite

V0=$CELLB/cohort/work/synthetic/cellB-pe-overlap/v0
mkdir -p "$V0" "$CELLB/cohort/results" "$OUT" "$SCRATCH"

echo "=== 1. simulate (seed 4242, deterministic) ==="
"$MIRROR/tools/fastsim/target/release/fastsim" \
    --pairs 400000 \
    --profile "$MIRROR/tools/fastsim/profiles/cellB-pe-overlap.json" \
    --out "$V0/sim.v0" \
    --seed 4242 \
    --force \
    --drop-rate 0 --dup-rate 0

echo "=== 2. winners.tsv (the published winner config) ==="
printf 'domain\trun_accession\tconfig_id\tmin_run\tmercy\tphred_sat\tgate_shape\tmin_len\tmean_q\tadapter_err\tmin_overlap\tr2_reverse\tinterior\tmin_len_effective\nsynthetic\tcellB-pe-overlap\tcellB\t12\t3\t20\t1.0\t0.9\t25.0\t0.2\t14\t0\t0\t0.9L\n' \
    > "$CELLB/cohort/results/winners.tsv"

echo "=== 3. audit the cohort (G-runs? readthrough?) ==="
python3 "$MIRROR/bench/work/cellB/audit_cohort.py" \
    "$V0/sim.v0.gt.tsv" "$V0/sim.v0_R1.fastq.gz" "$V0/sim.v0_R2.fastq.gz" || true

echo "=== 4a. report WITH --pe-overlap ==="
FASTSCRATCH="$SCRATCH/with-pe" \
FASTPREP_VARIANTS='[["fp-pe",["--pe-overlap"]]]' \
python3 "$MIRROR/bench/scripts/trim_report.py" --mode decomp \
    --work "$CELLB/cohort/work" --results "$CELLB/cohort/results" \
    --out "$OUT/cellB-with-pe.tsv" --tol 0 2>&1 | tail -40

echo "=== 4b. report WITHOUT (the baseline, no FASTPREP_VARIANTS) ==="
FASTSCRATCH="$SCRATCH/without-pe" \
FASTPREP_VARIANTS='' \
python3 "$MIRROR/bench/scripts/trim_report.py" --mode decomp \
    --work "$CELLB/cohort/work" --results "$CELLB/cohort/results" \
    --out "$OUT/cellB-without-pe.tsv" --tol 0 2>&1 | tail -40

echo "=== 5. outputs ==="
ls -la "$OUT"/cellB-*.tsv
