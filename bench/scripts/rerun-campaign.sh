#!/usr/bin/env bash
# Re-run the Phase D campaign with the fixed binaries.
#
# WHY A FULL RE-RUN, NOT A PATCH
#
# Three independent things invalidated every derived artefact on disk:
#
#   1. the polyG fitter fix (a933980) changes what `fastcheck --profile`
#      emits, so every fitted profile is stale;
#   2. the N-emission work (fc1bc78) adds draws from fastsim's seeded RNG
#      when a profile carries `n_run`, shifting the whole stream;
#   3. the v1/v2 duplication: sim.v1.gt.tsv and sim.v2.gt.tsv are
#      byte-identical in 100 of 102 datasets, so variant 2 was never an
#      independent draw. The recorded 167/306 is 306 rows over ~204 real
#      measurements, with all 18 duplicate-rate failures coming from the
#      duplicated pair.
#
# `chain.sh` guards each artefact with `if [[ ! -s file ]]`, so re-running it
# over an existing work tree silently reuses all of the above. This driver
# deletes the derived tree first. The raw units (units/) and the campaign
# results (results/) are never touched.
#
# Disk: work/ is the bulk of a 245G tree. Deleting it before regenerating is
# what keeps the run inside the 117G that is free.
#
# Usage:  rerun-campaign.sh [--jobs N] [--variants N] [--dry-run]

set -uo pipefail

S=/mnt/Fast/fastsuite-scratch
UNITS="$S/units"
WORK="$S/work"
RESULTS="$S/results"
SCRIPTS=/home/raphael/FastSuite/bench/scripts

JOBS=8
VARIANTS=3
DRY=0
THREADS="${THREADS:-4}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --jobs)     JOBS="$2"; shift 2 ;;
        --variants) VARIANTS="$2"; shift 2 ;;
        --dry-run)  DRY=1; shift ;;
        *) echo "rerun-campaign: unknown arg $1" >&2; exit 1 ;;
    esac
done

SCRATCH="$(mktemp -d "$S/rerun-XXXXXX")"
ROWS="$SCRATCH/rows.tsv"
LOGS="$SCRATCH/logs"
mkdir -p "$LOGS"

# Header must match sweep.py's ROW_FIELDS EXACTLY. This was the bug that made
# the first re-run useless: run_campaign.sh's header predates the head_dropout
# axis (commit c53dbd6), so it was one field short and every column from
# head_dropout onward was shifted — `score` read `should_keep` (a count), the
# 57 configs collapsed to 2-3 distinct values, and the winner came out as
# config 0 on all 102 datasets. sweep.py only writes the header when creating
# the file, so a wrong header here silently corrupts the whole run.
# Derive it from sweep.py rather than restating it.
ROWS_HEADER=$(python3 -c "
import sys; sys.path.insert(0, '$SCRIPTS')
import sweep; print('\t'.join(sweep.ROW_FIELDS))
") || { echo "rerun-campaign: cannot read sweep.ROW_FIELDS" >&2; exit 1; }
printf '%s\n' "$ROWS_HEADER" > "$ROWS"

units=()
for r1 in "$UNITS"/*/*.R1.fastq.gz; do
    [[ -e "$r1" ]] || continue
    u="${r1%.R1.fastq.gz}"
    [[ -s "$u.R2.fastq.gz" ]] || continue
    units+=("$u")
done
echo "rerun-campaign: ${#units[@]} datasets x $VARIANTS variants, jobs=$JOBS"
echo "rerun-campaign: scratch $SCRATCH"

# Space: work/ is 237G of a 245G tree with ~117G free, and the rebuild
# reproduces ~244G (306 variants x 798M). Wiping the stale tree up front is
# what makes that fit — the per-variant delete in run_one alone would leave
# the unprocessed variants holding their old data throughout.
#
# Nothing in work/ is worth keeping: every artefact there is stale for at
# least one of the three reasons in the header. The raw units and the
# results/ summaries are untouched.
wipe_tree() {
    for u in "${units[@]}"; do
        d="$(basename "$(dirname "$u")")/$(basename "$u")"
        rm -rf "$WORK/$d"
    done
}

if [[ $DRY -eq 1 ]]; then
    echo "--- dry run: would delete and rebuild"
    for u in "${units[@]}"; do
        d="$(basename "$(dirname "$u")")/$(basename "$u")"
        for v in $(seq 0 $(( VARIANTS - 1 ))); do
            echo "  rm -rf $WORK/$d/v$v"
        done
    done
    exit 0
fi

# Preserve the old tree's provenance before deleting it: the timestamps and
# the per-variant file hashes are the evidence for the v1/v2 finding.
if [[ ! -s "$RESULTS/STALE-BEFORE-RERUN.txt" ]]; then
    {
        echo "# campaign tree state before the fixed-binary re-run"
        echo "# generated $(date -Iseconds)"
        echo "# sim.v1.gt.tsv vs sim.v2.gt.tsv, per dataset:"
        for u in "${units[@]}"; do
            d="$(basename "$(dirname "$u")")/$(basename "$u")"
            f1="$WORK/$d/v1/sim.v1.gt.tsv"; f2="$WORK/$d/v2/sim.v2.gt.tsv"
            if [[ -s "$f1" && -s "$f2" ]]; then
                h1=$(md5sum "$f1" | cut -d' ' -f1); h2=$(md5sum "$f2" | cut -d' ' -f1)
                printf '%s/%s\tv1=%s\tv2=%s\t%s\n' "$d" "${h1:0:12}" "${h2:0:12}" \
                    "$([[ $h1 = $h2 ]] && echo IDENTICAL || echo differ)"
            fi
        done
    } > "$RESULTS/STALE-BEFORE-RERUN.txt"
    echo "rerun-campaign: wrote $RESULTS/STALE-BEFORE-RERUN.txt (the v1/v2 evidence)"
fi

echo "rerun-campaign: wiping the stale work tree..."
wipe_tree
df -h "$S" | tail -1

run_one() {
    local unit="$1" variant="$2"
    local domain acc vwork tag
    domain=$(basename "$(dirname "$unit")")
    acc=$(basename "$unit")
    tag="${domain}_${acc}.v${variant}"
    vwork="$WORK/$domain/$acc/v$variant"

    # THE POINT OF THIS DRIVER: rebuild from scratch, never reuse.
    rm -rf "$vwork"
    mkdir -p "$vwork"

    THREADS="$THREADS" VARIANT="$variant" \
        bash "$SCRIPTS/chain.sh" "$unit" "$vwork" "$ROWS" \
        > "$LOGS/$tag.log" 2>&1
    local rc=$?
    if [[ $rc -ne 0 ]]; then
        echo "  FAILED $tag (see $LOGS/$tag.log)" >&2
        return 1
    fi
    echo "  ok $tag"
}
export -f run_one
export WORK ROWS LOGS SCRIPTS THREADS

fail=0
for u in "${units[@]}"; do
    for v in $(seq 0 $(( VARIANTS - 1 ))); do
        run_one "$u" "$v" &
        while [[ "$(jobs -rp | wc -l)" -ge "$JOBS" ]]; do sleep 2; done
    done
done
wait

# ---- C1 aggregate -----------------------------------------------------------
cp "$ROWS" "$RESULTS/sweep_rows.tsv"
python3 "$SCRIPTS/aggregate.py" --rows "$RESULTS/sweep_rows.tsv" --outdir "$RESULTS" || fail=1

# ---- C2 the closed loop -----------------------------------------------------
python3 "$SCRIPTS/closed_loop.py" \
    --fastcheck /home/raphael/FastSuite/tools/fastcheck/target/release/fastcheck \
    --work "$WORK" --row "$RESULTS/closed_loop.tsv" --threads "$THREADS" || fail=1

echo
echo "rerun-campaign: done (scratch $SCRATCH, failure flag $fail)"
python3 - "$RESULTS/closed_loop.tsv" <<'PY'
import csv, sys, collections
rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
ok = lambda r: r['verdict'].strip().lower() == 'pass'
n = sum(1 for r in rows if ok(r))
print(f"  closed loop: {n}/{len(rows)} = {n/len(rows)*100:.1f}%")
by = collections.defaultdict(list)
for r in rows: by[r['variant']].append(r)
for v in sorted(by):
    s = by[v]; p = sum(1 for r in s if ok(r))
    print(f"    variant {v}: {p}/{len(s)}")
# the duplication check that motivated all this
g = collections.defaultdict(dict)
for r in rows: g[(r['domain'], r['run_accession'])][r['variant']] = r
cols = [c for c in rows[0] if c not in ('domain','run_accession','variant','verdict','failing')]
ident = sum(1 for k, v in g.items()
            if '1' in v and '2' in v and all(v['1'].get(c) == v['2'].get(c) for c in cols))
print(f"    v1==v2 on every metric: {ident}/{len(g)}  (must now be 0)")
PY
exit $fail
