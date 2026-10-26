#!/usr/bin/env bash
# Re-derive the closed loop (C2) on the post-2026-10-21 build.
#
# Why: the campaign's closed_loop.tsv was produced by binaries whose
#   * fastsim ignored profile.gc (emitted the --gc default 0.49),
#   * fastcheck measured insert length with a leftmost-anchor bias,
#   * fastcheck's profile fitter estimated polyG/quality from the per-position
#     G fraction instead of the exact histograms it held.
# So every number in that table describes the previous semantics. This driver
# regenerates, per (dataset, variant): the real QC + fitted profile, the
# simulation from that profile, and the simulation's QC — then lets
# closed_loop.py score them with the unchanged thresholds.
#
# Layout: $WORKROOT/<domain>/<acc>/v<variant>/ and $RESULTS/.
#
# Committed to the repo as bench/scripts/rerun-closed-loop.sh: a number whose
# producer is not in the repository cannot be re-derived (the same trap as
# doc/TRIM-ACCURACY.md citing a script that never existed).
# The old work/ and results/ trees are left untouched.

set -uo pipefail

# Defaults derived from this file's location, so a clone runs without editing:
# the binaries live at <repo>/tools/*/target/release/, the scorer at
# <repo>/bench/scripts/closed_loop.py.
SCRIPTS=${SCRIPTS:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
BIN=${BIN:-$(cd "$SCRIPTS/../.." && pwd)/tools}
FASTCHECK=$BIN/fastcheck/target/release/fastcheck
FASTSIM=$BIN/fastsim/target/release/fastsim
UNITS=${UNITS:-/mnt/Fast/fastsuite-scratch/units}
WORKROOT=${WORKROOT:-/mnt/Fast/fastsuite-scratch/work-gt2026}
RESULTS=${RESULTS:-/mnt/Fast/fastsuite-scratch/results-gt2026}
JOBS=${JOBS:-8}
THREADS=${THREADS:-2}
SIM_PAIRS=${SIM_PAIRS:-100000}

mkdir -p "$WORKROOT" "$RESULTS/logs"

one() {
    local unit=$1 v=$2
    local domain acc tag vwork
    domain=$(basename "$(dirname "$unit")")
    acc=$(basename "$unit")
    tag=${domain}_${acc}.v${v}
    vwork=$WORKROOT/$domain/$acc/v$v
    mkdir -p "$vwork" || return 1
    if [[ -s "$vwork/sim.v$v.report.json" && -s "$vwork/sim.v$v.gt.tsv" ]]; then
        echo "skip $tag"
        return 0
    fi
    cd "$vwork" || return 1

    # 1. real QC + the fitted profile, current build
    "$FASTCHECK" --json real.report.json --profile real.profile.json \
        --insert-hist --insert-pairs 20000 -j "$THREADS" -- \
        "$unit.R1.fastq.gz" "$unit.R2.fastq.gz" > real.check.log 2>&1 \
        || { echo "FAIL fastcheck-real $tag"; return 1; }

    # 2. fastsim knobs from the profile's own read length / insert floor
    local knobs READ_LEN INSERT_MIN
    knobs=$(python3 - <<'PY'
import json, sys
p = json.load(open("real.profile.json"))
rl = int(p.get("read_len") or 150)
ins = p.get("insert") or {}
lens = sorted(int(k) for k in ins)
imin = max(1, min(lens[0] if lens else 1, rl))
print(rl, imin)
PY
) || { echo "FAIL profile-knobs $tag"; return 1; }
    READ_LEN=$(echo "$knobs" | cut -d' ' -f1)
    INSERT_MIN=$(echo "$knobs" | cut -d' ' -f2)

    # 3. simulate: same seed scheme as the campaign (42 + 1000*variant),
    #    drop-rate 0 so the sim carries the profile's own insert distribution
    "$FASTSIM" --pairs "$SIM_PAIRS" --profile real.profile.json \
        --out "sim.v$v" --seed $(( 42 + v * 1000 )) \
        --read-len "$READ_LEN" --insert-min "$INSERT_MIN" --drop-rate 0 \
        --force > "sim.v$v.log" 2>&1 \
        || { echo "FAIL fastsim $tag"; return 1; }
    [[ -s "sim.v$v.gt.tsv" ]] || { echo "FAIL fastsim-gt $tag"; return 1; }

    # 4. the simulation's QC (closed_loop.py would otherwise do this serially)
    "$FASTCHECK" --json "sim.v$v.report.json" -q -j "$THREADS" -- \
        "sim.v${v}_R1.fastq.gz" "sim.v${v}_R2.fastq.gz" > "sim.v$v.check.log" 2>&1 \
        || { echo "FAIL fastcheck-sim $tag"; return 1; }

    echo "ok $tag"
}

start=$(date +%s)
for r1 in "$UNITS"/*/*.R1.fastq.gz; do
    [[ -e "$r1" ]] || continue
    unit=${r1%.R1.fastq.gz}
    [[ -s "$unit.R2.fastq.gz" ]] || continue
    for v in 0 1 2; do
        one "$unit" "$v" &
        while [[ $(jobs -rp | wc -l) -ge $JOBS ]]; do sleep 1; done
    done
done
wait
end=$(date +%s)
echo "driver: $(find "$WORKROOT" -name 'sim.v*.report.json' | wc -l) sim reports in $(( end - start ))s"

python3 "$SCRIPTS/closed_loop.py" --fastcheck "$FASTCHECK" \
    --work "$WORKROOT" --row "$RESULTS/closed_loop.tsv" --threads 4 --force \
    || { echo "driver: closed_loop.py failed"; exit 1; }

python3 - "$RESULTS/closed_loop.tsv" <<'PY'
import csv, sys, statistics
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
n = len(rows)
# closed_loop.py writes PASS in upper case; comparing against "pass" made this
# summary print "0/306 comparisons pass" for a run with 183. The TSV was always
# right — this line was not, and it is the line a reader quotes.
passed = [r for r in rows if r["verdict"].strip().upper() == "PASS"]
ds = {(r["domain"], r["run_accession"]) for r in rows}
ds_pass = {(r["domain"], r["run_accession"]) for r in passed}
g = [float(r["g_frac_max_diff"]) for r in rows if r["g_frac_max_diff"]]
print(f"closed_loop: {len(passed)}/{n} comparisons pass "
      f"({100.0*len(passed)/n:.1f}%), {len(ds_pass)}/{len(ds)} datasets pass all arms")
if g:
    print(f"closed_loop: median g_frac_max_diff {statistics.median(g):.4f} (tolerance 0.02)")
fails = {}
for r in rows:
    for tok in (r["failing"] or "").split(";"):
        tok = tok.strip().split(" ")[0]
        if tok:
            fails[tok] = fails.get(tok, 0) + 1
print("closed_loop: fail counts by gate:", fails)
PY
