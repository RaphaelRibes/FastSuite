#!/usr/bin/env bash
# Thread scaling of fastdedup, measured with the tool's own clock.
#
# WHY THIS EXISTS
#
# The 2026-10-26 benchmark of the fix wave needed to know whether the new
# ordering buffer (the dedup decision runs on results re-ordered by input
# ordinal) cost throughput. It does not — but the same harness showed that `-j`
# makes fastdedup SLOWER on every corpus tried, on the published 1.3.0 binary as
# well as on HEAD: 1.00x at -j 1 down to 0.09x at -j 16 on 4.8 M reads. See
# bench/results/wave2/thread-scaling.md for the numbers and the reproduction.
#
# This script is that measurement, kept so the finding can be re-checked instead
# of quoted.
#
# Usage:
#   thread-scaling.sh --head DIR --units DIR [--base DIR] [--out FILE]
#                     [--reps N] [--threads "1 2 4 8 16"]
#
#   --head   checkout whose tools/*/target/release/fastdedup is measured
#   --base   optional second checkout (e.g. a worktree at the published tag);
#            measured interleaved with --head so drift hits both equally
#   --units  a corpus directory holding *.R1.fastq.gz / *.R2.fastq.gz units
#
# HARNESS RULES (each learned the hard way — see thread-scaling.md):
#   * time in integer nanoseconds, reduce in one Python pass with LC_ALL=C
#     (`bc` is absent on the desktop and the locale uses a comma decimal);
#   * no /usr/bin/time (absent there: `-f '%M'` made every run fail in 0.7 ms);
#   * the summary line goes to STDOUT for file outputs (only `-` streams move it
#     to stderr) and `-q` suppresses it — never pass `-q` here.
set -uo pipefail
export LC_ALL=C

HEAD=""
BASE=""
UNITS=""
OUT=""
REPS=3
THREADS_LIST="1 2 4 8 16"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --head)    HEAD=$2; shift 2 ;;
        --base)    BASE=$2; shift 2 ;;
        --units)   UNITS=$2; shift 2 ;;
        --out)     OUT=$2; shift 2 ;;
        --reps)    REPS=$2; shift 2 ;;
        --threads) THREADS_LIST=$2; shift 2 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

[[ -n "$HEAD" && -n "$UNITS" ]] || { echo "need --head and --units" >&2; exit 2; }
FD_HEAD=$HEAD/tools/fastdedup/target/release/fastdedup
FD_BASE=${BASE:+$BASE/tools/fastdedup/target/release/fastdedup}
[[ -x "$FD_HEAD" ]] || { echo "no fastdedup at $FD_HEAD" >&2; exit 2; }
if [[ -n "$BASE" && ! -x "$FD_BASE" ]]; then
    echo "no fastdedup at $FD_BASE" >&2; exit 2
fi

WORK=$(mktemp -d "${TMPDIR:-/tmp}/thread-scaling.XXXXXX")
OUT=${OUT:-$WORK/thread-scaling.tsv}
RAW=$WORK/raw.tsv
: > "$RAW"
corpus=$WORK/corpus
mkdir -p "$corpus"

# One real unit, uncompressed, and eight copies of the same unit gzipped: the
# first is the clean case, the second is what a user actually has on disk.
first_r1=$(find "$UNITS" -name '*.R1.fastq.gz' | sort | head -1)
[[ -n "$first_r1" ]] || { echo "no *.R1.fastq.gz under $UNITS" >&2; exit 2; }
first_r2=${first_r1%.R1.fastq.gz}.R2.fastq.gz
zcat "$first_r1" > "$corpus/one.R1.fq"
zcat "$first_r2" > "$corpus/one.R2.fq"
for i in 1 2 3 4 5 6 7 8; do cat "$first_r1"; done > "$corpus/gz.R1.fastq.gz"
for i in 1 2 3 4 5 6 7 8; do cat "$first_r2"; done > "$corpus/gz.R2.fastq.gz"

elapsed_of() {  # "$summary line" -> seconds
    /usr/bin/python3 - "$1" <<'PY'
import re, sys
m = re.findall(r",\s*([0-9.]+)(ms|s|µs)\s*$", sys.argv[1].strip())
if not m:
    print(""); raise SystemExit
v, unit = float(m[0][0]), m[0][1]
print(f"{v/1000.0 if unit=='ms' else v/1e6 if unit=='µs' else v:.4f}")
PY
}

measure() {  # bin tag corpus j rep
    local bin=$1 tag=$2 c=$3 j=$4 rep=$5 line secs r1 r2
    case "$c" in
        one) r1=$corpus/one.R1.fq;    r2=$corpus/one.R2.fq ;;
        gz)  r1=$corpus/gz.R1.fastq.gz; r2=$corpus/gz.R2.fastq.gz ;;
        *)   echo "unknown corpus $c" >&2; return 1 ;;
    esac
    line=$("$bin" -1 "$r1" -2 "$r2" \
             -o "$WORK/o.fq" -p "$WORK/o2.fq" -j "$j" --force 2>&1 \
           | grep 'reads in' | tail -1)
    secs=$(elapsed_of "$line")
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$c" "$tag" "$j" "$rep" "$secs" "$line" >> "$RAW"
}

echo "== measured: $FD_HEAD${FD_BASE:+ and $FD_BASE}"
echo "== corpora: one=$(du -h "$corpus/one.R1.fq" | cut -f1) gz=$(du -h "$corpus/gz.R1.fastq.gz" | cut -f1)"
for c in one gz; do
    for j in $THREADS_LIST; do
        for rep in $(seq 1 "$REPS"); do
            [[ -n "$FD_BASE" ]] && measure "$FD_BASE" base "$c" "$j" "$rep"
            measure "$FD_HEAD" head "$c" "$j" "$rep"
        done
        echo "   $c j=$j done"
    done
done

cp "$RAW" "$OUT"
/usr/bin/python3 - "$OUT" <<'PY'
import sys, csv, statistics, collections
by = collections.defaultdict(list)
for r in csv.reader(open(sys.argv[1]), delimiter="\t"):
    if len(r) >= 6 and r[4]:
        by[(r[0], r[1], r[2])].append(float(r[4]))
for c in ("one", "gz"):
    tags = sorted({k[1] for k in by if k[0] == c})
    print(f"\n-- corpus {c}")
    print("   j  " + "  ".join(f"{t:>9}" for t in tags) + "   speedup vs j=1 (head)")
    base = statistics.median(by[(c, "head", "1")]) if (c, "head", "1") in by else None
    for j in sorted({k[2] for k in by if k[0] == c}, key=int):
        cells = []
        for t in tags:
            v = by.get((c, t, j))
            cells.append(f"{statistics.median(v):8.3f}s" if v else "        -")
        m = by.get((c, "head", j))
        sp = f"{base/statistics.median(m):5.2f}x" if (m and base) else "    -"
        print(f"  {j:>2}  " + "  ".join(cells) + f"   {sp}")
PY
echo "== tsv: $OUT"
