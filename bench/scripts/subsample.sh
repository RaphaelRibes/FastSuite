#!/usr/bin/env bash
# B2 — deterministic subsample of ENA runs to 100 k-pair units (PLAN-ML-LOCAL B2).
#
# For each manifest row: download the run's FASTQ(s) from the ENA ftp paths,
# md5-verify every file against the manifest, subsample deterministically to
# exactly PAIRS pairs, sha256 the unit, delete the raw. Resumable: a unit whose
# sha256 already exists is skipped.
#
# D-11: the raw NEVER persists beyond the verification. It lives in RAW_DIR and
# is rm -f'd by an EXIT trap in every worker, so a crash or a kill leaves no raw
# behind either. The derived units live in UNITS_DIR on /mnt/Fast/fastsuite-scratch
# and are gitignored; no data is ever committed to the repo.
#
# Determinism: stride = read_count / PAIRS and the selected records are
# 1, 1+stride, 1+2*stride, ... so the unit is a pure function of (run, PAIRS) —
# no RNG, no clock, byte-identical on every rerun.
#
# Usage:
#   subsample.sh --manifest real.tsv \
#                --units /mnt/Fast/fastsuite-scratch/units \
#                --raw   /mnt/Fast/fastsuite-scratch/raw \
#                --log   /mnt/Fast/fastsuite-scratch/logs/b2.tsv \
#                [--domains generalised,human_gut,virome] [--limit 10] \
#                [--jobs 4] [--pairs 100000]
#
# --jobs workers run concurrently, each holding at most one raw run (rolling
# window). A claim lock under RAW_DIR/../locks makes a concurrent rerun safe.

set -uo pipefail

PAIRS=100000
JOBS=4
DOMAINS=""
LIMIT=0
MANIFEST=""
UNITS=""
RAWDIR=""
LOG=""
FORCE=0

die() { printf 'subsample: %s\n' "$*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest) MANIFEST="$2"; shift 2 ;;
        --units)    UNITS="$2"; shift 2 ;;
        --raw)      RAWDIR="$2"; shift 2 ;;
        --log)      LOG="$2"; shift 2 ;;
        --domains)  DOMAINS="$2"; shift 2 ;;
        --limit)    LIMIT="$2"; shift 2 ;;
        --jobs)     JOBS="$2"; shift 2 ;;
        --pairs)    PAIRS="$2"; shift 2 ;;
        --force)    FORCE=1; shift ;;
        -h|--help)  sed -n '2,28p' "$0"; exit 0 ;;
        *)          die "unknown argument: $1" ;;
    esac
done

[[ -n "$MANIFEST" && -n "$UNITS" && -n "$RAWDIR" && -n "$LOG" ]] || die "need --manifest --units --raw --log"
[[ "$PAIRS" -gt 0 ]] || die "--pairs must be > 0"
[[ "$JOBS" -ge 1 ]] || die "--jobs must be >= 1"

mkdir -p "$UNITS" "$RAWDIR" "$(dirname "$LOG")" "$RAWDIR/../locks"

if [[ ! -s "$LOG" ]]; then
    printf 'domain\trun_accession\tstatus\tunit_sha256\tr1_reads\tr2_reads\traw_bytes\tunit_bytes\tseconds\tnote\n' > "$LOG"
fi

# --------------------------------------------------------------------------
# claim / release — a noclobber mkdir is the atomic test-and-set
# --------------------------------------------------------------------------
claim()   { mkdir "$RAWDIR/../locks/$1.claim" 2>/dev/null; }
unclaim() { rmdir "$RAWDIR/../locks/$1.claim" 2>/dev/null; }
claimed() { [[ -d "$RAWDIR/../locks/$1.claim" ]]; }

record() {  # record DOMAIN ACC STATUS SHA R1 R2 RAWBYTES UNITBYTES SECS NOTE
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$@" >> "$LOG"
}

# --------------------------------------------------------------------------
# subsample_stream IN OUT — every stride-th record, exactly PAIRS of them.
#
# awk emits the 4 lines of records 1, 1+stride, 1+2*stride, ... (a record is 4
# lines, so the stride in LINES is 4*stride). head -n then caps the output at
# PAIRS records and SIGPIPEs the upstream the moment the unit is complete, so
# a 60 M-read run stops being read at ~PAIRS*stride records rather than at the
# end of the file. The selection is a pure function of (stride, PAIRS).
# --------------------------------------------------------------------------
subsample_stream() {
    local in="$1" out="$2" stride="$3"
    local linestride=$(( 4 * stride ))
    zcat -f -- "$in" \
        | awk -v m="$linestride" '(NR-1) % m < 4' \
        | head -n $(( 4 * PAIRS )) \
        | gzip -1 > "$out"
}

count_records() { zcat -f -- "$1" | wc -l | awk '{print int($1/4)}'; }

# --------------------------------------------------------------------------
# do_run DOMAIN ACC FTP MD5 READS BYTES
# --------------------------------------------------------------------------
do_run() {
    local domain="$1" acc="$2" ftp="$3" md5="$4" reads="$5" bytes="$6"
    local udir="$UNITS/$domain/$acc"
    local t0 t1
    t0=$(date +%s)

    # The unit dir must exist before any writer touches it. It normally already
    # does (claim/skip-existing ran first), but the very first worker for a
    # domain creates it, and a concurrent `rm -rf` of a stale domain dir would
    # otherwise leave this worker writing into a path that no longer exists.
    mkdir -p "$UNITS/$domain" || {
        record "$domain" "$acc" "fail-mkdir" - - - - - - "$UNITS/$domain"
        unclaim "$acc"; return 1
    }

    # --- resumable: a sha256 beside the unit means this run is finished ---
    if [[ "$FORCE" -ne 1 && -s "$udir.sha256" ]]; then
        record "$domain" "$acc" "skip-existing" "$(cut -d' ' -f1 < "$udir.sha256")" - - - - - resumed
        return 0
    fi
    if claimed "$acc"; then
        record "$domain" "$acc" "skip-claimed" - - - - - - "another worker owns it"
        return 0
    fi
    claim "$acc" || { record "$domain" "$acc" "skip-claimed" - - - - - - "race"; return 0; }

    # D-11 backstop: however this worker ends, ITS OWN raws go away.
    # Scoped to this accession's files only: a bare *.fastq.gz glob would delete
    # another concurrent worker's in-flight download, which is exactly what
    # happened before this was scoped (a fail-md5 row whose 'got=' was empty
    # because the file had already been removed by a sibling's trap).
    MY_RAWS="$RAWDIR/${acc}_"*.fastq.gz
    trap 'rm -f "$RAWDIR"/"${acc}"_*.fastq.gz' EXIT INT TERM

    local -a urls md5s byte_list
    IFS=';' read -r -a urls <<< "$ftp"
    IFS=';' read -r -a md5s <<< "$md5"
    IFS=';' read -r -a byte_list <<< "$bytes"
    if [[ "${#urls[@]}" -ne "${#md5s[@]}" || "${#urls[@]}" -eq 0 ]]; then
        record "$domain" "$acc" "fail-manifest" - - - - - - "ftp/md5 count mismatch"
        unclaim "$acc"; return 1
    fi

    # ---- download + md5 verify, one file at a time (rolling window) ----
    local i url want raw got total_raw=0
    for i in "${!urls[@]}"; do
        url="${urls[$i]}"; want="${md5s[$i]}"
        [[ "$url" =~ ^ftp:// ]] || url="ftp://$url"
        raw="$RAWDIR/${acc}_${i}.fastq.gz"
        if ! curl -sS --fail --retry 4 --retry-delay 5 --connect-timeout 30 \
                  --max-time 5400 -o "$raw" "$url"; then
            record "$domain" "$acc" "fail-download" - - - "$total_raw" - - "$url"
            rm -f "$raw"; unclaim "$acc"; return 1
        fi
        # A partial file is a failed download, not a checksum failure: say so
        # distinctly, because the two have different fixes (retry vs bad manifest).
        local want_bytes="${byte_list[$i]}"
        if [[ -n "$want_bytes" && "$want_bytes" =~ ^[0-9]+$ ]]; then
            local got_bytes
            got_bytes=$(stat -c%s "$raw")
            if [[ "$got_bytes" -ne "$want_bytes" ]]; then
                record "$domain" "$acc" "fail-truncated" - - - "$total_raw" - - \
                    "want=${want_bytes}B got=${got_bytes}B"
                rm -f "$raw"; unclaim "$acc"; return 1
            fi
        fi
        got=$(md5sum "$raw" | cut -d' ' -f1)
        if [[ "$got" != "$want" ]]; then
            record "$domain" "$acc" "fail-md5" - - - "$total_raw" - - "want=$want got=$got"
            rm -f "$raw"; unclaim "$acc"; return 1
        fi
        total_raw=$(( total_raw + $(stat -c%s "$raw") ))
    done

    # ---- deterministic subsample ----
    local stride=$(( reads / PAIRS ))
    [[ "$stride" -lt 1 ]] && stride=1
    local r1=0 r2=0 ubytes sha note=""

    subsample_stream "$RAWDIR/${acc}_0.fastq.gz" "$udir.R1.fastq.gz" "$stride"
    r1=$(count_records "$udir.R1.fastq.gz")

    if [[ "${#urls[@]}" -ge 2 ]]; then
        subsample_stream "$RAWDIR/${acc}_1.fastq.gz" "$udir.R2.fastq.gz" "$stride"
        r2=$(count_records "$udir.R2.fastq.gz")
    else
        note="single-end"
    fi

    # ---- acceptance: exactly PAIRS on both mates, mates equal ----
    if [[ "$r1" -ne "$PAIRS" ]] || { [[ "${#urls[@]}" -ge 2 ]] && [[ "$r2" -ne "$PAIRS" ]]; }; then
        record "$domain" "$acc" "fail-count" - "$r1" "$r2" "$total_raw" - - "want=$PAIRS each"
        rm -f "$udir".R*.fastq.gz; unclaim "$acc"; return 1
    fi

    # ---- sha256 of the unit ----
    # R1-then-R2 concatenation with a mate label folded in, so a swapped pair
    # yields a different digest and the digest is order-sensitive.
    if [[ "${#urls[@]}" -ge 2 ]]; then
        sha=$( { printf 'R1\n'; cat "$udir.R1.fastq.gz"; printf 'R2\n'; cat "$udir.R2.fastq.gz"; } \
               | sha256sum | cut -d' ' -f1 )
        ubytes=$(( $(stat -c%s "$udir.R1.fastq.gz") + $(stat -c%s "$udir.R2.fastq.gz") ))
        printf '%s\n' "$sha" > "$udir.sha256"
        record "$domain" "$acc" "ok-pe" "$sha" "$r1" "$r2" "$total_raw" "$ubytes" $(( $(date +%s) - t0 )) "$note"
    else
        sha=$( { printf 'R1\n'; cat "$udir.R1.fastq.gz"; } | sha256sum | cut -d' ' -f1 )
        ubytes=$(stat -c%s "$udir.R1.fastq.gz")
        printf '%s\n' "$sha" > "$udir.sha256"
        record "$domain" "$acc" "ok-se" "$sha" "$r1" - "$total_raw" "$ubytes" $(( $(date +%s) - t0 )) "$note"
    fi

    rm -f "$RAWDIR/${acc}"_*.fastq.gz
    unclaim "$acc"
    trap - EXIT INT TERM
    return 0
}

# --------------------------------------------------------------------------
main() {
    # Parse the manifest in python (it is CRLF and carries ';'-joined
    # multi-file rows) and emit one TSV job line per selected run.
    local jobs_file
    jobs_file=$(mktemp)
    python3 - "$MANIFEST" "$DOMAINS" "$LIMIT" > "$jobs_file" <<'PY'
import csv, sys
manifest, domains, limit = sys.argv[1], sys.argv[2], int(sys.argv[3])
want = [d for d in domains.split(",") if d]
rows = list(csv.DictReader(open(manifest, newline=""), delimiter="\t"))
per = {}
for r in rows:
    d = r["domain"]
    if want and d not in want:
        continue
    if limit and per.get(d, 0) >= limit:
        continue
    per[d] = per.get(d, 0) + 1
    urls = [u for u in r["fastq_ftp"].split(";") if u]
    md5s = [m for m in r["fastq_md5"].split(";") if m]
    bytez = [b for b in r.get("fastq_bytes", "").split(";") if b]
    if len(urls) != len(md5s) or not urls:
        continue
    if len(bytez) != len(urls):
        bytez = [""] * len(urls)
    # fields: domain, acc, ftp, md5, read_count, per-file byte counts
    print("\t".join([d, r["run_accession"], ";".join(urls), ";".join(md5s),
                     r["read_count"], ";".join(bytez)]))
PY
    [[ -s "$jobs_file" ]] || { rm -f "$jobs_file"; die "no rows selected (domains='$DOMAINS' limit=$LIMIT)"; }

    printf 'subsample: %d run(s), pairs=%d, jobs=%d, units=%s\n' \
        "$(wc -l < "$jobs_file")" "$PAIRS" "$JOBS" "$UNITS"

    while IFS=$'\t' read -r domain acc ftp md5 reads bytes; do
        do_run "$domain" "$acc" "$ftp" "$md5" "$reads" "$bytes" >> "$LOG.err" 2>&1 &
        # job-slot throttler: at most JOBS live workers
        while [[ "$(jobs -rp | wc -l)" -ge "$JOBS" ]]; do sleep 2; done
    done < "$jobs_file"
    wait
    rm -f "$jobs_file"

    local ok fail
    ok=$(awk -F'\t' 'NR>1 && $3 ~ /^ok/'  "$LOG" | wc -l)
    fail=$(awk -F'\t' 'NR>1 && $3 ~ /^fail/' "$LOG" | wc -l)
    printf 'subsample: done — %d ok, %d failed (log: %s)\n' "$ok" "$fail" "$LOG"

    # D-11: assert no raw survives anywhere under scratch.
    local leftover
    leftover=$(find "$RAWDIR" -name '*.fastq.gz' 2>/dev/null | wc -l)
    if [[ "$leftover" -ne 0 ]]; then
        printf 'D-11 VIOLATION: %d raw file(s) under %s — removing\n' "$leftover" "$RAWDIR" >&2
        find "$RAWDIR" -name '*.fastq.gz' -delete
        return 1
    fi
    printf 'D-11: no raw file under %s\n' "$RAWDIR"
    return 0
}

main "$@"
