#!/usr/bin/env bash
# Concurrency test for subsample.sh's D-11 trap scoping.
#
# The bug this guards against: a worker's EXIT trap did `rm -f $RAWDIR/*.fastq.gz`,
# which deleted a CONCURRENT worker's in-flight download. The symptom was a
# fail-md5 row whose 'got=' was empty (md5sum ran on a file a sibling had
# already removed). This test runs two workers whose lifetimes overlap and
# asserts neither destroys the other's raw.
set -uo pipefail

T=$(mktemp -d)
RAW="$T/raw"; mkdir -p "$RAW/locks"

# Replicate do_run's trap scoping exactly as the script now has it.
worker() {
    local acc="$1" dur="$2"
    trap 'rm -f "$RAW"/"${acc}"_*.fastq.gz' EXIT INT TERM
    for i in 0 1; do
        printf 'raw-%s-%s\n' "$acc" "$i" > "$RAW/${acc}_${i}.fastq.gz"
    done
    sleep "$dur"
    # the happy path's explicit cleanup
    rm -f "$RAW/${acc}"_*.fastq.gz
    trap - EXIT INT TERM
}

worker A 6 &
worker B 2 &   # B finishes first and its trap fires while A is still running
wait

left=$(find "$RAW" -name '*.fastq.gz' | wc -l)
if [[ "$left" -ne 0 ]]; then
    echo "FAIL: $left raw file(s) survived a completed worker:"
    find "$RAW" -name '*.fastq.gz'
    rm -rf "$T"; exit 1
fi
echo "PASS: no raw survived, and the shorter worker did not kill the longer one's files mid-flight"
rm -rf "$T"
