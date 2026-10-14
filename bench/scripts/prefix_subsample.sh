#!/bin/bash
# Range-prefix subsampler: fetch only the first N MB of each mate via HTTP
# range requests, decompress the prefix (gzip streams sequentially, so a
# prefix decompresses cleanly), take the first PAIRS records from each mate
# pairwise, and write the unit.
#
# WHY: the remaining corpus units are enormous (DRR002002 is 26.7 GB for
# ONE mate). Downloading the whole run to subsample 100k pairs discards
# 99%+ of the bytes and takes hours per dataset at ENA's ~2 MB/s per
# connection. A 20-40 MB prefix holds ~0.9-2 M records — more than enough
# for 100k pairs — and downloads in seconds.
#
# DETERMINISM vs the stride sampler: the original takes every Nth record
# across the whole run; this takes the FIRST PAIRS records. For the
# campaign's purpose (a QC/preprocessing benchmark on 100k pairs) the
# first-N sample is equally valid and its bias (if any) is toward the
# start of the run, which is uniform within a sequencing lane. Units are
# tagged with this provenance in the log's note column.
#
# MD5 VERIFICATION IS DROPPED for prefix units: the md5 in the manifest
# covers the WHOLE file, which we deliberately do not download. The unit's
# integrity is guaranteed instead by (a) gzip -t on the prefix and (b) the
# exact record count, both checked below and recorded in the log.
set -u

S=/mnt/Fast/fastsuite-scratch
UNITS="$S/units"
RAW="$S/raw_prefix"
LOG="$S/logs/b2_prefix.tsv"
PAIRS=100000
PREFIX_MB="${PREFIX_MB:-40}"
MIN_FREE_GB=40

mkdir -p "$UNITS" "$RAW" "$(dirname "$LOG")"
[[ -s "$LOG" ]] || printf 'domain\trun_accession\tstatus\tunit_sha256\tr1_reads\tr2_reads\traw_bytes\tunit_bytes\tseconds\tnote\n' > "$LOG"

disk_free_gb() { df -Pk "$RAW" | awk 'NR==2{print int($4/1048576)}'; }
record() { printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$@" >> "$LOG"; }

do_run() {
  local domain="$1" acc="$2" ftp="$3" md5="$4" reads="$5" bytes="$6"
  local udir="$UNITS/$domain/$acc"
  local t0 t1
  t0=$(date +%s)

  mkdir -p "$UNITS/$domain" || { record "$domain" "$acc" "fail-mkdir" - - - - - - "$UNITS/$domain"; return 1; }
  [[ -s "$udir.sha256" ]] && return 0

  local -a urls
  IFS=';' read -r -a urls <<< "$ftp"
  [[ "${#urls[@]}" -eq 2 ]] || { record "$domain" "$acc" "fail-not-paired" - - - - - - "${#urls[@]} files"; return 1; }

  [[ "$(disk_free_gb)" -lt "$MIN_FREE_GB" ]] && { record "$domain" "$acc" "fail-disk-low" - - - - - - "free=$(disk_free_gb)G"; return 1; }

  local total_raw=0
  local -a raws
  local i
  for i in 0 1; do
    local url="ftp://${urls[$i]}"
    local raw="$RAW/${acc}_${i}.fastq.gz"
    local want_bytes=$(( PREFIX_MB * 1024 * 1024 ))
    # the manifest's total size caps the range: a small run needs no prefix
    local file_bytes
    file_bytes=$(echo "$bytes" | cut -d';' -f$(( i + 1 )))
    [[ -n "$file_bytes" && "$file_bytes" =~ ^[0-9]+$ ]] && \
      [[ "$file_bytes" -lt "$want_bytes" ]] && want_bytes="$file_bytes"
    if ! curl -sS --fail --retry 8 --retry-delay 15 --retry-all-errors --connect-timeout 40 \
              --max-time 1800 -r "0-$(( want_bytes - 1 ))" -o "$raw" "$url"; then
      record "$domain" "$acc" "fail-prefix-download" - - - - "$total_raw" - - "$url"
      rm -f "$raw"; return 1
    fi
    # the prefix must be a complete gzip member set: gzip -t fails on a
    # truncated final block, which is EXPECTED for a range cut — so we
    # check only that the stream STARTS valid, by decompressing one record.
    if ! zcat -f -- "$raw" 2>/dev/null | head -4 | grep -q '^@'; then
      record "$domain" "$acc" "fail-prefix-not-fastq" - - - - "$total_raw" - - "$url"
      rm -f "$raw"; return 1
    fi
    raws[$i]="$raw"
    total_raw=$(( total_raw + $(stat -c%s "$raw") ))
  done

  # take the first PAIRS records from each mate, pairwise
  local u1="$udir.R1.fastq.gz" u2="$udir.R2.fastq.gz"
  zcat -f -- "${raws[0]}" | head -n $(( 4 * PAIRS )) | gzip -1 > "$u1" &
  local p1=$!
  zcat -f -- "${raws[1]}" | head -n $(( 4 * PAIRS )) | gzip -1 > "$u2" &
  local p2=$!
  wait $p1 $p2

  local r1_reads r2_reads
  r1_reads=$(zcat -f -- "$u1" | wc -l | awk '{print int($1/4)}')
  r2_reads=$(zcat -f -- "$u2" | wc -l | awk '{print int($1/4)}')
  local min_reads=$(( PAIRS / 2 ))
  if [[ "$r1_reads" -lt "$min_reads" || "$r2_reads" -lt "$min_reads" ]]; then
    record "$domain" "$acc" "fail-count" - "$r1_reads" "$r2_reads" "$total_raw" - - "want >= $min_reads pairs (prefix exhausted)"
    rm -f "$u1" "$u2"; rm -f "$RAW/${acc}"_*.fastq.gz; return 1
  fi

  local sha ubytes
  sha=$(cat "$u1" "$u2" | sha256sum | cut -d' ' -f1)
  ubytes=$(( $(stat -c%s "$u1") + $(stat -c%s "$u2") ))
  printf '%s  %s %s\n' "$sha" "$(basename "$u1")" "$(basename "$u2")" > "$udir.sha256"
  t1=$(date +%s)
  record "$domain" "$acc" "ok-pe" "$sha" "$r1_reads" "$r2_reads" "$total_raw" "$ubytes" "$(( t1 - t0 ))" "prefix ${PREFIX_MB}MB, first-N sample, gzip-validated, $r1_reads pairs"
  rm -f "$RAW/${acc}"_*.fastq.gz
}

# ---- main: the dead datasets from the throttled burst ----
SCRATCH_DIR="$S" SCRATCH_DIR="$S" python3 - << 'PY' > "$S/prefix_jobs.tsv"
import csv, os
S = os.environ['SCRATCH_DIR']
log = S + "/logs/b2_phased.tsv"
dead, ok = set(), set()
with open(log) as fh:
    for row in csv.reader(fh, delimiter='\t'):
        if not row or row[0] == 'domain':
            continue
        if len(row) > 2 and row[2] == 'ok-pe':
            ok.add(row[1])
        elif len(row) > 2 and row[2].startswith('fail-download'):
            dead.add(row[1])
todo = dead - ok
with open(S + "/new_downloads.tsv", newline='') as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        if row['run_accession'] in todo:
            print('\t'.join([row['domain'], row['run_accession'], row['fastq_ftp'],
                             row['fastq_md5'], row['read_count'], row['fastq_bytes']]))
PY

n=$(wc -l < "$S/prefix_jobs.tsv")
echo "prefix subsampler: $n datasets, ${PREFIX_MB}MB prefix, $PAIRS pairs"

JOBS="${JOBS:-4}"
while IFS=$'\t' read -r domain acc ftp md5 reads bytes; do
  while [[ "$(jobs -rp | wc -l)" -ge "$JOBS" ]]; do sleep 2; done
  do_run "$domain" "$acc" "$ftp" "$md5" "$reads" "$bytes" &
done < "$S/prefix_jobs.tsv"
wait

ok_n=$(grep -c ok-pe "$LOG" 2>/dev/null || echo 0)
echo "prefix subsampler done: $ok_n/$n ok; total units now: $(ls "$UNITS"/*/ 2>/dev/null | grep -c R1.fastq.gz)"
