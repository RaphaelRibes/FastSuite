#!/bin/bash
#SBATCH --job-name=fs-acquire
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=48:00:00
#SBATCH --output=slurm/logs/acquire-%j.out
#SBATCH --error=slurm/logs/acquire-%j.err
#
# FastSuite PLAN-ML — phase 1: pull the real corpus from ENA and reduce each run
# to the 100k-pair sweep unit.
#
# Why a job and not a login-node loop: 700 runs × ~2 GB is ~1.4 TB, hours of
# transfer, and the rule is that nothing runs on the login node. A batch job is
# also resumable: curl -C - continues a partial file, and the per-run manifest
# makes an already-subsampled dataset a no-op, so relaunching after the 48 h
# wall clock costs nothing more than the queue wait.
#
# Egress must have been proven by the pilot (step 2). If it was NOT reachable,
# run this from the login node in a bounded loop instead — with your explicit
# go-ahead, since the standing rule puts nothing on the login node.
#
# 1 GB because curl buffers and python holds the manifest; disk is the real
# resource (each run is deleted after its subsample is verified).
set -euo pipefail

WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$WORKDIR"
mkdir -p slurm/logs data/raw data/work manifests

MANIFEST="manifests/real.tsv"          # accession, url_r1, url_r2, md5_r1, md5_r2, stratum…
UNITS=100000                           # pairs per sweep unit
STRIDE=$(( ${UNIT_STRIDE:-1} ))        # keep every STRIDE-th record (awk subsampler)

[[ -f "$MANIFEST" ]] || { echo "no $MANIFEST — generate it with scripts/ena_query.py first" >&2; exit 1; }

# Deterministic, dependency-free subsampling: seqkit is not installable on this
# cluster (pixi deadlocks), so every STRIDE-th record is kept, both mates with
# it. Same STRIDE ⇒ same reads, and R1/R2 stay in sync by construction.
subsample() {   # subsample <in.fastq.gz> <out.fastq.gz> <stride>
    zcat "$1" | awk -v s="$3" 'NR%4==1{n++} n%s==1' | gzip -1 > "$2"
}

n=0
while IFS=$'\t' read -r acc r1 r2 md5_1 md5_2 stratum; do
    [[ "$acc" =~ ^#.*$ || -z "$acc" ]] && continue
    out="data/work/$acc"
    if [[ -f "$out/done" ]]; then echo "skip $acc (done)"; continue; fi
    mkdir -p "$out"
    for mate in 1 2; do
        url=$([[ $mate == 1 ]] && echo "$r1" || echo "$r2")
        f="$out/raw_R$mate.fastq.gz"
        [[ -s "$f" ]] || curl -fL --retry 5 --retry-delay 10 -C - -o "$f" "$url"
        want=$([[ $mate == 1 ]] && echo "$md5_1" || echo "$md5_2")
        if [[ -n "$want" && "$want" != "-" ]]; then
            got=$(md5sum "$f" | awk '{print $1}')
            [[ "$got" == "$want" ]] || { echo "MD5 MISMATCH $acc R$mate: $got != $want" >&2; rm -f "$f"; continue 2; }
        fi
        subsample "$f" "$out/in_R$mate.fastq.gz" "$STRIDE"
    done
    # verify the unit, then delete the raw: /scratch peaks at ~60 GB, not 1.4 TB
    r1n=$(( $(zcat "$out/in_R1.fastq.gz" | wc -l) / 4 ))
    r2n=$(( $(zcat "$out/in_R2.fastq.gz" | wc -l) / 4 ))
    [[ "$r1n" == "$r2n" ]] || { echo "MATE DESYNC $acc: $r1n vs $r2n" >&2; continue; }
    sha1=$(sha256sum "$out/in_R1.fastq.gz" | awk '{print $1}')
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$acc" "$r1n" "$sha1" "$stratum" "$(hostname)" "$(date -Is)" \
        >> manifests/units.tsv
    rm -f "$out/raw_R1.fastq.gz" "$out/raw_R2.fastq.gz"
    touch "$out/done"
    echo "ok $acc: $r1n pairs, sha256 ${sha1:0:12}"
    n=$((n+1))
done < "$MANIFEST"

echo "acquired/updated $n unit(s); manifest: manifests/units.tsv"
