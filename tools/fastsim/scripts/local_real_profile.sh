#!/bin/bash
# Local parse check: the REAL measured mandril profile against the fixed binary.
set -e
cd /home/raphael/hermes_scratch/novaseq-sim

# pull the real profile from the cluster
timeout 60 scp -q -o BatchMode=yes Io:/scratch/users/ribesr/novaseq-sim/profiles/mandril_E.json profiles/mandril_E.json

echo "=== parse + small run with the real profile ==="
./target/release/novaseq-sim \
  --pairs 100000 \
  --profile profiles/mandril_E.json \
  --out /tmp/nvsq_real \
  --seed 7 \
  --no-gzip

echo "=== sanity: tail rate from GT ==="
python3 - <<'EOF'
import json
prof = json.load(open("/home/raphael/hermes_scratch/novaseq-sim/profiles/mandril_E.json"))
total = sum(prof["tail"]["R1"].values())
rate = sum(v for k, v in prof["tail"]["R1"].items() if int(k) >= 1) / total
print(f"profile R1 tail rate (>=1): {rate:.4f}")
# GT: count reads with tail >= 1
gt = open("/tmp/nvsq_real.gt.tsv")
next(gt)
n = nge1 = 0
for line in gt:
    p = line.split("\t")
    n += 1
    if int(p[3]) >= 1:  # r1_t (3' G-run)
        nge1 += 1
print(f"sim R1 tail rate (GT, >=1): {nge1/n:.4f}")

# sequence-level: count the actual 3' G-run in the emitted R1 reads
def seq_tail_run(seq):
    run = 0
    for b in reversed(seq):
        if b == "G":
            run += 1
        else:
            break
    return run

n2 = nge1_2 = 0
with open("/tmp/nvsq_real_R1.fastq") as f:
    for i, line in enumerate(f):
        if i % 4 == 1:
            n2 += 1
            if seq_tail_run(line.strip()) >= 1:
                nge1_2 += 1
print(f"sim R1 tail rate (SEQ, >=1): {nge1_2/n2:.4f}")
EOF
