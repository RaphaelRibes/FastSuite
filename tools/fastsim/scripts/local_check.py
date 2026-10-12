#!/usr/bin/env python3
"""Local sanity: compare GT and sequence-level tail rates vs profile."""
import json
import os

# The measured profile is not in the repo (it is built from the cluster
# stats), so resolve the checkout root from the env with the cluster
# literal as the fallback. Path indirection only — no logic change.
ROOT = os.environ.get("FASTSUITE_ROOT", "/scratch/users/ribesr/FastSuite")
prof = json.load(open(f"{ROOT}/tools/fastsim/profiles/mandril_E.json"))
total = sum(prof["tail"]["R1"].values())
rate = sum(v for k, v in prof["tail"]["R1"].items() if int(k) >= 1) / total
print(f"profile R1 tail rate (>=1): {rate:.4f}")

gt = open("/tmp/nvsq_real.gt.tsv")
next(gt)
n = nge1 = 0
for line in gt:
    p = line.split("\t")
    n += 1
    if int(p[3]) >= 1:  # r1_t (3' G-run)
        nge1 += 1
print(f"sim R1 tail rate (GT, >=1): {nge1/n:.4f}")


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
