#!/usr/bin/env python3
"""Where is the residual per-position Q diff? Mean vs max."""
import numpy as np
import json
import gzip


def load_fastq(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        while True:
            name = f.readline()
            if not name:
                return
            seq = f.readline().strip()
            f.readline()
            qual = f.readline().strip()
            yield seq, qual


prof = json.load(open("profiles/mandril_E.json"))
read_len = prof["read_len"]

for mate, fname in (("R1", "sim_R1.fastq"), ("R2", "sim_R2.fastq")):
    qsum = np.zeros(read_len)
    n = 0
    for seq, qual in load_fastq("local_work/sim/" + fname):
        for i, c in enumerate(qual):
            qsum[i] += ord(c) - 33
        n += 1
    sim_q = qsum / n
    profq = np.array(prof["pos_q_uncond"][mate][:read_len])
    diff = sim_q - profq
    mx = np.abs(diff).max()
    pos = int(np.argmax(np.abs(diff)))
    print(f"{mate}: max {mx:.2f} at pos {pos}, mean|diff| {np.abs(diff).mean():.3f}")
    print("   last 8 positions sim/prof/diff:")
    for i in range(read_len - 8, read_len):
        print(f"     {i}: {sim_q[i]:.2f} {profq[i]:.2f} {diff[i]:+.2f}")
    print("   first 4:")
    for i in range(4):
        print(f"     {i}: {sim_q[i]:.2f} {profq[i]:.2f} {diff[i]:+.2f}")
