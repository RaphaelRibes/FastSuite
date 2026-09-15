#!/usr/bin/env python3
"""Dump profile values that explain the R2 residual."""
import json
import sys

prof = json.load(open(sys.argv[1]))
read_len = prof["read_len"]

for mate in ("R1", "R2"):
    print("=== " + mate + " ===")
    print("run_q_mean: " + str(prof["run_q_mean"][mate]))
    print("tail_q_collapse: " + str(prof["tail_q_collapse"][mate]))
    uncond = prof["pos_q_uncond"][mate]
    out = prof["pos_q"][mate]
    t_hist = prof["tail"][mate]
    t_total = sum(t_hist.values())
    print("pos uncond out p_i")
    for i in list(range(0, 8)) + list(range(140, read_len)):
        p_i = sum(v for k, v in t_hist.items() if int(k) >= read_len - i) / t_total
        print("  " + str(i) + "  " + format(uncond[i], ".2f") + "  " + format(out[i], ".2f") + "  " + format(p_i, ".4f"))
    rmq = prof["read_meanq"][mate]
    total = sum(rmq.values())
    mean = sum(int(k) * v for k, v in rmq.items()) / total
    mean_uncond = sum(uncond) / len(uncond)
    print("mean(read_meanq)=" + format(mean, ".2f") + " mean(uncond)=" + format(mean_uncond, ".2f"))
    keys = sorted(int(k) for k in rmq)
    print("read_meanq bins: min=" + str(keys[0]) + " max=" + str(keys[-1]) + " n_bins=" + str(len(rmq)))
    for k in keys[:8]:
        print("  bin " + str(k) + " count " + str(rmq[str(k)]))
