#!/usr/bin/env python3
"""Confirm offset-centering bias: means of read_meanq vs pos_q curves.

Usage: check_offset_bias.py PROFILE.json
"""
import json
import sys

prof = json.load(open(sys.argv[1]))

for mate in ("R1", "R2"):
    rmq = prof["read_meanq"][mate]
    tot = sum(rmq.values())
    mean_rmq = sum(int(k) * v for k, v in rmq.items()) / tot
    mean_uncond = sum(prof["pos_q_uncond"][mate]) / len(prof["pos_q_uncond"][mate])
    mean_out = sum(prof["pos_q"][mate]) / len(prof["pos_q"][mate])
    print(f"{mate}: mean(read_meanq)={mean_rmq:.2f}  mean(pos_q_uncond)={mean_uncond:.2f}  mean(pos_q out-of-run)={mean_out:.2f}")
    print(f"   expected offset bias (rmq - out) = {mean_rmq - mean_out:+.2f} Phred")
    rq = prof["run_q_mean"][mate]
    print(f"   run_q_mean: head={rq['head']} tail={rq['tail']}  collapse={prof['tail_q_collapse'][mate]}")
