#!/usr/bin/env python3
"""Diagnose: joint_tail R1 marginal vs tail histogram in mandril_E.json."""
import json

prof = json.load(open("/scratch/users/ribesr/novaseq-sim/profiles/mandril_E.json"))
tail = prof["tail"]["R1"]
ttot = sum(tail.values())
jt = prof["joint_tail"]

# R1 marginal of joint: sum over R2
marg = {}
for k, v in jt.items():
    a, b = k.split(",")
    marg[int(a)] = marg.get(int(a), 0) + v
mtot = sum(marg.values())

print("len | tail_hist_frac | joint_marginal_frac")
for l in range(0, 70, 5):
    th = tail.get(str(l), 0) / ttot
    jm = marg.get(l, 0) / mtot
    flag = "  <-- MISMATCH" if abs(th - jm) > 0.005 else ""
    print(f"{l:3d} | {th:.5f} | {jm:.5f}{flag}")

# also check: is joint diagonal-heavy? (correlation)
diag = sum(v for k, v in jt.items() if k.split(",")[0] == k.split(",")[1]) / mtot
print(f"\njoint diagonal fraction: {diag:.4f}")

# check tail=1 marginal vs tail hist more precisely
print(f"\ntail_hist[1]={tail.get('1',0)/ttot:.5f}  joint_marg[1]={marg.get(1,0)/mtot:.5f}")
