#!/usr/bin/env python3
"""Compare mate-level read lengths and per-position depth, real vs sim, for one dataset."""
import json
import sys

root = "/mnt/Fast/fastsuite-scratch/work"
ds = sys.argv[1]
variant = sys.argv[2] if len(sys.argv) > 2 else "0"
d = f"{root}/{ds}/v{variant}"
real = json.load(open(f"{d}/real.report.json"))
sim = json.load(open(f"{d}/sim.v{variant}.report.json"))


def depths(rep):
    bp = rep["base_position"]
    return {int(k): sum(bp[k]) for k in bp}


print("=== ", ds, "v", variant)
print("real len hist:", real.get("read_length_hist"))
print("sim  len hist:", sim.get("read_length_hist"))
rd, sd = depths(real), depths(sim)
ks = sorted(rd)
print("real depth last8:", [(k, rd[k]) for k in ks[-8:]])
ks = sorted(sd)
print("sim  depth last8:", [(k, sd[k]) for k in ks[-8:]])
for rep, tag in ((real, "real"), (sim, "sim")):
    pf = rep.get("per_file") or []
    print(f"--- {tag} per_file ({len(pf)})")
    for f in pf:
        if isinstance(f, dict):
            print("   ", {k: f[k] for k in list(f)[:14]})
print("--- polyx blocks")
print("real polyx:", json.dumps(real.get("polyx"))[:1500])
print("sim  polyx:", json.dumps(sim.get("polyx"))[:1500])
