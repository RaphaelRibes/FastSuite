#!/usr/bin/env python3
"""Per-position curve + fitted head/base facts for the datasets whose g_frac
mismatch sits at position 1 (the 5' head)."""
import json
import sys

root = "/mnt/Fast/fastsuite-scratch/work"


def load(p):
    return json.load(open(p))


def gser(rep):
    bp = rep["base_position"]
    return {int(k): (bp[k][2] / sum(bp[k][:4]) if sum(bp[k][:4]) else None) for k in bp}


for ds in sys.argv[1:]:
    d = f"{root}/{ds}/v0"
    real = load(f"{d}/real.report.json")
    sim = load(f"{d}/sim.v0.report.json")
    prof = load(f"{d}/real.profile.json")
    rg, sg = gser(real), gser(sim)
    print(f"=== {ds}  read_len={prof['read_len']} gc={prof['gc']:.4f}")
    print("  head dist:", json.dumps(prof["head"]))
    print("  tail dist:", json.dumps(prof["tail"]))
    b = prof.get("base", {}).get("R1") or []
    print("  base.R1 G first8:", ["%.3f" % r[2] for r in b[:8]])
    print("  base.R1 first3 rows:", [[round(x, 3) for x in r] for r in b[:3]])
    print("  pos | realG  simG   diff")
    for i in range(1, 13):
        if i in rg and i in sg:
            print("   %3d | %.4f %.4f %+.4f" % (i, rg[i], sg[i], sg[i] - rg[i]))
    print("  metrics real vs sim: gc %.4f/%.4f  meanlen %s  reads %s" % (
        real["metrics"]["gc_fraction"], sim["metrics"]["gc_fraction"],
        real.get("read_length_hist"), sim.get("read_length_hist")))
    ra = real.get("base_composition") or {}
    sa = sim.get("base_composition") or {}
    print("  base_composition real:", json.dumps(ra)[:300])
    print("  base_composition sim :", json.dumps(sa)[:300])
