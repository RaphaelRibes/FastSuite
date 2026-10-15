#!/usr/bin/env python3
import json
real = json.load(open("/mnt/Fast/fastsuite-scratch/work/bacterial_isolate/DRR014283/v0/real.report.json"))
sim = json.load(open("/tmp/gfracverify2/sim.json"))
def g0(r):
    bp = r.get("base_position") or {}
    row = bp.get("1")
    if not row: return None
    t = sum(row)
    return row[2] / t if t else None
rg, sg = g0(real), g0(sim)
print(f"DRR014283  G@pos0 real={rg:.4f} sim={sg:.4f} diff={abs(rg-sg):.4f}")
print(f"  (was 0.1273 before the deconvolution, 0.0636 in the stale-output run)")
