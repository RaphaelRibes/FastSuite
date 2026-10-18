#!/usr/bin/env python3
"""Inspect report/profile structure for the g_frac diagnosis (scratch probe)."""
import json, sys

root = "/mnt/Fast/fastsuite-scratch/work"
ds = sys.argv[1] if len(sys.argv) > 1 else "metagenomics/DRR083326"
v = sys.argv[2] if len(sys.argv) > 2 else "0"
d = f"{root}/{ds}/v{v}"
r = json.load(open(f"{d}/real.report.json"))
print("report keys:", list(r.keys()))
bp = r["base_position"]
ks = sorted(bp, key=int)
print("npos", len(bp), "row len", len(bp[ks[0]]))
print("first", ks[0], bp[ks[0]], "last", ks[-1], bp[ks[-1]])
print("metrics keys:", list(r["metrics"].keys()))
for k in ("gc_fraction", "mean_read_mean_q", "duplication_level", "polyx5_frac_gt_1",
          "reads_with_non_acgt", "total_reads", "read_len", "ambiguous_bases"):
    print("  ", k, r["metrics"].get(k))
print("top-level:", {k: (v if not isinstance(v, (dict, list)) else type(v).__name__) for k, v in r.items()})
for other in ("polyx3", "polyx5", "per_base_quality", "quality_position"):
    if other in r:
        val = r[other]
        print(other, type(val).__name__, (list(val.keys())[:6] if isinstance(val, dict) else str(val)[:200]))
p = json.load(open(f"{d}/real.profile.json"))
print("profile keys:", list(p.keys()))
print("read_len", p["read_len"], "gc", p["gc"], "geometry", p.get("geometry"))
print("base lens:", {k: len(v) for k, v in p.get("base", {}).items()})
print("tail:", p["tail"])
print("head:", p["head"])
print("n_run:", p.get("n_run"))
print("adapter:", p["adapter"], "dup_rate:", p.get("dup_rate"), "full_length_rate:", p.get("full_length_rate"))
print("run_q_mean:", p.get("run_q_mean"), "deep_head_q:", p.get("deep_head_q"))
print("pos_q len:", {k: len(v) for k, v in p.get("pos_q", {}).items()})
