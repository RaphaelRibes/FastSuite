#!/usr/bin/env python3
"""Frame-shift test: is the simulation's per-position curve the real one
convolved with the head-run length distribution?

fastsim builds a read as [G; h] + insert[...], so a read whose head run is h
shows the fragment's base (p - h) at its own position p. The profile's `base`
plane is a MEASURED read-position composition, so if the simulator consumed it
as a fragment-position composition the realised read-position curve would be

    sim[p] = sum_h P(h) * real[p - h]

i.e. the measured curve smeared by the head-run distribution. This script tests
that prediction against the measured simulation curve, position by position, on
the pooled reports (no fastsim run needed), and reports the residual after the
prediction is subtracted.

Usage: probe_frame_shift.py <domain>/<acc> [...]
"""

import json
import sys

ROOT = "/mnt/Fast/fastsuite-scratch/work"


def load(p):
    return json.load(open(p))


def gser(rep):
    bp = rep["base_position"]
    return {int(k): (bp[k][2] / sum(bp[k][:4]) if sum(bp[k][:4]) else None) for k in bp}


def head_pmf(prof, mate):
    d = (prof["head"] or {}).get(mate) or {}
    tot = sum(d.values())
    return {int(k): v / tot for k, v in d.items()} if tot else {0: 1.0}


for ds in sys.argv[1:]:
    d = f"{ROOT}/{ds}/v0"
    real = gser(load(f"{d}/real.report.json"))
    sim = gser(load(f"{d}/sim.v0.report.json"))
    prof = load(f"{d}/real.profile.json")
    # pooled mates: the head distributions differ per mate, use the average pmf
    ph = head_pmf(prof, "R1")
    pr = head_pmf(prof, "R2")
    pmf = {}
    for k in set(ph) | set(pr):
        pmf[k] = 0.5 * ph.get(k, 0.0) + 0.5 * pr.get(k, 0.0)
    print(f"=== {ds}  head pmf {dict(sorted(pmf.items()))}")
    print("  pos | real   sim    rawdiff | pred(shift) pred-sim | residual_after_shift")
    worst_raw = worst_res = (0.0, -1)
    for p in range(4, 25):
        if p not in real or p not in sim or real[p] is None or sim[p] is None:
            continue
        pred = 0.0
        wsum = 0.0
        for h, w in pmf.items():
            q = p - h
            if q in real and real[q] is not None:
                pred += w * real[q]
                wsum += w
        if wsum <= 0:
            continue
        pred /= wsum
        raw = sim[p] - real[p]
        res = sim[p] - pred
        if abs(raw) > abs(worst_raw[0]):
            worst_raw = (raw, p)
        if abs(res) > abs(worst_res[0]):
            worst_res = (res, p)
        print("  %3d | %.4f %.4f %+.4f | %.4f %+.4f | %+.4f" % (
            p, real[p], sim[p], raw, pred, pred - sim[p], res))
    print("  worst raw diff %.4f @%d ;  worst residual after the shift model %.4f @%d"
          % (worst_raw[0], worst_raw[1], worst_res[0], worst_res[1]))
