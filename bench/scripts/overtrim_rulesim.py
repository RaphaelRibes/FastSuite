#!/usr/bin/env python3
"""Offline evaluation of boundary-cut rules on the run-CONTINUATION mechanism.

fastsim enforces a non-G base adjacent to every injected run (artefact.rs
"Boundary rule"), so its data never shows the run-continuation overtrim the
diagnosis describes. To evaluate the candidate rules on THAT mechanism we build
the generative model directly:

    artefact length   h ~ profile tail prior          (the true trim target)
    bio continuation  k ~ Geometric(p), P(k)=(1-p)p^k  (leading bio Gs after h)
    observed run      r = h + k                        (what the tool sees)

The tool sees only r (and that base r+1 is non-G). It must pick a cut c(r).
We score exact/over/under against h at tolerance 0 for each rule:

    whole-run   c = r                       (fastprep today)
    post-mean   c = round(E[h | r])         (candidate 1)
    post-MAP    c = argmax_h P(h | r)       (candidate 1, mode variant)
    bounded-h   c = min(r, h_max)           (candidate 4)

P(h | r) is Bayes with the profile's h prior and the geometric continuation:
  P(h|r) ∝ P(k = r-h) P(h) = (1-p) p^(r-h) P(h),  0 <= h <= r.
"""
from __future__ import annotations
import math, sys

# Profile tail priors (from bench/fixtures/overtrim/*.json), normalised.
MID_TAIL = {0: 9781, 1: 5, 3: 7, 4: 28, 5: 11, 6: 65, 7: 43, 8: 19, 9: 41}
HIGH_TAIL = {4: 300, 5: 700, 6: 900, 7: 1100, 8: 1100, 9: 900, 10: 700,
             11: 400, 12: 400, 14: 300, 16: 300, 18: 300, 20: 300, 0: 2300}

def geom(k, p):
    return (1 - p) * (p ** k)

def normalise(d):
    s = sum(d.values())
    return {k: v / s for k, v in d.items()}

def evaluate(prior, p, min_run=12, h_max_cap=12, label=""):
    prior = normalise(prior)
    hmax = max(prior) if h_max_cap is None else h_max_cap
    # Joint over (h, k): P(h)P(k). The tool sees r = h+k.
    # Build P(h, r) and per-rule cuts.
    rules = ["whole-run", "post-mean", "post-MAP", f"bounded-h({hmax})"]
    stat = {r_: {"exact": 0, "over": 0, "under": 0} for r_ in rules}
    # Precompute P(h|r) tables up to r_max.
    r_max = max(prior) + 40
    for r in range(0, r_max + 1):
        # P(h, r=h+k) = P(h) * P(k=r-h)
        joint = {}
        for h in range(0, r + 1):
            if h in prior and (r - h) >= 0:
                joint[h] = prior[h] * geom(r - h, p)
        mass = sum(joint.values())
        if mass <= 0:
            continue
        post = {h: w / mass for h, w in joint.items()}
        # Only reads where the scan fires (r >= min_run) get trimmed at all;
        # below min_run the tool leaves the run (cut 0). We score the boundary
        # decision for firing reads AND report the non-firing ones as cut 0.
        if r >= min_run:
            cuts = {
                "whole-run": r,
                "post-mean": int(round(sum(h * w for h, w in post.items()))),
                "post-MAP": max(post, key=lambda h: post[h]),
                f"bounded-h({hmax})": min(r, hmax),
            }
        else:
            cuts = {r_: 0 for r_ in rules}
        for name, c in cuts.items():
            for h, w in joint.items():
                # classify per unit mass
                if abs(c - h) <= 0:
                    stat[name]["exact"] += w
                elif c > h:
                    stat[name]["over"] += w
                else:
                    stat[name]["under"] += w
    print(f"\n===== {label}  (p_G={p}, min_run={min_run}) =====")
    print(f"{'rule':>16s} {'exact':>8s} {'over':>8s} {'under':>8s}")
    for name in rules:
        s = stat[name]
        tot = s["exact"] + s["over"] + s["under"]
        if tot <= 0:
            continue
        print(f"{name:>16s} {s['exact']/tot:8.4f} {s['over']/tot:8.4f} {s['under']/tot:8.4f}")

if __name__ == "__main__":
    p = float(sys.argv[1]) if len(sys.argv) > 1 else 0.235
    evaluate(MID_TAIL, p, 12, 12, "mid tail prior")
    evaluate(MID_TAIL, p, 12, 9, "mid tail prior (h_max=9)")
    evaluate(HIGH_TAIL, p, 12, 20, "high tail prior")
    evaluate(HIGH_TAIL, p, 12, 12, "high tail prior (h_max=12)")
