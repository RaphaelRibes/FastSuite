#!/usr/bin/env python3
"""Plot head and tail polyG CDFs: simulated (from GT) vs profile.

Usage: plot_heads_tails.py PROFILE_JSON SIM_GT OUT_PREFIX

Reads the ground-truth TSV (r1_h, r1_t, r2_h, r2_t per pair) and the
measured profile, then plots 2x2 CDF panels: head R1, head R2,
tail R1, tail R2. Simulated = solid, profile = dashed.

Requires: numpy, matplotlib.
"""

import argparse
import json
from collections import defaultdict

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def hist_to_cdf(h):
    """{len: count} (int or str keys) -> (lens, cumfrac)."""
    norm = {int(k): int(v) for k, v in h.items() if int(v) > 0}
    total = sum(norm.values())
    if total == 0:
        return np.array([]), np.array([])
    keys = sorted(norm)
    cum = np.cumsum([norm[k] for k in keys]) / total
    return np.array(keys), cum


def cdf_at(keys, cum, x):
    """CDF value at x (fraction of reads with run <= x)."""
    if len(keys) == 0:
        return 1.0
    idx = np.searchsorted(keys, x, side="right") - 1
    if idx < 0:
        return 0.0
    return cum[idx]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("profile_json")
    ap.add_argument("sim_gt")
    ap.add_argument("out_prefix")
    args = ap.parse_args()

    with open(args.profile_json) as f:
        prof = json.load(f)

    # ---- tally heads and tails from ground truth ----
    heads = {"R1": defaultdict(int), "R2": defaultdict(int)}
    tails = {"R1": defaultdict(int), "R2": defaultdict(int)}
    n = 0
    with open(args.sim_gt) as f:
        next(f)  # header
        for line in f:
            p = line.split("\t")
            if len(p) < 7:
                continue
            n += 1
            h1, t1, h2, t2 = (int(x) for x in p[2:6])
            heads["R1"][h1] += 1
            heads["R2"][h2] += 1
            tails["R1"][t1] += 1
            tails["R2"][t2] += 1
    print("pairs read: " + str(n))

    # ---- text metrics (conditional on run >= 1, mirroring validate_sim) ----
    for mate in ("R1", "R2"):
        hsim = {k: v for k, v in heads[mate].items() if k >= 1}
        tsim = {k: v for k, v in tails[mate].items() if k >= 1}
        hprof = {k: v for k, v in prof["head"][mate].items() if int(k) >= 1 and int(v) > 0}
        tprof = {k: v for k, v in prof["tail"][mate].items() if int(k) >= 1 and int(v) > 0}
        stot = float(sum(heads[mate].values()))
        hrate_sim = sum(hsim.values()) / stot
        hfull = prof["head"][mate]
        htot_all = sum(int(v) for v in hfull.values())
        hrate_prof = sum(int(v) for k, v in hfull.items() if int(k) >= 1) / htot_all if htot_all else 0.0
        k1, c1 = hist_to_cdf(hsim)
        k2, c2 = hist_to_cdf(hprof)
        d = np.abs(np.interp(np.arange(max(k1.max(), k2.max()) + 1), k1, c1) -
                   np.interp(np.arange(max(k1.max(), k2.max()) + 1), k2, c2)).max()
        print(mate + " head rate sim {:.4f} prof {:.4f}  head CDF KS D {:.4f}".format(hrate_sim, hrate_prof, d))

    # ---- figure: 2x2 ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    panels = []
    for ax, mate, kind in zip(axes.flat,
                              ("R1", "R2", "R1", "R2"),
                              ("head", "head", "tail", "tail")):
        sim_hist = {k: v for k, v in (heads if kind == "head" else tails)[mate].items() if k >= 1}
        prof_hist = {k: v for k, v in prof[kind][mate].items() if int(k) >= 1 and int(v) > 0}
        k1, c1 = hist_to_cdf(sim_hist)
        k2, c2 = hist_to_cdf(prof_hist)
        ax.plot(k1, c1, label="simulated")
        ax.plot(k2, c2, "--", label="profile")
        ax.set_xlabel(kind + " G-run length (bp)")
        ax.set_ylabel("CDF (reads with run)")
        ax.set_title(mate + " " + kind + " CDF", loc="left")
        ax.legend()
        if mate == "R1" and kind == "head":
            ax.set_xlim(0, 12)
        else:
            ax.set_xlim(0, 160)
        # annotate head/tail rate
        rate = sum(sim_hist.values()) / float(sum((heads if kind == "head" else tails)[mate].values()))
        ax.text(0.98, 0.05, "rate {:.3f}".format(rate), transform=ax.transAxes,
                ha="right", fontsize=8)
    fig.tight_layout()
    out = args.out_prefix + "_heads_tails.png"
    fig.savefig(out, dpi=130)
    print("figure: " + out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
