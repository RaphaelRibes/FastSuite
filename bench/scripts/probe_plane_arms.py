#!/usr/bin/env python3
"""Arms test: mutate ONE block of the fitted profile, re-simulate at the same
seed, and measure the per-position G curve. This isolates which fitted block
drives each position's G fraction -- no code is touched, so the experiment is
reproducible against any binary.

Why each arm exists
-------------------
The profile's `base` plane is the MEASURED per-position composition, so at the
positions an injected G-run reaches it already contains those Gs. fastsim
deconvolves the 5' artefact out of `base.R1` (the `head` block of its profile
loader: `p_mass = P(h > p)`), and does nothing else. Three contaminations are
therefore left in, one per (plane, artefact) pair:

  base.R1  is consumed as the fragment's 5' end  -> contaminated by R1's HEAD
           run (handled by fastsim) and by R1's TAIL run (NOT handled)
  base.R2  is consumed revcomp'd as the fragment's 3' end -> contaminated by
           R2's HEAD run (NOT handled) and by R2's TAIL run (NOT handled)

Arms:
  asis            the profile as fitted
  head_none       head.* = {0}          (no 5' run injected anywhere)
  tail_none       tail.* = {0}          (no 3' run injected anywhere)
  r1tail_bio      base.R1 deconvolved for R1's TAIL run only
  r2head_bio      base.R2 deconvolved for R2's HEAD run only
  r2tail_bio      base.R2 deconvolved for R2's TAIL run only
  fixed_planes    all three of the above at once (the minimal consumer-side fix)
  r1tail_exact    r1tail_bio plus the boundary-rule term P(t == L-1-p),
                  the base immediately before the run being forced non-G

A run from the 3' end covers 0-based position p iff t >= L - p; a run from the
5' end covers p iff h > p. Both arms are at the same seed and the same binary,
so the RNG stream is shared and only the mutated block differs.

Usage:
  probe_plane_arms.py <domain>/<acc> --fastsim P --fastcheck P --outdir D
                      [--pairs 50000] [--seed 42] [--arms a,b,c]
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

MATES = ("R1", "R2")


def load(p):
    with open(p) as fh:
        return json.load(fh)


def gser(rep):
    bp = rep["base_position"]
    return {int(k): (bp[k][2] / sum(bp[k][:4]) if sum(bp[k][:4]) else None) for k in bp}


def cover(dist, k):
    """P(run length >= k) for a weight map {str(len): weight}."""
    tot = sum(dist.values())
    return sum(v for key, v in dist.items() if int(key) >= k) / tot if tot else 0.0


def share(dist, k):
    """P(run length == k)."""
    tot = sum(dist.values())
    return dist.get(str(k), 0) / tot if tot else 0.0


def deconvolve(plane, dist, side, read_len, exact=False):
    """Carry the plane's G channel through the mirror deconvolution of one
    injected-Run distribution. Returns the number of positions changed."""
    n = 0
    for p in range(len(plane)):
        if side == "tail":
            pcov = cover(dist, read_len - p)            # t >= L - p covers p
            pforced = share(dist, read_len - 1 - p) if exact else 0.0
        else:
            pcov = cover(dist, p + 1)                   # h > p covers p
            pforced = share(dist, p) if exact else 0.0
        denom = 1.0 - pcov - pforced
        if pcov <= 1e-6 or denom <= 1e-6:
            continue
        row = list(plane[p])
        bio = min(max((row[2] - pcov) / denom, 0.0), 1.0)
        rest = sum(row[:2]) + row[3]
        if rest > 1e-9:
            sc = (1.0 - bio) / rest
            row[0] *= sc
            row[1] *= sc
            row[3] *= sc
        row[2] = bio
        s = sum(row)
        plane[p] = [x / s for x in row]
        n += 1
    return n


def mutate(prof, arm, read_len):
    p = json.loads(json.dumps(prof))
    if arm == "asis":
        return p, 0
    if arm == "head_none":
        for m in MATES:
            p["head"][m] = {"0": 10000}
        return p, 0
    if arm == "tail_none":
        for m in MATES:
            p["tail"][m] = {"0": 10000}
        return p, 0
    n = 0
    if arm in ("r1tail_bio", "r1tail_exact", "fixed_planes"):
        n += deconvolve(p["base"]["R1"], p["tail"].get("R1") or {}, "tail", read_len,
                        exact=(arm == "r1tail_exact"))
    if arm in ("r2head_bio", "fixed_planes"):
        n += deconvolve(p["base"]["R2"], p["head"].get("R2") or {}, "head", read_len)
    if arm in ("r2tail_bio", "fixed_planes"):
        n += deconvolve(p["base"]["R2"], p["tail"].get("R2") or {}, "tail", read_len)
    if n == 0 and arm not in ("asis", "head_none", "tail_none"):
        return None, 0
    return p, n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("--work", default="/mnt/Fast/fastsuite-scratch/work")
    ap.add_argument("--units", default="/mnt/Fast/fastsuite-scratch/units")
    ap.add_argument("--fastsim", required=True)
    ap.add_argument("--fastcheck", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--pairs", type=int, default=50000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--arms", default="asis,r1tail_bio,r2head_bio,r2tail_bio,fixed_planes")
    args = ap.parse_args()

    domain, acc = args.dataset.split("/")
    adir = os.path.join(args.work, domain, acc, "v0")
    prof = load(os.path.join(adir, "real.profile.json"))
    real = load(os.path.join(adir, "real.report.json"))
    rg = gser(real)
    read_len = int(prof["read_len"])
    ins = prof.get("insert") or {}
    ins_min = max(1, min(min(int(k) for k in ins), read_len)) if ins else 1
    os.makedirs(args.outdir, exist_ok=True)
    cache = os.path.join(args.outdir, "cache")
    os.makedirs(cache, exist_ok=True)

    def measure(fc, r1p, r2p, tag):
        outs, deps = [], []
        for path, m in ((r1p, "R1"), (r2p, "R2")):
            rep = os.path.join(cache, f"{tag}.{m}.json")
            subprocess.run([fc, "--json", rep, "-q", "-j", "4", "--", path],
                           capture_output=True)
            d = load(rep)
            outs.append(gser(d))
            bp = d["base_position"]
            deps.append({int(k): sum(bp[k][:4]) for k in bp})
        # the pooled curve exactly as fastcheck pools: depth-weighted over mates
        pool = {}
        for i in range(1, max(max(a) for a in outs) + 1):
            num = den = 0.0
            for a, dd in zip(outs, deps):
                if i in a and a[i] is not None and dd.get(i):
                    num += dd[i] * a[i]
                    den += dd[i]
            if den:
                pool[i] = num / den
        return pool, outs[0], outs[1]

    real_mate = measure(args.fastcheck,
                        os.path.join(args.units, domain, f"{acc}.R1.fastq.gz"),
                        os.path.join(args.units, domain, f"{acc}.R2.fastq.gz"), "real")
    # the stored v0 sim, if present, as a reference
    stored = os.path.join(adir, "sim.v0.report.json")
    if os.path.exists(stored):
        rs = load(stored)
        sg = {(k): None for k in []}  # noqa
    print(f"=== {args.dataset} read_len={read_len} ins_min={ins_min} pairs={args.pairs}")
    print("  real   p1=%.4f (R1 %.4f R2 %.4f) last=%.4f  head=%s  tail=%s" % (
        real_mate[0].get(1, -1), real_mate[1].get(1, -1), real_mate[2].get(1, -1),
        real_mate[0].get(read_len, -1), json.dumps(prof["head"]), json.dumps(prof["tail"])))

    for arm in args.arms.split(","):
        arm = arm.strip()
        p, nch = mutate(prof, arm, read_len)
        if p is None:
            print(f"  {arm:<14} skipped: that run-length block is empty for this library")
            continue
        pp = os.path.join(args.outdir, f"{acc}.{arm}.profile.json")
        with open(pp, "w") as fh:
            json.dump(p, fh, indent=1)
        out = os.path.join(args.outdir, f"{acc}.{arm}.sim")
        r = subprocess.run([args.fastsim, "--pairs", str(args.pairs), "--profile", pp,
                            "--out", out, "--seed", str(args.seed),
                            "--read-len", str(read_len), "--insert-min", str(ins_min),
                            "--drop-rate", "0", "--force"],
                           capture_output=True, cwd=args.outdir)
        if r.returncode != 0 or not os.path.exists(out + ".gt.tsv"):
            print(f"  {arm:<14} fastsim FAILED {r.stderr[-200:]!r}")
            continue
        pool, r1g, r2g = measure(args.fastcheck, out + "_R1.fastq.gz",
                                 out + "_R2.fastq.gz", arm)
        mx, at = 0.0, -1
        for i in sorted(set(rg) & set(pool)):
            d = abs(rg[i] - pool[i])
            if d > mx:
                mx, at = d, i
        # length-matched gate: pool the sim's mates with the REAL mate depths
        deep_r1 = load(os.path.join(cache, "real.R1.json"))["base_position"]
        deep_r2 = load(os.path.join(cache, "real.R2.json"))["base_position"]
        dr = [{int(k): sum(deep_r1[k]) for k in deep_r1},
              {int(k): sum(deep_r2[k]) for k in deep_r2}]
        lm, lat = 0.0, -1
        for i in sorted(set(rg) & set(r1g) & set(r2g)):
            num = den = 0.0
            for g, d in zip((r1g, r2g), dr):
                if i in g and g[i] is not None and d.get(i):
                    num += d[i] * g[i]
                    den += d[i]
            if not den:
                continue
            d0 = abs(rg[i] - num / den)
            if d0 > lm:
                lm, lat = d0, i
        print(f"  {arm:<14} n={nch:<3} pool@p1={pool.get(1, float('nan')):.4f} "
              f"(R1 {r1g.get(1, float('nan')):.4f}, R2 {r2g.get(1, float('nan')):.4f}) "
              f"poollast={pool.get(read_len, float('nan')):.4f} "
              f"(R1 {r1g.get(read_len, float('nan')):.4f}) "
              f"| GAte pooled={mx:.4f}@{at}  lengthmatched={lm:.4f}@{lat}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
