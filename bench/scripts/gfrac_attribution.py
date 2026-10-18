#!/usr/bin/env python3
"""g_frac attribution: WHERE the real-vs-simulated per-position G fraction
diverges, and WHICH mechanism carries the excess.

This replaces the single-dataset `verify_gfrac_deconv.py` (which only printed
one number for one dataset). Per (dataset, variant) it reports four things the
closed loop's `g_frac_max_diff` gate collapses into one scalar:

1. WHERE the mismatch lives: the real and simulated per-position G-fraction
   curves, their difference, the argmax position, and the mean |diff| in the
   head (1-10), body and tail (last 10) zones.

2. SAMPLE MIX. The report's `base_position` is POOLED over both mates, so a
   position that only one mate reaches in the real library is a two-mate
   average in the simulation -- and the simulation ALWAYS emits exactly
   `read_len` bases for both mates ("Sequencer output is always exactly
   read_len bases", fastsim main.rs). `gfrac_max_diff_lengthmatched` re-weights
   the simulation's per-mate curves by the REAL per-mate depths: the gate the
   closed loop would score if the simulated library had the real library's mate
   lengths. `scored - lengthmatched` is the pooling artefact.

3. THE 3' BOUNDARY PLANE. `tail_plane_pred_r1` is the model prediction of the
   simulation's R1 terminal G under "the measured base plane is used as a
   biology draw at the terminal position without removing the tail artefact":
       pred = P(t >= 1) + P(t = 0) * plane_G[last]
   (t = the fitted 3' G-run artefact length). With `--per-mate` the MEASURED
   per-mate terminal G of real and simulation sit next to it, so the model is
   checked against data rather than asserted.

4. N overlay and GC, as controls: `n_overlay_max_absd` is how far a position's
   G fraction moves when the denominator counts N as a base.

Per-mate curves need a fastcheck pass on each mate's FASTQ (`--per-mate`); the
stored reports only carry the pooled plane.

Usage:
  gfrac_attribution.py --work ROOT --datasets f.txt --out O.tsv [--per-mate]
                       [--fastcheck PATH] [--cache DIR] [--force]

The output is REFUSED if it already exists (docs/CLOSED-LOOP-RERUN.md §6):
pass --force or a new --out path.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys

ROW_FIELDS = [
    "domain", "run_accession", "variant",
    "real_len_r1", "real_len_r2", "sim_len_r1", "sim_len_r2", "real_asymmetric",
    "gfrac_max_diff_scored", "argmax_pos", "argmax_real", "argmax_sim", "argmax_zone",
    "head_mean_absd", "body_mean_absd", "tail_mean_absd",
    "gfrac_max_diff_lengthmatched", "lm_argmax_pos", "lm_argmax_real", "lm_argmax_sim",
    "realR1_G_last", "simR1_G_last", "realR2_G_last", "simR2_G_last",
    "tail_plane_pred_r1", "tail_plane_p_ge1_r1", "tail_plane_plane_G_r1",
    "real_tail_ge1", "sim_tail_ge1", "real_head_ge1", "sim_head_ge1",
    "n_overlay_max_absd", "gc_real", "gc_sim",
    "verdict_scored", "verdict_lengthmatched", "failing_scored",
]


def load(path):
    with open(path) as fh:
        return json.load(fh)


def g_series(report, n_aware=False):
    bp = report.get("base_position") or {}
    out = []
    for pos in sorted(bp, key=lambda k: int(k)):
        c = bp[pos]
        if len(c) < 4:
            out.append(None)
            continue
        den = sum(c[:5]) if (n_aware and len(c) > 4) else sum(c[:4])
        out.append(c[2] / den if den else None)
    return out


def depth_series(report):
    bp = report.get("base_position") or {}
    return [sum(bp[k]) for k in sorted(bp, key=lambda x: int(x))]


def polyx_runs(report):
    px = report.get("polyx") or {}
    out = {}
    for side, key in (("tail", "3p"), ("head", "5p")):
        d = (px.get(key) or {}).get("G") or {}
        out[side] = {int(k): int(v) for k, v in d.items()}
    return out


def run_ge1(dist):
    tot = sum(dist.values())
    return sum(v for k, v in dist.items() if k >= 1) / tot if tot else None


def mate_mean_len(report):
    pf = report.get("per_file") or []
    lens = [f.get("mean_len") for f in pf if isinstance(f, dict)]
    if len(lens) >= 2:
        return lens[0], lens[1]
    if lens:
        return lens[0], lens[0]
    return None, None


def fastcheck_curve(fc, path, cache, tag):
    """Per-position G curve + depth of one FASTQ, via fastcheck --json (cached)."""
    os.makedirs(cache, exist_ok=True)
    rep = os.path.join(cache, f"{tag}.report.json")
    if not os.path.exists(rep):
        r = subprocess.run([fc, "--json", rep, "-q", "-j", "4", "--", path],
                           capture_output=True)
        if r.returncode != 0:
            print(f"  fastcheck failed on {path}: {r.stderr[-300:]!r}", file=sys.stderr)
            return None, None
    d = load(rep)
    return g_series(d), depth_series(d)


def diff_stats(real, sim, n):
    best, bi = 0.0, -1
    for i in range(n):
        if real[i] is None or sim[i] is None:
            continue
        dd = abs(real[i] - sim[i])
        if dd > best:
            best, bi = dd, i
    zones = {"head(1-10)": range(0, min(10, n)),
             "body": range(min(10, n), max(min(10, n), n - 10)),
             "tail(last10)": range(max(0, n - 10), n)}
    means = {}
    for name, rg in zones.items():
        idx = [i for i in rg if i < n and real[i] is not None and sim[i] is not None]
        means[name] = sum(abs(real[i] - sim[i]) for i in idx) / len(idx) if idx else 0.0
    return best, bi, means


def length_matched(real_pool, sim_pool, sim_mate_g, real_mate_d):
    """Pool the sim's per-mate curves using the REAL per-mate depths."""
    n = min(len(real_pool), len(sim_pool))
    out = []
    for i in range(n):
        num = den = 0.0
        for mg, md in zip(sim_mate_g, real_mate_d):
            if mg and i < len(mg) and i < len(md) and mg[i] is not None and md[i] > 0:
                num += md[i] * mg[i]
                den += md[i]
        out.append(num / den if den else (sim_pool[i] if i < len(sim_pool) else None))
    return out


def last_g(g):
    if not g:
        return None
    for i in range(len(g) - 1, -1, -1):
        if g[i] is not None:
            return g[i]
    return None


def zone_of(pos1, n):
    if pos1 <= 10:
        return "head"
    if pos1 > n - 10:
        return "tail"
    return "body"


def tail_plane_model(prof, mate, read_len):
    """Terminal-G prediction under 'no 3' deconvolution of the base plane'.

    fastsim draws the R1 read's last base at fragment index `read_len - 1` from
    the head plane (the plane covers fragment positions 0..len-1). When the tail
    artefact length t is 0 no run is injected, so that base is the biology draw
    -- taken from the MEASURED plane value, which already contains the injected
    tail Gs. With per-read probability P(t>=1) the last base is forced G, so
        E[terminal G] = P(t>=1) + P(t=0) * plane_G[last].
    Returns (pred, p_ge1, plane_G) or (None, None, None).
    """
    t = (prof.get("tail") or {}).get(mate) or {}
    tot = sum(t.values())
    plane = (prof.get("base") or {}).get(mate) or []
    if not tot or not plane:
        return None, None, None
    p_ge1 = sum(v for k, v in t.items() if int(k) >= 1) / tot
    g_last = plane[min(read_len, len(plane)) - 1][2]
    return p_ge1 + (1.0 - p_ge1) * g_last, p_ge1, g_last


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--units", default=None)
    ap.add_argument("--datasets", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--variants", default="0,1,2")
    ap.add_argument("--per-mate", action="store_true")
    ap.add_argument("--fastcheck", default=None)
    ap.add_argument("--cache", default="/tmp/gfrac_attr_cache")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if os.path.exists(args.out) and not args.force:
        print(f"gfrac_attribution: {args.out} exists — refusing to overwrite. "
              f"Pass --force or a new path.", file=sys.stderr)
        return 2

    units_root = args.units or os.path.join(os.path.dirname(args.work.rstrip("/")), "units")
    variants = [v.strip() for v in args.variants.split(",") if v.strip() != ""]
    datasets = [l.strip() for l in open(args.datasets)
                if l.strip() and not l.startswith("#")]

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(ROW_FIELDS)
        for ds in datasets:
            domain, acc = ds.split("/")
            adir = os.path.join(args.work, domain, acc)
            real_r = os.path.join(adir, "v0", "real.report.json")
            if not os.path.exists(real_r):
                real_r = os.path.join(adir, "real.report.json")
            prof_p = os.path.join(adir, "v0", "real.profile.json")
            if not os.path.exists(real_r) or not os.path.exists(prof_p):
                print(f"skip {ds}: no real report/profile", file=sys.stderr)
                continue
            real = load(real_r)
            prof = load(prof_p)
            read_len = int(prof["read_len"])
            r1_real, r2_real = mate_mean_len(real)
            rg_pool = g_series(real)
            rg_pool_n = g_series(real, n_aware=True)
            rpx = polyx_runs(real)

            r_mate_g = r_mate_d = None
            if args.per_mate and args.fastcheck:
                r_mate_g, r_mate_d = [], []
                for m in ("R1", "R2"):
                    g, d = fastcheck_curve(
                        args.fastcheck,
                        os.path.join(units_root, domain, f"{acc}.{m}.fastq.gz"),
                        args.cache, f"{domain}_{acc}_real_{m}")
                    r_mate_g.append(g)
                    r_mate_d.append(d)
                if not all(g for g in r_mate_g):
                    r_mate_g = r_mate_d = None

            for variant in variants:
                vdir = os.path.join(adir, f"v{variant}")
                sim_r = os.path.join(vdir, f"sim.v{variant}.report.json")
                if not os.path.exists(sim_r):
                    continue
                sim = load(sim_r)
                s1, s2 = mate_mean_len(sim)
                sg_pool = g_series(sim)
                sg_pool_n = g_series(sim, n_aware=True)
                spx = polyx_runs(sim)
                n = min(len(rg_pool), len(sg_pool))
                scored, apos, zones = diff_stats(rg_pool, sg_pool, n)

                n_overlay = 0.0
                for i in range(min(n, len(rg_pool_n), len(sg_pool_n))):
                    if None in (rg_pool[i], sg_pool[i], rg_pool_n[i], sg_pool_n[i]):
                        continue
                    n_overlay = max(n_overlay, abs(rg_pool[i] - rg_pool_n[i]),
                                    abs(sg_pool[i] - sg_pool_n[i]))

                lm, lm_pos, lm_real, lm_sim = None, None, None, None
                s_mate_g = None
                if r_mate_g:
                    s_mate_g = []
                    for m in ("R1", "R2"):
                        g, _ = fastcheck_curve(
                            args.fastcheck,
                            os.path.join(vdir, f"sim.v{variant}_{m}.fastq.gz"),
                            args.cache, f"{domain}_{acc}_v{variant}_{m}")
                        s_mate_g.append(g)
                    if all(g for g in s_mate_g):
                        lm_curve = length_matched(rg_pool, sg_pool, s_mate_g, r_mate_d)
                        lm, lm_pos, _ = diff_stats(rg_pool, lm_curve, n)
                        if lm_pos >= 0:
                            lm_real, lm_sim = rg_pool[lm_pos], lm_curve[lm_pos]
                    else:
                        s_mate_g = None

                rR1 = last_g(r_mate_g[0]) if r_mate_g else None
                rR2 = last_g(r_mate_g[1]) if r_mate_g else None
                sR1 = last_g(s_mate_g[0]) if s_mate_g else None
                sR2 = last_g(s_mate_g[1]) if s_mate_g else None
                tp_r1, tp_ge1, tp_plane = tail_plane_model(prof, "R1", read_len)

                def f6(x):
                    return "" if x is None else "%.6f" % x

                asym = (r1_real is not None and r2_real is not None
                        and abs(r1_real - r2_real) > 0.5)
                row = [
                    domain, acc, variant,
                    f6(r1_real), f6(r2_real), f6(s1), f6(s2), int(asym),
                    "%.4f" % scored, apos + 1,
                    f6(rg_pool[apos]), f6(sg_pool[apos]), zone_of(apos + 1, n),
                    f6(zones["head(1-10)"]), f6(zones["body"]), f6(zones["tail(last10)"]),
                    f6(lm), (lm_pos + 1) if lm_pos is not None and lm_pos >= 0 else "",
                    f6(lm_real), f6(lm_sim),
                    f6(rR1), f6(sR1), f6(rR2), f6(sR2),
                    f6(tp_r1), f6(tp_ge1), f6(tp_plane),
                    f6(run_ge1(rpx["tail"])), f6(run_ge1(spx["tail"])),
                    f6(run_ge1(rpx["head"])), f6(run_ge1(spx["head"])),
                    "%.4f" % n_overlay,
                    f6(real["metrics"].get("gc_fraction")),
                    f6(sim["metrics"].get("gc_fraction")),
                    "PASS" if scored <= 0.02 else "FAIL",
                    "" if lm is None else ("PASS" if lm <= 0.02 else "FAIL"),
                    "" if scored <= 0.02 else "g_frac",
                ]
                w.writerow(row)
                fh.flush()
                print("%-16s %-12s v%s scored=%.4f @%d(%s) lm=%s  asym=%d" % (
                    domain, acc, variant, scored, apos + 1, zone_of(apos + 1, n),
                    ("%.4f" % lm) if lm is not None else "-", asym), flush=True)
    print(f"gfrac_attribution: wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
