#!/usr/bin/env python3
"""Summarise the g_frac attribution into one line per (dataset, variant), with
each candidate mechanism's measured contribution.

Reads the attribution TSV (gfrac_attribution.py) plus the corpus reports, and
computes:

  * the scored gate and the length-matched gate, and their difference = the
    SAMPLE-MIX contribution (the report's base_position is pooled over mates,
    the simulation always emits read_len for both, real libraries need not);
  * the 3'/5' G run-length agreement between the real and the simulated report
    (share of reads whose terminal/initial base is G, and the total-variation
    distance of the run-length distribution): if the fitter's run-length
    boundary deconvolution were the blocker this is where it would show;
  * the N-overlay and GC under-model columns, as controls.

Usage: gfrac_mechanism_summary.py --work ROOT --attr TSV --out TSV [--force]
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys


def load(p):
    with open(p) as fh:
        return json.load(fh)


def gser(rep):
    bp = rep["base_position"]
    return {int(k): (bp[k][2] / sum(bp[k][:4]) if sum(bp[k][:4]) else None) for k in bp}


def depths(rep):
    bp = rep["base_position"]
    return {int(k): sum(bp[k][:4]) for k in bp}


def run_dist(rep, side):
    px = (rep.get("polyx") or {}).get(side) or {}
    d = {int(k): int(v) for k, v in (px.get("G") or {}).items()}
    return d


def share_ge1(d, reads):
    return sum(d.values()) / reads if reads else None


def tv(a, b):
    """Total variation between two run-length distributions."""
    ta, tb = sum(a.values()), sum(b.values())
    if not ta or not tb:
        return None
    keys = set(a) | set(b)
    return 0.5 * sum(abs(a.get(k, 0) / ta - b.get(k, 0) / tb) for k in keys)


FIELDS = ["domain", "run_accession", "variant", "scored", "pos", "zone",
          "lengthmatched", "lm_pos", "sample_mix_component", "n_overlay",
          "body_mean_absd", "head_mean_absd", "tail_mean_absd",
          "real_3pG_share", "sim_3pG_share", "tail_run_TV",
          "real_5pG_share", "sim_5pG_share", "head_run_TV",
          "real_len_r1", "real_len_r2", "asym",
          "P_head_run_R1", "P_head_run_R2", "P_tail_run_R1", "P_tail_run_R2",
          "plane_G_last_R1"]


def run_ge1(dist):
    tot = sum(dist.values())
    return sum(v for k, v in dist.items() if int(k) >= 1) / tot if tot else 0.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--attr", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if os.path.exists(args.out) and not args.force:
        print(f"{args.out} exists — refusing to overwrite (pass --force)", file=sys.stderr)
        return 2

    rows = list(csv.DictReader(open(args.attr), delimiter="\t"))
    out = []
    for r in rows:
        d = os.path.join(args.work, r["domain"], r["run_accession"], "v0")
        rp = os.path.join(d, "real.report.json")
        sp = os.path.join(d, f"sim.v{r['variant']}.report.json")
        if not (os.path.exists(rp) and os.path.exists(sp)):
            continue
        real, sim = load(rp), load(sp)
        rr = real["metrics"]["reads"]
        sr = sim["metrics"]["reads"]
        r3, s3 = run_dist(real, "3p"), run_dist(sim, "3p")
        r5, s5 = run_dist(real, "5p"), run_dist(sim, "5p")
        s = float(r["gfrac_max_diff_scored"])
        lm = r["gfrac_max_diff_lengthmatched"]
        lm = float(lm) if lm else None
        prof = None
        pp = os.path.join(d, "real.profile.json")
        if os.path.exists(pp):
            prof = load(pp)
        ph = (prof or {}).get("head") or {}
        pt = (prof or {}).get("tail") or {}
        pl = ((prof or {}).get("base") or {}).get("R1") or []
        rl = int((prof or {}).get("read_len") or 0)
        out.append([
            r["domain"], r["run_accession"], r["variant"],
            "%.4f" % s, r["argmax_pos"], r["argmax_zone"],
            "" if lm is None else "%.4f" % lm, r["lm_argmax_pos"] or "",
            "" if lm is None else "%.4f" % (s - lm),
            r["n_overlay_max_absd"], r["body_mean_absd"], r["head_mean_absd"],
            r["tail_mean_absd"],
            "%.4f" % (share_ge1(r3, rr) or 0), "%.4f" % (share_ge1(s3, sr) or 0),
            "%.4f" % (tv(r3, s3) or 0),
            "%.4f" % (share_ge1(r5, rr) or 0), "%.4f" % (share_ge1(s5, sr) or 0),
            "%.4f" % (tv(r5, s5) or 0),
            r["real_len_r1"], r["real_len_r2"], r["real_asymmetric"],
            "%.4f" % run_ge1(ph.get("R1") or {}), "%.4f" % run_ge1(ph.get("R2") or {}),
            "%.4f" % run_ge1(pt.get("R1") or {}), "%.4f" % run_ge1(pt.get("R2") or {}),
            "%.4f" % (pl[rl - 1][2] if rl and len(pl) >= rl else 0.0),
        ])

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(FIELDS)
        w.writerows(out)
    print(f"wrote {args.out} ({len(out)} rows)")

    # a readable digest
    print(f"\n{'dataset':<13}{'v':<3}{'scored':<8}{'pos':<5}{'zone':<6}{'LM':<8}"
          f"{'mix':<8}{'3pG r/s':<18}{'tailTV':<8}{'5pG r/s':<18}{'headTV':<8}{'N':<8}")
    for row in out:
        print(f"{row[1]:<13}{row[2]:<3}{row[3]:<8}{row[4]:<5}{row[5]:<6}{row[6]:<8}"
              f"{row[8]:<8}{row[13]+'/'+row[14]:<18}{row[15]:<8}{row[16]+'/'+row[17]:<18}"
              f"{row[18]:<8}{row[9]:<8}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
