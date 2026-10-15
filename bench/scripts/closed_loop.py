#!/usr/bin/env python3
"""C2 — the closed loop on real data: does the simulated QC match the real QC?

This is GATE 2's input. The whole design rests on one premise: a profile fitted
from a real QC reproduces that dataset's QC closely enough that sweeping it
means something. A2/A3 tested that on synthetic data; C2 tests it on the real
thing, on all 30 pilot datasets, with the same thresholds as A3:

    per-position G fraction   within 0.02
    mean Phred Q              within 0.5
    GC fraction               within 0.01
    duplicate rate            within 2x

fastcheck's report schema (fastsuite.report.v1) is read at its real paths —
metrics.mean_read_mean_q, metrics.gc_fraction, metrics.duplication_level, and
base_position["<1-based pos>"] -> [A, C, G, T] counts (both mates pooled, which
is what the profile is fitted from).

Usage:
  closed_loop.py --fastcheck PATH --work ROOT --row closed_loop.tsv [--threads N]
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys

# A3's thresholds — the same numbers the recovery test was held to.
TOL_G_FRAC = 0.02
TOL_MEAN_Q = 0.5
TOL_GC = 0.01
TOL_DUP = 2.0  # ratio, both directions

ROW_FIELDS = ["domain", "run_accession", "variant",
              "g_frac_max_diff", "g_frac_pos_checked",
              "mean_q_real", "mean_q_sim", "mean_q_diff",
              "gc_real", "gc_sim", "gc_diff",
              "dup_rate_real", "dup_rate_sim", "dup_ratio",
              "polyx5_real", "polyx5_sim", "polyx5_diff",
              "verdict", "failing"]


def load_report(path: str) -> dict:
    with open(path) as fh:
        return json.load(fh)


def metrics(report: dict) -> dict:
    return report.get("metrics", {}) or {}


def g_frac_series(report: dict) -> list[float]:
    """Per-position G fraction over ACGT bases, 1-based positions in order.

    base_position maps "1".."N" -> [A, C, G, T] counts pooled over both mates.
    """
    bp = report.get("base_position", {}) or {}
    out = []
    for pos in sorted(bp.keys(), key=lambda k: int(k)):
        counts = bp[pos]
        if len(counts) < 4:
            continue
        acgt = counts[0] + counts[1] + counts[2] + counts[3]
        out.append(counts[2] / acgt if acgt else 0.0)
    return out


def mean_q(report: dict) -> float | None:
    v = metrics(report).get("mean_read_mean_q")
    return float(v) if isinstance(v, (int, float)) else None


def gc_frac(report: dict) -> float | None:
    v = metrics(report).get("gc_fraction")
    return float(v) if isinstance(v, (int, float)) else None


def dup_rate(report: dict) -> float | None:
    v = metrics(report).get("duplication_level")
    return float(v) if isinstance(v, (int, float)) else None


def polyx5_frac(report: dict) -> float | None:
    v = metrics(report).get("polyx5_frac_gt_1")
    return float(v) if isinstance(v, (int, float)) else None


def compare(real: dict, sim: dict) -> dict:
    """The four A3 thresholds, plus the polyX rate as a diagnostic."""
    fails = []

    # ---- per-position G fraction ----
    r, s = g_frac_series(real), g_frac_series(sim)
    max_diff, checked = 0.0, 0
    if r and s:
        n = min(len(r), len(s))
        for i in range(n):
            d = abs(r[i] - s[i])
            max_diff = max(max_diff, d)
        checked = n
        if max_diff > TOL_G_FRAC:
            fails.append(f"g_frac {max_diff:.4f}>{TOL_G_FRAC}")

    # ---- mean Q ----
    qr, qs = mean_q(real), mean_q(sim)
    qdiff = abs(qr - qs) if (qr is not None and qs is not None) else None
    if qdiff is not None and qdiff > TOL_MEAN_Q:
        fails.append(f"mean_q {qdiff:.3f}>{TOL_MEAN_Q}")

    # ---- GC ----
    gr, gs = gc_frac(real), gc_frac(sim)
    gdiff = abs(gr - gs) if (gr is not None and gs is not None) else None
    if gdiff is not None and gdiff > TOL_GC:
        fails.append(f"gc {gdiff:.4f}>{TOL_GC}")

    # ---- duplicate rate (a ratio, both directions) ----
    dr, ds = dup_rate(real), dup_rate(sim)
    ratio = None
    if dr is not None and ds is not None:
        hi, lo = max(dr, ds), min(dr, ds)
        ratio = (hi / lo) if lo > 0 else float("inf")
        if ratio > TOL_DUP:
            fails.append(f"dup_ratio {ratio:.2f}>{TOL_DUP}")

    # ---- polyX5 rate: a diagnostic, not a gate ----
    pr, ps = polyx5_frac(real), polyx5_frac(sim)
    pdiff = abs(pr - ps) if (pr is not None and ps is not None) else None

    # A metric that could not be read at all is a failure of the comparison,
    # not a pass: say so rather than reporting a green verdict on missing data.
    missing = [n for n, v in (("mean_q", qdiff), ("gc", gdiff), ("dup", ratio)) if v is None]
    if missing:
        fails.append("unreadable:" + ",".join(missing))

    return dict(g_frac_max_diff=max_diff, g_frac_pos_checked=checked,
                mean_q_real=qr, mean_q_sim=qs, mean_q_diff=qdiff,
                gc_real=gr, gc_sim=gs, gc_diff=gdiff,
                dup_rate_real=dr, dup_rate_sim=ds, dup_ratio=ratio,
                polyx5_real=pr, polyx5_sim=ps, polyx5_diff=pdiff,
                verdict="PASS" if not fails else "FAIL",
                failing=";".join(fails))


def fmt(v, spec="%.6g"):
    return spec % v if isinstance(v, (int, float)) else ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastcheck", required=True)
    ap.add_argument("--work", required=True, help="root holding <domain>/<acc>/")
    ap.add_argument("--row", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--force", action="store_true",
                    help="overwrite --row if it already exists (otherwise refuse)")
    args = ap.parse_args()

    # This script OVERWRITES the file it was given, it does not append.
    #
    # It used to append whenever --row existed, so every re-run silently doubled
    # the table: the 306-comparison truth became 612 rows after one repeat run,
    # and the pass rate was computed over rows that had been counted twice. A
    # generator is idempotent by construction (it re-reads the same reports and
    # recomputes the same verdicts), so re-running must reproduce the file, not
    # grow it.
    #
    # Refusing by default (rather than silently truncating) matters because the
    # old rows are the only record of a run that used different binaries.
    if os.path.exists(args.row) and not args.force:
        n_old = sum(1 for _ in open(args.row, errors="ignore"))
        print(f"closed_loop: {args.row} exists ({n_old} lines) — refusing to append. "
              f"Pass --force to overwrite it, or pick a new --row path.", file=sys.stderr)
        return 2

    os.makedirs(os.path.dirname(args.row) or ".", exist_ok=True)
    with open(args.row, "w", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(ROW_FIELDS)

    n = 0
    for domain in sorted(os.listdir(args.work)):
        dpath = os.path.join(args.work, domain)
        if not os.path.isdir(dpath):
            continue
        for acc in sorted(os.listdir(dpath)):
            adir = os.path.join(dpath, acc)
            # The REAL QC report sits at the dataset root only in the old
            # (shared-dir) layout. run_campaign.sh's per-variant layout puts
            # every artefact under v<variant>/, so the real report lives in
            # v0/ there. Without this fallback the loop silently compared only
            # the pilot-era datasets and reported 69 comparisons instead of 306.
            real_r = os.path.join(adir, "real.report.json")
            if not os.path.exists(real_r):
                alt = os.path.join(adir, "v0", "real.report.json")
                if os.path.exists(alt):
                    real_r = alt
            if not os.path.exists(real_r):
                continue
            for variant in range(3):
                # the sim outputs live either in the shared dir (old layout)
                # or under v<variant>/ (per-variant layout, current)
                vdir = os.path.join(adir, f"v{variant}")
                sdir = vdir if os.path.isdir(vdir) else adir
                sim_r = os.path.join(sdir, f"sim.v{variant}.report.json")
                if not os.path.exists(sim_r):
                    r1 = os.path.join(sdir, f"sim.v{variant}_R1.fastq.gz")
                    r2 = os.path.join(sdir, f"sim.v{variant}_R2.fastq.gz")
                    if not (os.path.exists(r1) and os.path.exists(r2)):
                        continue
                    subprocess.run([args.fastcheck, "--json", sim_r, "-q",
                                    "-j", str(args.threads), "--", r1, r2],
                                   capture_output=True)
                if not os.path.exists(sim_r):
                    continue
                m = compare(load_report(real_r), load_report(sim_r))
                row = [domain, acc, variant,
                       f"{m['g_frac_max_diff']:.4f}", m["g_frac_pos_checked"],
                       fmt(m["mean_q_real"]), fmt(m["mean_q_sim"]), fmt(m["mean_q_diff"]),
                       fmt(m["gc_real"]), fmt(m["gc_sim"]), fmt(m["gc_diff"]),
                       fmt(m["dup_rate_real"]), fmt(m["dup_rate_sim"]), fmt(m["dup_ratio"]),
                       fmt(m["polyx5_real"]), fmt(m["polyx5_sim"]), fmt(m["polyx5_diff"]),
                       m["verdict"], m["failing"]]
                with open(args.row, "a", newline="") as fh:
                    csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(row)
                n += 1
    print(f"closed_loop: {n} (dataset, variant) comparisons written to {args.row}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
