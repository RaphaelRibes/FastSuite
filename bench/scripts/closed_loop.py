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

POOLING vs LENGTH-MATCHED  (docs/GFRAC-MECHANISM.md §2, §7.1)
-------------------------------------------------------------
`base_position` is POOLED over both mates, and the profile carries ONE scalar
`read_len` that fastsim emits for BOTH mates ("Sequencer output is always
exactly read_len bases", fastsim main.rs). Real libraries are not symmetric:
the failing metagenomics set is R1 101 / R2 99, and DRR083073 is 101 / 91. At
a position only R1 reaches, the POOLED simulated value is a two-mate average,
so the gate scores the pooling convention, not the composition error. Measured:
DRR083073 is 0.1029 pooled and 0.0358 length-matched; for nine group-A datasets
the sign flips (0.028-0.041 pooled but 0.098-0.112 length-matched — the pooling
HIDES a 2.5-3x larger per-mate error).

With `--per-mate` the harness runs fastcheck on each mate's FASTQ (real and
simulated, cached), and scores two extra statistics:

  * `g_frac_lm_max_diff` — the LENGTH-MATCHED gate: the sim's two per-mate
    curves are pooled using the REAL per-mate depths, so the simulated library
    is compared to the real one at the real one's mate lengths. This is the
    primary statistic (docs/GFRAC-MECHANISM.md §7.1).
  * `g_frac_permate_max_diff` — max over both mates of the per-mate |diff|, the
    strictest reading: it never averages two mates into one number at all.

The per-mate reports are cached under --cache with the same names
gfrac_attribution.py uses ({domain}_{acc}_real_{m}.report.json /
{domain}_{acc}_v{variant}_{m}.report.json), so the two scripts share a cache.

Usage:
  closed_loop.py --fastcheck PATH --work ROOT --row closed_loop.tsv [--threads N]
                 [--per-mate --units UNITS_DIR --cache DIR] [--domains a,b]
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
              "verdict", "failing",
              # --- length-matched / per-mate block (empty unless --per-mate) ---
              "g_frac_lm_max_diff", "g_frac_lm_pos_checked",
              "g_frac_permate_max_diff", "g_frac_permate_pos",
              "verdict_lm", "failing_lm"]


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
            out.append(None)
            continue
        acgt = counts[0] + counts[1] + counts[2] + counts[3]
        out.append(counts[2] / acgt if acgt else None)
    return out


def depth_series(report: dict) -> list[int]:
    """Per-position total base count (pooled over both mates)."""
    bp = report.get("base_position", {}) or {}
    return [sum(bp[k]) for k in sorted(bp.keys(), key=lambda k: int(k))]


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


def max_diff(real: list, sim: list) -> tuple[float, int]:
    """max |real[i]-sim[i]| over positions where both are non-None (0-based idx)."""
    best, bi = 0.0, -1
    n = min(len(real), len(sim))
    for i in range(n):
        if real[i] is None or sim[i] is None:
            continue
        d = abs(real[i] - sim[i])
        if d > best:
            best, bi = d, i
    return best, bi


def length_matched_curve(sim_mate_g: list[list], real_mate_d: list[list]) -> list:
    """Pool the sim's per-mate curves using the REAL per-mate depths.

    The simulated library only ever has two equal-length mates; the real one
    may not. Weighting the sim's flexible per-mate curves by the real depths
    answers: what would the pooled simulation look like if it had the real
    library's mate lengths? (docs/GFRAC-MECHANISM.md §2.)
    """
    n = max((len(g) for g in sim_mate_g if g), default=0)
    out = [None] * n
    for i in range(n):
        num = den = 0.0
        for mg, md in zip(sim_mate_g, real_mate_d):
            if mg and md and i < len(mg) and i < len(md) and mg[i] is not None and md[i] > 0:
                num += md[i] * mg[i]
                den += md[i]
        if den:
            out[i] = num / den
    return out


def per_mate_max(real_mate_g: list[list], sim_mate_g: list[list]) -> tuple[float, int, int]:
    """max over (mate, position) of |real - sim|; returns (diff, 0-based pos, mate)."""
    best, bp, bm = 0.0, -1, -1
    for m, (rg, sg) in enumerate(zip(real_mate_g, sim_mate_g)):
        d, i = max_diff(rg, sg)
        if d > best:
            best, bp, bm = d, i, m
    return best, bp, bm


def compare(real: dict, sim: dict) -> dict:
    """The four A3 thresholds, plus the polyX rate as a diagnostic."""
    fails = []

    # ---- per-position G fraction ----
    r, s = g_frac_series(real), g_frac_series(sim)
    max_d, _ = max_diff(r, s)
    checked = sum(1 for i in range(min(len(r), len(s))) if r[i] is not None and s[i] is not None)
    if checked and max_d > TOL_G_FRAC:
        fails.append(f"g_frac {max_d:.4f}>{TOL_G_FRAC}")

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

    # Fail reasons OTHER than the pooled g_frac — reused verbatim by the
    # length-matched verdict, which swaps the g_frac term for the LM one.
    non_g = [f for f in fails if not f.startswith("g_frac ")]

    return dict(g_frac_max_diff=max_d, g_frac_pos_checked=checked,
                mean_q_real=qr, mean_q_sim=qs, mean_q_diff=qdiff,
                gc_real=gr, gc_sim=gs, gc_diff=gdiff,
                dup_rate_real=dr, dup_rate_sim=ds, dup_ratio=ratio,
                polyx5_real=pr, polyx5_sim=ps, polyx5_diff=pdiff,
                non_g_fails=non_g,
                verdict="PASS" if not fails else "FAIL",
                failing=";".join(fails))


def verdict_with_g(max_d: float | None, non_g_fails: list[str]) -> tuple[str, str]:
    """The verdict when the g_frac term is `max_d` (LM or per-mate)."""
    fails = list(non_g_fails)
    if max_d is not None and max_d > TOL_G_FRAC:
        fails.insert(0, f"g_frac {max_d:.4f}>{TOL_G_FRAC}")
    return ("PASS" if not fails else "FAIL"), ";".join(fails)


def fastcheck_mate_report(fc: str, path: str, cache: str, tag: str, threads: int) -> dict | None:
    """Per-position G curve + depth of ONE FASTQ, via fastcheck --json (cached).

    Cache names match gfrac_attribution.py so the two scripts share a cache.
    """
    os.makedirs(cache, exist_ok=True)
    rep = os.path.join(cache, f"{tag}.report.json")
    if not os.path.exists(rep):
        r = subprocess.run([fc, "--json", rep, "-q", "-j", str(threads), "--", path],
                           capture_output=True)
        if r.returncode != 0 or not os.path.exists(rep):
            print(f"  fastcheck failed on {path}: {r.stderr[-300:]!r}", file=sys.stderr)
            return None
    return load_report(rep)


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
    # --- length-matched / per-mate comparison (docs/GFRAC-MECHANISM.md §7.1) ---
    ap.add_argument("--per-mate", action="store_true",
                    help="also score the length-matched and per-mate g_frac gates")
    ap.add_argument("--units", default=None,
                    help="units root (<domain>/<acc>.R1.fastq.gz); required by --per-mate")
    ap.add_argument("--cache", default="/tmp/gfrac_attr_cache",
                    help="per-mate fastcheck report cache (shared with gfrac_attribution.py)")
    ap.add_argument("--domains", default=None,
                    help="comma-separated domain filter (default: all under --work)")
    args = ap.parse_args()

    if args.per_mate and not args.units:
        print("closed_loop: --per-mate needs --units", file=sys.stderr)
        return 2

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

    want = None
    if args.domains:
        want = {d.strip() for d in args.domains.split(",") if d.strip()}

    os.makedirs(os.path.dirname(args.row) or ".", exist_ok=True)
    with open(args.row, "w", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(ROW_FIELDS)

    n = 0
    for domain in sorted(os.listdir(args.work)):
        if want is not None and domain not in want:
            continue
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

            # ---- per-mate real curves/depths (computed once per dataset) ----
            real_mate_g = real_mate_d = None
            if args.per_mate:
                real_mate_g, real_mate_d = [], []
                ok = True
                for m in ("R1", "R2"):
                    rep = fastcheck_mate_report(
                        args.fastcheck,
                        os.path.join(args.units, domain, f"{acc}.{m}.fastq.gz"),
                        args.cache, f"{domain}_{acc}_real_{m}", args.threads)
                    if rep is None:
                        ok = False
                        break
                    real_mate_g.append(g_frac_series(rep))
                    real_mate_d.append(depth_series(rep))
                if not ok:
                    real_mate_g = real_mate_d = None

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
                real, sim = load_report(real_r), load_report(sim_r)
                m = compare(real, sim)

                # ---- length-matched and per-mate gates (spec §7.1) ----
                lm_diff = pm_diff = None
                lm_checked = 0
                pm_pos = ""
                if args.per_mate and real_mate_g is not None:
                    sim_mate_g = []
                    ok = True
                    for mm in ("R1", "R2"):
                        rep = fastcheck_mate_report(
                            args.fastcheck,
                            os.path.join(sdir, f"sim.v{variant}_{mm}.fastq.gz"),
                            args.cache, f"{domain}_{acc}_v{variant}_{mm}", args.threads)
                        if rep is None:
                            ok = False
                            break
                        sim_mate_g.append(g_frac_series(rep))
                    if ok:
                        rg_pool = g_frac_series(real)
                        lm_curve = length_matched_curve(sim_mate_g, real_mate_d)
                        lm_diff, _ = max_diff(rg_pool, lm_curve)
                        lm_checked = sum(1 for i in range(min(len(rg_pool), len(lm_curve)))
                                         if rg_pool[i] is not None and lm_curve[i] is not None)
                        pm_diff, pm_i, pm_m = per_mate_max(real_mate_g, sim_mate_g)
                        if pm_i >= 0:
                            pm_pos = f"{pm_i + 1}/{('R1', 'R2')[pm_m]}"

                v_lm, f_lm = verdict_with_g(lm_diff, m["non_g_fails"])
                row = [domain, acc, variant,
                       f"{m['g_frac_max_diff']:.4f}", m["g_frac_pos_checked"],
                       fmt(m["mean_q_real"]), fmt(m["mean_q_sim"]), fmt(m["mean_q_diff"]),
                       fmt(m["gc_real"]), fmt(m["gc_sim"]), fmt(m["gc_diff"]),
                       fmt(m["dup_rate_real"]), fmt(m["dup_rate_sim"]), fmt(m["dup_ratio"]),
                       fmt(m["polyx5_real"]), fmt(m["polyx5_sim"]), fmt(m["polyx5_diff"]),
                       m["verdict"], m["failing"],
                       fmt(lm_diff, "%.4f"), lm_checked,
                       fmt(pm_diff, "%.4f"), pm_pos,
                       (v_lm if lm_diff is not None else ""), (f_lm if lm_diff is not None else "")]
                with open(args.row, "a", newline="") as fh:
                    csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(row)
                n += 1
    print(f"closed_loop: {n} (dataset, variant) comparisons written to {args.row}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
