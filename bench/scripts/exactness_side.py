#!/usr/bin/env python3
"""exactness_side.py — where the overtrim actually is, per dataset.

The arm sweep (`exactness_sweep.py`) says HOW MUCH each rule moves the counts.
This answers WHERE the remaining error sits, which is what decides whether a
boundary rule can touch it at all:

  * 5' vs 3' vs both, for every kept mate whose cut overshoots its ground-truth
    boundary;
  * the overshoot distance histogram on the 5' side (1, 2, 3, 4, 5, 6+ bp) —
    a 1..=4 bp structure is the signature of `--head-dropout`, whose window is
    exactly that, and of nothing else in the tool;
  * how many of those 5' cuts land on a mate the ground truth calls clean
    (h = 0): those are false trims by construction.

Reads the same winner config as the sweep (phase-D `winners.tsv` through
`overtrim_vs_fastp`'s helpers), one fastprep run per dataset, PXMT v3 unchanged.

    exactness_side.py --jobs 12 --out bench/results/exactness/side.tsv
"""

import argparse
import collections
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import score_weighted as sw            # noqa: E402
import exactness_sweep as es          # noqa: E402

BUCKETS = (1, 2, 3, 4, 5, 6)          # 6 means "6 or more"


def run_one(task):
    fp, ds, r1, r2, gt, cfg = task
    d = tempfile.mkdtemp(dir=task[6], prefix="side-")
    try:
        px = os.path.join(d, "cuts.pxmt")
        argv = [fp, "-1", r1, "-2", r2, "-o", os.path.join(d, "o1.fq"),
                "-p", os.path.join(d, "o2.fq"), "--metadata", px, "-j", "1"] + cfg
        r = subprocess.run(argv, capture_output=True, text=True)
        if r.returncode != 0:
            return {"ds": ds, "error": (r.stderr or r.stdout)[-300:]}
        dec = sw.parse_pxmt(px)
        grow = sw.read_gt(gt)
        c = collections.Counter()
        hist = collections.Counter()
        for (pid, mate), (status, t5, t3, _rl) in dec.items():
            row = grow.get(pid)
            if row is None:
                continue
            c["mates"] += 1
            if status != 0:
                c["dropped"] += 1
                continue
            c["kept"] += 1
            h, t, _full = row[mate]
            a, b = t5 > h, t3 > t
            if a and b:
                c["over_both"] += 1
            elif a:
                c["over_5"] += 1
            elif b:
                c["over_3"] += 1
            else:
                c["clean_cut"] += 1
            if a:
                hist[min(t5 - h, 6)] += 1
                # a 5' cut on a mate the GT calls clean is a false trim
                c["over_5_gt_clean" if h == 0 else "over_5_gt_head"] += 1
        out = {"ds": ds, **{k: c[k] for k in
                            ("mates", "kept", "dropped", "over_5", "over_3",
                             "over_both", "clean_cut")},
               "over_5_gt_clean": c["over_5_gt_clean"],
               "over_5_gt_head": c["over_5_gt_head"]}
        for b in BUCKETS:
            out[f"overshoot_{b}"] = hist[b]
        return out
    except Exception as exc:
        return {"ds": ds, "error": f"{type(exc).__name__}: {exc}"}
    finally:
        shutil.rmtree(d, ignore_errors=True)


FIELDS = ["ds", "mates", "kept", "dropped", "over_5", "over_3", "over_both",
          "clean_cut", "over_5_gt_clean", "over_5_gt_head"] + \
         [f"overshoot_{b}" for b in BUCKETS]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus", default="/mnt/Fast/fastsuite-scratch/work-gt2026")
    ap.add_argument("--winners", default="bench/results/phase-d/winners.tsv")
    ap.add_argument("--fastprep", default="tools/fastprep/target/release/fastprep")
    ap.add_argument("--work", default="/mnt/Fast/fastsuite-scratch/work-exactness")
    ap.add_argument("--out", default="bench/results/exactness/side.tsv")
    ap.add_argument("--variants", default="v0")
    ap.add_argument("--jobs", type=int, default=12)
    a = ap.parse_args(argv)

    variants = [] if a.variants == "all" else a.variants.split(",")
    datasets = es.discover(a.corpus, variants)
    winners = es.load_winners(a.winners)
    fp = os.path.abspath(a.fastprep)
    os.makedirs(a.work, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    tasks = []
    for d in datasets:
        w = winners.get((d["dom"], d["acc"]))
        if w is None:
            continue
        tasks.append((fp, d["ds"], d["r1"], d["r2"], d["gt"],
                      es.winner_flags(w), a.work))
    print(f"{len(tasks)} datasets")

    rows, errors = [], 0
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        futs = [pool.submit(run_one, t) for t in tasks]
        for i, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            if "error" in rec:
                errors += 1
                print(f"  FAIL {rec['ds']}: {rec['error'][:200]}")
            else:
                rows.append(rec)
            if i % 25 == 0:
                print(f"  {i}/{len(tasks)}", flush=True)

    rows.sort(key=lambda r: r["ds"])
    with open(a.out, "w") as fh:
        fh.write("\t".join(FIELDS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k, 0)) for k in FIELDS) + "\n")
        tot = {k: sum(r.get(k, 0) for r in rows) for k in FIELDS[1:]}
        fh.write("\t".join(["TOTAL"] + [str(tot[k]) for k in FIELDS[1:]]) + "\n")

    mate = tot["mates"] or 1
    print(f"\n{len(rows)} datasets, {errors} failed -> {a.out}")
    print(f"mates {tot['mates']}, kept {tot['kept']}")
    print(f"overtrim: 5' {tot['over_5']}  3' {tot['over_3']}  "
          f"both {tot['over_both']}  (of which 5' on GT-clean mates: "
          f"{tot['over_5_gt_clean']})")
    print("5' overshoot bp: " + ", ".join(
        f"{b if b < 6 else '6+'}: {tot[f'overshoot_{b}']}" for b in BUCKETS))


if __name__ == "__main__":
    main()
