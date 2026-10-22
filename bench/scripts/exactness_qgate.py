#!/usr/bin/env python3
"""exactness_qgate.py — can quality tell a dark-cycle 5' G-run from biology?

The exactness campaign found the whole tol-0 overtrim is the 5' G-run class, and
`--head-dropout` cuts it by POLICY (it ignores quality by design). The obvious
rescue would be a quality gate: only trim a leading G-run the sequencer was
unconfident about. This script tests whether that gate has anything to see, by
comparing the mean quality of the leading 1..4 bp G-run between the mates whose
ground-truth window is the whole read (biology) and the mates inside an artefact
window (where the run really is a dark cycle).

Output: one row per dataset x class x run length with count and Q quantiles, plus
a TOTAL row per class. Emits `bench/results/exactness/qgate.tsv`; a gate is
useful only if the two classes separate, and a dataset where they do not is the
finding, not a failure.

    exactness_qgate.py --jobs 12
"""

import argparse
import csv
import glob
import gzip
import os
import statistics
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

FIELDS = ["dataset", "gt_class", "run_len", "reads", "mean_q", "p10", "p50", "p90"]


def quantiles(vals):
    v = sorted(vals)
    n = len(v)
    return (v[n // 10], v[n // 2], v[9 * n // 10]) if n else (0, 0, 0)


def scan(task):
    gt_path, fq_path = task
    gt = {}
    with open(gt_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gt[int(r["pair_id"])] = int(r["r1_h"])
    buckets = defaultdict(list)
    with gzip.open(fq_path, "rt") as fh:
        while True:
            head = fh.readline()
            if not head:
                break
            seq = fh.readline().strip()
            fh.readline()
            qual = fh.readline().strip()
            fields = head.split()[0].split(":")
            if len(fields) < 7:
                continue
            pid = int(fields[5]) + 10000 * int(fields[6])
            h = gt.get(pid)
            if h is None:
                continue
            run = 0
            while run < len(seq) and seq[run] == "G":
                run += 1
            if not 1 <= run <= 4:
                continue
            # a dark cycle and a sequencing error are both low Q; biology is not
            buckets[("gt_h=0" if h == 0 else "gt_h>0", run)].append(
                statistics.mean(ord(c) - 33 for c in qual[:run]))
    out = []
    for (cls, run), vals in sorted(buckets.items()):
        p10, p50, p90 = quantiles(vals)
        out.append({"dataset": os.path.basename(os.path.dirname(gt_path)),
                    "gt_class": cls, "run_len": run, "reads": len(vals),
                    "mean_q": round(statistics.mean(vals), 3),
                    "p10": p10, "p50": p50, "p90": p90})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus", default="/mnt/Fast/fastsuite-scratch/work-gt2026")
    ap.add_argument("--out", default="bench/results/exactness/qgate.tsv")
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--variants", default="v0")
    a = ap.parse_args(argv)

    tasks = []
    for gt in sorted(glob.glob(os.path.join(a.corpus, "**", "*.gt.tsv"), recursive=True)):
        stem = gt[: -len(".gt.tsv")]
        if a.variants != "all" and not stem.endswith(tuple(a.variants.split(","))):
            continue
        fq = stem + "_R1.fastq.gz"
        if os.path.exists(fq):
            tasks.append((gt, fq))
    print(f"{len(tasks)} datasets")
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    rows = []
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        for i, part in enumerate(pool.map(scan, tasks), 1):
            rows.extend(part)
            if i % 20 == 0:
                print(f"  {i}/{len(tasks)}", flush=True)

    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"\nwrote {a.out} ({len(rows)} rows)")
    for cls in ("gt_h=0", "gt_h>0"):
        sel = [r for r in rows if r["gt_class"] == cls]
        reads = sum(r["reads"] for r in sel)
        mq = sum(r["mean_q"] * r["reads"] for r in sel) / reads
        print(f"  {cls}: {reads:,} reads, weighted mean Q {mq:.2f}, "
              f"median-of-datasets p50 {statistics.median([r['p50'] for r in sel]):.1f}")


if __name__ == "__main__":
    main()
