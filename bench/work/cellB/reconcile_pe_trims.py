#!/usr/bin/env python3
"""Cell B: reconcile pe_trims with the scored-class movement it produces.

`pe_trims` (trim_report.py) counts (pair, mate) whose (status, t5, t3) differ
from the fastprep-hd0 baseline. The scored counters only move when the CLASS
also changes. This script, on the full 400k-pair cohort, joins both PXMTs to
the GT and reports:

  * how many changed decisions moved a scored class, and which way;
  * the net change in each class, which must reconcile against the report's
    own `cellB-with-pe.tsv` minus `cellB-without-pe.tsv` deltas.

The 20k-cohort probe could not answer this: the report deletes its scratch, so
the PXMTs it used are gone, and re-running the probe standalone reproduces only
a subset of the decisions.

Usage:
    python3 reconcile_pe_trims.py <gt.tsv> <R1.fastq.gz> <R2.fastq.gz> <outdir>
"""
from __future__ import annotations

import collections
import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from explain_pe_effect import parse_pxmt, _key, _classify, BASE_ARGS, FP  # noqa: E402


def run(fastq1: str, fastq2: str, outdir: str, tag: str, extra: list[str]):
    px = os.path.join(outdir, f"{tag}.pxmt")
    subprocess.run([FP, "-1", fastq1, "-2", fastq2,
                    "-o", os.path.join(outdir, f"{tag}_1.fq.gz"),
                    "-p", os.path.join(outdir, f"{tag}_2.fq.gz"),
                    "--metadata", px] + BASE_ARGS + extra, check=True,
                   capture_output=True)
    return parse_pxmt(px)


def main() -> None:
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    gt_path, r1, r2, outdir = sys.argv[1:5]
    os.makedirs(outdir, exist_ok=True)
    gt = {int(r["pair_id"]): r
          for r in csv.DictReader(open(gt_path), delimiter="\t")}
    base = run(r1, r2, outdir, "base", [])
    pe = run(r1, r2, outdir, "pe", ["--pe-overlap"])

    c = collections.Counter()
    for (rid, mate), b in base.items():
        p = pe.get((rid, mate))
        if p is None or b == p:
            continue
        row = gt.get(_key(rid))
        if row is None:
            c["changed (no GT row)"] += 1
            continue
        m = "r1" if mate == 1 else "r2"
        h, t = int(row[f"{m}_h"]), int(row[f"{m}_t"])
        c[f"changed (mate{mate})"] += 1
        cb = "dropped" if b[0] != 0 else _classify(b[1], b[2], h, t)
        cp = "dropped" if p[0] != 0 else _classify(p[1], p[2], h, t)
        if cb == cp:
            c[f"  class-neutral ({cb})"] += 1
        else:
            c[f"  moved {cb} -> {cp}"] += 1
    print("=== Cell B (400k): changed decisions vs scored-class movement ===")
    for k in sorted(c):
        print(f"  {k:36} {c[k]}")
    print()
    print("  Reconciliation: the sum of 'moved X -> Y' per class must equal the")
    print("  report's cellB-with-pe.tsv minus cellB-without-pe.tsv deltas.")


if __name__ == "__main__":
    main()
