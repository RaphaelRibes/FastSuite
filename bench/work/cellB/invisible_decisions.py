#!/usr/bin/env python3
"""Cell B: what is the 4866-counted 'dropped -> dropped' change?

verify_deltas.py found 4866 changed decisions whose class did not move, all
"dropped -> dropped" -- i.e. the tool dropped the mate in both arms but with
different trim positions recorded on the drop record. That is a real decision
change that the scorer cannot see, and it is worth naming.

This script breaks those 4866 down: did the drop's recorded t3 shrink (the
mercy-walk boundary moving), and is the mate an adapter mate?
"""
from __future__ import annotations

import collections
import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from explain_pe_effect import parse_pxmt, _key, BASE_ARGS, FP  # noqa: E402


def run(fastq1, fastq2, outdir, tag, extra):
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
            continue
        m = "r1" if mate == 1 else "r2"
        ad = int(row[f"{m}_adapter"])
        c["changed"] += 1
        if ad > 0:
            c["changed on an adapter mate"] += 1
        if b[0] != 0 and p[0] != 0:
            c["dropped in both arms"] += 1
            if p[2] < b[2]:
                c["  drop record t3 shrank"] += 1
            elif p[2] > b[2]:
                c["  drop record t3 grew"] += 1
            else:
                c["  drop record t3 unchanged (t5/status only)"] += 1
        elif p[2] < b[2]:
            c["t3 shrank"] += 1
    print("=== Cell B: the invisible decisions ===")
    for k in sorted(c):
        print(f"  {k:46} {c[k]}")


if __name__ == "__main__":
    main()
