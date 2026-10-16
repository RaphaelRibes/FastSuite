#!/usr/bin/env python3
"""Is there ANY point in (min_run, head-dropout, mercy, min_len) that fixes the
5' overtrim? Exhaustive 2-D sweep on the repo's fitted profile.

The 5' overtrim has one structural cause: the leading-G rule fires on a 1 bp
G-run. The useful question is whether ANY configuration of the existing knobs
avoids that, not whether a few hand-picked points do.
"""
from __future__ import annotations

import csv
import itertools
import os
import subprocess
import sys

FASTPREP = "/home/raphael/hermes_scratch/FastSuite/tools/fastprep/target/release/fastprep"
WORK = "/home/raphael/.hermes/cache/scratch/overtrim/mandril"


def parse_pxmt(path: str):
    import struct
    out = {}
    buf = open(path, "rb").read()
    off, n = 5, len(buf)
    while off < n:
        (id_len,) = struct.unpack_from("<H", buf, off); off += 2
        rid = buf[off:off + id_len]; off += id_len
        mate = buf[off]; off += 1
        status = buf[off]; off += 1
        t5a, t3a, t5p, t3p, _rl = struct.unpack_from("<IIIII", buf, off); off += 20
        nr = buf[off]; off += 1
        off += 8 * nr
        parts = rid.split(b" ", 1)[0].decode("ascii").lstrip("@").split(":")
        out[(int(parts[-2]) + 10_000 * int(parts[-1]), mate)] = (status, t5a + t5p, t3a + t3p)
    return out


def parse_gt(path: str):
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[int(r["pair_id"])] = (int(r["r1_h"]), int(r["r1_t"]),
                                      int(r["r2_h"]), int(r["r2_t"]))
    return out


def score(px, gt, tol):
    c = {"exact": 0, "overtrim": 0, "undertrim": 0, "o5": 0, "o3": 0}
    for pid, (h1, t1, h2, t2) in gt.items():
        for mate, (h, t) in ((0, (h1, t1)), (1, (h2, t2))):
            s = px.get((pid, mate))
            if s is None:
                continue
            st, t5, t3 = s
            if st != 0:
                continue
            if abs(t5 - h) <= tol and abs(t3 - t) <= tol:
                c["exact"] += 1
            elif t5 > h or t3 > t:
                c["overtrim"] += 1
            else:
                c["undertrim"] += 1
            if t5 > h:
                c["o5"] += 1
            if t3 > t:
                c["o3"] += 1
    n = max(1, sum(c[k] for k in ("exact", "overtrim", "undertrim")))
    return c["exact"] / n, c["overtrim"] / n, c["undertrim"] / n, c["o5"] / n, c["o3"] / n


def main() -> int:
    gt = parse_gt(os.path.join(WORK, "sim.v0.gt.tsv"))
    results = []
    # min_run 4..30, head_dropout 1..(min_run-1), mercy {3}, min_len 0.9
    for mr in list(range(4, 31)):
        for hd in range(1, min(mr, 12)):
            argv = [FASTPREP, "-1", f"{WORK}/sim.v0_R1.fastq.gz",
                    "-2", f"{WORK}/sim.v0_R2.fastq.gz",
                    "-o", f"{WORK}/g1.fq.gz", "-p", f"{WORK}/g2.fq.gz",
                    "--metadata", f"{WORK}/g.pxmt", "--force", "-q", "-j", "4",
                    "--min-run", str(mr), "--head-dropout", str(hd),
                    "--mercy", "3", "--min-len", "0.9L", "--mean-q", "25",
                    "--phred-sat", "20", "--gz-level", "6"]
            p = subprocess.run(argv, capture_output=True, text=True)
            if p.returncode != 0:
                continue
            e0, o0, u0, o5, o3 = score(parse_pxmt(f"{WORK}/g.pxmt"), gt, 0)
            e2, o2, u2, _, _ = score(parse_pxmt(f"{WORK}/g.pxmt"), gt, 2)
            results.append(dict(min_run=mr, head_dropout=hd, exact0=e0, over0=o0,
                                under0=u0, exact2=e2, over2=o2, o5=o5, o3=o3))
    with open(f"{WORK}/grid2d.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(results)
    print(f"{len(results)} points")
    print("\nBest 12 by exact0 (any min_run x head_dropout):")
    for r in sorted(results, key=lambda r: -r["exact0"])[:12]:
        print(f"  min_run={r['min_run']:>2d} head_dropout={r['head_dropout']:>2d} | "
              f"exact0={r['exact0']:.4f} over0={r['over0']:.4f} und0={r['under0']:.4f} | "
              f"exact2={r['exact2']:.4f} over2={r['over2']:.4f} | 5'over={r['o5']:.4f} 3'over={r['o3']:.4f}")
    print("\nBest 8 by over0 then exact0:")
    for r in sorted(results, key=lambda r: (r["over0"], -r["exact0"]))[:8]:
        print(f"  min_run={r['min_run']:>2d} head_dropout={r['head_dropout']:>2d} | "
              f"exact0={r['exact0']:.4f} over0={r['over0']:.4f} und0={r['under0']:.4f} | "
              f"exact2={r['exact2']:.4f} over2={r['over2']:.4f} | 5'over={r['o5']:.4f} 3'over={r['o3']:.4f}")
    # The theoretical ceiling: never trim 5'
    print("\nCeiling check — the BEST head_dropout is the smallest one:")
    for hd in (1, 2, 3):
        rs = [r for r in results if r["head_dropout"] == hd]
        best = max(rs, key=lambda r: r["exact0"])
        print(f"  hd={hd}: best exact0={best['exact0']:.4f} over0={best['over0']:.4f} (min_run={best['min_run']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
