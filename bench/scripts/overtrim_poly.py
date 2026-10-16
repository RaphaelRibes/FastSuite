#!/usr/bin/env python3
"""Clean poly-only overtrim measurement + head-dropout policy trade.

Separates polyX trim (t5p/t3p) from adapter trim (t5a/t3a). The correct poly
answer is (h, t); adapter is a separate, correct trim that the older scorer
folded into the comparison (inflating 3' 'overtrim'). We score POLY ONLY.

Also, to size the algorithmic room on the dominant mechanism (a 1 bp leading
biological G trimmed by the head-dropout rule on reads whose GT h = 0), we run
fastprep under several head-dropout settings and report the exact/over/under
trade at tolerance 0 and 2.
"""
from __future__ import annotations
import csv, gzip, struct, subprocess, sys, os, collections

FASTPREP = "/home/raphael/hermes_scratch/FastSuite/tools/fastprep/target/release/fastprep"
WORKROOT = "/home/raphael/.hermes/cache/scratch/overtrim"

def parse_pxmt(path):
    """-> {(pid, mate): (status, t5a, t3a, t5p, t3p)}"""
    out = {}
    buf = open(path, "rb").read()
    if buf[:4] != b"PXMT":
        return out
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
        pid = int(parts[-1]) * 10_000 + int(parts[-2])
        out[(pid, mate)] = (status, t5a, t3a, t5p, t3p)
    return out

def parse_gt(path):
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[int(r["pair_id"])] = (int(r["r1_h"]), int(r["r1_t"]),
                                      int(r["r2_h"]), int(r["r2_t"]))
    return out

def classify(p5, p3, h, t, tol):
    if abs(p5 - h) <= tol and abs(p3 - t) <= tol:
        return "exact"
    if p5 > h or p3 > t:
        return "over"
    return "under"

def score(px, gt, tol):
    c = collections.Counter()
    dropped = 0
    for pid, (h1, t1, h2, t2) in gt.items():
        for mate, (h, t) in ((0, (h1, t1)), (1, (h2, t2))):
            s = px.get((pid, mate))
            if s is None or s[0] != 0:
                dropped += 1
                continue
            _, _t5a, _t3a, p5, p3 = s
            c[classify(p5, p3, h, t, tol)] += 1
    n = max(1, c["exact"] + c["over"] + c["under"])
    return c["exact"] / n, c["over"] / n, c["under"] / n, n, dropped

def run(work, hd, tag):
    out = f"{work}/_pol_{tag}.pxmt"
    argv = [FASTPREP, "-1", f"{work}/sim.v0_R1.fastq.gz", "-2", f"{work}/sim.v0_R2.fastq.gz",
            "-o", f"{work}/_p1.fastq.gz", "-p", f"{work}/_p2.fastq.gz",
            "--metadata", out, "--json", f"{work}/_p.json", "--force", "-q", "-j", "4",
            "--min-run", "12", "--mercy", "3", "--phred-sat", "20", "--gate-shape", "1.0",
            "--min-len", "0.9L", "--mean-q", "25.0", "--adapter-error-rate", "0.2",
            "--min-adapter-overlap", "14", "--gz-level", "6",
            "--match-leniency", "0.1", "--mismatch-strength", "0.5",
            "--head-dropout", str(hd)]
    p = subprocess.run(argv, capture_output=True, text=True)
    if p.returncode != 0:
        print(f"  hd={hd}: rc={p.returncode} {p.stderr[-160:]}")
        return None
    return parse_pxmt(out)

def main():
    regimes = sys.argv[1:] or ["low", "mid", "high"]
    for regime in regimes:
        work = f"{WORKROOT}/{regime}"
        gt = parse_gt(f"{work}/sim.v0.gt.tsv")
        print(f"\n===== {regime} ===== (POLY-only scoring)")
        print(f"{'head-dropout':>13s} | {'exact0':>7s} {'over0':>7s} {'und0':>7s} | {'exact2':>7s} {'over2':>7s} {'und2':>7s}")
        for hd in (0, 1, 2, 4):
            px = run(work, hd, f"hd{hd}")
            if px is None:
                continue
            e0, o0, u0, n, dr = score(px, gt, 0)
            e2, o2, u2, _, _ = score(px, gt, 2)
            print(f"hd={hd:>11d} | {e0:7.4f} {o0:7.4f} {u0:7.4f} | {e2:7.4f} {o2:7.4f} {u2:7.4f}   (n={n} dropped={dr})")

if __name__ == "__main__":
    main()
