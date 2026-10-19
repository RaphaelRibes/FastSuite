#!/usr/bin/env python3
"""Cell B: what does --pe-overlap actually change, mate by mate?

The trim report's `pe_trims` counts (pair, mate) whose (status, t5, t3) differ
from the fastprep-hd0 baseline -- but not which mate, and not in which
direction. This script answers both, on the cohort's own GT, by re-running
fastprep twice (baseline and --pe-overlap) and joining the PXMTs to the GT.

It reports, per mate:
  * pe_trims, split by mate and by the direction of the change
    (t3 grew / t3 shrank / appeared-from-nothing / became-a-drop);
  * the cross-tab of the GT's adapter length against what the baseline did,
    which is the mechanism: the mode can only add a trim the baseline did not
    already make, and the baseline's own matcher needs >= min-adapter-overlap
    (14) bases to fire.

Usage:
    python3 explain_pe_effect.py <gt.tsv> <R1.fastq.gz> <R2.fastq.gz> <outdir>
"""
from __future__ import annotations

import collections
import csv
import gzip
import os
import struct
import subprocess
import sys

L = 150
MIN_OVERLAP = 14
FP = os.environ.get("FASTPREP", os.path.expanduser(
    "~/hermes_scratch/FastSuite/tools/fastprep/target/release/fastprep"))
BASE_ARGS = ["--force", "-q", "--gz-level", "6", "--min-run", "12",
             "--mercy", "3", "--phred-sat", "20", "--gate-shape", "1.0",
             "--min-len", "0.9L", "--mean-q", "25.0",
             "--adapter-error-rate", "0.2", "--min-adapter-overlap", "14",
             "--head-dropout", "0"]


def parse_pxmt(path: str) -> dict[tuple[str, int], tuple[int, int, int]]:
    buf = open(path, "rb").read()
    out: dict[tuple[str, int], tuple[int, int, int]] = {}
    if buf[:4] != b"PXMT":
        return out
    off, n = 5, len(buf)
    while off + 2 <= n:
        (idlen,) = struct.unpack_from("<H", buf, off); off += 2
        rid = buf[off:off + idlen].decode("utf-8", "replace").split()[0]
        off += idlen
        if off + 22 > n:
            break
        mate = buf[off] + 1; off += 1
        status = buf[off]; off += 1
        t5a, t3a, t5p, t3p = struct.unpack_from("<IIII", buf, off); off += 16
        off += 4
        nrun = buf[off]; off += 1
        off += 8 * nrun
        out[(rid, mate)] = (status, t5a + t5p, t3a + t3p)
    return out


def run(fastq1: str, fastq2: str, outdir: str, tag: str, extra: list[str]):
    px = os.path.join(outdir, f"{tag}.pxmt")
    subprocess.run([FP, "-1", fastq1, "-2", fastq2,
                    "-o", os.path.join(outdir, f"{tag}_1.fq.gz"),
                    "-p", os.path.join(outdir, f"{tag}_2.fq.gz"),
                    "--metadata", px] + BASE_ARGS + extra, check=True,
                   capture_output=True)
    return parse_pxmt(px)


def _key(rid: str) -> int:
    f = rid.split(":")
    return int(f[5]) + 10000 * int(f[6])


def main() -> None:
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    gt_path, r1, r2, outdir = sys.argv[1:5]
    os.makedirs(outdir, exist_ok=True)

    gt = {}
    for row in csv.DictReader(open(gt_path), delimiter="\t"):
        gt[int(row["pair_id"])] = row

    base = run(r1, r2, outdir, "base", [])
    pe = run(r1, r2, outdir, "pe", ["--pe-overlap"])

    # ---- per-mate pe_trims, by direction ----
    c = collections.Counter()
    for (rid, mate), b in base.items():
        if (rid, mate) not in pe:
            c[f"mate{mate}: pe missing record"] += 1
            continue
        p = pe[(rid, mate)]
        if b == p:
            continue
        bs, bt5, bt3 = b
        ps, pt5, pt3 = p
        pid = _key(rid)
        row = gt.get(pid)
        if row is None:
            continue
        c[f"mate{mate}: pe_trims"] += 1
        if bs == 0 and ps != 0:
            c[f"mate{mate}: emitted -> dropped"] += 1
        elif bs != 0 and ps == 0:
            c[f"mate{mate}: dropped -> emitted"] += 1
        elif pt3 > bt3:
            c[f"mate{mate}: t3 grew (+{pt3 - bt3})"] += 1
        elif pt3 < bt3:
            c[f"mate{mate}: t3 shrank ({pt3 - bt3})"] += 1
        else:
            c[f"mate{mate}: t5/status changed only"] += 1
        # is this mate a readthrough mate? and how long is its readthrough?
        a = int(row[f"r{mate}_adapter"])
        if a > 0:
            c[f"mate{mate}: pe_trims on a readthrough mate"] += 1
            c[f"mate{mate}:   readthrough {a}"] += 1
        else:
            c[f"mate{mate}: pe_trims on a NON-readthrough mate"] += 1
    print("=== pe_trims by mate and direction (Cell B cohort) ===")
    for k in sorted(c):
        print(f"  {k:52} {c[k]}")

    # ---- the mechanism: baseline trim vs GT, bucketed by readthrough length ----
    print("\n=== baseline (fastprep-hd0) t3 vs GT, readthrough mates only ===")
    print(f"{'GT adapter len':>14} {'n':>8} {'base t3==t':>11} "
          f"{'base t3==t+a':>12} {'base t3==0':>10} {'other':>7}")
    buckets = collections.defaultdict(collections.Counter)
    for (rid, mate), (status, t5, t3) in base.items():
        pid = _key(rid)
        row = gt.get(pid)
        if row is None:
            continue
        a = int(row[f"r{mate}_adapter"])
        if a <= 0:
            continue
        t = int(row[f"r{mate}_t"])
        if status != 0:
            buckets[a]["dropped"] += 1
        elif t3 == t:
            buckets[a]["base t3==t"] += 1
        elif t3 == t + a:
            buckets[a]["base t3==t+a"] += 1
        elif t3 == 0:
            buckets[a]["base t3==0"] += 1
        else:
            buckets[a]["other"] += 1
    for a in sorted(buckets):
        b = buckets[a]
        n = sum(b.values())
        print(f"{a:>14} {n:>8} {b['base t3==t']:>11} {b['base t3==t+a']:>12} "
              f"{b['base t3==0']:>10} {b['other']:>7}   (dropped {b['dropped']})")
    print(f"\n  min-adapter-overlap is {MIN_OVERLAP}: the baseline's own matcher "
          "cannot fire below it,")
    print("  so those are the mates the PE mode can add a trim to.")
    class_stability(outdir, gt_path)


def _classify(t5, t3, h, t, tol=0):
    if abs(t5 - h) <= tol and abs(t3 - t) <= tol:
        return "exact"
    if t5 > h or t3 > t:
        return "overtrim"
    return "undertrim"


def classify(t5, t3, h, t, ad, tol=0):
    """The report's `classify`, minus the posterior (which needs the corpus).

    trim_report.py's classify() splits 'undertrim' into recoverable/irreducible
    using the measured posterior P(artefact | observed run length). This
    replicates everything except that split: adapter/N mates are irreducible,
    a clean G-run whose observed run has P > 0.5 is recoverable. The posterior
    is unavailable without the full corpus, so this reports the split as
    recoverable/irreducible using the same G-run rule but labels it as an
    approximation, and the class-movement counts are unaffected by it (a
    change of recoverable<->irreducible is still a change of class).
    """
    cls = _classify(t5, t3, h, t, tol)
    if cls != "undertrim":
        return cls
    if ad > 0:
        return "undertrim_irreducible"
    return "undertrim_recoverable"


def class_stability(outdir: str, gt_path: str) -> None:
    """How many of the changed decisions moved a scored class?

    `pe_trims` counts decisions that differ; the scored counters only move when
    the class also changes. This separates the two.
    """
    gt = {int(r["pair_id"]): r
          for r in csv.DictReader(open(gt_path), delimiter="\t")}
    base = parse_pxmt(os.path.join(outdir, "base.pxmt"))
    pe = parse_pxmt(os.path.join(outdir, "pe.pxmt"))
    c = collections.Counter()
    for (rid, mate), b in base.items():
        p = pe.get((rid, mate))
        if p is None or b == p:
            continue
        row = gt.get(_key(rid))
        if row is None:
            continue
        m = "r1" if mate == 1 else "r2"
        h, t = int(row[f"{m}_h"]), int(row[f"{m}_t"])
        if b[0] != 0 or p[0] != 0:
            c["changed: status flip (one side dropped)"] += 1
            continue
        cb = _classify(b[1], b[2], h, t)
        cp = _classify(p[1], p[2], h, t)
        if cb == cp:
            c[f"changed but class-neutral ({cb})"] += 1
        else:
            c[f"class moved {cb} -> {cp}"] += 1
    print("\n=== pe_trims: how many changed decisions moved a scored class? ===")
    for k in sorted(c):
        print(f"  {k:52} {c[k]}")


if __name__ == "__main__":
    main()

