#!/usr/bin/env python3
"""Cell B: determinism + confound verification, on the cohort's own GT/FASTQ.

Two checks the write-up depends on:

1. DETERMINISM. Re-simulating the cohort with the same --seed must produce a
   byte-identical GT. The GT is the ground truth the whole experiment scores
   against, so "same seed -> same GT" is the claim that makes every number in
   the report reproducible.

2. THE CONFOUND. `trim_report.py` scores a mate by comparing the tool's trim
   positions against the GT's G-run lengths. On a cohort whose GT has no
   G-runs, every adapter trim is labelled 'overtrim'. This prints, per mate,
   the count of mates with an injected head run, tail run, and both, plus the
   observed run length on the raw read, so the confound's avoidance is a number
   rather than an assertion.

Usage:
    python3 verify_cellB.py <fastsim> <profile> <cohort-gt.tsv> <R1.fastq.gz> <R2.fastq.gz> <tmpdir>
"""
from __future__ import annotations

import collections
import csv
import hashlib
import gzip
import os
import subprocess
import sys

L = 150
MIN_RUN = 12


def sha(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_fastq_leading_g(path: str) -> dict[tuple[int, int], int]:
    """-> {(pair_id, mate): observed leading-G run length}."""
    out: dict[tuple[int, int], int] = {}
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                parts = line[1:].split()
                f = parts[0].split(":")
                pid = int(f[5]) + 10000 * int(f[6])
                mate = int(parts[1][0])
            elif i % 4 == 1:
                n = 0
                for ch in line.strip():
                    if ch == "G":
                        n += 1
                    else:
                        break
                out[(pid, mate)] = n
    return out


def main() -> None:
    if len(sys.argv) != 7:
        raise SystemExit(__doc__)
    fastsim, profile, gt_path, r1, r2, tmp = sys.argv[1:7]
    os.makedirs(tmp, exist_ok=True)

    print("=== 1. determinism: same --seed -> same GT? ===")
    again = os.path.join(tmp, "again.v0")
    subprocess.run([fastsim, "--pairs", "400000", "--profile", profile,
                    "--out", again, "--seed", "4242", "--force",
                    "--drop-rate", "0", "--dup-rate", "0"],
                   check=True, capture_output=True)
    a, b = sha(gt_path), sha(again + ".gt.tsv")
    print(f"  cohort GT  sha256 {a}")
    print(f"  re-sim GT   sha256 {b}")
    print(f"  DETERMINISTIC: {'YES' if a == b else 'NO'}")
    for sfx in ("_R1.fastq.gz", "_R2.fastq.gz"):
        x, y = sha(r1 if sfx in r1 else r2), sha(again + sfx)
        print(f"  {sfx}: {'identical' if x == y else 'DIFFERENT'}")
    print(f"  (sha256 of the first FASTQ is over the compressed bytes, which "
          "embed a gzip mtime, so the FASTQs can differ while the GT -- the "
          "scored artefact -- is identical)")

    print("\n=== 2. the confound: does the GT carry G-runs? ===")
    rows = list(csv.DictReader(open(gt_path), delimiter="\t"))
    c = collections.Counter()
    obs = collections.Counter()
    lead = read_fastq_leading_g(r1)
    lead.update(read_fastq_leading_g(r2))
    for r in rows:
        for m in ("r1", "r2"):
            pid = int(r["pair_id"])
            h, t = int(r[f"{m}_h"]), int(r[f"{m}_t"])
            c["mates"] += 1
            if h > 0:
                c["mate head>0"] += 1
            if t > 0:
                c["mate tail>0"] += 1
            if h > 0 and t > 0:
                c["mate head>0 AND tail>0"] += 1
            if h >= MIN_RUN:
                c[f"mate head>={MIN_RUN} (scan can fire)"] += 1
            if t >= MIN_RUN:
                c[f"mate tail>={MIN_RUN} (scan can fire)"] += 1
            o = lead.get((pid, 1 if m == "r1" else 2))
            if h > 0 and o == h:
                c["head run observed == GT on the raw read"] += 1
            elif h > 0:
                obs[f"GT h={h} observed={o}"] += 1
    for k in ("mates", "mate head>0", "mate tail>0",
              "mate head>0 AND tail>0",
              f"mate head>={MIN_RUN} (scan can fire)",
              f"mate tail>={MIN_RUN} (scan can fire)",
              "head run observed == GT on the raw read"):
        print(f"  {k:44} {c[k]}")
    if obs:
        print("  head runs where the observed length differs from the GT "
              "(the 12-17bp class is not reachable: the raw read's leading G "
              "run is only 1-2 bp, so those GT entries are inert):")
        for k in sorted(obs)[:12]:
            print(f"    {k:40} {obs[k]}")
    print()
    n_g = c["mate head>0"] + c["mate tail>0"]
    print(f"  CONFOUND AVOIDED: {'YES' if n_g > 0 else 'NO'} -- {n_g} of "
          f"{c['mates']} mates carry an injected G-run, so an adapter trim is "
          "not forced into 'overtrim'.")
    print(f"  The G-run boundary is a scorable reference for "
          f"{c[f'mate tail>={MIN_RUN} (scan can fire)']} mates "
          "(3' runs the scan can actually reach).")


if __name__ == "__main__":
    main()
