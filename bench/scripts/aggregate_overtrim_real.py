#!/usr/bin/env python3
"""Aggregate the real-corpus head_dropout / boundary-anchored comparison."""
import csv, collections, sys

import sys
R = sys.argv[1] if len(sys.argv) > 1 else "/mnt/Fast/fastsuite-scratch/results/overtrim-real2.tsv"
rows = list(csv.DictReader(open(R), delimiter="\t"))
print(f"rows={len(rows)}")

ds = {(r["domain"], r["acc"]) for r in rows}
print(f"datasets with data: {len(ds)}")

VARIANTS = ["baseline", "hd1", "hd0", "ba"]

for tol in (0, 2):
    print(f"\n=== tolerance {tol}")
    print(f"{'variant':10s} {'mates':>12s} {'exact':>8s} {'overtrim':>9s} "
          f"{'undertrim':>10s} {'dropped':>9s}")
    for v in VARIANTS:
        sub = [r for r in rows if r["variant"] == v and int(r["tol"]) == tol]
        if not sub:
            continue
        s = lambda k: sum(int(r[k]) for r in sub)
        tot = s("exact") + s("overtrim") + s("undertrim")
        print(f"{v:10s} {tot:12d} {s('exact')/tot*100:7.2f}% "
              f"{s('overtrim')/tot*100:8.2f}% {s('undertrim')/tot*100:9.2f}% "
              f"{s('dropped'):9d}")

print("\n=== per-variant, mates scored (tol 0), for reference")
for v in VARIANTS:
    sub = [r for r in rows if r["variant"] == v and int(r["tol"]) == 0]
    print(f"  {v:10s} n={len(sub)} datasets")

# which dataset is missing?
allw = list(csv.DictReader(
    open("/mnt/Fast/fastsuite-scratch/results/winners.tsv"), delimiter="\t"))
want = {(w["domain"], w["run_accession"]) for w in allw}
missing = want - ds
print(f"\nmissing datasets ({len(missing)}): {sorted(missing)}")
