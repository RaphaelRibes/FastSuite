#!/usr/bin/env python3
"""Aggregate the trim-accuracy benchmark: overtrim / undertrim / exact.

SUPERSEDED — do not cite the tables this prints as current. This reads
bench/results/phase-d/trim-accuracy.tsv (606 rows), the original 3-tool run
whose fastprep baseline (fastprep-winner, 66.5%) is under open investigation
and does not reproduce (75.59% fresh). The current trim-accuracy measurement
is bench/results/phase-d/overtrim-vs-fastp.tsv (all four tools under one
harness), produced by bench/scripts/overtrim_vs_fastp.py and documented in
docs/TRIM-ACCURACY.md. This script is kept only to reproduce the historical
table.

Reads bench/results/phase-d/trim-accuracy.tsv (or $1) and prints the
tolerance-0 / tolerance-2 tables for every tool, plus per-domain exact
rates. Deterministic: same input, same output.
"""
import csv
import sys

DEFAULT = "bench/results/phase-d/trim-accuracy.tsv"
TOOLS = ("fastprep-winner", "fastp-default", "fastp-tuned")


def load(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def totals(rows):
    keys = ("exact", "overtrim", "undertrim", "unaligned", "dropped_or_missing")
    s = {k: 0 for k in keys}
    for r in rows:
        for k in keys:
            s[k] += int(r[k])
    denom = s["exact"] + s["overtrim"] + s["undertrim"] + s["unaligned"]
    return s, denom


def pct(n, d):
    return f"{n / d * 100:5.1f}%" if d else "  n/a"


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    rows = load(path)
    print(f"rows={len(rows)}  datasets={len(set(r['run_accession'] for r in rows))}")

    for tol in (0, 2):
        print()
        print(f"=== tolerance {tol} ({'strict' if tol == 0 else 'forgiving'})")
        print(f"{'tool':18s} {'mates':>9s} {'exact':>8s} {'overtrim':>9s} "
              f"{'undertrim':>10s} {'unaligned':>10s} {'dropped':>8s}")
        for tool in TOOLS:
            sub = [r for r in rows if r["tool"] == tool and int(r["tol"]) == tol]
            if not sub:
                continue
            s, denom = totals(sub)
            print(f"{tool:18s} {denom:9d} {s['exact']:8d} {s['overtrim']:9d} "
                  f"{s['undertrim']:10d} {s['unaligned']:10d} "
                  f"{s['dropped_or_missing']:8d}")
            print(f"{'':18s} {'100.0%':>9s} {pct(s['exact'], denom):>8s} "
                  f"{pct(s['overtrim'], denom):>9s} "
                  f"{pct(s['undertrim'], denom):>10s} "
                  f"{pct(s['unaligned'], denom):>10s}")

    print()
    print("=== per-domain exact rate (tol 0 / tol 2)")
    domains = sorted(set(r["domain"] for r in rows))
    for tool in TOOLS:
        print(f"  {tool}")
        for dom in domains:
            line = f"    {dom:20s}"
            for tol in (0, 2):
                sub = [r for r in rows if r["tool"] == tool
                       and r["domain"] == dom and int(r["tol"]) == tol]
                s, denom = totals(sub)
                line += f"  tol{tol}={pct(s['exact'], denom):>7s}"
            print(line)


if __name__ == "__main__":
    main()
