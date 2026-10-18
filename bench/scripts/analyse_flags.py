#!/usr/bin/env python3
"""Aggregate the flag-measurement decomposition TSV.

Reads the per-tool decomposition written by trim_report.py and prints, per
config and tolerance: exact / overtrim / undertrim-recoverable /
undertrim-irreducible / dropped / scored, the exact rate on the KEPT set, and
retention (kept / total) -- plus the deltas against the fastprep-hd0 baseline.

Usage: analyse_flags.py <decomp.tsv>
"""
import collections, csv, sys

path = sys.argv[1]
rows = list(csv.DictReader(open(path), delimiter="\t"))
per = collections.defaultdict(dict)          # (tool,tol) -> row
for r in rows:
    per[(r["tool"], int(r["tol"]))] = {k: int(r[k]) for k in
        ("exact", "overtrim", "undertrim_recoverable",
         "undertrim_irreducible", "dropped", "unaligned", "scored")}

ORDER = ["fastprep-hd0", "fp-3p6", "fp-n", "fp-3p6-n", "fastprep-hd4"]
LABEL = {"fastprep-hd0": "a) baseline hd0", "fp-3p6": "b) + --min-run-3p 6",
         "fp-n": "c) + --trim-terminal-n", "fp-3p6-n": "d) + both",
         "fastprep-hd4": "(control) hd4"}

for tol in (0, 2):
    print(f"\n===== tolerance {tol} "
          f"{'=' * 58}")
    print(f"{'config':22} {'exact':>12} {'overtrim':>10} {'und-rec':>10} "
          f"{'und-irr':>10} {'dropped':>10} {'scored':>12} {'ret%':>7} "
          f"{'exact/kept':>10}")
    base = per[("fastprep-hd0", tol)]
    for tag in ORDER:
        s = per.get((tag, tol))
        if not s:
            continue
        tot = s["scored"] + s["dropped"]
        ret = 100 * s["scored"] / tot
        rate = 100 * s["exact"] / s["scored"]
        print(f"{LABEL[tag]:22} {s['exact']:12d} {s['overtrim']:10d} "
              f"{s['undertrim_recoverable']:10d} "
              f"{s['undertrim_irreducible']:10d} {s['dropped']:10d} "
              f"{s['scored']:12d} {ret:7.2f} {rate:9.3f}%")
    print(f"\n  deltas vs baseline (absolute mate counts):")
    print(f"{'config':22} {'d-exact':>12} {'d-overtrim':>12} {'d-und-rec':>12} "
          f"{'d-und-irr':>12} {'d-dropped':>12} {'d-total_err':>12}")
    for tag in ORDER:
        s = per.get((tag, tol))
        if not s or tag == "fastprep-hd0":
            continue
        e = s["exact"] - base["exact"]
        o = s["overtrim"] - base["overtrim"]
        rc = s["undertrim_recoverable"] - base["undertrim_recoverable"]
        ir = s["undertrim_irreducible"] - base["undertrim_irreducible"]
        dr = s["dropped"] - base["dropped"]
        derr = (s["overtrim"] + s["undertrim_recoverable"]
                + s["undertrim_irreducible"]) - \
               (base["overtrim"] + base["undertrim_recoverable"]
                + base["undertrim_irreducible"])
        print(f"{LABEL[tag]:22} {e:+12d} {o:+12d} {rc:+12d} {ir:+12d} "
              f"{dr:+12d} {derr:+12d}")
        kept_pp = 100 * e / base["scored"]
        print(f"{'':22} exact-rate delta on the SAME kept-set basis: "
              f"{e / base['scored'] * 100:+.4f} pp")
