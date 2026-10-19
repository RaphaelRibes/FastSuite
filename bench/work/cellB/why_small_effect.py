#!/usr/bin/env python3
"""Cell B: where could --pe-overlap possibly add value, and where is it blind?

Two mechanisms decide whether the mode can help a mate:

  (A) the 3' polyG scan. The scorer compares the tool's 3' trim against the
      GT's 3' G-run length `t`. A run shorter than `--min-run` (12) cannot be
      trimmed by the scan at all, so there is no run for a PE-derived boundary
      to interact with.
  (B) the baseline's own adapter matcher. It needs
      `--min-adapter-overlap` (14) bases of adapter at the read's 3' end before
      it fires. Below that it is blind, and the PE mode is the only thing that
      could trim there.

This script buckets the cohort's readthrough mates on (A) x (B) and reports the
count in each cell -- the size of the addressable set, and the size of the set
the mode is structurally unable to reach.

Usage:
    python3 why_small_effect.py <gt.tsv>
"""
from __future__ import annotations

import collections
import csv
import sys

L = 150
MIN_RUN = 12
MIN_ADAPTER_OVERLAP = 14


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
    c = collections.Counter()
    pairs = collections.Counter()
    for r in rows:
        for m in ("r1", "r2"):
            a = int(r[f"{m}_adapter"])
            if a == 0:
                continue
            t = int(r[f"{m}_t"])
            scan = "t>=min_run" if t >= MIN_RUN else "t<min_run"
            matcher = "a>=14" if a >= MIN_ADAPTER_OVERLAP else "a<14"
            c[f"{scan} / {matcher}"] += 1
            pairs[f"pairs {scan} / {matcher}"] += 1
    print("=== Cell B: readthrough mates by (polyG scan active?) x (baseline matcher can see?) ===")
    print(f"{'bucket':28} {'mates':>9} {'pairs':>9}")
    for k in sorted(c):
        print(f"  {k:26} {c[k]:>9} {pairs['pairs ' + k]:>9}")
    tot = sum(c.values())
    addable = c[f"t>=min_run / a<14"]
    blind = c[f"t<min_run / a<14"]
    already = c[f"t>=min_run / a>=14"]
    print()
    print(f"  total readthrough mates                    {tot}")
    print(f"  scan active + matcher blind  (addressable) {addable}")
    print(f"  scan active + matcher sees  (already)     {already}")
    print(f"  scan blind   + matcher blind (unreachable) {blind}")


if __name__ == "__main__":
    main()
