#!/usr/bin/env python3
"""Cell B: is the in_window=1 slice EMPTY of readthrough pairs a bug?

The trim report's `in_window` slice (L <= insert <= 2L-12) shows `rt=1` with
`in_window=1` as 0 rows, and `readthrough pairs in_window` as 0. This script
shows why, from the GT alone, and whether it is structural.

    readthrough on a mate  <=>  adapter_len = L - h - t - min(I, L-h-t) > 0
                           <=>  I < L - t
    in_window               <=>  I >= L

    =>  in_window AND readthrough is EMPTY for every t >= 0. QED.

So on a cohort whose readthrough comes only from short inserts (the mechanism
fastsim actually models), the two slices cannot overlap. It is not a bug in the
cohort and not a bug in the report: it is the definitions.

The PE mode's real operating region is therefore the BELOW-window readthrough
mates (I < L), and that is where its effect must be read. This script reports
the counts so the write-up can state the limit instead of leaving a suspicious
empty cell.

Usage:
    python3 window_readthrough_check.py <gt.tsv>
"""
from __future__ import annotations

import collections
import csv
import sys

L = 150
CEILING = 2 * L - 12


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
    c = collections.Counter()
    for r in rows:
        i = int(r["insert"])
        win = L <= i <= CEILING
        for m in ("r1", "r2"):
            a = int(r[f"{m}_adapter"])
            t = int(r[f"{m}_t"])
            c["mates"] += 1
            if a > 0:
                c["readthrough mates"] += 1
                c["readthrough & in_window" if win else "readthrough & below_window"] += 1
                # the mechanism, checked directly:
                c["  matches I < L - t" if i < L - t else "  CONTRADICTION"] += 1
            else:
                c["non-readthrough mates"] += 1
        if win:
            c["pairs in_window"] += 1
        else:
            c["pairs below_window"] += 1
    print("=== Cell B: readthrough x window ===")
    for k in ("mates", "pairs in_window", "pairs below_window",
              "readthrough mates", "readthrough & in_window",
              "readthrough & below_window", "  matches I < L - t",
              "  CONTRADICTION", "non-readthrough mates"):
        print(f"  {k:34} {c[k]}")
    print()
    print("  readthrough & in_window == 0 is STRUCTURAL: readthrough needs "
          "I < L - t, in_window needs I >= L.")
    print("  The two slices are disjoint by definition, so the PE mode's effect "
          "is read on the below-window mates.")


if __name__ == "__main__":
    main()
