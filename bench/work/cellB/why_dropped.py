#!/usr/bin/env python3
"""Cell B: why are so many readthrough mates DROPPED?

The with/without report shows `dropped` as the second-largest bucket
(fastprep-hd0: 131273 of 800000 mates, 16.4%). A drop is not a trim decision the
scorer can score, so if the mode changes the drop rate it changes the kept-set
the exact rate is earned on. This script answers the mechanism question directly
from the GT + the PXMT:

  * what does the GT say the correct answer is for a dropped mate? A mate is
    `full_flag=1` (entirely artefact) when the pair is all artefact -- those
    SHOULD be dropped;
  * how long is the GT's keep window `[h, min(I, L-t))`? Below `--min-len`
    (0.9L = 135) a correct trimmer must drop it;
  * so: of the mates the baseline dropped, how many are a CORRECT drop (keep
    window < min_len, or full_flag) and how many are a drop the GT would not
    require?

Usage:
    python3 why_dropped.py <gt.tsv> <pxmt>
"""
from __future__ import annotations

import collections
import csv
import struct
import sys

L = 150
MIN_LEN = int(0.9 * L)   # --min-len 0.9L


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


def key(rid: str) -> int:
    f = rid.split(":")
    return int(f[5]) + 10000 * int(f[6])


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    gt_path, pxmt = sys.argv[1], sys.argv[2]
    gt = {int(r["pair_id"]): r for r in csv.DictReader(open(gt_path), delimiter="\t")}
    p = parse_pxmt(pxmt)
    c = collections.Counter()
    for (rid, mate), (status, t5, t3) in p.items():
        row = gt.get(key(rid))
        if row is None:
            continue
        m = "r1" if mate == 1 else "r2"
        h, t = int(row[f"{m}_h"]), int(row[f"{m}_t"])
        i = int(row["insert"])
        keep = max(0, min(i, L - t) - h)
        full = row["full_flag"] == "1"
        if status == 0:
            c["emitted"] += 1
            continue
        c["dropped"] += 1
        c[f"dropped mate{mate}"] += 1
        if full:
            c["dropped: GT says full (correct drop)"] += 1
        elif keep < MIN_LEN:
            c["dropped: keep window < min_len (correct drop)"] += 1
            c[f"  keep window {keep}"] += 1 if keep < 20 else 0
        else:
            c["dropped: GT says KEEP, tool dropped anyway"] += 1
            c[f"  dropped with keep window {keep}"] += 1 if keep < 145 else 0
    print("=== Cell B: the 'dropped' bucket, from the GT ===")
    for k in ("emitted", "dropped", "dropped mate1", "dropped mate2",
              "dropped: GT says full (correct drop)",
              "dropped: keep window < min_len (correct drop)",
              "dropped: GT says KEEP, tool dropped anyway"):
        print(f"  {k:52} {c[k]}")


if __name__ == "__main__":
    main()
