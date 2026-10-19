#!/usr/bin/env python3
"""Cell B probe: what did each fastprep config actually decide, per mate?

Reads the PXMT files the trim report left behind (it deletes its scratch, so
run it once with the scratch dir kept) and reports, per (config, mate), the
counts by decision: emit / drop-with-trim / drop-without-trim, plus how often
the drop carried a non-zero trim. That last number is the mechanism question
behind the 'dropped' bucket: a read the tool drops for having too little left
after a trim is not the same event as a read it dropped for quality, and the
trim report's 'dropped' column cannot tell them apart.

Usage: probe_pxmt.py <pxmt> [<pxmt> ...]
"""
from __future__ import annotations

import collections
import struct
import sys


def parse_pxmt(path: str) -> dict[tuple[str, int], tuple[int, int, int]]:
    """-> {(read_id, mate): (status, t5, t3)}. Mirrors the harness's parser."""
    buf = open(path, "rb").read()
    if buf[:4] != b"PXMT":
        return {}
    out = {}
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


def main() -> None:
    for path in sys.argv[1:]:
        p = parse_pxmt(path)
        c = collections.Counter()
        for (rid, mate), (status, t5, t3) in p.items():
            if status != 0:
                c[f"mate{mate} dropped (t5={t5},t3={t3})"] += 1
                if t5 or t3:
                    c[f"mate{mate} dropped carrying a trim"] += 1
            else:
                c[f"mate{mate} emitted"] += 1
        print(f"--- {path} ({len(p)} records) ---")
        for k in sorted(c):
            print(f"  {k:44} {c[k]}")


if __name__ == "__main__":
    main()
