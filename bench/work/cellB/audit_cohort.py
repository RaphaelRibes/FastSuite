#!/usr/bin/env python3
"""Cell B audit: does the cohort carry BOTH G-runs and adapter readthrough?

The crux of the whole experiment. `bench/scripts/trim_report.py` classifies a
mate by comparing the tool's trim positions against the GROUND TRUTH's G-run
lengths (`r1_h/r1_t/r2_h/r2_t`):

    abs(t5-h) <= tol and abs(t3-t) <= tol  ->  'exact'
    t5 > h or t3 > t                      ->  'overtrim'

An ADAPTER trim is not a G-run. So on a cohort whose GT has NO G-runs every
adapter trim is labelled 'overtrim' -- the false verdict that made Cell A
unreadable. This script states, from the GT itself, whether that confound is
avoided here:

  * the G-run side: how many mates carry an injected head and/or tail run, and
    of those, how many are LONG enough (>= min_run) for the quality-gated polyX
    scan to fire at all -- a run shorter than min_run cannot be trimmed by the
    scan, so it is not a reference the scorer can score against;
  * the readthrough side: how many pairs carry adapter on at least one mate,
    how many of those are inside the overlap window, and how many carry >= 14
    adapter bases (fastprep's `--min-adapter-overlap`, the floor for its own
    sequence matcher to fire);
  * the crossing: readthrough pairs that ALSO carry a G-run. Those are the pairs
    where the two questions are answered on the same mate.

It also replicates fastprep's PE-overlap detector in pure Python on the raw
reads (same anchor, same mismatch budgets, same uniqueness rule, same
"the readthrough must really be adapter" guard) and reports how often it fires
on this cohort -- i.e. the size of the cohort the mode can act on.

Usage:
    python3 bench/work/cellB/audit_cohort.py <gt.tsv> <R1.fastq[.gz]> <R2.fastq[.gz]>
"""
from __future__ import annotations

import collections
import csv
import gzip
import sys

L = 150
ANCHOR = 12
CEILING = 2 * L - ANCHOR
MIN_ADAPTER_OVERLAP = 14
ERROR_RATE = 0.2
# fastprep's built-in catalogue (adapter.rs): TruSeq R1, TruSeq R2, Nextera.
CATALOGUE = (
    b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    b"CTGTCTCTTATACACATCT",
)

_COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def revcomp(s: bytes) -> bytes:
    return s.translate(_COMP)[::-1]


def read_fastq(path: str) -> dict[tuple[int, int], bytes]:
    """-> {(pair_id, mate): seq}, keyed by the sim header's x + 10000*y."""
    out: dict[tuple[int, int], bytes] = {}
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                parts = line[1:].split()
                f = parts[0].split(":")
                pid = int(f[5]) + 10000 * int(f[6])
                mate = int(parts[1][0])
            elif i % 4 == 1:
                out[(pid, mate)] = line.strip().encode()
    return out


def match_adapters(seq: bytes, min_overlap: int = MIN_ADAPTER_OVERLAP):
    """fastprep's match_adapters, restricted to the case this cohort produces:
    an adapter that starts at p and runs to the read's 3' end. Returns the
    leftmost (p, k) that matches, else None."""
    best = None
    for ad in CATALOGUE:
        a = len(ad)
        max_p = len(seq) - min_overlap
        for p in range(0, max_p + 1):
            k = min(a, len(seq) - p)
            if k < min_overlap:
                break
            max_mm = int(k * ERROR_RATE + 0.999999)
            mm = 0
            for x, y in zip(seq[p:p + k], ad[:k]):
                if x != y:
                    mm += 1
                    if mm > max_mm:
                        break
            if mm <= max_mm:
                if best is None or p < best[0]:
                    best = (p, k)
                break
    return best


def detect_pe_overlap(r1: bytes, r2: bytes):
    """Replicates fastprep's pe_overlap::detect_pe_overlap_trim."""
    l1, l2 = len(r1), len(r2)
    if l1 < ANCHOR or l2 < ANCHOR:
        return None
    max_insert = min(l1, l2)
    if max_insert < ANCHOR:
        return None
    rc2 = revcomp(r2)
    admissible = []
    for o in range(ANCHOR, max_insert + 1):
        a = r1[o - ANCHOR:o]
        b = rc2[l2 - ANCHOR:]
        if sum(1 for x, y in zip(a, b) if x != y) > 1:
            continue
        mm = 0
        ok = True
        for j in range(o):
            if r1[j] != rc2[l2 - o + j]:
                mm += 1
                if mm > 3:
                    ok = False
                    break
        if ok:
            admissible.append(o)
    admissible = sorted(set(admissible))
    if len(admissible) != 1:
        return None
    o = admissible[0]
    if o >= l1:
        return None
    readthrough = l1 - o
    # readthrough_is_adapter: the tail must be adapter starting at the boundary
    tail = r1[l1 - readthrough:]
    hit = False
    for ad in CATALOGUE:
        a = len(ad)
        if a < MIN_ADAPTER_OVERLAP:
            continue
        k = min(a, len(tail))
        p = len(tail) - k
        if p != 0 or k != readthrough or k < MIN_ADAPTER_OVERLAP:
            continue
        max_mm = int(k * ERROR_RATE + 0.999999)
        mm = sum(1 for x, y in zip(tail, ad[:k]) if x != y)
        if mm <= max_mm:
            hit = True
            break
    return readthrough if hit else None


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    gt_path, r1_path, r2_path = sys.argv[1:4]
    rows = list(csv.DictReader(open(gt_path), delimiter="\t"))
    reads = read_fastq(r1_path)
    reads.update(read_fastq(r2_path))

    c: collections.Counter = collections.Counter()
    for r in rows:
        pid = int(r["pair_id"])
        i = int(r["insert"])
        win = L <= i <= CEILING
        rt = int(r["r1_adapter"]) > 0 or int(r["r2_adapter"]) > 0
        mates = (("r1", 1), ("r2", 2))
        any_g = False
        for key, mate in mates:
            h, t = int(r[f"{key}_h"]), int(r[f"{key}_t"])
            a = int(r[f"{key}_adapter"])
            c["mates"] += 1
            if h > 0:
                c["mate head>0"] += 1
            if t > 0:
                c["mate tail>0"] += 1
            if t > 0 and t >= 12:
                c["mate tail>=min_run (12)"] += 1
            if h > 0 and h >= 12:
                c["mate head>=min_run (12)"] += 1
            if a > 0:
                c["mate adapter>0"] += 1
            if a >= MIN_ADAPTER_OVERLAP:
                c[f"mate adapter>={MIN_ADAPTER_OVERLAP}"] += 1
            if h > 0 or t > 0:
                any_g = True
        c["pairs"] += 1
        if win:
            c["pairs in_window"] += 1
        if rt:
            c["pairs readthrough"] += 1
            if win:
                c["readthrough pairs in_window"] += 1
            if any_g:
                c["readthrough pairs ALSO carrying a G-run"] += 1
            r1, r2 = reads.get((pid, 1)), reads.get((pid, 2))
            if r1 and r2:
                d = detect_pe_overlap(r1, r2)
                if d is None:
                    c["readthrough pairs: PE detector refused"] += 1
                else:
                    c["readthrough pairs: PE detector FIRED"] += 1
                    if d == max(int(r["r1_adapter"]), int(r["r2_adapter"])):
                        c["  fired with readthrough == GT"] += 1
                    else:
                        c["  fired with readthrough != GT"] += 1

    print("=== Cell B cohort audit ===")
    for k in ("pairs", "pairs in_window", "pairs readthrough",
              "readthrough pairs in_window",
              "readthrough pairs ALSO carrying a G-run",
              "readthrough pairs: PE detector FIRED",
              "  fired with readthrough == GT",
              "  fired with readthrough != GT",
              "readthrough pairs: PE detector refused",
              "mates", "mate head>0", "mate tail>0",
              "mate head>=min_run (12)", "mate tail>=min_run (12)",
              "mate adapter>0", f"mate adapter>={MIN_ADAPTER_OVERLAP}"):
        print(f"  {k:52} {c[k]}")
    print()
    print("CONFOUND CHECK: the GT carries G-runs "
          f"({'YES' if c['mate head>0'] + c['mate tail>0'] > 0 else 'NO'}) -- "
          f"{c['mate head>0'] + c['mate tail>0']} of {c['mates']} mates carry an "
          "injected G-run, so an adapter trim is no longer forced into "
          "'overtrim'.")
    print(f"VALUE CHECK: {c['pairs readthrough']} of {c['pairs']} pairs carry "
          "adapter readthrough, so the PE mode has a boundary to find.")


if __name__ == "__main__":
    main()
