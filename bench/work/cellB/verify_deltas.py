#!/usr/bin/env python3
"""Cell B: independent per-mate counter, recomputed from PXMT + GT alone.

The report's `cellB-with-pe.tsv` vs `cellB-without-pe.tsv` deltas are the
headline. This script recomputes the SAME decomposition from the two PXMTs and
the cohort GT, with no reference to the report, and asserts it matches. It is
the verification that the deltas are not an artefact of the harness's own
accumulation.

It also answers the question the report cannot: which class did each changed
decision move, and does the net movement reconcile with the report's deltas?

Usage:
    python3 verify_deltas.py <gt.tsv> <base.pxmt> <pe.pxmt> <report-with.tsv> <report-without.tsv>
"""
from __future__ import annotations

import collections
import csv
import sys

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from explain_pe_effect import parse_pxmt, _key, _classify, classify  # noqa: E402

MIN_RUN = 12


def main() -> None:
    if len(sys.argv) != 6:
        raise SystemExit(__doc__)
    gt_path, base_px, pe_px, rep_with, rep_without = sys.argv[1:6]
    gt = {int(r["pair_id"]): r
          for r in csv.DictReader(open(gt_path), delimiter="\t")}
    base, pe = parse_pxmt(base_px), parse_pxmt(pe_px)

    # ---- independent decomposition, both arms, tolerance 0 ----
    def decomp(px):
        c = collections.Counter()
        for (rid, mate), (status, t5, t3) in px.items():
            row = gt.get(_key(rid))
            if row is None:
                continue
            m = "r1" if mate == 1 else "r2"
            h, t = int(row[f"{m}_h"]), int(row[f"{m}_t"])
            ad = int(row[f"{m}_adapter"])
            if status != 0:
                c["dropped"] += 1
                continue
            c[classify(t5, t3, h, t, ad)] += 1
        return c

    cb, cp = decomp(base), decomp(pe)
    print("=== independent decomposition (from PXMT + GT, tolerance 0) ===")
    print(f"{'class':24} {'baseline':>10} {'fp-pe':>10} {'delta':>8}")
    for k in ("exact", "overtrim", "undertrim_recoverable",
              "undertrim_irreducible", "dropped"):
        print(f"  {k:22} {cb[k]:>10} {cp[k]:>10} {cp[k] - cb[k]:>+8}")

    # ---- the report's own deltas ----
    def row(path, tool):
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r["tool"] == tool and int(r["tol"]) == 0:
                return r
        raise SystemExit(f"{tool} not in {path}")

    rw, rwo = row(rep_with, "fp-pe"), row(rep_without, "fastprep-hd0")
    print("\n=== report deltas (cellB-with-pe minus cellB-without-pe, fp-pe vs hd0) ===")
    ok = True
    for k in ("exact", "overtrim", "undertrim_recoverable",
              "undertrim_irreducible", "dropped"):
        got = int(rw[k]) - int(rwo[k])
        want = cp[k] - cb[k]
        flag = "OK" if got == want else "MISMATCH"
        if got != want:
            ok = False
        print(f"  {k:22} report {got:>+6}   recomputed {want:>+6}   {flag}")
    print(f"\n  RECONCILED: {'YES' if ok else 'NO'}")

    # ---- class movement, to name what moved ----
    print("\n=== which class did each changed decision move? ===")
    mv = collections.Counter()
    for (rid, mate), b in base.items():
        p = pe.get((rid, mate))
        if p is None or b == p:
            continue
        row = gt.get(_key(rid))
        if row is None:
            continue
        m = "r1" if mate == 1 else "r2"
        h, t = int(row[f"{m}_h"]), int(row[f"{m}_t"])
        ad = int(row[f"{m}_adapter"])
        if b[0] != 0:
            cb = "dropped"
        else:
            cb = classify(b[1], b[2], h, t, ad)
        if p[0] != 0:
            cq = "dropped"
        else:
            cq = classify(p[1], p[2], h, t, ad)
        mv[f"{cb} -> {cq}"] += 1
    for k in sorted(mv):
        print(f"  {k:46} {mv[k]}")


if __name__ == "__main__":
    main()
