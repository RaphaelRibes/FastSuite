#!/usr/bin/env python3
"""Targeted comparison of the configs that decide fastprep's 5' overtrim.

The Phase D bench scored the sweep winner at strict tolerance 0 by reading its
PXMT metadata. This script re-scores a short list of candidate configs on the
same synthetic profiles the sensitivity analysis uses, with the head-dropout
value pinned EXPLICITLY on every run: fastprep's CLI default is 4, and pinning
it matters because (a) the default is a real decision, not a zero, and (b)
fastprep refuses head-dropout >= min-run, so a silent default changes the grid
underneath you when min_run drops.

Candidate list (all scored at tol 0 and tol 2, plus the sweep's own score):

  base                   the shipped default (min_run 12, mercy 3, min_len 0.9)
  base+hd=1              the shipped default with head-dropout pinned to 1
  hd=2/3/4/6             the head-dropout axis alone
  min_run 9..16          the axis the Phase D winners actually sit on
  min_len 0.95..1.0      the length floor
  best-combo              min_run 14 + min_len 0.95 + head-dropout 3
  mercy low/high         the gap-mass axis
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys

FASTPREP = "/home/raphael/hermes_scratch/FastSuite/tools/fastprep/target/release/fastprep"
WORK = "/home/raphael/.hermes/cache/scratch/overtrim"

BASE = dict(min_run=12, mercy=3, phred_sat=20, gate_shape=1.0, min_len=0.9,
            mean_q=25.0, adapter_err=0.2, min_overlap=14)


def build(cfg: dict) -> list[str]:
    ml = cfg["min_len"]
    ml = str(int(ml)) if float(ml).is_integer() else f"{ml}L"
    a = ["--min-run", str(cfg["min_run"]), "--mercy", str(cfg["mercy"]),
         "--phred-sat", str(cfg["phred_sat"]), "--gate-shape", str(cfg["gate_shape"]),
         "--min-len", ml, "--mean-q", str(cfg["mean_q"]),
         "--adapter-error-rate", str(cfg["adapter_err"]),
         "--min-adapter-overlap", str(cfg["min_overlap"]), "--gz-level", "6"]
    if cfg.get("head_dropout") is not None:
        a += ["--head-dropout", str(cfg["head_dropout"])]
    else:
        a += ["--head-dropout", "4"]  # the shipped default, stated explicitly
    if cfg.get("r2_reverse"):
        a.append("--r2-scan-reverse")
    if cfg.get("interior"):
        a.append("--interior")
    return a


CANDIDATES = [
    ("base(default hd=4)", dict()),
    ("hd=1", dict(head_dropout=1)),
    ("hd=2", dict(head_dropout=2)),
    ("hd=3", dict(head_dropout=3)),
    ("hd=6 min_run=8", dict(head_dropout=6, min_run=8)),
    ("mr=9", dict(min_run=9)),
    ("mr=10", dict(min_run=10)),
    ("mr=12", dict(min_run=12)),
    ("mr=14", dict(min_run=14)),
    ("mr=16", dict(min_run=16)),
    ("ml=0.95", dict(min_len=0.95)),
    ("ml=1.0", dict(min_len=1.0)),
    ("mercy=1", dict(mercy=1)),
    ("mercy=2", dict(mercy=2)),
    ("mercy=8", dict(mercy=8)),
    ("mr=14 + ml=0.95", dict(min_run=14, min_len=0.95)),
    ("mr=16 + ml=0.95", dict(min_run=16, min_len=0.95)),
    ("mr=14 + hd=3", dict(min_run=14, head_dropout=3)),
    ("mr=16 + hd=3 + ml=0.95", dict(min_run=16, head_dropout=3, min_len=0.95)),
    ("mr=16 + hd=3 + ml=0.95 + m=2", dict(min_run=16, head_dropout=3, min_len=0.95, mercy=2)),
]


def parse_pxmt(path: str):
    import struct
    out = {}
    buf = open(path, "rb").read()
    off, n = 5, len(buf)
    while off < n:
        (id_len,) = struct.unpack_from("<H", buf, off); off += 2
        rid = buf[off:off + id_len]; off += id_len
        mate = buf[off]; off += 1
        status = buf[off]; off += 1
        t5a, t3a, t5p, t3p, _rl = struct.unpack_from("<IIIII", buf, off); off += 20
        nr = buf[off]; off += 1
        off += 8 * nr
        parts = rid.split(b" ", 1)[0].decode("ascii").lstrip("@").split(":")
        pid = int(parts[-2]) + 10_000 * int(parts[-1])
        out[(pid, mate)] = (status, t5a + t5p, t3a + t3p)
    return out


def parse_gt(path: str):
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[int(r["pair_id"])] = (int(r["r1_h"]), int(r["r1_t"]),
                                      int(r["r2_h"]), int(r["r2_t"]))
    return out


def classify(t5, t3, h, t, tol):
    if abs(t5 - h) <= tol and abs(t3 - t) <= tol:
        return "exact"
    if t5 > h or t3 > t:
        return "overtrim"
    return "undertrim"


def score(px, gt, tol):
    c = {"exact": 0, "overtrim": 0, "undertrim": 0}
    e5 = e3 = 0.0
    for pid, (h1, t1, h2, t2) in gt.items():
        for mate, (h, t) in ((0, (h1, t1)), (1, (h2, t2))):
            s = px.get((pid, mate))
            if s is None:
                continue
            st, t5, t3 = s
            if st != 0:
                continue
            c[classify(t5, t3, h, t, tol)] += 1
            e5 += t5 - h; e3 += t3 - t
    n = max(1, sum(c.values()))
    return dict(exact=c["exact"] / n, over=c["overtrim"] / n,
                under=c["undertrim"] / n, e5=e5 / n, e3=e3 / n, n=n)


def main() -> int:
    out = {}
    for regime in ("low", "mid", "high"):
        work = os.path.join(WORK, regime)
        gt = parse_gt(os.path.join(work, "sim.v0.gt.tsv"))
        rows = []
        for name, over in CANDIDATES:
            cfg = dict(BASE); cfg.update(over)
            tag = name.replace(" ", "_").replace("=", "").replace("+", "_")
            pxf = os.path.join(work, f"cand_{tag}.pxmt")
            argv = [FASTPREP, "-1", os.path.join(work, "sim.v0_R1.fastq.gz"),
                    "-2", os.path.join(work, "sim.v0_R2.fastq.gz"),
                    "-o", os.path.join(work, "z1.fastq.gz"),
                    "-p", os.path.join(work, "z2.fastq.gz"),
                    "--metadata", pxf, "--json", os.path.join(work, "cand.json"),
                    "--force", "-q", "-j", "4"] + build(cfg)
            p = subprocess.run(argv, capture_output=True, text=True)
            if p.returncode != 0:
                print(f"{regime}/{name}: rc={p.returncode} {p.stderr[-160:]}",
                      file=sys.stderr)
                continue
            r0 = score(parse_pxmt(pxf), gt, 0)
            r2 = score(parse_pxmt(pxf), gt, 2)
            rows.append((name, r0, r2))
        out[regime] = rows
        print(f"===== {regime} =====")
        hdr = f"{'config':30s} {'exact0':>7s} {'over0':>7s} {'und0':>6s} | {'exact2':>7s} {'over2':>6s} | {'e5':>6s} {'e3':>6s}"
        print(hdr)
        for name, r0, r2 in rows:
            print(f"{name:30s} {r0['exact']:7.4f} {r0['over']:7.4f} {r0['under']:6.4f} | "
                  f"{r2['exact']:7.4f} {r2['over']:6.4f} | {r0['e5']:+6.3f} {r0['e3']:+6.3f}")
    with open(os.path.join(WORK, "candidates.json"), "w") as fh:
        json.dump({k: [(n, a, b) for n, a, b in v] for k, v in out.items()}, fh, indent=1)
    print("wrote candidates.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
