#!/usr/bin/env python3
"""Sensitivity of fastprep's strict-tolerance overtrim to its knobs.

For each dataset (a fastsim profile simulated with ground truth) and each
config, run fastprep with --metadata, parse PXMT, and classify every mate
against the GT's exact artefact window:

    exact     |t5-h| <= tol and |t3-t| <= tol
    overtrim  t5 > h or t3 > t
    undertrim t5 < h or t3 < t

t5/t3 are the adapter+polyX totals recorded per read, which is what the
"sweep.py in-depth" scoring does; the classification is the trim_accuracy.py
one. Both tolerances are reported.

cfg is the same dict bench/scripts/sweep.py builds, extended with the knobs
that GRID does not reach (match_leniency, mismatch_strength, min-len
fractions below 0.8, head_dropout 1..8, min_run 3..10, mercy 0.1..2.0).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import struct
import subprocess
import sys
import time

FASTPREP = "/home/raphael/hermes_scratch/FastSuite/tools/fastprep/target/release/fastprep"
MAGIC = b"PXMT"
TOLS = (0, 2)


# ---------------------------------------------------------------- PXMT parse
def parse_pxmt(path: str) -> dict[tuple[int, int], tuple[int, int, int]]:
    """-> {(pair_id, mate): (status, t5_total, t3_total)}; status 0 = kept.

    pair_id is recovered from the Illumina coordinate header exactly as
    sweep.py does: x = pid % 10000, y = pid // 10000.
    """
    out: dict[tuple[int, int], tuple[int, int, int]] = {}
    buf = open(path, "rb").read()
    if buf[:4] != MAGIC:
        raise ValueError(f"{path}: not PXMT")
    off, n = 5, len(buf)
    while off < n:
        (id_len,) = struct.unpack_from("<H", buf, off); off += 2
        rid = buf[off:off + id_len]; off += id_len
        mate = buf[off]; off += 1
        status = buf[off]; off += 1
        t5a, t3a, t5p, t3p, _rl = struct.unpack_from("<IIIII", buf, off); off += 20
        n_runs = buf[off]; off += 1
        off += n_runs * 8
        out[coord_to_pid(rid, mate)] = (status, t5a + t5p, t3a + t3p)
    return out


def coord_to_pid(rid: bytes, mate: int) -> tuple[int, int]:
    try:
        parts = rid.split(b" ", 1)[0].decode("ascii").lstrip("@").split(":")
        x, y = int(parts[-2]), int(parts[-1])
        return (y * 10_000 + x, mate)
    except (ValueError, UnicodeDecodeError, IndexError):
        return (-1, mate)


def parse_gt(path: str) -> dict[int, tuple[int, int, int, int]]:
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


def score(px: dict, gt: dict, tol: int, read_len: int = 150) -> dict:
    """Per-mate classification, and the mean signed error on each boundary.

    err5 = t5 - h and err3 = t3 - t are the signed boundary errors; their
    means say WHERE the overtrim sits, which is the part a knob can move.
    """
    c = {"exact": 0, "overtrim": 0, "undertrim": 0, "dropped": 0}
    e5 = e3 = 0.0
    n5 = n3 = 0
    over5 = over3 = 0
    d5 = d3 = 0
    for pid, (h1, t1, h2, t2) in gt.items():
        for mate, (h, t) in ((0, (h1, t1)), (1, (h2, t2))):
            s = px.get((pid, mate))
            if s is None:
                c["dropped"] += 1
                continue
            status, t5, t3 = s
            if status != 0:
                c["dropped"] += 1
                continue
            c[classify(t5, t3, h, t, tol)] += 1
            e5 += t5 - h; n5 += 1
            e3 += t3 - t; n3 += 1
            if t5 > h:
                over5 += 1; d5 += t5 - h
            if t3 > t:
                over3 += 1; d3 += t3 - t
    n = max(1, c["exact"] + c["overtrim"] + c["undertrim"])
    return dict(exact=c["exact"] / n, over=c["overtrim"] / n, under=c["undertrim"] / n,
                dropped=c["dropped"], scored=n,
                mean_err5=e5 / max(1, n5), mean_err3=e3 / max(1, n3),
                over5_frac=over5 / n, over3_frac=over3 / n,
                over5_bp=d5 / max(1, n5), over3_bp=d3 / max(1, n3))


# ---------------------------------------------------------------- the grid
BASE = dict(min_run=12, mercy=3, phred_sat=20, gate_shape=1.0, min_len=0.9,
            mean_q=25.0, adapter_err=0.2, min_overlap=14,
            r2_reverse=False, interior=False, head_dropout=None,
            match_leniency=0.1, mismatch_strength=0.5)


def argv(cfg: dict, work: str, cid: int) -> list[str]:
    ml = cfg["min_len"]
    ml = str(int(ml)) if float(ml).is_integer() else f"{ml}L"
    a = ["-1", f"{work}/sim.v0_R1.fastq.gz", "-2", f"{work}/sim.v0_R2.fastq.gz",
         "-o", f"{work}/c{cid}_R1.fastq.gz", "-p", f"{work}/c{cid}_R2.fastq.gz",
         "--metadata", f"{work}/c{cid}.pxmt", "--json", f"{work}/c{cid}.json",
         "--force", "-q", "-j", "4",
         "--min-run", str(cfg["min_run"]), "--mercy", str(cfg["mercy"]),
         "--phred-sat", str(cfg["phred_sat"]), "--gate-shape", str(cfg["gate_shape"]),
         "--min-len", ml, "--mean-q", str(cfg["mean_q"]),
         "--adapter-error-rate", str(cfg["adapter_err"]),
         "--min-adapter-overlap", str(cfg["min_overlap"]), "--gz-level", "6",
         "--match-leniency", str(cfg["match_leniency"]),
         "--mismatch-strength", str(cfg["mismatch_strength"])]
    if cfg.get("r2_reverse"):
        a.append("--r2-scan-reverse")
    if cfg.get("interior"):
        a.append("--interior")
    if cfg.get("head_dropout") is None:
        # fastprep refuses to combine a head-dropout >= min_run, which the old
        # CLI default of 4 hit on --min-run < 5 (the default is now 0, commit
        # 3d76a31). Always pin the value explicitly so every config in the grid
        # is runnable and the value is a visible choice, not a silent clash:
        # 1 means "only trim a 1 bp leading G-run", which for min_run <= 5 is
        # what the scan would own anyway.
        a += ["--head-dropout", "1"]
    else:
        n = min(cfg["head_dropout"], max(1, cfg["min_run"] - 1))
        a += ["--head-dropout", str(n)]
    return a


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--regimes", nargs="+", default=["low", "mid", "high"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--grid", default="base", help="base | extended | full")
    args = ap.parse_args()

    configs: list[tuple[str, dict]] = []

    def add(name, **over):
        c = dict(BASE)
        c.update(over)
        configs.append((name, c))

    add("base")
    # ---- axis 1: min_run
    for mr in (3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20):
        add(f"min_run={mr}", min_run=mr)
    # ---- axis 2: mercy
    for m in (0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 15.0):
        add(f"mercy={m}", mercy=m)
    # ---- axis 3: quality windows
    for q in (10, 15, 20, 25, 30, 35):
        add(f"phred_sat={q}", phred_sat=q)
    for x in (0.25, 0.5, 1.5, 2.0, 3.0):
        add(f"gate_shape={x}", gate_shape=x)
    # ---- axis 4: min-len floor
    for ml in (0.3, 0.5, 0.6, 0.7, 0.8, 0.85, 0.95, 1.0, 30, 50):
        add(f"min_len={ml}", min_len=ml)
    # ---- axis 5: mean_q
    for q in (15.0, 20.0, 25.0, 30.0):
        add(f"mean_q={q}", mean_q=q)
    # ---- axis 6: head dropout
    for n in (1, 2, 3, 4, 5, 6, 8):
        add(f"head_dropout={n}", head_dropout=n)
    # ---- axis 7: match leniency / mismatch strength
    for mu in (0.0, 0.05, 0.2, 0.4, 0.8):
        add(f"match_leniency={mu}", match_leniency=mu)
    for nu in (0.0, 0.25, 0.75, 1.0):
        add(f"mismatch_strength={nu}", mismatch_strength=nu)
    # ---- axis 8: scan structure
    add("interior", interior=True)
    add("r2_reverse", r2_reverse=True)

    fields = ["regime", "config"]
    rows = []
    for regime in args.regimes:
        work = f"/home/raphael/.hermes/cache/scratch/overtrim/{regime}"
        gt = parse_gt(f"{work}/sim.v0.gt.tsv")
        for name, cfg in configs:
            cid = len(rows)
            t0 = time.time()
            p = subprocess.run([FASTPREP] + argv(cfg, work, cid),
                               capture_output=True, text=True)
            secs = time.time() - t0
            if p.returncode != 0:
                print(f"{regime}/{name}: rc={p.returncode} {p.stderr[-200:]}",
                      file=sys.stderr)
                continue
            px = parse_pxmt(f"{work}/c{cid}.pxmt")
            r0 = score(px, gt, 0)
            r2 = score(px, gt, 2)
            rows.append(dict(regime=regime, config=name, secs=f"{secs:.3f}",
                             **{f"{k}0": f"{v:.5f}" for k, v in r0.items()
                                if isinstance(v, float)},
                             **{f"{k}2": f"{v:.5f}" for k, v in r2.items()
                                if isinstance(v, float)}, **cfg))
            for f in (f"c{cid}_R1.fastq.gz", f"c{cid}_R2.fastq.gz", f"c{cid}.json"):
                try:
                    os.unlink(f"{work}/{f}")
                except FileNotFoundError:
                    pass
        print(f"{regime}: {len(configs)} configs", file=sys.stderr)

    keys = list(rows[0].keys())
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.out} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
