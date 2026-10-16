#!/usr/bin/env python3
"""Diagnose WHERE fastprep's strict-tolerance overtrim lives on synthetic data.

fastsim forces a non-G base adjacent to every injected run (artefact.rs
"Boundary rule"), so the run-length ambiguity the probe (/tmp/probe_overtrim.py)
measures is REMOVED here. Whatever overtrim remains is a fastprep algorithm
effect, not boundary ambiguity. This script classifies every kept mate against
the GT and, for the overtrims, reports the geometry: which end, the error
magnitude, the GT artefact length, and the actual leading/trailing G-run.

Reuses the scoring convention of bench/scripts/overtrim_sensitivity.py.
"""
from __future__ import annotations
import csv, gzip, struct, subprocess, sys, os, collections

FASTPREP = "/home/raphael/hermes_scratch/FastSuite/tools/fastprep/target/release/fastprep"
WORKROOT = "/home/raphael/.hermes/cache/scratch/overtrim"

# ---- PXMT parse (overtrim_sensitivity convention) ----
def parse_pxmt(path):
    out = {}
    buf = open(path, "rb").read()
    if buf[:4] != b"PXMT":
        return out
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
        pid = int(parts[-1]) * 10_000 + int(parts[-2])
        out[(pid, mate)] = (status, t5a + t5p, t3a + t3p)
    return out

def parse_gt(path):
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[int(r["pair_id"])] = (int(r["r1_h"]), int(r["r1_t"]),
                                      int(r["r2_h"]), int(r["r2_t"]))
    return out

def load_seqs(path):
    seqs = {}
    with gzip.open(path, "rt") as fh:
        i = 0; k = None
        for l in fh:
            if i % 4 == 0:
                c = l[1:].split()[0]
                parts = c.split(":")
                pid = int(parts[-1]) * 10_000 + int(parts[-2])
                k = pid
            elif i % 4 == 1:
                seqs[k] = l.rstrip()
            i += 1
    return seqs

def argv(work, out):
    return [FASTPREP, "-1", f"{work}/sim.v0_R1.fastq.gz", "-2", f"{work}/sim.v0_R2.fastq.gz",
            "-o", f"{work}/_d_R1.fastq.gz", "-p", f"{work}/_d_R2.fastq.gz",
            "--metadata", out, "--json", f"{work}/_d.json", "--force", "-q", "-j", "4",
            "--min-run", "12", "--mercy", "3", "--phred-sat", "20", "--gate-shape", "1.0",
            "--min-len", "0.9L", "--mean-q", "25.0", "--adapter-error-rate", "0.2",
            "--min-adapter-overlap", "14", "--gz-level", "6",
            "--match-leniency", "0.1", "--mismatch-strength", "0.5",
            "--head-dropout", "1"]

def main():
    regimes = sys.argv[1:] or ["low", "mid", "high"]
    for regime in regimes:
        work = f"{WORKROOT}/{regime}"
        gt = parse_gt(f"{work}/sim.v0.gt.tsv")
        pxout = f"{work}/_diag.pxmt"
        p = subprocess.run(argv(work, pxout), capture_output=True, text=True)
        if p.returncode != 0:
            print(f"[{regime}] rc={p.returncode} {p.stderr[-200:]}")
            continue
        px = parse_pxmt(pxout)
        s1 = load_seqs(f"{work}/sim.v0_R1.fastq.gz")
        s2 = load_seqs(f"{work}/sim.v0_R2.fastq.gz")
        seqs = {**{(k, 0): v for k, v in s1.items()}, **{(k, 1): v for k, v in s2.items()}}

        cnt = {0: collections.Counter(), 2: collections.Counter()}
        over5_mag = collections.Counter(); over3_mag = collections.Counter()
        over5_ctx = collections.Counter()  # (h, leading-run-len) for 5' overtrims
        over3_ctx = collections.Counter()  # (t, trailing-run-len) for 3' overtrims
        dropped = 0
        for pid, (h1, t1, h2, t2) in gt.items():
            for mate, (h, t) in ((0, (h1, t1)), (1, (h2, t2))):
                s = px.get((pid, mate))
                if s is None or s[0] != 0:
                    dropped += 1
                    continue
                _, t5, t3 = s
                for tol in (0, 2):
                    if abs(t5 - h) <= tol and abs(t3 - t) <= tol:
                        cnt[tol]["exact"] += 1
                    elif t5 > h or t3 > t:
                        cnt[tol]["over"] += 1
                    else:
                        cnt[tol]["under"] += 1
                if t5 > h:
                    over5_mag[t5 - h] += 1
                    seq = seqs.get((pid, mate), "")
                    run = 0
                    while run < len(seq) and seq[run] == 'G':
                        run += 1
                    over5_ctx[(h, run)] += 1
                if t3 > t:
                    over3_mag[t3 - t] += 1
                    seq = seqs.get((pid, mate), "")
                    run = 0
                    while run < len(seq) and seq[len(seq) - 1 - run] == 'G':
                        run += 1
                    over3_ctx[(t, run)] += 1

        print(f"\n===== {regime} =====  kept_scored={cnt[0]['exact']+cnt[0]['over']+cnt[0]['under']} dropped={dropped}")
        for tol in (0, 2):
            c = cnt[tol]; n = max(1, sum(c.values()))
            print(f"  tol={tol}: exact={100*c['exact']/n:.2f}%  over={100*c['over']/n:.2f}%  under={100*c['under']/n:.2f}%")
        print(f"  5' overtrim error magnitudes: {dict(sorted(over5_mag.items()))}")
        print(f"  3' overtrim error magnitudes: {dict(sorted(over3_mag.items()))}")
        print(f"  5' overtrim (gt_h, leading_Grun) -> count: {dict(sorted(over5_ctx.items()))}")
        print(f"  3' overtrim (gt_t, trailing_Grun) -> count: {dict(sorted(over3_ctx.items()))}")

if __name__ == "__main__":
    main()
