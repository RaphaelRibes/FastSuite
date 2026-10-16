#!/usr/bin/env python3
"""Measure the boundary-anchored cut against the shipped whole-run cut.

Two modes:

  accuracy --regimes low mid high
      Runs fastprep on each fastsim profile (exact GT) with and without
      --boundary-anchored and scores POLY-only exact/over/under at tolerance
      0 and 2 (adapter trims are excluded — see overtrim_poly.py). Also
      reports the read-drop count so a change in the length budget's drops
      is visible.

  speed --pairs N
      Generates an N-pair paired dataset from the mid profile and runs both
      cuts on it, reporting wall time, throughput (pairs/s), peak RSS
      (/usr/bin/time -v "Maximum resident set size"), and drops.

Ground truth scoring reuses the bench/scripts/overtrim_sensitivity.py PXMT
convention (pair_id = parts[-1]*10000 + parts[-2]).
"""
from __future__ import annotations
import argparse, csv, gzip, os, re, struct, subprocess, sys, time

FASTPREP = "/home/raphael/hermes_scratch/FastSuite/tools/fastprep/target/release/fastprep"
FASTSIM = "/home/raphael/hermes_scratch/FastSuite/tools/fastsim/target/release/fastsim"
WORKROOT = "/home/raphael/.hermes/cache/scratch/overtrim"
PROFILE_DIR = "/home/raphael/hermes_scratch/FastSuite/bench/fixtures/overtrim"

# ---------------- PXMT / GT ----------------
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
        out[(pid, mate)] = (status, t5p, t3p)
    return out

def parse_gt(path):
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[int(r["pair_id"])] = (int(r["r1_h"]), int(r["r1_t"]),
                                      int(r["r2_h"]), int(r["r2_t"]))
    return out

def classify(p5, p3, h, t, tol):
    if abs(p5 - h) <= tol and abs(p3 - t) <= tol:
        return "exact"
    if p5 > h or p3 > t:
        return "over"
    return "under"

def score(px, gt, tol):
    c = {"exact": 0, "over": 0, "under": 0}
    dropped = 0
    for pid, (h1, t1, h2, t2) in gt.items():
        for mate, (h, t) in ((0, (h1, t1)), (1, (h2, t2))):
            s = px.get((pid, mate))
            if s is None or s[0] != 0:
                dropped += 1
                continue
            c[classify(s[1], s[2], h, t, tol)] += 1
    n = max(1, c["exact"] + c["over"] + c["under"])
    return {"exact": c["exact"] / n, "over": c["over"] / n, "under": c["under"] / n,
            "dropped": dropped, "scored": n}

# ---------------- fastprep invocation ----------------
def base_argv(work, r1, r2, out_prefix, pxmt, extra):
    return [FASTPREP, "-1", r1, "-2", r2,
            "-o", f"{out_prefix}_R1.fastq.gz", "-p", f"{out_prefix}_R2.fastq.gz",
            "--metadata", pxmt, "--json", f"{out_prefix}.json", "--force", "-q", "-j", "4",
            "--min-run", "12", "--mercy", "3", "--phred-sat", "20", "--gate-shape", "1.0",
            "--min-len", "0.9L", "--mean-q", "25.0", "--adapter-error-rate", "0.2",
            "--min-adapter-overlap", "14", "--gz-level", "6",
            "--match-leniency", "0.1", "--mismatch-strength", "0.5",
            "--head-dropout", "1"] + extra

def _child_peak_rss_mb(pid):
    """Peak RSS (VmHWM) of a live child in MB; 0.0 once gone.

    ru_maxrss (wait4/getrusage) is a running high-water mark over every child
    the process ever spawned and is contaminated by fork-time inheritance of
    the parent's RSS — it cannot measure one child. Poll the child's own
    /proc VmHWM instead (same method as bench/scripts/fastp_compare.py).
    """
    try:
        with open(f"/proc/{pid}/status") as fh:
            for line in fh:
                if line.startswith("VmHWM:"):
                    return int(line.split()[1]) / 1024.0
    except (FileNotFoundError, ProcessLookupError, ValueError, IndexError):
        pass
    return 0.0

def run_measured(argv):
    """Run and return (returncode, wall_s, peak_rss_mb, stderr).

    Peak RSS is the child's own VmHWM polled every 20 ms for its whole life.
    """
    err_path = argv[argv.index("--metadata") + 1] + ".speed.err"
    t0 = time.time()
    with open(err_path, "w") as errf:
        p = subprocess.Popen(argv, stdout=subprocess.DEVNULL, stderr=errf)
        peak = 0.0
        while p.poll() is None:
            peak = max(peak, _child_peak_rss_mb(p.pid))
            time.sleep(0.02)
        peak = max(peak, _child_peak_rss_mb(p.pid))
        rc = p.returncode
    wall = time.time() - t0
    err = open(err_path).read()
    return rc, wall, peak, err

# ---------------- accuracy mode ----------------
def accuracy(regimes):
    print(f"{'regime':6s} {'cut':14s} | {'exact0':>7s} {'over0':>7s} {'und0':>7s} | {'exact2':>7s} {'over2':>7s} {'und2':>7s} | {'drop':>5s}")
    for regime in regimes:
        work = f"{WORKROOT}/{regime}"
        gt = parse_gt(f"{work}/sim.v0.gt.tsv")
        for tag, extra in [("whole-run", []), ("boundary-anchored", ["--boundary-anchored"])]:
            stem = tag.replace('-', '_')
            pxmt = f"{work}/_impl_{stem}.pxmt"
            argv = base_argv(work, f"{work}/sim.v0_R1.fastq.gz", f"{work}/sim.v0_R2.fastq.gz",
                             f"{work}/_impl_{stem}", pxmt, extra)
            rc = subprocess.run(argv, capture_output=True, text=True).returncode
            if rc != 0:
                print(f"{regime:6s} {tag:14s} rc={rc}")
                continue
            px = parse_pxmt(pxmt)
            r0 = score(px, gt, 0); r2 = score(px, gt, 2)
            print(f"{regime:6s} {tag:14s} | {r0['exact']:7.4f} {r0['over']:7.4f} {r0['under']:7.4f} | "
                  f"{r2['exact']:7.4f} {r2['over']:7.4f} {r2['under']:7.4f} | {r0['dropped']:5d}")

# ---------------- speed mode ----------------
def gen_pairs(n, outdir):
    os.makedirs(outdir, exist_ok=True)
    prof = f"{PROFILE_DIR}/mid.json"
    argv = [FASTSIM, "--pairs", str(n), "--profile", prof, "-o", f"{outdir}/speed",
            "--seed", "7", "--dup-rate", "0.02", "--drop-rate", "0.05", "-f", "-q"]
    p = subprocess.run(argv, capture_output=True, text=True)
    if p.returncode != 0:
        print("fastsim failed:", p.stderr[-300:]); sys.exit(1)

def speed(pairs):
    outdir = "/home/raphael/.hermes/cache/scratch/overtrim/speed"
    r1, r2 = f"{outdir}/speed_R1.fastq.gz", f"{outdir}/speed_R2.fastq.gz"
    if not (os.path.exists(r1) and os.path.exists(r2)):
        print(f"generating {pairs} pairs ...")
        gen_pairs(pairs, outdir)
    print(f"{'cut':20s} {'wall_s':>8s} {'pairs/s':>10s} {'Mbp/s':>8s} {'peakRSS':>10s} {'drop':>6s}")
    for tag, extra in [("whole-run", []), ("boundary-anchored", ["--boundary-anchored"])]:
        stem = tag.replace('-', '_')
        pxmt = f"{outdir}/_speed_{stem}.pxmt"
        argv = base_argv(outdir, r1, r2, f"{outdir}/_speed_{stem}", pxmt, extra)
        rc, wall, rss, err = run_measured(argv)
        if rc != 0:
            print(f"{tag:20s} rc={rc} {err[-200:]}"); continue
        gt = parse_gt(f"{outdir}/speed.gt.tsv")
        px = parse_pxmt(pxmt)
        drops = sum(1 for pid, (h1,t1,h2,t2) in gt.items()
                    for mate in (0,1)
                    if (s := px.get((pid, mate))) is None or s[0] != 0)
        pairs_per_s = pairs / wall if wall else 0
        mbp_per_s = (pairs * 2 * 150 / 1e6) / wall if wall else 0
        print(f"{tag:20s} {wall:8.3f} {pairs_per_s:10.0f} {mbp_per_s:8.1f} {rss:8.1f}MB {drops:6d}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["accuracy", "speed"])
    ap.add_argument("--regimes", nargs="+", default=["low", "mid", "high"])
    ap.add_argument("--pairs", type=int, default=100_000)
    a = ap.parse_args()
    if a.mode == "accuracy":
        accuracy(a.regimes)
    else:
        speed(a.pairs)

if __name__ == "__main__":
    main()
