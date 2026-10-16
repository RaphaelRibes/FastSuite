#!/usr/bin/env python3
"""Timing + peak RSS for the head-dropout settings and boundary-anchored."""
import os, subprocess, tempfile, time, resource, json

FS = "tools/fastsim/target/release/fastsim"
FP = "tools/fastprep/target/release/fastprep"

work = tempfile.mkdtemp(prefix="spd_")
prof = "bench/fixtures/overtrim/mid.json"
out = f"{work}/sim"
subprocess.run([FS, "--pairs", "100000", "--profile", prof, "--out", out,
                "--seed", "42", "--read-len", "150", "--insert-min", "1",
                "--drop-rate", "0", "--force"], check=True, capture_output=True)

variants = [
    ("hd=4 default", []),
    ("hd=1", ["--head-dropout", "1"]),
    ("hd=0", ["--head-dropout", "0"]),
    ("hd=4 + boundary-anchored", ["--boundary-anchored"]),
    ("hd=0 + boundary-anchored", ["--head-dropout", "0", "--boundary-anchored"]),
]
print(f"{'variant':30s} {'median_s':>9s} {'pairs/s':>10s} {'peak_MB':>9s}")
for name, extra in variants:
    times, rss = [], []
    for rep in range(3):
        d = f"{work}/t{len(extra)}_{rep}"
        os.makedirs(d, exist_ok=True)
        usage0 = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        t0 = time.monotonic()
        r = subprocess.run([FP, "-1", f"{out}_R1.fastq.gz", "-2", f"{out}_R2.fastq.gz",
                            "-o", f"{d}/r1.fq.gz", "-p", f"{d}/r2.fq.gz",
                            "--no-dedup", "--no-mean-q", "-q", "-j", "4"] + extra,
                           capture_output=True)
        dt = time.monotonic() - t0
        if r.returncode != 0:
            print(f"{name:30s} FAILED rc={r.returncode}")
            break
        times.append(dt)
        rss.append(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    if times:
        times.sort()
        print(f"{name:30s} {times[len(times)//2]:9.3f} "
              f"{100000/times[len(times)//2]:10.0f} {max(rss)/1024:9.1f}")
