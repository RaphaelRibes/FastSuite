#!/usr/bin/env python3
"""Does `fasthelper concatenate` parse each report once, or once per format?

Task 2.20 claims `check -> aggregate -> concatenate -> html` re-reads every
report. The source says otherwise: `collect()` runs once and `write_reports`
renders every format from the same in-memory Vec. This probe measures it
rather than trusting the read.

Method: make the PARSE dominate. Each report carries a large heatmap-like
`counts` array, so a re-parse would multiply that cost. Then time the same
input with one output format against four:

  * parse once  -> t(4 formats) = t(1) + 3 * render, and render is cheap
  * parse per format -> t(4 formats) ~= 4 * t(1)

Prints the ratio. A ratio near 1.0-1.5 means parsed once; near 4.0 means
parsed four times.
"""
import json
import os
import shutil
import statistics
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "tools", "fasthelper", "target", "release", "fasthelper")
N_SAMPLES = 40
GRID = 200  # 200x200 cells -> 40 000 counts per report


def make_corpus(d):
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d, exist_ok=True)
    counts = [((i * 7) % 5000) for i in range(GRID * GRID)]
    for s in range(N_SAMPLES):
        sd = os.path.join(d, f"sample{s:03d}")
        os.makedirs(sd, exist_ok=True)
        rep = {
            "schema": "fastsuite.report.v1",
            "tool": "fastmaps",
            "tool_version": "0.1.0",
            "inputs": {"files": [f"s{s}_R1.fq.gz", f"s{s}_R2.fq.gz"]},
            "metrics": {"total_reads": 1000000 + s, "size": GRID},
            "counts": counts,
        }
        with open(os.path.join(sd, "fastmaps.json"), "w") as fh:
            json.dump(rep, fh)
    total = sum(
        os.path.getsize(os.path.join(d, f"sample{s:03d}", "fastmaps.json"))
        for s in range(N_SAMPLES)
    )
    print(f"corpus: {N_SAMPLES} reports, {total / 1e6:.1f} MB total")


def time_once(d, out, fmt, reps):
    times = []
    for _ in range(reps):
        shutil.rmtree(out, ignore_errors=True)
        t0 = time.perf_counter()
        r = subprocess.run(
            [BIN, "concatenate", d, "-o", out, "--format", fmt],
            capture_output=True, text=True)
        times.append(time.perf_counter() - t0)
        if r.returncode != 0:
            print(r.stdout, r.stderr, file=sys.stderr)
            raise SystemExit(f"concatenate failed: {r.returncode}")
    return statistics.median(times)


def main():
    if not os.path.exists(BIN):
        raise SystemExit(f"build fasthelper first: {BIN} not found")
    d = "/tmp/fh-scaling"
    make_corpus(d)
    reps = 5

    one = time_once(d, "/tmp/fh-out1", "json", reps)
    four = time_once(d, "/tmp/fh-out4", "json,tsv,md,html", reps)
    print(f"1 format  (json)          : {one:.3f} s")
    print(f"4 formats (json,tsv,md,html): {four:.3f} s")
    print(f"ratio 4/1                 : {four / one:.2f}x")
    print()
    if four / one < 2.0:
        print("VERDICT: parsed ONCE — adding three output formats costs render only.")
    elif four / one > 3.0:
        print("VERDICT: parsed PER FORMAT — every format re-reads the input.")
    else:
        print("VERDICT: ambiguous; the corpus is not parse-dominated enough.")
    shutil.rmtree("/tmp/fh-out1", ignore_errors=True)
    shutil.rmtree("/tmp/fh-out4", ignore_errors=True)


if __name__ == "__main__":
    main()
