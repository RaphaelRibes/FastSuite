#!/usr/bin/env python3
"""Task 2.12: is the dup module's memory LINEAR in reads, or does it double?

`--only dup` peaks at 688 MB while its fingerprint table projects at 272 MB.
If the peak is linear in read count and ~2.5x the table everywhere, something
holds a parallel copy of the table. If it jumps at a capacity boundary, it is
rehash. If it is linear with a slope of ~68 B/read, there is a per-read
structure nobody accounted for.

Sweeps `--reads` and prints peak RSS, the projected table size parsed from the
tool's own stderr, and the ratio between them. Peak RSS is from wait4.
"""
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FC = os.path.join(REPO, "tools", "fastcheck", "target", "release", "fastcheck")
CORPUS = "/mnt/Fast/fastsuite-scratch/fm-ab5m/prep"
OUT = "/mnt/Fast/fastsuite-scratch/fc-scale"


def run(cmd, logpath):
    pid = os.fork()
    if pid == 0:
        try:
            fd = os.open(logpath, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
            os.dup2(fd, 1)
            os.dup2(fd, 2)
            os.execvp(cmd[0], cmd)
        finally:
            os._exit(127)
    _, status, ru = os.wait4(pid, 0)
    return ru.ru_maxrss / 1024.0, os.waitstatus_to_exitcode(status)


def main():
    if not os.path.exists(FC):
        raise SystemExit(f"build fastcheck first: {FC}")
    os.makedirs(OUT, exist_ok=True)
    r1 = os.path.join(CORPUS, "clean_R1.fastq.gz")
    r2 = os.path.join(CORPUS, "clean_R2.fastq.gz")
    print(f"{'reads':>10} {'peak MB':>9} {'table MB':>9} {'peak/table':>11}")
    rows = []
    for reads in (250_000, 500_000, 1_000_000, 2_000_000, 5_000_000):
        log = os.path.join(OUT, f"r{reads}.log")
        cmd = [FC, r1, r2, "-j", "8", "-q", "--only", "dup",
               "--reads", str(reads)]
        peak, rc = run(cmd, log)
        if rc != 0:
            print(f"{reads}: exit {rc}", file=sys.stderr)
            print(open(log).read()[-400:], file=sys.stderr)
            continue
        text = open(log).read()
        m = re.search(r"projected at ([\d.]+) MB", text)
        table = float(m.group(1)) if m else 0.0
        ratio = peak / table if table else 0.0
        print(f"{reads:>10,} {peak:>9.1f} {table:>9.1f} {ratio:>11.2f}")
        rows.append((reads, peak, table))

    print()
    if len(rows) >= 2:
        (r_a, p_a, _), (r_b, p_b, _) = rows[0], rows[-1]
        slope = (p_b - p_a) / ((r_b - r_a) / 1e6)
        print(f"slope over the sweep: {slope:.1f} MB per million reads")
        print(f"  a 12 B/entry table would be ~12 MB/M-read; 17 B/bucket at "
              f"0.875 load is ~19.4 MB/M-read")
        if slope > 40:
            print("  -> far above BOTH: something per-read is being kept that is")
            print("     not the fingerprint table.")
        else:
            print("  -> consistent with the fingerprint table dominating.")


if __name__ == "__main__":
    main()
