#!/usr/bin/env python3
"""Same-machine, same-input profile of the whole FastSuite pipeline.

Every stage runs on the output of the previous one, on one machine, so the
per-stage costs are comparable to each other (they are NOT comparable to the
numbers in the individual papers, which were measured on older builds).

Peak RSS is the reaped child's ru_maxrss from wait4 — not a /proc poll, which
under-reports anything that finishes between polls.
"""
import csv
import os
import sys
import time

# Repo root from this file's location (bench/profile/ -> repo), so the script
# works from whichever checkout it lives in: both machines keep two, and a
# hardcoded home path silently ran the wrong binaries.
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TOOLS = os.path.join(REPO, "tools")
BIN = {t: f"{TOOLS}/{t}/target/release/{t}" for t in
       ("fastsim", "fastdedup", "fastprep", "fastmaps", "fastcheck")}
PROFILE = f"{TOOLS}/fastsim/profiles/cellB-pe-overlap.json"


def run(cmd, logpath):
    """Run cmd, return (wall_s, peak_mb, exit_code) for THAT child."""
    t0 = time.perf_counter()
    pid = os.fork()
    if pid == 0:  # child
        try:
            fd = os.open(logpath, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
            os.dup2(fd, 1)
            os.dup2(fd, 2)
            os.execvp(cmd[0], cmd)
        finally:
            os._exit(127)
    _, status, ru = os.wait4(pid, 0)
    return (time.perf_counter() - t0, ru.ru_maxrss / 1024.0,
            os.waitstatus_to_exitcode(status), ru.ru_utime + ru.ru_stime)


def main():
    pairs = int(sys.argv[1])
    j = sys.argv[2]
    out = sys.argv[3]
    logs = os.path.join(out, "logs")
    os.makedirs(logs, exist_ok=True)
    s = os.path.join(out, "sim")
    os.makedirs(s, exist_ok=True)
    rows = []

    stages = [
        ("fastsim", [BIN["fastsim"], "--profile", PROFILE, "--pairs", str(pairs),
                     "--seed", "1", "--out", f"{s}/sim"]),
        ("fastdedup", [BIN["fastdedup"], "-1", f"{s}/sim_R1.fastq.gz",
                       "-2", f"{s}/sim_R2.fastq.gz", "-o", f"{s}/dedup_R1.fastq.gz",
                       "-p", f"{s}/dedup_R2.fastq.gz", "-j", j, "--gz-level", "6",
                       "--json", f"{s}/fastdedup.json", "--force", "-q"]),
        ("fastprep", [BIN["fastprep"], "-1", f"{s}/dedup_R1.fastq.gz",
                      "-2", f"{s}/dedup_R2.fastq.gz", "-o", f"{s}/clean_R1.fastq.gz",
                      "-p", f"{s}/clean_R2.fastq.gz", "-j", j, "--gz-level", "6",
                      "--json", f"{s}/fastprep.json", "--force", "-q"]),
        ("fastmaps", [BIN["fastmaps"], f"{s}/clean_R1.fastq.gz",
                      f"{s}/clean_R2.fastq.gz", "-o", f"{s}/maps", "-j", j,
                      "--json", f"{s}/maps/fastmaps.json", "--force", "-q"]),
        ("fastcheck", [BIN["fastcheck"], f"{s}/clean_R1.fastq.gz",
                       f"{s}/clean_R2.fastq.gz", "-j", j,
                       "--json", f"{s}/fastcheck.json",
                       "--html", f"{s}/fastcheck.html", "-q"]),
    ]

    os.makedirs(f"{s}/maps", exist_ok=True)  # fastmaps needs --json's parent

    for name, cmd in stages:
        wall, peak, rc, cpu = run(cmd, os.path.join(logs, f"{name}.log"))
        rows.append(dict(stage=name, pairs=pairs, threads=j, wall_s=round(wall, 3),
                         cpu_s=round(cpu, 3), peak_mb=round(peak, 1), exit=rc))
        print(f"{name:10} {wall:8.2f} s  {peak:8.1f} MB  cpu {cpu:7.2f} s  exit {rc}")
        if rc != 0:
            print(f"  stage failed, see {logs}/{name}.log")

    for f in ("sim_R1.fastq.gz", "dedup_R1.fastq.gz", "clean_R1.fastq.gz"):
        p = os.path.join(s, f)
        if os.path.exists(p):
            print(f"  {f}: {os.path.getsize(p) / 1e6:.1f} MB")

    with open(os.path.join(out, "timing.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
