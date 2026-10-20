#!/usr/bin/env python3
"""Does fastdedup's peak memory scale with -j?

Same 5 M-pair input, same machine, only the thread count changes. Peak RSS is the
reaped child's ru_maxrss from wait4. Each result is appended to threads.csv
immediately, so an interrupted run still leaves the measurements it completed.
"""
import csv
import os
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "tools/fastdedup/target/release/fastdedup")
SRC = os.environ.get("FS_SIM_DIR",
                     "/home/raphael/.hermes/cache/scratch/fsprobe/pipe5m/sim")
OUT = os.environ.get("FS_OUT_DIR",
                     "/home/raphael/.hermes/cache/scratch/fsprobe/dedup_threads")
CSV_PATH = os.path.join(OUT, "threads.csv")


def run(cmd, log):
    t0 = time.perf_counter()
    pid = os.fork()
    if pid == 0:
        try:
            fd = os.open(log, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
            os.dup2(fd, 1)
            os.dup2(fd, 2)
            os.execvp(cmd[0], cmd)
        finally:
            os._exit(127)
    _, status, ru = os.wait4(pid, 0)
    return (time.perf_counter() - t0, ru.ru_maxrss / 1024.0,
            os.waitstatus_to_exitcode(status), ru.ru_utime + ru.ru_stime)


def main():
    os.makedirs(OUT, exist_ok=True)
    threads = [int(x) for x in (sys.argv[1:] or ["1", "2", "4", "8"])]
    if not os.path.exists(CSV_PATH):
        with open(CSV_PATH, "w", newline="") as fh:
            csv.writer(fh).writerow(
                ["threads", "pairs", "wall_s", "cpu_s", "peak_mb", "exit"])

    for j in threads:
        cmd = [BIN, "-1", f"{SRC}/sim_R1.fastq.gz", "-2", f"{SRC}/sim_R2.fastq.gz",
               "-o", f"{OUT}/t{j}_R1.fastq.gz", "-p", f"{OUT}/t{j}_R2.fastq.gz",
               "-j", str(j), "--gz-level", "6", "--force", "-q"]
        wall, peak, rc, cpu = run(cmd, f"{OUT}/t{j}.log")
        row = [j, 5_000_000, round(wall, 2), round(cpu, 2), round(peak, 1), rc]
        with open(CSV_PATH, "a", newline="") as fh:
            csv.writer(fh).writerow(row)
        print(f"-j {j:<2} wall {wall:8.2f}s cpu {cpu:7.2f}s peak {peak:8.1f} MB exit {rc}",
              flush=True)
        for f in (f"{OUT}/t{j}_R1.fastq.gz", f"{OUT}/t{j}_R2.fastq.gz"):
            if os.path.exists(f):
                os.remove(f)


if __name__ == "__main__":
    main()
