#!/usr/bin/env python3
"""Task 2.6 (FD2-4): diagnose the fastdedup -j 1..4 flatline.

The recorded measurement (bench/profile/dedup_threads.csv, post-queue-fix):
j1 94.84 s, j2 92.13 s, j4 92.2 s, but j8 22.73 s. Four threads doing no
better than one, then one thread more than tripling the rate, is not a
parallel-scaling curve — something switches on at 8.

Two candidate explanations, and this probe separates them:

  1. the writer queue bound (53d6686), so workers stall on compression;
  2. single-threaded deflate, so the wall time is compression, not hashing.

The separating experiment is the plan's own Step 1: run with PLAIN output,
where there is no compression at all. If every -j lands on the same short wall
time, the 82 s gap is entirely in the compressed-output path.

`parallel.rs` uses gzippy (parallel compression) only when
`num_threads > 1 && compression_level >= 6`, so gz6 at j2/j4 SHOULD take the
parallel path. If it does and is still slow, the level floor is not the
explanation either.

Peak RSS is the reaped child's ru_maxrss from wait4.
"""
import os
import shutil
import statistics
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FD = os.path.join(REPO, "tools", "fastdedup", "target", "release", "fastdedup")
CORPUS = "/mnt/Fast/fastsuite-scratch/fm-ab5m/prep"
OUT = "/mnt/Fast/fastsuite-scratch/dedup-flatline"


def run(cmd, logpath):
    t0 = time.perf_counter()
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
    return (time.perf_counter() - t0, ru.ru_utime + ru.ru_stime,
            ru.ru_maxrss / 1024.0, os.waitstatus_to_exitcode(status))


def arm(label, level, threads, gz, reps=2):
    r1 = os.path.join(CORPUS, "clean_R1.fastq.gz")
    r2 = os.path.join(CORPUS, "clean_R2.fastq.gz")
    walls, cpus, peaks = [], [], []
    for _ in range(reps):
        d = os.path.join(OUT, label)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d, exist_ok=True)
        ext = "fq.gz" if gz else "fq"
        o1 = os.path.join(d, f"o1.{ext}")
        o2 = os.path.join(d, f"o2.{ext}")
        cmd = [FD, "-1", r1, "-2", r2, "-o", o1, "-p", o2,
               "-j", str(threads), "--gz-level", str(level), "--force", "-q"]
        wall, cpu, peak, rc = run(cmd, os.path.join(OUT, f"{label}.log"))
        if rc != 0:
            print(f"{label}: exit {rc}", file=sys.stderr)
            print(open(os.path.join(OUT, f"{label}.log")).read()[-500:], file=sys.stderr)
            return None
        walls.append(wall); cpus.append(cpu); peaks.append(peak)
    w = statistics.median(walls)
    c = statistics.median(cpus)
    print(f"{label:14} wall {w:6.2f}s  cpu {c:7.2f}s  busy {c/w:5.2f}  peak {max(peaks):7.1f} MB")
    return dict(label=label, wall=w, cpu=c, busy=c / w, peak=max(peaks))


def main():
    if not os.path.exists(FD):
        raise SystemExit(f"build fastdedup first: {FD}")
    if not os.path.isdir(CORPUS):
        raise SystemExit(f"corpus missing: {CORPUS}")
    os.makedirs(OUT, exist_ok=True)
    print("5M pairs, median of 2 reps\n")

    print("=== gz level 6, the recorded flatline ===")
    g = {}
    for j in (1, 2, 4, 8):
        g[j] = arm(f"gz6_j{j}", 6, j, True)

    print("\n=== PLAIN output: no compression at all (plan Step 1) ===")
    p = {}
    for j in (1, 4, 8):
        p[j] = arm(f"plain_j{j}", 1, j, False)

    if any(v is None for v in list(g.values()) + list(p.values())):
        raise SystemExit("an arm failed; see the logs")

    print("\n=== diagnosis ===")
    plain_flat = max(p[j]["wall"] for j in p) / min(p[j]["wall"] for j in p)
    print(f"  plain output, j1/j4/j8 = {p[1]['wall']:.2f} / {p[4]['wall']:.2f} / "
          f"{p[8]['wall']:.2f} s  (spread {plain_flat:.2f}x)")
    if plain_flat < 1.3:
        print("  -> the hash/parse pipeline is FLAT in -j, as the plan expected.")
        print("     The whole gap is therefore in the COMPRESSED-output path.")
    else:
        print("  -> even with plain output the cost scales with -j, so parsing")
        print("     is not flat and the queue is not the whole story either.")

    print(f"\n  gz6 j1/j2/j4/j8 = {g[1]['wall']:.2f} / {g[2]['wall']:.2f} / "
          f"{g[4]['wall']:.2f} / {g[8]['wall']:.2f} s")
    print(f"  busy cores      = {g[1]['busy']:.2f} / {g[2]['busy']:.2f} / "
          f"{g[4]['busy']:.2f} / {g[8]['busy']:.2f}")
    print(f"  peak RSS        = {g[1]['peak']:.0f} / {g[2]['peak']:.0f} / "
          f"{g[4]['peak']:.0f} / {g[8]['peak']:.0f} MB")

    j4_vs_j1 = g[1]["wall"] / g[4]["wall"]
    print(f"\n  j4 is {j4_vs_j1:.2f}x j1's speed; j8 is {g[1]['wall']/g[8]['wall']:.2f}x j1's")
    if g[4]["busy"] < 2.0:
        print("  -> j4 uses under 2 cores' worth of CPU: compression is NOT running")
        print("     in parallel there, despite threads>1 and level>=6.")
    else:
        print("  -> j4 does use multiple cores, so the parallel path is engaged;")
        print("     the wall time is then dominated by something other than deflate.")


if __name__ == "__main__":
    main()
