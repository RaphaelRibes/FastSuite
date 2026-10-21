#!/usr/bin/env python3
"""Task 2.4a: would flipping fastprep's default gzip level to 1 pay?

Task 2.4b measured fastdedup and found level 1 is 1.16x faster at -j 8 with
25x less peak RSS, at 1.87x the bytes. fastprep is a different tool with a
different writer (BGZF via gzp, not gzippy), so it needs its own measurement
rather than an extrapolation.

Arms, on one corpus, at the shipped thread count:
  * --gz-level 6 (today's default)
  * --gz-level 1 (proposed)
and the same pair at -j 1, so the thread effect is separable.

Reports wall, CPU, busy cores, peak RSS, output bytes, and whether `gzip -t`
accepts the output. Peak RSS is the reaped child's ru_maxrss from wait4.
"""
import os
import shutil
import statistics
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FP = os.path.join(REPO, "tools", "fastprep", "target", "release", "fastprep")
CORPUS = "/mnt/Fast/fastsuite-scratch/fm-ab5m/prep"
OUT = "/mnt/Fast/fastsuite-scratch/prep-gz"


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


def arm(label, level, threads, reps=3):
    r1 = os.path.join(CORPUS, "clean_R1.fastq.gz")
    r2 = os.path.join(CORPUS, "clean_R2.fastq.gz")
    walls, cpus, peaks, sizes = [], [], [], []
    for _ in range(reps):
        d = os.path.join(OUT, label)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d, exist_ok=True)
        o1 = os.path.join(d, "o1.fq.gz")
        o2 = os.path.join(d, "o2.fq.gz")
        cmd = [FP, "-1", r1, "-2", r2, "-o", o1, "-p", o2,
               "-j", str(threads), "--gz-level", str(level), "--force", "-q"]
        wall, cpu, peak, rc = run(cmd, os.path.join(OUT, f"{label}.log"))
        if rc != 0:
            print(f"{label}: exit {rc}", file=sys.stderr)
            print(open(os.path.join(OUT, f"{label}.log")).read()[-600:], file=sys.stderr)
            return None
        walls.append(wall); cpus.append(cpu); peaks.append(peak)
        sizes.append(os.path.getsize(o1))
    w = statistics.median(walls)
    c = statistics.median(cpus)
    res = dict(label=label, level=level, threads=threads, wall=w, cpu=c,
               busy=c / w, peak=max(peaks), size=max(sizes))
    print(f"{label:12} wall {w:6.2f}s  cpu {c:7.2f}s  busy {c/w:4.2f}  "
          f"peak {res['peak']:7.1f} MB  out {res['size']/1e6:6.1f} MB")
    return res


def gzip_ok(path):
    return subprocess.run(["gzip", "-t", path], capture_output=True).returncode == 0


def main():
    if not os.path.exists(FP):
        raise SystemExit(f"build fastprep first: {FP}")
    if not os.path.isdir(CORPUS):
        raise SystemExit(f"corpus missing: {CORPUS}")
    os.makedirs(OUT, exist_ok=True)
    print(f"in: {os.path.getsize(os.path.join(CORPUS,'clean_R1.fastq.gz'))/1e6:.0f} MB per mate\n")

    print("=== -j 8 ===")
    a = arm("gz6_j8", 6, 8)
    b = arm("gz1_j8", 1, 8)
    print("\n=== -j 1 ===")
    c = arm("gz6_j1", 6, 1)
    d = arm("gz1_j1", 1, 1)
    for name, res in (("gz6_j8", a), ("gz1_j8", b), ("gz6_j1", c), ("gz1_j1", d)):
        if res is None:
            raise SystemExit(f"arm {name} failed; see the logs")
    assert a is not None and b is not None and c is not None and d is not None

    print("\n=== standard gzip? ===")
    for lab in ("gz6_j8", "gz1_j8"):
        p = os.path.join(OUT, lab, "o1.fq.gz")
        print(f"  {lab:8} gzip -t: {'OK' if gzip_ok(p) else 'FAILED'}")

    print("\n=== verdict ===")
    sp = a["wall"] / b["wall"]
    print(f"  level 1 is {sp:.2f}x {'faster' if sp > 1 else 'SLOWER'} at -j 8 "
          f"({b['wall']:.2f}s vs {a['wall']:.2f}s)")
    print(f"  peak RSS {a['peak']:.0f} -> {b['peak']:.0f} MB")
    print(f"  bytes {a['size']/1e6:.0f} -> {b['size']/1e6:.0f} MB ({b['size']/a['size']:.2f}x)")


if __name__ == "__main__":
    main()
