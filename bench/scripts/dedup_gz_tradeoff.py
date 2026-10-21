#!/usr/bin/env python3
"""Task 2.4b: is flipping fastdedup's default gzip level to 1 a net loss?

The plan's own warning, verbatim: "This interaction is the single most likely
way the gzip lever backfires." The interaction is that `parallel.rs` uses
gzippy — the parallel compressor that makes `-j 8` fast — ONLY when
`compression_level >= 6`, because below 6 gzippy emits a proprietary "GZ"
multi-block format that standard tools cannot read.

So a level-1 default does not merely trade bytes for speed. It drops the
output onto the single-threaded flate2 path, which is what the flat
`-j 1..4` measurements (94.84 / 92.13 / 92.2 s) were showing.

Measures, on one corpus, at -j 8 and -j 1:
  * wall, CPU, busy cores, peak RSS
  * output bytes
  * whether the output is standard gzip (readable by the system `gzip -t`)

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
OUT = "/mnt/Fast/fastsuite-scratch/dedup-gz"


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
    wall = time.perf_counter() - t0
    cpu = ru.ru_utime + ru.ru_stime
    return wall, cpu, ru.ru_maxrss / 1024.0, os.waitstatus_to_exitcode(status)


def arm(label, level, threads, reps=3):
    walls, cpus, peaks, sizes = [], [], [], []
    r1 = os.path.join(CORPUS, "clean_R1.fastq.gz")
    r2 = os.path.join(CORPUS, "clean_R2.fastq.gz")
    for _ in range(reps):
        d = os.path.join(OUT, label)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d, exist_ok=True)
        o1 = os.path.join(d, "o1.fq.gz")
        o2 = os.path.join(d, "o2.fq.gz")
        cmd = [FD, "-1", r1, "-2", r2, "-o", o1, "-p", o2,
               "-j", str(threads), "--gz-level", str(level), "--force", "-q"]
        wall, cpu, peak, rc = run(cmd, os.path.join(OUT, f"{label}.log"))
        if rc != 0:
            print(f"{label}: exit {rc}", file=sys.stderr)
            print(open(os.path.join(OUT, f"{label}.log")).read()[-600:], file=sys.stderr)
            return None
        walls.append(wall); cpus.append(cpu); peaks.append(peak)
        sizes.append(os.path.getsize(o1))
    res = dict(label=label, level=level, threads=threads,
               wall=statistics.median(walls), cpu=statistics.median(cpus),
               busy=statistics.median(cpus) / statistics.median(walls),
               peak=max(peaks), size=max(sizes))
    print(f"{label:22} wall {res['wall']:7.2f}s  cpu {res['cpu']:8.2f}s  "
          f"busy {res['busy']:4.2f}  peak {res['peak']:7.1f} MB  "
          f"out {res['size']/1e6:7.1f} MB")
    return res


def gzip_ok(path):
    r = subprocess.run(["gzip", "-t", path], capture_output=True)
    return r.returncode == 0


def main():
    if not os.path.exists(FD):
        raise SystemExit(f"build fastdedup first: {FD}")
    if not os.path.isdir(CORPUS):
        raise SystemExit(f"corpus missing: {CORPUS}")
    os.makedirs(OUT, exist_ok=True)
    print(f"corpus: {CORPUS}")
    print(f"in: {os.path.getsize(os.path.join(CORPUS,'clean_R1.fastq.gz'))/1e6:.1f} MB per mate\n")

    print("=== the shipped default, and the proposed one, at -j 8 ===")
    a = arm("gz6_j8", 6, 8)
    b = arm("gz1_j8", 1, 8)

    print("\n=== the same two at -j 1, to separate parse from compression ===")
    c = arm("gz6_j1", 6, 1)
    d = arm("gz1_j1", 1, 1)

    for name, arm_res in (("gz6_j8", a), ("gz1_j8", b), ("gz6_j1", c), ("gz1_j1", d)):
        if arm_res is None:
            raise SystemExit(f"arm {name} failed; see the logs")
    # Narrow the Optionals for the reader and the type checker.
    assert a is not None and b is not None and c is not None and d is not None

    print("\n=== is the output standard gzip? ===")
    for lab in ("gz6_j8", "gz1_j8"):
        p = os.path.join(OUT, lab, "o1.fq.gz")
        print(f"  {lab:8} gzip -t: {'OK (standard)' if gzip_ok(p) else 'FAILED (non-standard format)'}")

    print("\n=== verdict ===")
    if b["wall"] > a["wall"]:
        print(f"  FLIPPING IS A NET LOSS: level 1 is {b['wall']/a['wall']:.2f}x SLOWER "
              f"at -j 8 ({b['wall']:.2f}s vs {a['wall']:.2f}s).")
        print(f"  busy cores {a['busy']:.2f} -> {b['busy']:.2f}: the parallel compressor "
              f"stops engaging below level 6.")
    else:
        print(f"  level 1 is {a['wall']/b['wall']:.2f}x faster at -j 8; the flip is viable "
              f"if the byte growth is acceptable.")
    print(f"  bytes: {a['size']/1e6:.0f} MB -> {b['size']/1e6:.0f} MB "
          f"({b['size']/a['size']:.2f}x)")
    print(f"  peak RSS: {a['peak']:.0f} MB -> {b['peak']:.0f} MB")


if __name__ == "__main__":
    main()
