#!/usr/bin/env python3
"""Task 2.12 Step 1: where does fastcheck's 785 MB peak actually live?

The tool reports its own accounting ("fingerprint table projected at 272.0 MB",
"accumulator state 64.3 MB"), which sums to 336 MB -- less than half the
measured 785 MB peak. So the plan's premise ("a single 32 MB allocation
dominating VmHWM", "start with the fingerprint table's pre-allocation") does
not fit the numbers.

`--only MODULE,...` materialises only what the enabled modules need, so
running each module in isolation localises the footprint instead of guessing.

Peak RSS is the reaped child's ru_maxrss from wait4.
"""
import json
import os
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FC = os.path.join(REPO, "tools", "fastcheck", "target", "release", "fastcheck")
CORPUS = "/mnt/Fast/fastsuite-scratch/fm-ab5m/prep"
OUT = "/mnt/Fast/fastsuite-scratch/fc-only"


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
    return (time.perf_counter() - t0, ru.ru_maxrss / 1024.0,
            os.waitstatus_to_exitcode(status), ru.ru_utime + ru.ru_stime)


def arm(label, only=None):
    r1 = os.path.join(CORPUS, "clean_R1.fastq.gz")
    r2 = os.path.join(CORPUS, "clean_R2.fastq.gz")
    os.makedirs(OUT, exist_ok=True)
    log = os.path.join(OUT, f"{label}.log")
    cmd = [FC, r1, r2, "-j", "8", "-q"]
    if only:
        cmd += ["--only", only]
    wall, peak, rc, cpu = run(cmd, log)
    if rc != 0:
        print(f"{label}: exit {rc}", file=sys.stderr)
        print(open(log).read()[-500:], file=sys.stderr)
        return None
    # the tool's own accounting, from its stderr
    acc = ""
    for line in open(log):
        if "projected at" in line:
            acc = line.strip().split("projected at")[1].split(" if")[0].strip()
    print(f"{label:34} peak {peak:7.1f} MB  cpu {cpu:6.1f}s  wall {wall:5.1f}s"
          + (f"   [table {acc}]" if acc else ""))
    return dict(label=label, peak=peak, cpu=cpu)


def main():
    if not os.path.exists(FC):
        raise SystemExit(f"build fastcheck first: {FC}")
    if not os.path.isdir(CORPUS):
        raise SystemExit(f"corpus missing: {CORPUS}")
    print("5M pairs, -j 8\n")
    print("=== each module alone ===")
    arms = {}
    for m in ("dup", "bases", "quality", "length", "gc", "heatmap", "polyx", "kmer"):
        arms[m] = arm(f"only {m}", m)

    print("\n=== everything (the default) ===")
    everything = arm("all modules")
    print("\n=== the plan's suggested starting point, unchecked ===")
    nodup = arm("all except dup", "bases,quality,length,gc,heatmap,polyx,kmer")

    if everything is None:
        raise SystemExit("the full run failed")

    print("\n=== attribution ===")
    ranked = sorted(
        ((m, r["peak"]) for m, r in arms.items() if r),
        key=lambda kv: -kv[1])
    for m, p in ranked:
        print(f"  {m:10} {p:7.1f} MB")
    print()
    if nodup:
        drop = everything["peak"] - nodup["peak"]
        print(f"  full {everything['peak']:.1f} MB - without dup {nodup['peak']:.1f} MB "
              f"= {drop:.1f} MB attributable to dup")
    biggest = ranked[0]
    print(f"  largest single module: {biggest[0]} at {biggest[1]:.1f} MB "
          f"({biggest[1]/everything['peak']*100:.0f}% of the full peak)")
    other = everything["peak"] - sum(p for _, p in ranked[:1])
    print(f"  remaining {other:.1f} MB is the shared read path + accumulators, not one module")


if __name__ == "__main__":
    main()
