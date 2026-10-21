#!/usr/bin/env python3
"""fastmaps-only A/B harness (Task 2.14, 2.13, 2.15, 2.16).

The full `profile_pipeline.py` chain spends ~280 s in fastsim to produce an
input that fastmaps then reads in ~4 s. For a fastmaps-only comparison that
ratio is absurd, so this harness generates the corpus ONCE and then measures
the binary under test against it, many times.

Peak RSS is the reaped child's ru_maxrss from wait4 (not a /proc poll, which
under-reports a child that finishes between polls).

Usage:
    fastmaps_ab.py prepare <dir> <pairs>     # generate the corpus
    fastmaps_ab.py measure <dir> <label> <binary> [reps]
"""
import json
import os
import statistics
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TOOLS = os.path.join(REPO, "tools")
PROFILE = f"{TOOLS}/fastsim/profiles/cellB-pe-overlap.json"


def run(cmd, logpath):
    """Run cmd, return (wall_s, peak_mb, exit_code, cpu_s) for THAT child."""
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


def prepare(d, pairs):
    os.makedirs(d, exist_ok=True)
    sim = os.path.join(d, "sim")
    os.makedirs(sim, exist_ok=True)
    print(f"simulating {pairs} pairs -> {sim} (this is the slow part)")
    wall, peak, rc, cpu = run(
        [f"{TOOLS}/fastsim/target/release/fastsim", "--profile", PROFILE,
         "--pairs", str(pairs), "--seed", "1", "--out", f"{sim}/sim"],
        os.path.join(d, "sim.log"))
    print(f"  fastsim: {wall:.1f} s wall, {cpu:.1f} s cpu, {peak:.0f} MB, exit {rc}")
    if rc != 0:
        sys.exit(rc)
    # Trim the simulated reads the way the pipeline does, so fastmaps sees
    # the same input it sees in the real chain (adapter-free, deduplicated).
    prep = os.path.join(d, "prep")
    os.makedirs(prep, exist_ok=True)
    for stage, cmd in [
        ("dedup", [f"{TOOLS}/fastdedup/target/release/fastdedup",
                   "-1", f"{sim}/sim_R1.fastq.gz", "-2", f"{sim}/sim_R2.fastq.gz",
                   "-o", f"{prep}/dedup_R1.fastq.gz", "-p", f"{prep}/dedup_R2.fastq.gz",
                   "-j", "8", "--gz-level", "1", "--json", f"{prep}/fastdedup.json",
                   "--force", "-q"]),
        ("prep", [f"{TOOLS}/fastprep/target/release/fastprep",
                  "-1", f"{prep}/dedup_R1.fastq.gz", "-2", f"{prep}/dedup_R2.fastq.gz",
                  "-o", f"{prep}/clean_R1.fastq.gz", "-p", f"{prep}/clean_R2.fastq.gz",
                  "-j", "8", "--gz-level", "1", "--json", f"{prep}/fastprep.json",
                  "--force", "-q"]),
    ]:
        wall, peak, rc, cpu = run(cmd, os.path.join(d, f"{stage}.log"))
        print(f"  {stage}: {wall:.1f} s wall, {peak:.0f} MB, exit {rc}")
        if rc != 0:
            sys.exit(rc)
    for f in ("clean_R1.fastq.gz", "clean_R2.fastq.gz"):
        print(f"  {f}: {os.path.getsize(os.path.join(prep, f)) / 1e6:.1f} MB")


def measure(d, label, binary, reps):
    prep = os.path.join(d, "prep")
    out = os.path.join(d, f"maps_{label}")
    os.makedirs(out, exist_ok=True)
    walls, cpus, peaks = [], [], []
    json_path = os.path.join(out, "fastmaps.json")
    for i in range(reps):
        cmd = [binary, f"{prep}/clean_R1.fastq.gz", f"{prep}/clean_R2.fastq.gz",
               "-o", f"{out}/maps", "-j", "8", "--json", json_path, "--force", "-q"]
        wall, peak, rc, cpu = run(cmd, os.path.join(d, f"maps_{label}.log"))
        if rc != 0:
            print(f"  fastmaps {label} rep {i}: exit {rc}, see {d}/maps_{label}.log")
            sys.exit(rc)
        walls.append(wall)
        cpus.append(cpu)
        peaks.append(peak)
    res = {
        "label": label,
        "binary": binary,
        "reps": reps,
        "wall_s_min": round(min(walls), 3),
        "wall_s_median": round(statistics.median(walls), 3),
        "wall_s_mean": round(statistics.mean(walls), 3),
        "cpu_s_min": round(min(cpus), 3),
        "cpu_s_median": round(statistics.median(cpus), 3),
        "peak_mb": round(max(peaks), 1),
    }
    with open(os.path.join(d, f"result_{label}.json"), "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"{label:10} wall min {res['wall_s_min']:7.3f} s  median {res['wall_s_median']:7.3f} s"
          f"  cpu min {res['cpu_s_min']:7.3f} s  peak {res['peak_mb']:6.1f} MB")
    return res


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "prepare":
        prepare(sys.argv[2], int(sys.argv[3]))
    elif what == "measure":
        measure(sys.argv[2], sys.argv[3], sys.argv[4],
                int(sys.argv[5]) if len(sys.argv) > 5 else 5)
