#!/usr/bin/env python3
"""Task 2.8: what the metadata writer costs, and what the codec fix buys.

The defect this measures around: `open_output(tp, ...)` inferred the codec from
the ".part" TEMP path, so `--metadata x.pxmt.gz` was written UNCOMPRESSED (the
temp path's extension is ".part", which falls through to plain). Fixed by
inferring from the final path via `open_output_as`, which is what the FASTQ
writer already did.

Reports wall, CPU and peak RSS for a run with and without `--metadata`, and the
metadata file size. Peak RSS is the reaped child's ru_maxrss from wait4.
"""
import os
import statistics
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FP = os.path.join(REPO, "tools", "fastprep", "target", "release", "fastprep")
D = "/tmp/meta-cost"
N = 200_000


def make_corpus():
    os.makedirs(D, exist_ok=True)
    p = os.path.join(D, "in.fq")
    if os.path.exists(p) and os.path.getsize(p) > 0:
        return p
    # Deterministic, stdlib only: a xorshift, no random module state.
    x = 12345
    with open(p, "w") as fh:
        for i in range(N):
            seq = []
            for _ in range(150):
                x ^= (x << 13) & 0xFFFFFFFFFFFFFFFF
                x ^= x >> 7
                x ^= (x << 17) & 0xFFFFFFFFFFFFFFFF
                seq.append("ACGT"[(x >> 3) & 3])
            q = "".join(chr(33 + 20 + ((x >> (j % 17)) & 7)) for j in range(150))
            fh.write(f"@read{i} extra\n{''.join(seq)}\n+\n{q}\n")
    return p


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


def arm(label, meta, reps=5):
    walls, cpus, peaks, sizes = [], [], [], []
    for _ in range(reps):
        out = os.path.join(D, label)
        subprocess.run(["rm", "-rf", out], check=False)
        os.makedirs(out, exist_ok=True)
        cmd = [FP, "-1", os.path.join(D, "in.fq"), "-o", os.path.join(out, "o.fq.gz"),
               "-j", "8", "-q", "--force"]
        if meta:
            cmd += ["--metadata", os.path.join(out, "m.pxmt.gz")]
        wall, peak, rc, cpu = run(cmd, os.path.join(D, f"{label}.log"))
        if rc != 0:
            print(f"{label}: exit {rc}", file=sys.stderr)
            print(open(os.path.join(D, f"{label}.log")).read()[-800:], file=sys.stderr)
            raise SystemExit(1)
        walls.append(wall); cpus.append(cpu); peaks.append(peak)
        mp = os.path.join(out, "m.pxmt.gz")
        if os.path.exists(mp):
            sizes.append(os.path.getsize(mp))
    res = dict(label=label, wall=statistics.median(walls), cpu=statistics.median(cpus),
               peak=max(peaks), size=max(sizes) if sizes else 0)
    print(f"{label:12} wall {res['wall']:6.3f}s  cpu {res['cpu']:6.3f}s  "
          f"peak {res['peak']:7.1f} MB  metadata {res['size']:>9,} B")
    return res


def main():
    if not os.path.exists(FP):
        raise SystemExit(f"build fastprep first: {FP}")
    make_corpus()
    print(f"corpus: {N:,} reads, {os.path.getsize(os.path.join(D,'in.fq'))/1e6:.1f} MB")
    a = arm("nometa", False)
    b = arm("withmeta", True)
    dw = b["wall"] - a["wall"]
    dc = b["cpu"] - a["cpu"]
    print()
    print(f"metadata cost: +{dw:.3f} s wall, +{dc:.3f} s cpu "
          f"({b['size']/N:.2f} bytes/read on disk)")
    if b["size"]:
        print(f"codec: {'gzip (fixed)' if b['size'] < N * 30 else 'UNCOMPRESSED (defect present)'}")


if __name__ == "__main__":
    main()
