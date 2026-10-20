"""One benchmark cell = one tool invocation, measured under taskset.

Why this script exists rather than a `shell:` line: a bench row must be
impossible to fake. This one refuses to emit a row unless

  * the child exited 0,
  * its output exists and is non-empty,
  * stderr contains nothing except the known citation line,

and it measures the child's own peak RSS with `os.wait4` (which returns the
rusage OF THAT CHILD) rather than `getrusage(RUSAGE_CHILDREN)` — the latter is a
never-resetting high-water mark across every reaped child, a defect that once
made three different tools report the same 1209.8 MB.

Invoked once per cell by the Snakefile; prints exactly one TSV row on stdout.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import time


def die(msg):
    print(f"measure: FAILED: {msg}", file=sys.stderr)
    raise SystemExit(1)


def tool_argv(args):
    """Command line for one measured run. Every tool gets the same discipline:
    explicit thread count, explicit gzip level where it has one, --force."""
    reads = args.reads
    r1 = os.path.join(args.reads_dir, f"{reads}_R1.fastq.gz")
    r2 = os.path.join(args.reads_dir, f"{reads}_R2.fastq.gz")
    out = args.out_dir
    t = str(args.threads)

    if args.tool == "fastprep":
        return [args.bin, "-1", r1, "-2", r2, "-j", t,
                "-o", os.path.join(out, "trim_R1.fastq.gz"),
                "-p", os.path.join(out, "trim_R2.fastq.gz"),
                "--json", os.path.join(out, "fastprep.json"),
                "--gz-level", "6", "--force", "-q"]
    if args.tool == "fastcheck":
        # fastcheck takes positional FILEs (single file or mate pair), not -1/-2.
        return [args.bin, r1, r2, "-j", t,
                "--json", os.path.join(out, "fastcheck.json"), "-q"]
    if args.tool == "fastmaps":
        return [args.bin, r1, r2, "-j", t, "-o", out,
                "--json", os.path.join(out, "heatmap.json"), "--force", "-q"]
    if args.tool == "fastdedup":
        return [args.bin, "-1", r1, "-2", r2, "-j", t,
                "-o", os.path.join(out, "dedup_R1.fastq.gz"),
                "-p", os.path.join(out, "dedup_R2.fastq.gz"),
                "--gz-level", "6", "--json", os.path.join(out, "fastdedup.json"),
                "--force", "-q"]
    if args.tool == "fastsim":
        # fastsim is single-threaded by construction; --threads is accepted and
        # ignored, so the threads axis is pinned to 1 for it in the config.
        argv = [args.bin, "--pairs", str(reads),
                "--profile", args.profile, "--seed", "42",
                "-o", os.path.join(out, "sim"), "--force", "-q"]
        if args.gz_level:
            # Fail rather than silently measure the default: a cell labelled
            # gz1 that really ran gz6 would read as "the change did nothing".
            if not supports_gz_level(args.bin):
                die(f"{args.bin} has no --gz-level; the gz axis is meaningless for it")
            argv += ["--gz-level", str(args.gz_level)]
        return argv
    if args.tool == "fasthelper":
        # Aggregation over a pre-built report tree; `reads` carries the sample
        # count, which is the axis that matters for this tool.
        return [args.bin, "concatenate", args.tree, "-o", out, "--format", "json"]
    die(f"unknown tool {args.tool}")


_GZ_LEVEL_CACHE = {}


def supports_gz_level(bin_path):
    """`--gz-level` was added to fastsim in this round; a baseline binary does
    not have it, and passing it would abort the cell."""
    if bin_path not in _GZ_LEVEL_CACHE:
        help_txt = subprocess.run([bin_path, "--help"], capture_output=True, text=True).stdout
        _GZ_LEVEL_CACHE[bin_path] = "--gz-level" in help_txt
    return _GZ_LEVEL_CACHE[bin_path]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--tool", required=True)
    ap.add_argument("--reads", type=int, required=True)
    ap.add_argument("--threads", type=int, required=True)
    ap.add_argument("--gz-level", type=int, default=0)
    ap.add_argument("--profile", default="")
    ap.add_argument("--reads-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--kind", required=True)
    ap.add_argument("--tree", default="")
    args = ap.parse_args()

    if not os.path.exists(args.bin):
        die(f"binary not found: {args.bin}")
    os.makedirs(args.out_dir, exist_ok=True)

    # taskset pins the child to exactly `threads` cpus, so `-j N` cannot win a
    # wall-clock cell by spending cores it was not given.
    cpus = ",".join(str(i) for i in range(args.threads))
    argv = ["taskset", "-c", cpus] + tool_argv(args)

    # Do NOT use communicate()/wait(): they reap the child, and the child's own
    # rusage (ru_maxrss, cpu time) is only available from wait4 on the pid that
    # is still a child of this process.
    stdout_path = os.path.join(args.out_dir, ".cell_stdout")
    stderr_path = os.path.join(args.out_dir, ".cell_stderr")
    peak_kb = 0
    with open(stdout_path, "wb") as so, open(stderr_path, "wb") as se:
        t0 = time.monotonic()
        proc = subprocess.Popen(argv, stdout=so, stderr=se)
        # Poll the child's own VmHWM while it lives: the only per-run peak POSIX
        # offers, and the reason the memory column of a bench is trustworthy.
        while True:
            try:
                with open(f"/proc/{proc.pid}/status") as f:
                    for line in f:
                        if line.startswith("VmHWM:"):
                            peak_kb = max(peak_kb, int(line.split()[1]))
                            break
            except FileNotFoundError:
                pass
            pid, status, ru = os.wait4(proc.pid, os.WNOHANG)
            if pid:
                break
            time.sleep(0.02)
        wall = time.monotonic() - t0

    rc = os.waitstatus_to_exitcode(status)
    with open(stderr_path, "rb") as f:
        err = f.read()
    for tmp in (stdout_path, stderr_path):
        try:
            os.unlink(tmp)
        except OSError:
            pass

    # Benign progress noise, matched explicitly. Everything else on stderr is a
    # failure: that rule is what caught a 0.02 s / 0.0 MB phantom row when a tool
    # could not create its output directory.
    benign = re.compile(
        r"^(If you use|\[\d+/\d+\] |Classified files: "
        r"|(fastcheck|fastdedup|fastmaps|fastprep|fastsim|fasthelper): )"
    )
    fatal = re.compile(r"(?i)\b(error|panicked|failed to|cannot|no such file)\b")
    stderr_lines = [l for l in err.decode("utf-8", "replace").splitlines() if l.strip()]
    unexpected = [l for l in stderr_lines if fatal.search(l) or not benign.match(l)]
    if rc != 0:
        die(f"{args.tool} exited {rc} on {args.reads} reads at -j{args.threads}: "
            f"{' | '.join(unexpected[-3:]) or 'no stderr'}")
    if unexpected:
        die(f"{args.tool} wrote unexpected stderr: {unexpected[:3]}")

    produced = []
    for root, _dirs, files in os.walk(args.out_dir):
        for f in files:
            if f.startswith(".cell_"):
                continue
            produced.append(os.path.join(root, f))
    if not produced:
        die(f"{args.tool} produced no output files in {args.out_dir}")
    empty = [p for p in produced if os.path.getsize(p) == 0]
    if empty:
        die(f"{args.tool} produced empty output: {empty[:3]}")
    out_bytes = sum(os.path.getsize(p) for p in produced)

    cpu_s = ru.ru_utime + ru.ru_stime
    cols = ["tool", "kind", "reads", "threads", "gz_level", "wall_s", "cpu_s",
            "cpu_over_wall", "peak_rss_kb", "out_bytes", "files", "rc"]
    row = [args.tool, args.kind, str(args.reads), str(args.threads),
           str(args.gz_level or ""), f"{wall:.3f}", f"{cpu_s:.3f}",
           f"{(cpu_s / wall if wall > 0 else 0):.3f}", str(peak_kb),
           str(out_bytes), str(len(produced)), str(rc)]
    print("\t".join(cols))
    print("\t".join(row))


if __name__ == "__main__":
    main()
