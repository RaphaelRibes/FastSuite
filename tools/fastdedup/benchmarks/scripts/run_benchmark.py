#!/usr/bin/env python3
"""
Run a single fdedup benchmark and capture performance metrics.

Metrics captured via resource.getrusage():
  - User CPU time
  - System CPU time
  - Peak RSS (KB)
  - Filesystem I/O blocks in/out

Also captures:
  - Wall clock time (time.monotonic)
  - Read counts from --verbose output

Output: single TSV line written to the output file.
"""

import argparse
import os
import re
import resource
import subprocess
import sys
import time


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run one fdedup benchmark and emit metrics."
    )
    parser.add_argument("--fdedup", required=True, help="Path to fdedup binary")
    parser.add_argument("--dataset", required=True, help="Dataset label (e.g. 100K)")
    parser.add_argument("--io-mode", required=True, help="I/O mode (u2u/c2u/c2c)")
    parser.add_argument("--se-pe", required=True, choices=["SE", "PE"])
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--variant", required=True,
                        choices=["normal", "abundance", "abundance_tempdir"])
    parser.add_argument("--input-r1", required=True)
    parser.add_argument("--input-r2", default=None)
    parser.add_argument("--output-r1", required=True)
    parser.add_argument("--output-r2", default=None)
    parser.add_argument("--output-file", required=True,
                        help="Path to write the TSV result line")
    return parser.parse_args()


def run_fdedup(args: argparse.Namespace) -> dict:
    """Run fdedup and return metrics dict."""
    cmd = [args.fdedup]

    # Input
    cmd.extend(["-1", args.input_r1])

    if args.se_pe == "PE":
        assert args.input_r2, "PE mode requires --input-r2"
        cmd.extend(["-2", args.input_r2])

    # Output
    cmd.extend(["-o", args.output_r1])
    if args.se_pe == "PE":
        assert args.output_r2, "PE mode requires --output-r2"
        cmd.extend(["-p", args.output_r2])

    # Force overwrite (clean slate each run)
    cmd.append("--force")

    # Verbose for read count parsing
    cmd.append("--verbose")

    # Threads
    cmd.extend(["-j", str(args.threads)])

    # Variant
    if args.variant == "abundance":
        cmd.append("--abundance")
    elif args.variant == "abundance_tempdir":
        cmd.append("--abundance")
        cmd.extend(["--temp-dir", "/tmp"])

    # Wall clock start
    t0 = time.monotonic()

    # Run process
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )

    t1 = time.monotonic()
    wall_sec = round(t1 - t0, 3)

    # Get resource usage of child process
    rusage = resource.getrusage(resource.RUSAGE_CHILDREN)

    # Parse verbose output for read counts
    # Typical output:
    #   Processed fragments: 1000
    #   Duplicates removed: 15.00%
    #   Total execution time: 0.05s
    stdout = proc.stdout
    reads_processed = 0
    duplicates_pct = 0.0
    for line in stdout.splitlines():
        m = re.match(r"Processed fragments:\s+(\d+)", line)
        if m:
            reads_processed = int(m.group(1))
            continue
        m = re.match(r"Duplicates removed:\s+([0-9.]+)%", line)
        if m:
            duplicates_pct = float(m.group(1))
            continue

    return {
        "wall_sec": wall_sec,
        "cpu_user_sec": round(rusage.ru_utime, 3),
        "cpu_sys_sec": round(rusage.ru_stime, 3),
        "max_rss_kb": rusage.ru_maxrss,
        "io_in_blocks": rusage.ru_inblock,
        "io_out_blocks": rusage.ru_oublock,
        "reads_processed": reads_processed,
        "duplicates_pct": round(duplicates_pct, 2),
        "returncode": proc.returncode,
    }


def main():
    args = parse_args()

    # Ensure output dir exists
    os.makedirs(os.path.dirname(args.output_file), exist_ok=True)

    # Run
    metrics = run_fdedup(args)

    # Build TSV line
    fields = [
        args.dataset,
        args.io_mode,
        args.se_pe,
        str(args.threads),
        args.variant,
        str(metrics["wall_sec"]),
        str(metrics["cpu_user_sec"]),
        str(metrics["cpu_sys_sec"]),
        str(metrics["max_rss_kb"]),
        str(metrics["io_in_blocks"]),
        str(metrics["io_out_blocks"]),
        str(metrics["reads_processed"]),
        str(metrics["duplicates_pct"]),
        str(metrics["returncode"]),
    ]

    with open(args.output_file, "w") as f:
        f.write("\t".join(fields) + "\n")

    print(f"Wrote {args.output_file}")


if __name__ == "__main__":
    main()
