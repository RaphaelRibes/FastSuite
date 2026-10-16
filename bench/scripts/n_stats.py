#!/usr/bin/env python3
"""Measure the ambiguous-base (N) structure of real FASTQ libraries.

The closed loop's duplicate-rate gate fails because fastsim emits **no** N at
all: `random_dna_bias` draws only ACGT, so padding-driven duplicates -- reads
whose identity comes from a poly-N run rather than from biology -- have no
counterpart in the simulation. This script is the measurement that a profile
must reproduce, so the emission is calibrated from data rather than guessed.

It is deliberately independent of fastcheck: fastcheck counts *that* a read
carries a non-ACGT base (`metrics.reads_with_non_acgt`) but says nothing about
run length, run position, or how the N mass concentrates in duplicate reads,
and those are the three things the emission model needs.

Per input FASTQ (one mate) it reports:

  reads               reads scanned
  reads_with_n        reads carrying >= 1 N
  reads_with_n_frac   reads_with_n / reads
  n_bases             total N bases
  n_base_frac         n_bases / total bases
  runs                number of maximal N-runs
  runs_per_n_read     mean runs among reads that carry one
  run_len_*           the run-length histogram (as `<len>:count` in the TSV's
                      `run_len_hist` column, so the whole distribution travels)
  run_len_p50/p90/max
  pos_leading_frac    runs touching position 0
  pos_trailing_frac   runs touching the last base
  pos_interior_frac   every other run
  dup_reads           reads whose exact sequence occurs more than once
  dup_reads_with_n    of those, how many carry an N
  dup_reads_with_n_frac
  dup_runs_frac       share of ALL N-runs that sit in duplicate reads
  approval            share of reads whose N mass is a run of >= 10 (padding)

Usage:
  n_stats.py --out n_stats.tsv SAMPLE_R1.fastq.gz [SAMPLE_R1b.fastq.gz ...]
  n_stats.py --out n_stats.tsv --label human_gut/DRR425489 -- \
      DRR425489.R1.fastq.gz DRR425489.R2.fastq.gz

Each positional argument is ONE mate of ONE library; a library measured on both
mates contributes two rows, which is what the profile (per-mate CDFs) needs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import Counter

ROW_FIELDS = [
    "sample", "label", "reads", "reads_with_n", "reads_with_n_frac",
    "n_bases", "n_base_frac", "runs", "runs_per_n_read",
    "run_len_p50", "run_len_p90", "run_len_max", "run_len_hist",
    "pos_leading_frac", "pos_trailing_frac", "pos_interior_frac",
    "dup_reads", "dup_reads_with_n", "dup_reads_with_n_frac", "dup_runs_frac",
    "padding_frac",
]

PADDING_MIN_RUN = 10  # a run this long is instrumentation padding, not a call


def open_fastq(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def n_runs(seq: bytes) -> list[tuple[int, int]]:
    """Maximal N-runs as (start, length), in read coordinates."""
    out = []
    i, n = 0, len(seq)
    while i < n:
        if seq[i] in (ord("N"), ord("n")):
            j = i
            while j < n and seq[j] in (ord("N"), ord("n")):
                j += 1
            out.append((i, j - i))
            i = j
        else:
            i += 1
    return out


def percentile(sorted_vals: list[int], q: float) -> float:
    if not sorted_vals:
        return 0.0
    idx = min(len(sorted_vals) - 1, int(q * len(sorted_vals)))
    return float(sorted_vals[idx])


def measure(path: str, sample: str, label: str) -> dict:
    seqs: list[bytes] = []
    reads = 0
    reads_with_n = 0
    n_bases = 0
    total_bases = 0
    run_lens: list[int] = []
    leading = trailing = interior = 0
    per_read_runs: list[int] = []
    n_read_run_indices: list[tuple[int, int, int]] = []  # (read_idx, start, len)
    padding_reads = 0

    with open_fastq(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:  # sequence line of each 4-line record
                continue
            seq = line.strip()
            seqs.append(seq)
            reads += 1
            total_bases += len(seq)
            runs = n_runs(seq)
            if not runs:
                continue
            reads_with_n += 1
            per_read_runs.append(len(runs))
            has_padding = False
            for (start, ln) in runs:
                n_bases += ln
                run_lens.append(ln)
                run_ends = start + ln >= len(seq)
                if start == 0:
                    leading += 1
                elif run_ends:
                    trailing += 1
                else:
                    interior += 1
                if ln >= PADDING_MIN_RUN:
                    has_padding = True
                n_read_run_indices.append((reads - 1, start, ln))
            if has_padding:
                padding_reads += 1

    # Duplicate concentration: exact-sequence multiplicity over this one mate.
    counts = Counter(seqs)
    dup_reads = sum(c for c in counts.values() if c > 1)
    dup_reads_with_n = 0
    dup_run_count = 0
    dup_read_idx = set()
    idx = 0
    for seq in seqs:
        if counts[seq] > 1:
            dup_read_idx.add(idx)
            if n_runs(seq):
                dup_reads_with_n += 1
        idx += 1
    for (read_idx, _start, _ln) in n_read_run_indices:
        if read_idx in dup_read_idx:
            dup_run_count += 1

    total_runs = len(run_lens)
    run_lens_sorted = sorted(run_lens)
    f = lambda a, b: (a / b) if b else 0.0
    return {
        "sample": sample,
        "label": label,
        "reads": reads,
        "reads_with_n": reads_with_n,
        "reads_with_n_frac": round(f(reads_with_n, reads), 6),
        "n_bases": n_bases,
        "n_base_frac": round(f(n_bases, total_bases), 6),
        "runs": total_runs,
        "runs_per_n_read": round(f(total_runs, reads_with_n), 4),
        "run_len_p50": percentile(run_lens_sorted, 0.50),
        "run_len_p90": percentile(run_lens_sorted, 0.90),
        "run_len_max": run_lens_sorted[-1] if run_lens_sorted else 0,
        "run_len_hist": ";".join(f"{k}:{v}" for k, v in sorted(Counter(run_lens).items())),
        "pos_leading_frac": round(f(leading, total_runs), 6),
        "pos_trailing_frac": round(f(trailing, total_runs), 6),
        "pos_interior_frac": round(f(interior, total_runs), 6),
        "dup_reads": dup_reads,
        "dup_reads_with_n": dup_reads_with_n,
        "dup_reads_with_n_frac": round(f(dup_reads_with_n, dup_reads), 6),
        "dup_runs_frac": round(f(dup_run_count, total_runs), 6),
        "padding_frac": round(f(padding_reads, reads_with_n), 6),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="one FASTQ(.gz) per measurement")
    ap.add_argument("--out", required=True, help="output TSV")
    ap.add_argument("--label", default="", help="library label (domain/accession)")
    ap.add_argument("--force", action="store_true", help="overwrite --out")
    args = ap.parse_args()

    import os
    if os.path.exists(args.out) and not args.force:
        print(f"n_stats: {args.out} exists; pass --force", file=sys.stderr)
        return 2

    rows = []
    for path in args.inputs:
        sample = os.path.basename(path)
        rows.append(measure(path, sample, args.label))
        r = rows[-1]
        print(f"{sample}: {r['reads_with_n_frac']:.4%} of reads carry N "
              f"({r['n_bases']} N bases), dup-N {r['dup_reads_with_n_frac']:.4f}",
              file=sys.stderr)

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ROW_FIELDS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"n_stats: {len(rows)} row(s) -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
