#!/usr/bin/env python3
"""Independent reference implementation of fastcheck's core statistics.

Not a copy of the Rust code: this walks records with a different parser, uses
Python integers, and derives the numbers the obvious way. The point is to catch
a systematic mistake in the Rust accumulators (a wrong denominator, a dropped
position, an off-by-one in a histogram), which unit tests written by the same
hand tend to miss.

Usage:
    reference_qc.py FILE [--max-pos N] [--offset 33|64]
Prints one JSON object with the same metric names as
`fastsuite.report.v1`'s `metrics` block plus the per-position tables, so a
caller can diff the two.
"""

import argparse
import gzip
import json
import sys
from collections import Counter

QBINS = 256


def open_maybe_gz(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def quantile_nearest_rank(hist: Counter, q: float) -> float:
    """Smallest value whose cumulative count reaches ceil(q * n)."""
    n = sum(hist.values())
    if n == 0:
        return 0.0
    rank = max(1, int(-(-(q * n) // 1)))  # ceil(q*n)
    seen = 0
    for value in sorted(hist):
        seen += hist[value]
        if seen >= rank:
            return float(value)
    return float(max(hist))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--max-pos", type=int, default=1000)
    ap.add_argument("--offset", type=int, default=33)
    args = ap.parse_args()

    reads = 0
    bases = 0
    reads_with_non_acgt = 0
    len_hist = Counter()
    acgt_hist = Counter()
    gc_hist = Counter()
    qual_byte = Counter()
    meanq_raw = Counter()
    base_pos = {}  # pos -> Counter over b"A",b"C",b"G",b"T",b"other"
    qual_pos = {}  # pos -> Counter over raw byte

    with open_maybe_gz(args.path) as fh:
        header = b""
        while True:
            header = fh.readline()
            if not header:
                break
            if header.strip() == b"":
                continue
            if not header.startswith(b"@"):
                raise SystemExit(f"not FASTQ at {header[:40]!r}")
            seq = fh.readline().rstrip(b"\r\n")
            plus = fh.readline()
            qual = fh.readline().rstrip(b"\r\n")
            if not plus.startswith(b"+"):
                raise SystemExit("missing + line")
            if len(qual) != len(seq):
                raise SystemExit("wrapped FASTQ not supported")

            reads += 1
            bases += len(seq)
            len_hist[len(seq)] += 1

            acgt = 0
            gc = 0
            has_non = False
            for i, raw in enumerate(seq):
                b = bytes([raw]).upper()
                if b == b"A":
                    code = 0
                elif b == b"C":
                    code = 1
                elif b == b"G":
                    code = 2
                elif b == b"T":
                    code = 3
                else:
                    code = 4
                if code == 4:
                    has_non = True
                else:
                    acgt += 1
                    if code in (1, 2):
                        gc += 1
                if i < args.max_pos:
                    slot = base_pos.setdefault(i, Counter())
                    slot[code] += 1
            if has_non:
                reads_with_non_acgt += 1
            acgt_hist[acgt] += 1
            gc_hist[gc] += 1

            qsum = 0
            for i, q in enumerate(qual):
                qual_byte[q] += 1
                qsum += q
                if i < args.max_pos:
                    qual_pos.setdefault(i, Counter())[q] += 1
            meanq_raw[qsum // len(qual) if qual else 0] += 1

    def hist_to_json(h):
        return {str(k): v for k, v in sorted(h.items()) if v}

    positions = []
    for i in range(max(base_pos) + 1 if base_pos else 0):
        bp = base_pos.get(i, Counter())
        qp = qual_pos.get(i, Counter())
        depth = sum(qp.values())
        if depth:
            shifted = Counter()
            total = 0.0
            for byte, c in qp.items():
                phred = max(0, byte - args.offset)
                shifted[phred] += c
                total += phred * c
            mean_q = total / depth
            q10 = quantile_nearest_rank(shifted, 0.10)
            q50 = quantile_nearest_rank(shifted, 0.50)
            q90 = quantile_nearest_rank(shifted, 0.90)
        else:
            mean_q = q10 = q50 = q90 = 0.0
        positions.append(
            {
                "pos": i + 1,
                "a": bp.get(0, 0),
                "c": bp.get(1, 0),
                "g": bp.get(2, 0),
                "t": bp.get(3, 0),
                "n": bp.get(4, 0),
                "depth": depth,
                "mean_q": mean_q,
                "q10": q10,
                "q50": q50,
                "q90": q90,
            }
        )

    acgt_total = sum(k * v for k, v in acgt_hist.items())
    gc_total = sum(k * v for k, v in gc_hist.items())
    mean_q_reads = (
        sum(max(0, b - args.offset) * c for b, c in meanq_raw.items()) / reads
        if reads
        else 0.0
    )
    q20 = sum(c for b, c in meanq_raw.items() if b - args.offset >= 20) / reads if reads else 0.0
    q30 = sum(c for b, c in meanq_raw.items() if b - args.offset >= 30) / reads if reads else 0.0

    out = {
        "metrics": {
            "reads": reads,
            "bases": bases,
            "bases_acgt": acgt_total,
            "bases_non_acgt": bases - acgt_total,
            "reads_with_non_acgt": reads_with_non_acgt,
            "min_len": float(min(len_hist)) if len_hist else 0.0,
            "max_len": float(max(len_hist)) if len_hist else 0.0,
            "median_len": quantile_nearest_rank(len_hist, 0.5),
            "mean_len": (bases / reads) if reads else 0.0,
            "gc_fraction": (gc_total / acgt_total) if acgt_total else 0.0,
            "mean_read_mean_q": mean_q_reads,
            "reads_q20_frac": q20,
            "reads_q30_frac": q30,
            "quality_byte_min": float(min(qual_byte)) if qual_byte else 0.0,
            "quality_byte_max": float(max(qual_byte)) if qual_byte else 0.0,
            "partial": 0,
        },
        "read_length_hist": hist_to_json(len_hist),
        "read_mean_q_hist": {
            str(max(0, b - args.offset)): c for b, c in sorted(meanq_raw.items()) if c
        },
        "positions": positions,
    }
    json.dump(out, sys.stdout)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
