#!/usr/bin/env python3
"""
Generate synthetic FASTQ reads for FastDedup benchmarking.

Produces R1 and R2 FASTQ files for each requested dataset size, with
~15% duplication rate for realistic benchmarking. 150 bp reads,
Phred 'F' (score 37). Deterministic (seed=42).
"""

import argparse
import gzip
import os
import random
import sys

SEED = 42
READ_LEN = 150
QUAL = "F" * READ_LEN
DUP_RATE = 0.15  # ~15% of sequences will be duplicates

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choices("ACGT", k=length))


def generate_dataset(n_reads: int, output_dir: str) -> None:
    """Generate R1/R2 FASTQ files for a given number of reads."""
    rng = random.Random(SEED + hash(str(n_reads)))

    # Generate unique sequences (85% of total)
    n_unique = int(n_reads * (1 - DUP_RATE))
    n_dups = n_reads - n_unique

    print(f"  Generating {n_unique} unique + {n_dups} duplicates = {n_reads} total reads")

    # Pre-generate all unique sequences for R1 and R2
    unique_r1 = [random_dna(READ_LEN, rng) for _ in range(n_unique)]
    unique_r2 = [random_dna(READ_LEN, rng) for _ in range(n_unique)]

    # Assign which indices will be duplicated (choose randomly from unique pool)
    dup_indices = rng.choices(range(n_unique), k=n_dups)

    # Build the read list: [unique_seqs, dup_seqs]
    reads_r1 = list(unique_r1)
    reads_r2 = list(unique_r2)
    for idx in dup_indices:
        reads_r1.append(unique_r1[idx])
        reads_r2.append(unique_r2[idx])

    # Shuffle so duplicates aren't all at the end
    combined = list(zip(reads_r1, reads_r2))
    rng.shuffle(combined)
    reads_r1, reads_r2 = zip(*combined) if combined else ([], [])
    reads_r1 = list(reads_r1)
    reads_r2 = list(reads_r2)

    r1_path = os.path.join(output_dir, f"{n_reads}_R1.fastq")
    r2_path = os.path.join(output_dir, f"{n_reads}_R2.fastq")

    with open(r1_path, "w") as f1, open(r2_path, "w") as f2:
        for i, (seq1, seq2) in enumerate(zip(reads_r1, reads_r2), start=1):
            f1.write(f"@BENCH_SEQ_{i} 1:N:0:REF\n{seq1}\n+\n{QUAL}\n")
            f2.write(f"@BENCH_SEQ_{i} 2:N:0:REF\n{seq2}\n+\n{QUAL}\n")

    print(f"  Wrote {r1_path}  ({os.path.getsize(r1_path)} bytes)")
    print(f"  Wrote {r2_path}  ({os.path.getsize(r2_path)} bytes)")

    # Also create gzipped versions
    for path in [r1_path, r2_path]:
        gz_path = path + ".gz"
        with open(path, "rb") as src, gzip.open(gz_path, "wb", compresslevel=6) as dst:
            dst.writelines(src)
        print(f"  Wrote {gz_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic FASTQ reads for FastDedup benchmarking")
    parser.add_argument(
        "sizes", nargs="*", type=int, default=[100_000, 1_000_000, 10_000_000],
        help="Dataset sizes in reads (default: 100000 1000000 10000000)"
    )
    parser.add_argument(
        "--output-dir", default=None,
        help="Output directory (default: benchmarks/data/)"
    )
    args = parser.parse_args()

    output_dir = args.output_dir or os.path.join(BASE_DIR, "data")
    os.makedirs(output_dir, exist_ok=True)

    print(f"Generating synthetic FASTQ datasets in {output_dir}")
    print(f"  Read length: {READ_LEN} bp")
    print(f"  Duplication rate: {DUP_RATE*100:.0f}%")
    print(f"  Seed: {SEED}")

    for n in args.sizes:
        print(f"\n--- {n:,} reads ---")
        generate_dataset(n, output_dir)

    print(f"\n✅ All datasets generated in {output_dir}")


if __name__ == "__main__":
    main()
