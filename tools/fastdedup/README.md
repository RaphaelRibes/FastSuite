![Version](https://img.shields.io/badge/Version-1.3.0-blue)
[![Crates.io](https://img.shields.io/crates/v/fastdedup?color=orange&label=cargo)](https://crates.io/crates/fastdedup)
[![Conda Receipt](https://img.shields.io/badge/install%20with-bioconda-brightgreen)](https://anaconda.org/bioconda/fdedup)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19630793.svg)](https://doi.org/10.5281/zenodo.19630793)

# FastDedup

FastDedup (FDedup) is a fast and memory-efficient FASTX PCR deduplication tool written in Rust.
It utilizes [needletail](https://github.com/onecodex/needletail) for high-performance sequence parsing, [xxh3](https://github.com/DoumanAsh/xxhash-rust) for rapid hashing, and [rustc-hash](https://github.com/rust-lang/rustc-hash) for a low-overhead memory cache.

![benchmark](imgs/PE.svg)

Read the paper on [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.04.29.721745v1).

## Features

- **Fast & Memory Efficient**: Uses zero-allocation sequence parsing and a non-cryptographic high-speed hashing cache, which automatically scales based on the estimated input file size.
- **Supports Compressed Formats**: Transparently reads both uncompressed and GZIP compressed (`.gz`) FASTQ/FASTA files. Writes to both uncompressed and GZIP compressed formats.
- **OOM-Safe Allocations**: All internal hash tables use `try_reserve` to detect out-of-memory conditions and return a clear error instead of panicking silently.
- **Abundance Tracking** (`--abundance`): Counts occurrences of each unique sequence and annotates the output header with `;size=N`. The second pass can be accelerated with `--temp-dir` to avoid re-decompressing the input gzip.
- **Multi-Threaded Mode** (`--threads` / `-j`): Parallel hashing and writing pipeline using a producer/consumer pattern with crossbeam channels, plus in-process gzip compression via gzippy (no subprocess dependency).
- **Pre-Counting** (`--pre-count`): Scans the input first to obtain the exact sequence count, giving the hash table a perfect capacity estimate — avoids resize-induced memory spikes and reduces peak RAM by up to 32%.
- **Incremental Deduplication & Auto-Recovery**: By default, FDedup appends new sequences to an existing uncompressed output file. It safely pre-loads existing hashes to prevent duplicates. If an uncompressed output file is corrupted due to a previous crash, FDedup automatically truncates it to the last valid sequence and resumes safely.

## Installation

### From pre-compiled binaries or container images
You can download the latest pre-compiled binaries from the [releases page](https://github.com/RaphaelRibes/FastDedup/releases)

### From bioconda
The recommended way to install FastDedup is with pixi through bioconda:

```bash
pixi add bioconda::fdedup
```

### From Cargo
You can install FastDedup directly from Cargo:

```bash
cargo install fastdedup
```

## CLI Usage

```bash
fdedup [OPTIONS] --input <INPUT>
```

- `-1, --input <INPUT>`: Path to the input FASTA/FASTQ/GZ file (R1 or Single-End).
- `-2, --input-r2 <INPUT_R2>`: Path to the input R2 file (Optional, enables Paired-End mode).
- `-o, --output <OUTPUT>`: Path to the output file (R1 or Single-End). Required unless `--dry-run` is set.
- `-p, --output-r2 <OUTPUT_R2>`: Path to the output R2 file (Required if `-2` is provided).
- `-f, --force`: Overwrite the output file if it exists (instead of pre-loading hashes and appending).
- `-v, --verbose`: Print processing stats, such as execution time, number of sequences, and duplication rates.
- `-s, --dry-run`: Calculate duplication rate without creating an output file.
- `-t, --threshold <THRESHOLD>`: Threshold for automatic hash size selection$^1$ (default: 0.001).
- `-H, --hash <HASH>`: Manually specify hash size (64 or 128 bits).
- `-c, --compression <LEVEL>`: GZIP compression level, 1–9 (default: 6).
- `-P, --read-length <LENGTH>`: Expected read length in base pairs, used to tune I/O buffers (default: 150).
- `-j, --threads <N>`: Number of threads for parallel processing (0 = auto-detect, default: 1).
- `--abundance`: Enable abundance tracking mode — annotate output headers with `;size=N`.
- `--abundance-sep <SEP>`: Separator for abundance annotation (default: `;size=`).
- `--pre-count`: Pre-count sequences for exact hash table sizing (reduces peak RAM).
- `--temp-dir <DIR>`: Directory for temp file in abundance mode (avoids re-reading gzip).
- `-` in place of `-1/--input`, `-2/--input-r2`, `-o/--output`, `-p/--output-r2`: read from stdin / write to stdout (D-4, standard streams).
- `--out-format <fastq|fasta>`: record format when writing to a stream. Named paths keep inferring it from the extension. Default: `fastq`.
- `--out-compression <none|gz>`: compression when writing to a stream. Named paths keep inferring it from the extension. Default: `none`.

1: The probability $p$ of having at least one collision is calculated as $p= \frac{x^2}{2 \cdot 2^{64}}$ where $x$ is the estimated number of hashes.
If the probability is higher than the specified threshold, FDedup will automatically switch to 128-bit hashing to nullify the risk of collisions.

### Streaming (zero intermediate files)

```bash
zcat in.fastq.gz | fastdedup -1 - -o - | fastprep -1 - -o clean.fastq.gz
```

- Input compression is detected from **magic bytes**, not the extension, so a stream and a mislabelled file behave identically.
- Output format and compression are explicit (`--out-format`, `--out-compression`): a stream has no name to infer from, and the pipe's bytes must not silently depend on the producer's format.
- Streams are written fresh: no preload, no resume, no truncation. Resume is a file-only feature.
- With `-o -` (or `--json -`) the summary line goes to **stderr**, so stdout carries data only.
- Refusals (exit 1, reason on stderr): both inputs as `-`; both outputs as `-`; `--json -` together with `-o -`; `--pre-count` with stdin; `--abundance` with stdin but no `--temp-dir`.

### Run it from Cargo

You can run it directly from Cargo:

```bash
cargo run --release -- --input <INPUT> [OPTIONS]
```

### Run with Pixi

You can also rely on Pixi to run:

```bash
pixi run cargo build --release
pixi run fdedup --input <INPUT> [OPTIONS]
```

### Run with Singularity / Apptainer

You can download the [latest release](https://github.com/RaphaelRibes/FastDedup/releases/tag/v1.0.0) and run the containerized version of FDedup:

Using Apptainer:
```bash
apptainer run fdedup.sif fdedup --input <INPUT> [OPTIONS]
```

Using Singularity:
```bash
singularity run fdedup.sif fdedup --input <INPUT> [OPTIONS]
```

> Note: `--force` is very slow when used in a Singularity container. We recommend just deleting the output file before running the container if you want to start from scratch.

### Run directly with the binary

After downloading the latest release, you can run the binary directly:

```bash
./fdedup --input <INPUT> [OPTIONS]
```

## Build it yourself
### Requirements
If you want to build it from source, you need to have the following dependencies installed:
- [Rust](https://rustup.rs/) (>= 1.85)
- [Pixi](https://pixi.sh) and [pixitainer](https://github.com/RaphaelRibes/pixitainer) to build the container

### Build and run with Cargo
You can build and run FastDedup directly with Cargo:

```bash
cargo build --release
cargo run --release -- --input <INPUT> [OPTIONS]
```

### Build and run with pixitainer
You can also build and run FastDedup with pixitainer:

```bash
pixi containerize
apptainer run fdedup.sif fdedup --input <INPUT> [OPTIONS]
```

## Why `-j 1..4` are all the same speed

With paired input and gzip output, `-j 1`, `-j 2` and `-j 4` all take about the
same wall time, and `-j 8` is several times faster. It looks like a bug; it is
the thread budget.

Measured on 5 M pairs, median of 2 reps:

| arm | wall | CPU | busy cores | peak RSS |
|---|---:|---:|---:|---:|
| `--gz-level 6 -j 1` | 84.20 s | 84.03 s | 1.00 | 79.5 MB |
| `--gz-level 6 -j 2` | 80.67 s | 87.10 s | 1.08 | 101.0 MB |
| `--gz-level 6 -j 4` | 81.80 s | 92.62 s | 1.13 | 104.4 MB |
| `--gz-level 6 -j 8` | **20.44 s** | 125.48 s | **6.14** | 3 776 MB |
| uncompressed output, `-j 1` | 5.83 s | 5.77 s | 0.99 | 76.9 MB |
| uncompressed output, `-j 4` | 5.71 s | 10.73 s | 1.88 | 100.9 MB |
| uncompressed output, `-j 8` | 8.18 s | 31.21 s | 3.82 | 149.6 MB |

Two things follow from that table.

**The hashing is not the problem.** With uncompressed output the whole run is
~6 s at any `-j`, so the 75 s gap between `-j 1` and `-j 8` is entirely in the
compressed-output path.

**Parallel compression only switches on at `-j` around 6.** The busy-core count
shows it directly: 1.00, 1.08 and 1.13 cores at `-j 1/2/4`, then 6.14 at `-j 8`.
Compression is single-threaded for the first three. The cause is the thread
budget in `parallel.rs`:

```rust
// paired: two output writers, so the spare threads are split between them
let gzip_threads = (num_workers - 2) / 2;   // j2->1, j4->1, j8->3
```

and the parallel compressor requires `num_threads > 1`. One spare thread per
writer does not engage it, so `-j 2` and `-j 4` fall back to single-threaded
flate2 — which is why they match `-j 1` rather than beating it. The single-end
path uses a more generous budget (`num_workers - 2`, no halving) and engages
parallel compression already at `-j 4`.

**So `-j 8` (or higher) is what makes gzip output pay on paired input.** Below
that, spending more threads buys nothing for compression; it only helps the
hashing, which is already close to I/O-bound.

**And `--gz-level 1` sidesteps the whole thing.** On the same corpus, level 1
single-threaded runs in 12.9 s — faster than level 6 with parallel compression
at `-j 8` (20.4 s) — for 1.87× more bytes and 25× less peak memory. When speed
is what you want, change the level before adding threads.

## Output order is not part of the contract

The **kept set** of reads is stable across `-j`; the **order** in which they are
written is not. Measured on a 200 000-pair corpus at `-j 1/2/4/8`: every run
produced the same 200 000 records, and one single digest of the sorted records
matched across all four — but three distinct start-to-finish orderings appeared
(j1 and j2 happened to agree; j4 and j8 each differed).

So a check that compares output byte-for-byte between two runs is comparing
something the tool never promised. Compare **sorted** records instead. If you
need a deterministic byte-for-byte output, sort the file once at the end; the
tool will not do it for you, because serialising the emit order costs wall time.

## Recommendations

If you are using FDedup in a pre-processing step, we recommend you to not export your file to a `.gz` format.
Incremental/resumable deduplication, only works with uncompressed output files.
If you output to a compressed format, FDedup requires `--force` to restart from scratch on any subsequent run.
However, if you output to an uncompressed format, FDedup will automatically detect any crash-induced corruption, safely truncate the file to the last valid sequence, and seamlessly resume deduplication.

## To-Do List

- [x] Support for **Paired-End read deduplication**.
- [x] Add **Multithreading** to parallelize sequence hashing and processing.
- [x] Support tracking sequence **abundances** (counts) instead of naive discarding.
- [x] Add a possibility for exporting sequences as **FASTA**.
- [x] Improve error handling.

## Citation

If you use FastDedup in your work, please cite:

> Ribes R., Mandier C., Baniel A. FastDedup: a fast and memory-efficient tool for read deduplication. *bioRxiv* 2026.04.29.721745 (2026). doi:[10.64898/2026.04.29.721745](https://www.biorxiv.org/content/10.64898/2026.04.29.721745v1)

Additionally, the specific software version can be cited via Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19630793.svg)](https://doi.org/10.5281/zenodo.19630793)

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Author

[Raphaël Ribes](https://www.raphaelrib.es) - ISDM, Univ Montpellier, CIRAD, INSERM, Montpellier, France (coding and design)

[Céline Mandier](https://gitlab.in2p3.fr/celine.mandier1) - ISDM, Univ Montpellier, CIRAD, INSERM, Montpellier, France (design)

[Alice Baniel](mailto:alice.baniel@umontpellier.fr) - ISEM, Univ Montpellier, CNRS, IRD, Montpellier, France (funding)

# Acknowledgements

Computations were performed on the ISDM-MESO HPC platform, funded in the framework of State-region planning contracts (Contrat de plan État-région – CPER) by the French Government, the Occitanie/Pyrénées-Méditerranée Region, Montpellier Méditerranée Métropole, and the University of Montpellier.

Funding came from a grant awarded to Alice Baniel by the University of Montpellier (I-SITE Excellence Programme 2022 - SOFIT).
