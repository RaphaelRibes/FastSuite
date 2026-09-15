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
- `-o, --output <OUTPUT>`: Path to the output file (R1 or Single-End). Defaults to `output_R1.fastq.gz`.
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

1: The probability $p$ of having at least one collision is calculated as $p= \frac{x^2}{2 \cdot 2^{64}}$ where $x$ is the estimated number of hashes.
If the probability is higher than the specified threshold, FDedup will automatically switch to 128-bit hashing to nullify the risk of collisions.

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
