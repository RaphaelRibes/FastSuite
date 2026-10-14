# Changelog

## Unreleased

### ✨ New Features

- **Standard streams (D-4)**: `-` is accepted for `-1/--input`, `-2/--input-r2`, `-o/--output` and `-p/--output-r2`. `zcat in.fastq.gz | fastdedup -1 - -o - | ...` runs with zero intermediate files.
- `--out-format <fastq|fasta>` and `--out-compression <none|gz>`: explicit format and compression when writing to a stream; named paths keep inferring both from the extension.
- `--json -` writes the report to stdout; the summary line then moves to stderr.
- Input compression is detected from magic bytes, so gzip on a stream is decoded like gzip in a file.

## v1.3.0 — 2025-05-18

### ✨ New Features

- **Multi-threaded pipeline** (`--threads` / `-j`): Producer/consumer pattern with crossbeam channels for parallel hashing and writing. Auto-detect with `-j 0`.
- **In-process gzip** (flate2+gzippy): Replaced external `pigz` subprocess with in-memory `GzDecoder` for input and `gzippy` multi-threaded compression for output. Up to 3.5× faster at high thread counts vs pigz.
- **Abundance tracking** (`--abundance`): Count occurrences per unique sequence and annotate output headers with `;size=N`.
- **Temp-file acceleration** (`--temp-dir`): Avoid re-decompressing gzip input on the abundance second pass by writing a raw temp file.
- **OOM-safe allocations**: All hash tables use `try_reserve` instead of panicking on out-of-memory conditions.
- **Pre-counting** (`--pre-count`): Scan input first for exact hash table sizing, reducing peak RAM by up to 32%.
- **Paired-end dry-run desync fix**: Corrected desynchronization detection in paired-end `--dry-run` mode.
- **Dry-run without output**: `--output` is now optional when `--dry-run` is set. Useful for a quick stats check before the real run.

### 🔧 Changes

- `src/writer.rs` extracted from `utils.rs` for cleaner module separation.
- README overhaul with comprehensive CLI docs, abundance/MT/OOM/pre-count flags.
- Added `Alice Baniel` to authors.

### 🐛 Bug Fixes

- Fixed PE dry-run desync detection (task2).
- General error handling improvements.
- **Fix `-j 0` auto-detect:** `cfg.threads == 0` now correctly dispatches to the multi-threaded pipeline instead of silently falling through to single-threaded.
- **Fix hardcoded duplicate count:** `--abundance --temp-dir` (SE and PE) now correctly computes `duplicates = total_reads - tracker.len()` instead of returning `0`.
- **Fix PE desync detection in abundance mode:** Added `base_read_id` comparison + trailing record check to PE abundance (both temp-dir and classic paths).
- **Fix classic abundance PASS 2:** Replaced `tracker.track()` (which corrupted counts) with `tracker.get_count()` + `written` set, so unique reads are correctly written to output instead of all reads.

### 📝 Docs

- Restored I-SITE / SOFIT funding line in README Acknowledgements.
- Fixed stale `pigz` docstrings in `src/parallel.rs` — input is now streamed via flate2, output via gzippy.
- `demo.sh` reworked: real dedup demo (no dry-run) + separate dry-run section.
- `--abundance` help now notes "incompatible with --threads (always runs single-threaded)"; runtime warning emitted when combined with `-j > 1`.
- `--compression` help now notes "levels < 6 disable multi-threaded gzip output".

### 🧪 Tests

- Added 7 integration tests for `--abundance`, `--pre-count`, `--temp-dir`, and PE desync detection in abundance mode.
- Added `flate2` to `[dev-dependencies]` explicitly.
- Removed `__echo` leftover task in `pixi.toml`.

### ⚙️ Build

- pixi.toml: full build pipeline (conda publish, archive, Apptainer container, Docker).
- pixi-build preview with `RUSTFLAGS="-C link-self-contained=yes"`.