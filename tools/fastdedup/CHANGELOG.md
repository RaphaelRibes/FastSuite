# Changelog

## v1.4.0 — 2026-10-26

### ✨ New Features

- **Standard streams (D-4)**: `-` is accepted for `-1/--input`, `-2/--input-r2`, `-o/--output` and `-p/--output-r2`. `zcat in.fastq.gz | fastdedup -1 - -o - | ...` runs with zero intermediate files.
- `--out-format <fastq|fasta>` and `--out-compression <none|gz>`: explicit format and compression when writing to a stream; named paths keep inferring both from the extension.
- `--json -` writes the report to stdout; the summary line then moves to stderr.
- Input compression is detected from magic bytes, so gzip on a stream is decoded like gzip in a file.

### 🐛 Fixes

- **Which copy of a duplicated record survives no longer depends on thread
  scheduling.** At `-j 8` under load, six runs kept six different header
  multisets (sequence multiset identical); on a quiet box the same binary was
  24/24 identical. The writer now re-orders worker results by input ordinal
  before the decision, so the retained copy is the **first occurrence in input
  order** — a documented contract, not an accident.
- **A refused FASTA destination is no longer destroyed.** Writing FASTA to a
  `.fq` path truncated the file and then refused it; with `--force` it deleted
  it first. The format check now runs before any output is opened, so the
  destination survives intact and no zero-byte file is left behind.
- **`-p -` no longer writes the summary line into the R2 stream.** The summary
  moves to stderr whenever stdout carries data, and `-p - --json -` is refused
  instead of interleaved.
- **The paired-end summary counts reads**, matching its own `--json` report; it
  used to count pairs and disagree with itself by a factor of two.
- **A missing output parent directory is created**, as the other five tools do,
  instead of failing with "Failed to create file".

### ⚠️ Behaviour changes

- The report publishes the resume accounting: `reads_preloaded`,
  `reads_already_present`, `reads_new`, with `reads_out == reads_preloaded +
  reads_new` and `reads_in == reads_new + reads_already_present +
  reads_duplicate`. A resumed run says so on stderr unless `-q`.
- `-p -` routes the summary to stderr; `--json -` alongside any `-` output is
  refused.
- The PE summary's unit changes from pairs to reads.

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