# fastdedup: stdin/stdout streaming

> **For Hermes:** Use subagent-driven-development skill to implement this plan task-by-task.

**Goal:** Make fastdedup accept `-` for `-1/-2/-o/-p` (decision D-4), so `zcat in.fastq.gz | fastdedup -1 - -o - | fastprep -1 - -o clean.fastq.gz` runs with zero intermediate files.

**Architecture:** A new `src/stream.rs` owns the `-` convention: `is_stream(path)` and `open_fastx(path)` (stdin via `needletail::parse_fastx_reader`, which sniffs gzip/bzip2/xz/zstd from magic bytes — exactly the D-4 rule that compression is decided from bytes, never from a filename). Output to a stream gets an explicit decision via two new flags, `--out-format fastq|fasta` and `--out-compression none|gz` (defaults: fastq, none); inheritance from the input is rejected because a pipe whose bytes silently depend on the producer's format is a trap for the next command. Streams are written fresh: no preload, no resume, no truncation. Everything that cannot work on a stream refuses with a reason on stderr (exit 1): two inputs as `-`, two outputs as `-`, `--json -` with `-o -`, `--pre-count` with stdin, `--abundance` with stdin but no `--temp-dir`.

**Tech Stack:** Rust (edition 2024, rustc 1.85+), needletail 0.6.3, flate2, clap 4 derive, anyhow, assert_cmd/assert_fs/predicates (dev).

---

## Context (verified 2026-10-12)

* Current failure, confirmed:
  ```
  $ ./tools/fastdedup/target/release/fastdedup -1 - -o /tmp/o.fastq < /dev/null
  Error: Input file not found: -
  ```
  The existence check in `src/main.rs:180-188` rejects `-` before anything else runs.
* Test suite today: **87 tests**, all green (`cargo test --offline` from `tools/fastdedup`):
  29 unit (in-bin), 7 in `tests/abundance_tests.rs`, 51 in `tests/cli.rs`.
  This plan adds 2 unit + 16 integration tests → **105 expected** at the end.
* Reference implementation for the reader: `tools/fastcheck/src/reader.rs` (`open()` peeks 18 bytes, detects BGZF/gzip from magic bytes, chains the peeked bytes back, labels the mode `stdin-bgzf-parallel` / `stdin-gzip-serial` / `stdin-plain`).
* needletail already does the D-4 input rule for us: `parse_fastx_reader` reads the first two bytes and picks `MultiGzDecoder` on gzip magic, so handing it `io::stdin()` gives plain + gzip on a stream for free. Serial gzip on stdin is accepted (a stream cannot be indexed for parallel BGZF; that stays file-only and is not a regression).
* D-4 (`docs/DECISIONS.md`): all tools except `fastmaps` take `-`; magic-byte detection; streams cannot be seeked; the output-compression question was left open and this plan **settles it** as `--out-compression none|gz`.
* CLI-SPEC §4 gap 7 lists `-`/stdin-stdout as fastdedup's remaining gap.

**Known pre-existing issue, out of scope (do not fix):** `fastdedup -s -j 4` (dry-run with threads > 1) fails with `Failed to create file: ""` because `execute_deduplication_parallel` has no dry-run branch. Observed while grounding this plan; unrelated to streaming. All tests here use the default `-j 1` for dry-runs.

## Files

| file | role |
|---|---|
| `tools/fastdedup/src/stream.rs` | **create** — `STREAM`, `is_stream`, `open_fastx` + unit tests |
| `tools/fastdedup/src/main.rs` | modify — validation exemptions, stream refusals, `--json -`, summary-to-stderr |
| `tools/fastdedup/src/cli.rs` | modify — `--out-format`, `--out-compression`, help text for `-` |
| `tools/fastdedup/src/writer.rs` | modify — `StreamCompression` enum, `prepare_stream_writer` |
| `tools/fastdedup/src/processor.rs` | modify — `DedupConfig` fields, `open_output` helper, `open_fastx` at every input site |
| `tools/fastdedup/src/parallel.rs` | modify — `open_input_reader` + `create_output_writer` stream branches |
| `tools/fastdedup/tests/streaming_tests.rs` | **create** — 16 integration tests |
| `tools/fastdedup/README.md` | modify — CLI docs + streaming section |
| `tools/fastdedup/CHANGELOG.md` | modify — Unreleased entry |
| `docs/DECISIONS.md` | modify — settle D-4's output-compression bullet |
| `docs/CLI-SPEC.md` | modify — gap 7 + two common-core flags |

## Conventions for the implementer

* All `cargo` commands run from `tools/fastdedup/`. Binary smoke tests run from the repo root.
* `tests/cli.rs` is French-named; new tests are English-named. Match the existing test style (one behavior per test, `predicate::str::contains` for stderr messages).
* Existing error phrasings to preserve (tests depend on them): `Input file not found`, `Failed to read input file`, `FASTA → FASTQ conversion not supported`.
* Commit after every task, from the repo root: `feat:`/`fix:`/`test:`/`docs:` prefixes.

---

### Task 1: Create `stream.rs` with the `-` primitive and unit tests

**Objective:** One module owns the stream convention so no other file re-implements `== "-"`.

**Files:**
- Create: `tools/fastdedup/src/stream.rs`
- Modify: `tools/fastdedup/src/main.rs:1-6` (add `mod stream;`)

**Step 1: Write the failing unit tests**

`tools/fastdedup/src/stream.rs`:

```rust
//! D-4 standard streams: `-` means stdin (input) or stdout (output).
//!
//! A stream has no filename, so nothing may be inferred from one:
//!   * compression on INPUT is sniffed from magic bytes (needletail's
//!     auto-detection: gzip/bzip2/xz/zstd, plain as fallback),
//!   * format and compression on OUTPUT are explicit flags
//!     (`--out-format`, `--out-compression`), defaulting to plain FASTQ.

/// The path that means "the stream" (docs/DECISIONS.md, D-4).
pub const STREAM: &str = "-";

/// True when `path` designates stdin/stdout rather than a file.
pub fn is_stream(path: &str) -> bool {
    path == STREAM
}

/// Open a FASTX input: stdin for `-`, the file otherwise.
///
/// `parse_fastx_reader` sniffs the first two bytes and selects the decoder,
/// so gzip on stdin behaves exactly like gzip in a mislabelled file.
pub fn open_fastx(path: &str) -> Result<Box<dyn needletail::FastxReader>> {
    if is_stream(path) {
        needletail::parse_fastx_reader(std::io::stdin())
            .map_err(|e| anyhow::anyhow!("reading FASTX from stdin: {e}"))
    } else {
        needletail::parse_fastx_file(path)
            .map_err(|e| anyhow::anyhow!("Failed to read input file {path}: {e}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dash_is_the_stream() {
        assert!(is_stream("-"));
        assert!(!is_stream("--"));
        assert!(!is_stream("./-"));
        assert!(!is_stream("in.fastq"));
    }

    #[test]
    fn a_stream_is_just_a_reader() {
        // A stream is a `Read`; the same sniffing path must apply to it.
        let body = b"@r1\nACGT\n+\nIIII\n";
        let mut rdr = needletail::parse_fastx_reader(&body[..]).unwrap();
        let rec = rdr.next().expect("one record").unwrap();
        assert_eq!(rec.id(), b"r1");
        assert!(rdr.next().is_none());
    }
}
```

In `src/main.rs`, insert after line 5 (`mod processor;`):

```rust
mod stream;
```

**Step 2: Run the tests to verify failure**

Run: `cargo test --offline --bin fastdedup dash_is_the_stream a_stream_is_just_a_reader`
Expected: FAIL — compile error `cannot find module 'stream'` (the file does not exist yet — write it first, then run; on the very first run the module is missing).

**Step 3: (already written above)** — create the file, add the `mod` line.

**Step 4: Run the tests to verify pass**

Run: `cargo test --offline --bin fastdedup dash_is_the_stream`
Expected: `test result: ok. 2 passed; 0 failed` (the second assertion lives in the same module run).

**Step 5: Commit**

```bash
git add tools/fastdedup/src/stream.rs tools/fastdedup/src/main.rs
git commit -m "feat: add stream module owning the '-' stdin/stdout convention (D-4)"
```

---

### Task 2: Exempt `-` from file-existence checks and add the stream refusals

**Objective:** `-1 -` no longer dies at `Input file not found`, and the combinations that cannot work on a stream refuse with a reason.

**Files:**
- Modify: `tools/fastdedup/src/main.rs:178-215` (input validation block)

**Step 1: Write the failing tests**

Create `tools/fastdedup/tests/streaming_tests.rs`:

```rust
use assert_cmd::Command;
use assert_cmd::cargo;
use assert_fs::prelude::*;
use predicates::prelude::*;
use std::fs;

fn cmd() -> Command {
    cargo::cargo_bin_cmd!("fastdedup")
}

/// Two reads, second identical to the first: 2 in, 1 out.
const FASTQ: &[u8] = b"@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nACGTACGT\n+\n!!!!!!!!\n";

#[test]
fn pre_count_with_stdin_is_refused() {
    cmd()
        .arg("-1")
        .arg("-")
        .arg("-s")
        .arg("--pre-count")
        .write_stdin(FASTQ)
        .assert()
        .failure()
        .stderr(predicate::str::contains("--pre-count"));
}

#[test]
fn both_inputs_as_stdin_is_refused() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    cmd()
        .arg("-1")
        .arg("-")
        .arg("-2")
        .arg("-")
        .arg("-o")
        .arg(temp_dir.child("out1.fastq").path())
        .arg("-p")
        .arg(temp_dir.child("out2.fastq").path())
        .write_stdin(FASTQ)
        .assert()
        .failure()
        .stderr(predicate::str::contains("cannot be stdin"));
}

#[test]
fn abundance_from_stdin_requires_temp_dir() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    cmd()
        .arg("-1")
        .arg("-")
        .arg("--abundance")
        .arg("-o")
        .arg(temp_dir.child("out.fastq").path())
        .write_stdin(FASTQ)
        .assert()
        .failure()
        .stderr(predicate::str::contains("--temp-dir"));
}
```

**Step 2: Run the tests to verify failure**

Run: `cargo test --offline --test streaming_tests`
Expected: 3 failed — stderr is `Input file not found: -`, not the expected messages.

**Step 3: Implement**

Replace `src/main.rs:180-188` (the two existence checks) with:

```rust
    let streaming_input = crate::stream::is_stream(&args.input);
    if !streaming_input && !Path::new(&args.input).exists() {
        bail!("Input file not found: {}", args.input);
    }

    if let Some(ref r2) = args.input_r2 {
        if !crate::stream::is_stream(r2) && !Path::new(r2).exists() {
            bail!("R2 input file not found: {r2}");
        }
    }
```

Then insert the refusal block right after the output-path resolution (after the `let output = match ...` block ending at `src/main.rs:197`):

```rust
    // --- D-4: what a single stream cannot carry ---
    let out1_stream = args.output.as_deref().is_some_and(crate::stream::is_stream);
    let out2_stream = args.output_r2.as_deref().is_some_and(crate::stream::is_stream);

    if streaming_input && args.input_r2.as_deref().is_some_and(crate::stream::is_stream) {
        bail!(
            "Both --input-r1 (-1) and --input-r2 (-2) cannot be stdin: one stream cannot carry two files. Pipe one mate and name the other."
        );
    }
    if streaming_input && args.pre_count {
        bail!(
            "--pre-count needs a second pass over the input, impossible on a stream. Drop --pre-count or pass a file path."
        );
    }
    if streaming_input && args.abundance && args.temp_dir.is_none() {
        bail!(
            "--abundance needs a second pass over the input. With stdin, pass --temp-dir so unique reads are spilled to disk and the stream is read once."
        );
    }
    if out1_stream && out2_stream {
        bail!(
            "--out1 (-o) and --out2 (-p) cannot both be stdout: R1 and R2 would interleave into one corrupted stream. Send one mate to a file."
        );
    }
    if out1_stream && args.json.as_deref() == Some("-") {
        bail!(
            "--out1 (-) and --json (-) cannot both write to stdout. Send the report to a file, or drop one of them."
        );
    }
```

Also update the verbose input echo (`src/main.rs:218`) so a stream reads sensibly:

```rust
    if args.verbose {
        println!(
            "Primary input file (R1): {}",
            if streaming_input { "(stdin)" } else { &args.input }
        );
```

**Step 4: Run the tests to verify pass**

Run: `cargo test --offline --test streaming_tests`
Expected: `test result: ok. 3 passed; 0 failed`.

**Step 5: Commit**

```bash
git add tools/fastdedup/src/main.rs tools/fastdedup/tests/streaming_tests.rs
git commit -m "feat: accept '-' as input, refuse the stream combinations that cannot work (D-4)"
```

---

### Task 3: Read stdin in the single-end serial path

**Objective:** `fastdedup -1 - -s` computes real statistics from a pipe.

**Files:**
- Modify: `tools/fastdedup/src/processor.rs:60` (dry-run reader) and `:126` (main reader)

**Step 1: Write the failing tests**

Append to `tests/streaming_tests.rs`:

```rust
#[test]
fn dry_run_from_stdin() {
    cmd()
        .arg("-1")
        .arg("-")
        .arg("-s")
        .arg("-v")
        .write_stdin(FASTQ)
        .assert()
        .success()
        .stdout(predicate::str::contains("Processed fragments: 2"))
        .stdout(predicate::str::contains("Duplicates removed: 50.00%"));
}

#[test]
fn empty_stdin_is_a_clean_refusal() {
    cmd()
        .arg("-1")
        .arg("-")
        .arg("-s")
        .write_stdin(b"")
        .assert()
        .failure()
        .stderr(predicate::str::contains("empty"));
}
```

**Step 2: Run the tests to verify failure**

Run: `cargo test --offline --test streaming_tests dry_run_from_stdin empty_stdin`
Expected: both fail — the reader still calls `parse_fastx_file("-")` → `Failed to read input file -`.

**Step 3: Implement**

In `src/processor.rs`, replace both occurrences of

```rust
    let mut reader = parse_fastx_file(input_path).context("Failed to read input file")?;
```

(dry-run branch, line 60; main branch, line 126 — both are byte-identical, use the surrounding lines to disambiguate) with:

```rust
    let mut reader = crate::stream::open_fastx(input_path)?;
```

**Step 4: Run the tests to verify pass**

Run: `cargo test --offline --test streaming_tests dry_run_from_stdin empty_stdin`
Expected: `test result: ok. 2 passed; 0 failed`.

**Step 5: Commit**

```bash
git add tools/fastdedup/src/processor.rs tools/fastdedup/tests/streaming_tests.rs
git commit -m "feat: read single-end input from stdin in the serial path"
```

---

### Task 4: Stream output — flags, writer, and the summary-to-stderr rule

**Objective:** `fastdedup -1 - -o -` writes deduplicated FASTQ to stdout; the summary line moves to stderr so stdout carries only data.

**Files:**
- Modify: `tools/fastdedup/src/cli.rs` (two new flags + `-` in help text)
- Modify: `tools/fastdedup/src/writer.rs` (`StreamCompression`, `prepare_stream_writer`)
- Modify: `tools/fastdedup/src/processor.rs` (`DedupConfig` fields, `open_output` helper, SE writer + preload guard)
- Modify: `tools/fastdedup/src/main.rs` (`cfg` construction, summary destination)

**Step 1: Write the failing test**

Append to `tests/streaming_tests.rs`:

```rust
#[test]
fn stdin_to_stdout_plain() {
    let assert = cmd()
        .arg("-1")
        .arg("-")
        .arg("-o")
        .arg("-")
        .write_stdin(FASTQ)
        .assert()
        .success();
    let out = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    assert_eq!(
        out.lines().filter(|l| l.starts_with('@')).count(),
        1,
        "2 reads in, 1 unique out"
    );
    assert!(out.starts_with('@'), "stdout must be FASTQ only, no summary");
    let err = String::from_utf8(assert.get_output().stderr.clone()).unwrap();
    assert!(
        err.contains("fastdedup: 2 reads in, 1 out"),
        "summary belongs on stderr when stdout carries data: {err}"
    );
}
```

**Step 2: Run the test to verify failure**

Run: `cargo test --offline --test streaming_tests stdin_to_stdout_plain`
Expected: FAIL — `-o -` is treated as a filename (`Failed to create file: "-"`).

**Step 3: Implement**

`src/cli.rs` — add the two enums after `HashMode` (end of file):

```rust
/// Record format for stream output (`-o -`): D-4 gives a stream no
/// extension to infer from, so the choice is explicit. Default: fastq.
#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutFormat {
    Fastq,
    Fasta,
}

/// Compression for stream output (`-o -`). Default: none — plain bytes in
/// the pipe; ask for `gz` when the consumer wants compressed bytes.
#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutCompression {
    None,
    Gz,
}
```

`src/cli.rs` — add the two args after the `output_r2` field (after line 47), and amend the help text of the I/O args to mention `-`:

```rust
    /// Output record format when the output is a stream (`-`).
    /// Named paths keep inferring the format from their extension.
    #[arg(long, value_enum, default_value_t = OutFormat::Fastq)]
    pub out_format: OutFormat,

    /// Compression when the output is a stream (`-`).
    /// Named paths keep inferring compression from their extension.
    #[arg(long, value_enum, default_value_t = OutCompression::None)]
    pub out_compression: OutCompression,
```

Also change the doc comment of `input` (line 29) to:
`/// Path to the input FASTX file (R1 or Single-End). \`-\` reads stdin (D-4).`
and of `output` (line 38) to:
`/// Path to the output file (R1 or Single-End). \`-\` writes stdout (D-4).`

`src/writer.rs` — add after `OutputFormat::is_gz` (after line 38):

```rust
/// Compression decision for stream output (D-4): a stream has no
/// filename, so `-o -` cannot infer gzip from an extension.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StreamCompression {
    None,
    Gz,
}

/// Open stdout as an output stream. The record format is chosen per
/// record by the caller from `--out-format`; this only decides the bytes'
/// compression. Gzip finishes on drop, exactly like the file path.
pub fn prepare_stream_writer(
    compression: StreamCompression,
    compression_level: u32,
) -> Box<dyn Write> {
    const STREAM_BUFFER: usize = 1 << 20; // 1 MiB, same as the file path
    match compression {
        StreamCompression::None => {
            Box::new(io::BufWriter::with_capacity(STREAM_BUFFER, io::stdout()))
        }
        StreamCompression::Gz => Box::new(GzEncoder::new(
            io::BufWriter::with_capacity(STREAM_BUFFER, io::stdout()),
            Compression::new(compression_level),
        )),
    }
}
```

`src/processor.rs` — add the two fields to `DedupConfig` (after `temp_dir`, line 33):

```rust
    /// Output record format, used when the output is a stream; named paths
    /// infer it from their extension.
    pub out_format: OutputFormat,
    /// Output compression, used when the output is a stream; named paths
    /// infer it from their extension.
    pub out_compression: crate::writer::StreamCompression,
```

`src/processor.rs` — add the `open_output` helper right above `execute_deduplication` (before line 44's `// ====` banner):

```rust
/// Open the output for a run: stdout when the path is `-` (D-4), the file
/// otherwise. The stream branch carries the explicit format/compression
/// decision from `--out-format` / `--out-compression`.
fn open_output(path: &str, cfg: &DedupConfig) -> Result<(Box<dyn Write>, OutputFormat)> {
    if crate::stream::is_stream(path) {
        Ok((
            crate::writer::prepare_stream_writer(cfg.out_compression, cfg.compression_level),
            cfg.out_format,
        ))
    } else {
        crate::writer::prepare_writer(Path::new(path), cfg.force, cfg.compression_level)
    }
}
```

`src/processor.rs` — in `execute_deduplication` only (lines 79-117), guard the preload/resume block and route the writer:

```rust
    // --- PRELOADING ---
    let output_path_bound = Path::new(output_path);
    let is_gz = OutputFormat::is_gz(output_path_bound);

    if crate::stream::is_stream(output_path) {
        // D-4: a stream is written fresh — nothing exists to preload hashes
        // from and nothing to truncate. Resume is a file-only feature.
        if cfg.verbose {
            println!("Output is a stream: written fresh (no preload, no resume).");
        }
    } else if cfg.force {
        if cfg.verbose {
            println!("--force enabled: overwriting output.");
        }
    } else {
        let (preloads, valid_bytes) =
            preload_existing_hashes::<T>(output_path, &mut verifier, cfg.verbose)?;
        // ... unchanged body of the existing else branch ...
    }
```

(replace only the `if cfg.force {` line at 83 with the three-branch form above; the existing `else` body stays as-is), and replace the writer preparation (lines 116-117):

```rust
    // --- WRITER PREPARATION ---
    let (writer, output_format) = open_output(output_path, cfg)?;
```

`src/main.rs` — build the two new config fields (in the `let cfg = DedupConfig { ... }` block, after `temp_dir: args.temp_dir.clone(),`):

```rust
        out_format: match args.out_format {
            OutFormat::Fastq => OutputFormat::Fastq,
            OutFormat::Fasta => OutputFormat::Fasta,
        },
        out_compression: match args.out_compression {
            OutCompression::None => StreamCompression::None,
            OutCompression::Gz => StreamCompression::Gz,
        },
```

with these imports added at the top of `main.rs`:

```rust
use crate::cli::{Cli, OutCompression, OutFormat};
use crate::writer::{OutputFormat, StreamCompression};
```

`src/main.rs` — send the summary line to stderr when stdout carries data (replace the `if !args.quiet {` block at lines 368-378):

```rust
    if !args.quiet {
        // D-4 / CLI-SPEC §0: with `-o -` or `--json -` stdout carries data,
        // so the summary goes to stderr and never pollutes the pipe.
        let macro_ = if out1_stream || args.json.as_deref() == Some("-") {
            eprintln!
        } else {
            println!
        };
        macro_!(
            "fastdedup: {} reads in, {} out, {} duplicates ({:.2}%), {} threads, {:.2?}",
            processed,
            reads_out,
            duplicates,
            dup_pct,
            args.threads,
            elapsed
        );
    }
```

(`out1_stream` was defined in Task 2; `macro_!` works because both `println!` and `eprintln!` have identical format-string syntax.)

**Step 4: Run the tests to verify pass**

Run: `cargo test --offline --test streaming_tests`
Expected: `test result: ok. 4 passed; 0 failed`.

**Step 5: Commit**

```bash
git add tools/fastdedup/src/cli.rs tools/fastdedup/src/writer.rs tools/fastdedup/src/processor.rs tools/fastdedup/src/main.rs tools/fastdedup/tests/streaming_tests.rs
git commit -m "feat: stream output to stdout with explicit --out-format/--out-compression (D-4)"
```

---

### Task 5: Gzip on stdin is decoded by magic bytes

**Objective:** Prove the D-4 input rule: a gzip *stream* decodes like a gzip file, with no filename involved. (No production code — needletail already sniffs the magic bytes.)

**Files:**
- Modify: `tools/fastdedup/tests/streaming_tests.rs`

**Step 1: Write the failing test**

Add the fixture helper and the test at the top of `tests/streaming_tests.rs` (after `const FASTQ`):

```rust
fn gz(body: &[u8]) -> Vec<u8> {
    use flate2::Compression;
    use flate2::write::GzEncoder;
    let mut enc = GzEncoder::new(Vec::new(), Compression::default());
    enc.write_all(body).unwrap();
    enc.finish().unwrap()
}
```

(with `use std::io::Write;` added to the imports) and the test:

```rust
#[test]
fn gzipped_stdin_is_decoded_from_magic_bytes() {
    let assert = cmd()
        .arg("-1")
        .arg("-")
        .arg("-o")
        .arg("-")
        .write_stdin(gz(FASTQ))
        .assert()
        .success();
    let out = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    assert_eq!(
        out.lines().filter(|l| l.starts_with('@')).count(),
        1,
        "the gzip stream must be transparently decoded"
    );
}
```

**Step 2: Run the test to verify failure**

Run: `cargo test --offline --test streaming_tests gzipped_stdin_is_decoded_from_magic_bytes`
Expected: should already pass after Task 3-4 (needletail sniffs `0x1f 0x8b` on any reader, including stdin). If it fails, the regression is in `open_fastx` — fix that, not this test.

**Step 3: No implementation needed.** This task exists to lock the behavior with a regression test.

**Step 4: Run the test to verify pass**

Run: `cargo test --offline --test streaming_tests`
Expected: `test result: ok. 5 passed; 0 failed`.

**Step 5: Commit**

```bash
git add tools/fastdedup/tests/streaming_tests.rs
git commit -m "test: gzip on stdin is decoded from magic bytes (D-4)"
```

---

### Task 6: `--out-compression gz` writes gzip to stdout

**Objective:** A consumer that wants compressed bytes in the pipe asks for them explicitly.

**Files:** none new (flag wired in Task 4; this task only tests it)

**Step 1: Write the failing test**

Append to `tests/streaming_tests.rs`:

```rust
#[test]
fn stdout_can_be_gzip() {
    let assert = cmd()
        .arg("-1")
        .arg("-")
        .arg("-o")
        .arg("-")
        .arg("--out-compression")
        .arg("gz")
        .write_stdin(FASTQ)
        .assert()
        .success();
    let bytes = assert.get_output().stdout.clone();
    assert_eq!(bytes[0], 0x1f, "gzip magic byte");
    assert_eq!(bytes[1], 0x8b, "gzip magic byte");
}
```

**Step 2: Run the test to verify failure**

Run: `cargo test --offline --test streaming_tests stdout_can_be_gzip`
Expected: passes immediately (Task 4 wired `prepare_stream_writer`). It exists as the regression guard for the gz branch of that function — if it fails, Task 4's writer plumbing regressed.

**Step 3: No implementation needed.**

**Step 4: Run the test to verify pass**

Run: `cargo test --offline --test streaming_tests`
Expected: `test result: ok. 6 passed; 0 failed`.

**Step 5: Commit**

```bash
git add tools/fastdedup/tests/streaming_tests.rs
git commit -m "test: --out-compression gz emits gzip magic bytes on stdout"
```

---

### Task 7: `--out-format fasta` on stdout, and FASTA-in → FASTQ-out refusal

**Objective:** Stream output format is explicit; FASTA input can never be silently turned into quality-less FASTQ, even through a pipe.

**Files:** none new

**Step 1: Write the failing tests**

Append to `tests/streaming_tests.rs`:

```rust
#[test]
fn stdout_can_be_fasta() {
    let assert = cmd()
        .arg("-1")
        .arg("-")
        .arg("-o")
        .arg("-")
        .arg("--out-format")
        .arg("fasta")
        .write_stdin(FASTQ)
        .assert()
        .success();
    let out = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    assert!(out.starts_with('>'), "FASTA output must start with '>'");
    assert!(!out.contains('+'), "FASTA has no '+' separator line");
}

#[test]
fn fasta_stdin_to_fastq_stdout_is_refused() {
    cmd()
        .arg("-1")
        .arg("-")
        .arg("-o")
        .arg("-")
        .write_stdin(b">r1\nACGTACGT\n")
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "FASTA → FASTQ conversion not supported",
        ));
}
```

**Step 2: Run the tests to verify failure**

Run: `cargo test --offline --test streaming_tests stdout_can_be_fasta fasta_stdin_to_fastq`
Expected: both should pass after Task 4 (`output_format` now comes from `cfg.out_format`, and `validate_format_compatibility` runs on the first streamed record). A failure means Task 4 did not route `output_format` through `open_output` — fix that.

**Step 3: No implementation needed.**

**Step 4: Run the tests to verify pass**

Run: `cargo test --offline --test streaming_tests`
Expected: `test result: ok. 8 passed; 0 failed`.

**Step 5: Commit**

```bash
git add tools/fastdedup/tests/streaming_tests.rs
git commit -m "test: --out-format fasta on stdout, FASTA to FASTQ refusal over a stream"
```

---

### Task 8: `--json -` writes the report to stdout, summary to stderr

**Objective:** CLI-SPEC §0 — the report may go to the stream; when it does, the summary must move to stderr so a consumer parsing stdout sees only JSON.

**Files:**
- Modify: `tools/fastdedup/src/main.rs:380-409` (report writing)

**Step 1: Write the failing tests**

Append to `tests/streaming_tests.rs`:

```rust
#[test]
fn json_report_to_stdout_moves_summary_to_stderr() {
    let assert = cmd()
        .arg("-1")
        .arg("-")
        .arg("-s")
        .arg("--json")
        .arg("-")
        .write_stdin(FASTQ)
        .assert()
        .success();
    let out = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    assert!(
        out.contains("fastsuite.report.v1"),
        "stdout must carry the report: {out}"
    );
    let err = String::from_utf8(assert.get_output().stderr.clone()).unwrap();
    assert!(
        err.contains("fastdedup: 2 reads in"),
        "summary must move to stderr: {err}"
    );
}

#[test]
fn json_and_fastq_on_stdout_is_refused() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    cmd()
        .arg("-1")
        .arg("-")
        .arg("-o")
        .arg("-")
        .arg("--json")
        .arg("-")
        .write_stdin(FASTQ)
        .assert()
        .failure()
        .stderr(predicate::str::contains("--json (-)"));
}
```

**Step 2: Run the tests to verify failure**

Run: `cargo test --offline --test streaming_tests json_report_to_stdout json_and_fastq`
Expected: `json_and_fastq_on_stdout_is_refused` passes (Task 2 bail); `json_report_to_stdout...` fails — `std::fs::write("-", ...)` creates a file literally named `-`.

**Step 3: Implement**

Replace the report-writing block at `src/main.rs:407-408`:

```rust
        std::fs::write(path, serde_json::to_string_pretty(&report)? + "\n")
            .with_context(|| format!("writing JSON report to {path}"))?;
```

with:

```rust
        let body = serde_json::to_string_pretty(&report)? + "\n";
        if path == "-" {
            // D-4: `-` means the stream, for the report too.
            use std::io::Write;
            let mut out = std::io::stdout();
            out.write_all(body.as_bytes())
                .and_then(|_| out.flush())
                .context("writing the JSON report to stdout")?;
        } else {
            std::fs::write(path, body)
                .with_context(|| format!("writing JSON report to {path}"))?;
        }
```

**Step 4: Run the tests to verify pass**

Run: `cargo test --offline --test streaming_tests json_report_to_stdout_moves_summary_to_stderr`
Expected: PASS.

**Step 5: Commit**

```bash
git add tools/fastdedup/src/main.rs tools/fastdedup/tests/streaming_tests.rs
git commit -m "feat: --json - writes the report to stdout, summary moves to stderr (D-4)"
```

---

### Task 9: Paired-end serial path streams (stdin R1, stream outputs)

**Objective:** `-1 - -2 r2.fastq -o - -p out2.fastq` works; both PE writers accept `-`.

**Files:**
- Modify: `tools/fastdedup/src/processor.rs:182-183` (readers), `:239-304` (preload guard + writers)

**Step 1: Write the failing test**

Append to `tests/streaming_tests.rs`:

```rust
#[test]
fn paired_stdin_r1_with_file_r2() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let r2 = temp_dir.child("R2.fastq");
    r2.write_str("@p1 2:N:0:ATGC\nTTTTTTTT\n+\nJJJJJJJJ\n").unwrap();
    let out2 = temp_dir.child("out2.fastq");

    let assert = cmd()
        .arg("-1")
        .arg("-")
        .arg("-2")
        .arg(r2.path())
        .arg("-o")
        .arg("-")
        .arg("-p")
        .arg(out2.path())
        .write_stdin(b"@p1 1:N:0:ATGC\nACGTACGT\n+\nIIIIIIII\n")
        .assert()
        .success();

    let out1 = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    assert_eq!(
        out1.lines().filter(|l| l.starts_with('@')).count(),
        1,
        "one unique pair survives"
    );
    let out2_body = fs::read_to_string(out2.path()).unwrap();
    assert_eq!(
        out2_body.lines().filter(|l| l.starts_with('@')).count(),
        1,
        "R2 file must receive the same single pair"
    );
}
```

**Step 2: Run the test to verify failure**

Run: `cargo test --offline --test streaming_tests paired_stdin_r1_with_file_r2`
Expected: FAIL — `parse_fastx_file("-")` and `prepare_writer(Path::new("-"))`.

**Step 3: Implement**

In `src/processor.rs` `execute_paired_deduplication`:

* Replace lines 182-183:
  ```rust
      let mut reader_r1 = parse_fastx_file(input_r1_path).context("Failed to read R1")?;
      let mut reader_r2 = parse_fastx_file(input_r2_path).context("Failed to read R2")?;
  ```
  with:
  ```rust
      let mut reader_r1 = crate::stream::open_fastx(input_r1_path)?;
      let mut reader_r2 = crate::stream::open_fastx(input_r2_path)?;
  ```

* Guard the preload block: replace the `if cfg.force {` at line 245 with the same three-branch form as Task 4 (stream → fresh, no preload; then `else if cfg.force`; then the existing `else` with `preload_existing_paired_hashes` and the R1/R2 truncation logic unchanged).

* Replace the writer preparation at lines 297-300:
  ```rust
      let (writer_r1, output_format_r1) =
          prepare_writer(path_r1_bound, cfg.force, cfg.compression_level)?;
      let (writer_r2, output_format_r2) =
          prepare_writer(path_r2_bound, cfg.force, cfg.compression_level)?;
  ```
  with:
  ```rust
      let (writer_r1, output_format_r1) = open_output(output_r1_path, cfg)?;
      let (writer_r2, output_format_r2) = open_output(output_r2_path, cfg)?;
  ```

(The existing `if output_format_r1 != output_format_r2 { bail!(...) }` check stays: with `--out-format fasta -o - -p out.fastq` it correctly refuses the mixed pair. `path_r1_bound`/`path_r2_bound` are still used by the truncation branch and stay defined.)

**Step 4: Run the tests to verify pass**

Run: `cargo test --offline --test streaming_tests paired_stdin_r1_with_file_r2`
Expected: PASS.

**Step 5: Commit**

```bash
git add tools/fastdedup/src/processor.rs tools/fastdedup/tests/streaming_tests.rs
git commit -m "feat: paired-end serial path reads stdin R1 and writes stream outputs"
```

---

### Task 10: Abundance mode streams (with `--temp-dir`)

**Objective:** `--abundance` needs two passes; with stdin the second pass reads the temp spill, so the stream is consumed exactly once.

**Files:**
- Modify: `tools/fastdedup/src/processor.rs:419` (SE temp-dir pass 1), `:471` (SE temp-dir pass 2 writer)

**Step 1: Write the failing test**

Append to `tests/streaming_tests.rs`:

```rust
#[test]
fn abundance_from_stdin_with_temp_dir_works() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let spill = temp_dir.child("spill");
    spill.create_dir_all().unwrap();

    let assert = cmd()
        .arg("-1")
        .arg("-")
        .arg("--abundance")
        .arg("--temp-dir")
        .arg(spill.path())
        .arg("-o")
        .arg("-")
        .write_stdin(FASTQ)
        .assert()
        .success();

    let out = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    assert!(
        out.contains(";size=2"),
        "the duplicated read must be annotated with its abundance: {out}"
    );
}
```

**Step 2: Run the test to verify failure**

Run: `cargo test --offline --test streaming_tests abundance_from_stdin_with_temp_dir_works`
Expected: FAIL — pass 1 calls `parse_fastx_file("-")`.

**Step 3: Implement**

In `execute_deduplication_with_abundance`, temp-dir branch:

* Line 419: `let mut reader = parse_fastx_file(input_path).context("Failed to read input file")?;` →
  `let mut reader = crate::stream::open_fastx(input_path)?;`
* Line 471 (pass 2 output): replace
  ```rust
          let output_path_bound = Path::new(output_path);
          let (writer, output_format) = prepare_writer(output_path_bound, true, cfg.compression_level)?;
  ```
  with
  ```rust
          let (writer, output_format) = open_output(output_path, cfg)?;
  ```

Leave the classic (no `--temp-dir`) abundance branches on `parse_fastx_file` — Task 2 refuses that combination with stdin, so they are unreachable with `-`.

**Step 4: Run the tests to verify pass**

Run: `cargo test --offline --test streaming_tests abundance_from_stdin_with_temp_dir_works`
Expected: PASS.

**Step 5: Commit**

```bash
git add tools/fastdedup/src/processor.rs tools/fastdedup/tests/streaming_tests.rs
git commit -m "feat: abundance mode reads stdin once when --temp-dir is set"
```

---

### Task 11: Multi-threaded pipeline streams (stdin in, stdout out)

**Objective:** The `-j > 1` producer/consumer path accepts `-` on both ends. Note its output is serial flate2 gzip when `--out-compression gz` (no gzippy thread — stdout cannot be reopened) and no resume.

**Files:**
- Modify: `tools/fastdedup/src/parallel.rs:137-145` (`open_input_reader`), `:175-215` (`create_output_writer`), call sites at `:310-311` and `:519-522`

**Step 1: Write the failing test**

Append to `tests/streaming_tests.rs`:

```rust
#[test]
fn threaded_stdin_to_stdout() {
    let assert = cmd()
        .arg("-1")
        .arg("-")
        .arg("-o")
        .arg("-")
        .arg("-j")
        .arg("4")
        .write_stdin(FASTQ)
        .assert()
        .success();
    let out = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    assert_eq!(
        out.lines().filter(|l| l.starts_with('@')).count(),
        1,
        "the threaded pipeline must read stdin and write stdout identically"
    );
}
```

**Step 2: Run the test to verify failure**

Run: `cargo test --offline --test streaming_tests threaded_stdin_to_stdout`
Expected: FAIL — `open_input_reader` tries `File::open("-")`.

**Step 3: Implement**

In `src/parallel.rs`:

* `open_input_reader` (line 137): add the stream branch first:
  ```rust
  fn open_input_reader(path: &str) -> Result<Box<dyn Read + Send + 'static>> {
      if crate::stream::is_stream(path) {
          // needletail sniffs gzip from magic bytes on any reader, so a
          // stdin stream decodes exactly like a mislabelled file.
          return Ok(Box::new(std::io::stdin()));
      }
      // ... existing body unchanged ...
  }
  ```

* `create_output_writer`: change the signature and add the stream branch:
  ```rust
  fn create_output_writer(
      path: &str,
      compression_level: u32,
      num_threads: usize,
      out_format: OutputFormat,
      out_compression: StreamCompression,
  ) -> Result<(
      Box<dyn Write>,
      OutputFormat,
      Option<std::thread::JoinHandle<Result<()>>>,
  )> {
      if crate::stream::is_stream(path) {
          // D-4: stdout cannot be seeked or reopened, so there is no gzippy
          // thread and no resume — serial flate2 gzip when asked for.
          let writer = crate::writer::prepare_stream_writer(out_compression, compression_level);
          return Ok((writer, out_format, None));
      }
      // ... existing body unchanged ...
  }
  ```
  and add at the top of the file: `use crate::writer::StreamCompression;`

* Update the two call sites:
  ```rust
      let (writer, output_format, gzip_handle) = create_output_writer(
          output_path_str,
          cfg.compression_level,
          gzip_threads,
          cfg.out_format,
          cfg.out_compression,
      )?;
  ```
  (in `execute_deduplication_parallel`) and in `execute_paired_deduplication_parallel`:
  ```rust
      let (writer_r1, output_format, gzip_handle1) = create_output_writer(
          output_r1_path,
          cfg.compression_level,
          gzip_threads,
          cfg.out_format,
          cfg.out_compression,
      )?;
      let (writer_r2, _output_format_r2, gzip_handle2) = create_output_writer(
          output_r2_path,
          cfg.compression_level,
          gzip_threads,
          cfg.out_format,
          cfg.out_compression,
      )?;
  ```

**Step 4: Run the tests to verify pass**

Run: `cargo test --offline --test streaming_tests`
Expected: `test result: ok. 13 passed; 0 failed`.

**Step 5: Commit**

```bash
git add tools/fastdedup/src/parallel.rs tools/fastdedup/tests/streaming_tests.rs
git commit -m "feat: threaded pipeline accepts '-' for stdin and stdout (D-4)"
```

---

### Task 12: Paired-end abundance paths accept `-` too

**Objective:** Close the remaining `parse_fastx_file`/`prepare_writer` call sites so PE abundance streaming does not hit a `-` filename. (Tested via the existing abundance suite staying green plus one PE smoke in Task 15; the unit surface is identical to Tasks 9-10.)

**Files:**
- Modify: `tools/fastdedup/src/processor.rs:638-639` (PE abundance temp-dir pass 1 readers), `:696-697` (pass 2 writers), `:803-804` and `:857-858` (classic path readers), `:838-839` (classic path writers)

**Step 1: Write the failing test**

Run the existing suite first — there is no dedicated new test here; the guard is the full suite plus the manual PE-abundance pipe in Task 15. The failing signal during implementation is a compile/logic error on `-` in these branches.

**Step 2: Verify the gap exists**

Run: `grep -n "parse_fastx_file\|prepare_writer" tools/fastdedup/src/processor.rs`
Expected: the six sites listed in Files above still call the path-only API.

**Step 3: Implement**

Mechanical replacements in `src/processor.rs`:

* Every remaining `let mut reader... = parse_fastx_file(...).context(...)?;` in the **abundance functions** becomes `crate::stream::open_fastx(...)?`:
  - line 638-639 (`reader_r1`, `reader_r2`, temp-dir pass 1)
  - line 709 reads the temp file by path — leave it (it is a real file)
  - lines 803-804, 857-858 (classic path passes 1 and 2)
* Every `prepare_writer(path_..._bound, true, cfg.compression_level)?` in the abundance functions becomes `open_output(output_..._path, cfg)?`:
  - lines 696-697 (PE temp-dir pass 2)
  - lines 838-839 (PE classic pass 2)
  (the single-end abundance writer was converted in Task 10)

**Step 4: Run the tests to verify pass**

Run: `cargo test --offline`
Expected: `test result: ok. 31 passed` (bin), `7 passed` (abundance_tests), `67 passed` (cli + streaming) — 0 failed. If a count differs, an earlier task's expectation was off; the invariant is **0 failed** and the streaming tests all present.

**Step 5: Commit**

```bash
git add tools/fastdedup/src/processor.rs
git commit -m "feat: paired-end abundance paths accept '-' for input and output"
```

---

### Task 13: Update the fastdedup README

**Objective:** The tool's own docs describe the new flags and the streaming rules.

**Files:**
- Modify: `tools/fastdedup/README.md:44-70` (CLI Usage list)

**Step 1: Write the change**

In the `## CLI Usage` bullet list, after the `--output-r2` bullet, add:

```markdown
- `-` in place of `-1/--input`, `-2/--input-r2`, `-o/--output`, `-p/--output-r2`: read from stdin / write to stdout (D-4, standard streams).
- `--out-format <fastq|fasta>`: record format when writing to a stream. Named paths keep inferring it from the extension. Default: `fastq`.
- `--out-compression <none|gz>`: compression when writing to a stream. Named paths keep inferring it from the extension. Default: `none`.
```

And after the flag list, add:

```markdown
### Streaming (zero intermediate files)

```bash
zcat in.fastq.gz | fastdedup -1 - -o - | fastprep -1 - -o clean.fastq.gz
```

- Input compression is detected from **magic bytes**, not the extension, so a stream and a mislabelled file behave identically.
- Output format and compression are explicit (`--out-format`, `--out-compression`): a stream has no name to infer from, and the pipe's bytes must not silently depend on the producer's format.
- Streams are written fresh: no preload, no resume, no truncation. Resume is a file-only feature.
- With `-o -` (or `--json -`) the summary line goes to **stderr**, so stdout carries data only.
- Refusals (exit 1, reason on stderr): both inputs as `-`; both outputs as `-`; `--json -` together with `-o -`; `--pre-count` with stdin; `--abundance` with stdin but no `--temp-dir`.
```

**Step 2: Verify**

Run: `sed -n 44,110p tools/fastdedup/README.md`
Expected: the new bullets and section are present and the markdown fences are balanced.

**Step 3: Commit**

```bash
git add tools/fastdedup/README.md
git commit -m "docs: document streaming mode and the new output flags in the fastdedup README"
```

---

### Task 14: Changelog + settle D-4 in the suite docs

**Objective:** Record the feature and the settled decision so it is not re-litigated.

**Files:**
- Modify: `tools/fastdedup/CHANGELOG.md` (prepend an Unreleased section)
- Modify: `docs/DECISIONS.md:64-69` (D-4 output-compression bullet)
- Modify: `docs/CLI-SPEC.md` (§1 common-core table, §4 gap 7)

**Step 1: Write the changes**

`tools/fastdedup/CHANGELOG.md` — insert above `## v1.3.0`:

```markdown
## Unreleased

### ✨ New Features

- **Standard streams (D-4)**: `-` is accepted for `-1/--input`, `-2/--input-r2`, `-o/--output` and `-p/--output-r2`. `zcat in.fastq.gz | fastdedup -1 - -o - | ...` runs with zero intermediate files.
- `--out-format <fastq|fasta>` and `--out-compression <none|gz>`: explicit format and compression when writing to a stream; named paths keep inferring both from the extension.
- `--json -` writes the report to stdout; the summary line then moves to stderr.
- Input compression is detected from magic bytes, so gzip on a stream is decoded like gzip in a file.
```

`docs/DECISIONS.md` — replace the bullet at lines 64-69:

```markdown
* Output to a stream needs an explicit format/compression decision, because
  there is no filename to infer from: `--out-format fastq|fasta` and
  `--out-compression none|gz` (defaults: fastq, none) — or `-o -` inherits the
  input's compression. To be settled when the first tool implements output
  streaming.
```

with:

```markdown
* Settled 2026-10-12 (fastdedup, the first tool with output streaming): output
  to a stream carries an explicit decision — `--out-format fastq|fasta` and
  `--out-compression none|gz` (defaults: fastq, none). `-o -` does **not**
  inherit the input's compression: a pipe whose bytes silently depend on the
  producer's format is a trap for the next command in the chain. Named paths
  keep inferring both from their extension.
```

`docs/CLI-SPEC.md` — in the §1 common-core table, after the `--gz-level` row:

```markdown
| output format on a stream | `--out-format fastq\|fasta` | enum | `-o -` only; named paths use the extension (D-4) |
| output compression on a stream | `--out-compression none\|gz` | enum | `-o -` only; default `none` (D-4) |
```

and replace gap 7 in §4:

```markdown
7. `-`/stdin-stdout support (streaming mode) — needed to compete with
   `fastp`'s pipe usage; blocked on the reader architecture in fastprep.
```

with:

```markdown
7. `-`/stdin-stdout support (streaming mode) — **fastdedup done** (2026-10-12:
   `-` on `-1/-2/-o/-p`, plus `--out-format`/`--out-compression`); fastprep and
   fastmaps remain, needed to compete with `fastp`'s pipe usage.
```

**Step 2: Verify**

Run: `git diff --stat docs/ tools/fastdedup/CHANGELOG.md`
Expected: three files changed, no source changes in this task.

**Step 3: Commit**

```bash
git add docs/DECISIONS.md docs/CLI-SPEC.md tools/fastdedup/CHANGELOG.md
git commit -m "docs: settle D-4 output compression (explicit flags), mark fastdedup streaming done"
```

---

### Task 15: Final verification — full suite + pipe smoke tests

**Objective:** Prove the whole suite is green and the pipes work end to end with real bytes, from a directory that stays empty.

**Files:** none

**Step 1: Build release**

Run (from `tools/fastdedup/`):
```bash
cargo build --release --offline
```
Expected: `Compiling fastdedup` … `Finished release [optimized] target(s)` , exit 0.

**Step 2: Full test suite**

Run (from `tools/fastdedup/`):
```bash
cargo test --offline 2>&1 | grep "test result"
```
Expected — three lines, zero failed:
```
test result: ok. 31 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out
test result: ok. 7 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out
test result: ok. 67 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out
```
(29+2 unit, 7 abundance, 51+16 integration = 105 total. Any red means stop and fix before continuing.)

**Step 3: The canonical pipe smoke test**

Run (from the repo root):
```bash
printf '@r\nACGT\n+\nIIII\n' | ./tools/fastdedup/target/release/fastdedup -1 - -o - | wc -l
```
Expected: `4` (one record = 4 lines; the summary line is on stderr and does not reach `wc`).

**Step 4: Gzip in, gzip out**

Run (from the repo root):
```bash
printf '@r\nACGT\n+\nIIII\n@r\nACGT\n+\nIIII\n' \
  | gzip -c \
  | ./tools/fastdedup/target/release/fastdedup -1 - -o - --out-compression gz 2>/dev/null \
  | zcat
```
Expected:
```
@r
ACGT
+
IIII
```

**Step 5: Two fastdedups chained, zero intermediate files**

Run (from the repo root):
```bash
mkdir -p /tmp/fd-pipe && cd /tmp/fd-pipe && \
printf '@r\nACGT\n+\nIIII\n@r\nACGT\n+\nIIII\n' \
  | /home/raphael/hermes_scratch/FastSuite/tools/fastdedup/target/release/fastdedup -1 - -o - 2>/dev/null \
  | /home/raphael/hermes_scratch/FastSuite/tools/fastdedup/target/release/fastdedup -1 - -s -v 2>&1 | grep Processed; \
ls /tmp/fd-pipe
```
Expected:
```
Processed fragments: 1
```
and `ls` prints nothing — the directory is empty, proving no intermediate file was written.

**Step 6: Paired-end + abundance manual checks**

Run (from the repo root):
```bash
printf '@p 1:N:0:ATGC\nACGTACGT\n+\nIIIIIIII\n' \
  | ./tools/fastdedup/target/release/fastdedup -1 - -2 R2.fastq -o - -p /tmp/fd-pipe/out2.fastq 2>/dev/null \
  | head -1
printf '@r\nACGTACGT\n+\nIIIIIIII\n@r2\nACGTACGT\n+\nJJJJJJJJ\n' \
  | ./tools/fastdedup/target/release/fastdedup -1 - --abundance --temp-dir /tmp/fd-pipe -o - 2>/dev/null
```
Expected: `@p 1:N:0:ATGC` then `@r;size=2` (with a `R2.fastq` containing a `@p 2:N:0:ATGC` record created first with `printf`).

**Step 7: Final commit if any drift appeared**

```bash
git status --short
git add -A tools/fastdedup/tests/streaming_tests.rs
git commit -m "test: streaming suite complete" # only if Step 1-6 left edits
```

---

## Definition of done

* [ ] `cargo test --offline` from `tools/fastdedup/`: 105 tests, 0 failed.
* [ ] `printf '@r\nACGT\n+\nIIII\n' | ./tools/fastdedup/target/release/fastdedup -1 - -o - | wc -l` prints `4`.
* [ ] The refusals all exit 1 with a reason on stderr: both inputs `-`, both outputs `-`, `--json -` + `-o -`, `--pre-count` + stdin, `--abundance` + stdin without `--temp-dir`.
* [ ] `docs/DECISIONS.md` D-4 records the settled output-compression decision; `docs/CLI-SPEC.md` marks fastdedup's gap 7 done.
* [ ] README documents the flags, the magic-byte rule, the fresh-write rule, and the summary-to-stderr rule.
