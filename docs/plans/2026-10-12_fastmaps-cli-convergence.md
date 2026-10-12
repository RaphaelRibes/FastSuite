# fastmaps: CLI convergence (commit -j, v1 json wrapper, force/quiet)

> **For Hermes:** Use subagent-driven-development skill to implement this plan task-by-task.

**Goal:** Close every remaining `fastmaps` gap in `docs/CLI-SPEC.md` §2/§4 (gaps 3, 4, 5, 6): commit the already-written `-j/--threads`, add `--json <PATH>`, emit the `fastsuite.report.v1` wrapper with the corner metrics `fasthelper` currently derives, add `--force`, `-q/--quiet`, and `FASTSUITE_THREADS` env default, and document `--threads 0 = all cores`.

**Architecture:** Everything is local to `tools/fastmaps` (D-1: no shared crate). CLI surface changes go in `src/cli.rs`; the v1 wrapper is a new pure module `src/report.rs` that takes the already-computed count grids (from `SpillerState`/`write_format`) and emits the documented schema by hand (no serde — `src/encode.rs` already has `json_escape`, and `format.rs` hand-writes JSON today); `src/main.rs` wires the flags. The corner metrics (`polyT_frac`, `polyG_frac`, `r2_polyG_frac`, off-diagonal mass) are computed by a new pure function `corner_metrics()` that reproduces the consumer-side convention in `tools/fasthelper/src/report.rs:197-249` *verbatim* — the tool now owns that maths; `fasthelper` keeps its heatmap.json fallback untouched for one release.

**Tech Stack:** Rust (clap derive, hand-written JSON, no new deps), cargo test, git.

---

## Facts the implementer needs (verified on this checkout)

* 27 tests pass today: `cd tools/fastmaps && cargo test --offline` → `test result: ok. 27 passed; 0 failed`.
* `tools/fastmaps/src/cli.rs` and `src/main.rs` have **uncommitted** modifications adding `-j/--threads`; `tools/fasthelper/src/compare.rs` is also dirty but **out of scope — do not stage it**.
* `--threads 0` currently means auto (all cores, bounded by file count) in `main.rs:115-122`; that matches fastdedup's `resolve_thread_count` (`tools/fastdedup/src/parallel.rs:601`) and the CLI-SPEC §1 rule 3. `fastprep` is the other tool still treating 0 as 1 — not our scope.
* `FASTSUITE_THREADS` is a doc-only target today; no binary reads it yet (grep: 0 hits in `tools/*/src`). fastdedup's pattern to copy is `resolve_thread_count` in `tools/fastdedup/src/parallel.rs:601-609`.
* The v1 wrapper shape (normative copy): `docs/REPORT-SCHEMA.md` and `fasthelper schema`. Reference implementation: `tools/fastcheck/src/report.rs:284-295` (keys `schema/tool/tool_version/command/mode/inputs/elapsed_s/metrics/...`). fasthelper's parser: `tools/fasthelper/src/report.rs:147-183`.
* The corner maths fasthelper derives from `heatmap.json` today, to move into fastmaps (`tools/fasthelper/src/report.rs:197-249`):
  ```rust
  let c = (size / 50).max(2); // corner window
  // per (x,y) with n>0, over the flat row-major grid index y*size+x:
  //   low_x = x <= c; high_x = x >= size-1-c; low_y = y <= c; high_y = y >= size-1-c
  //   low_x && low_y  → polyA; high_x && high_y → polyT
  //   low_x && high_y → polyG; high_x && low_y → polyC
  //   (x - y).abs() > size / 10 → off_diagonal
  // metrics: {pre}polyA_frac/{pre}polyT_frac/{pre}polyG_frac/{pre}polyC_frac/
  //          {pre}off_diagonal_frac with pre "" (all), "r1_", "r2_"; plus total_reads.
  ```
* Verified expected numbers on a hand-written fixture (2 gz files, 4 reads of length 10, size=11, c=2):
  | window | polyA_frac | polyT_frac | polyG_frac | polyC_frac | off_diagonal_frac | total |
  |---|---|---|---|---|---|---|
  | all | 0.25 | 0.25 | 0.5 | 0.0 | 0.5 | 4 |
  | r1 | 0.0 | 0.5 | 0.5 | 0.0 | 0.5 | 2 |
  | r2 | 0.5 | 0.0 | 0.5 | 0.0 | 0.5 | 2 |
  Use exactly this fixture (`/tmp/mapsfix/...`) in the verification steps.
* Output lines today are all `println!` (main.rs:200/233, format.rs:268, output.rs:147/174/246); progress is `eprintln!`. "Quiet the summary line" = suppress the `println!` summary block, keep stderr progress/warnings. Follow `tools/fastcheck`'s minimal reading of `-q` (`tools/fastcheck/src/cli.rs:83-85`).
* fastmaps writes `heatmap.json` unconditionally into `-o` (format.rs:264). `--json` must not change that default: with no `--json`, `heatmap.json` still lands in the output dir; `--json PATH` controls where **both** the grid json and the v1 wrapper go. `--json -` (stdout) is out of scope (D-4, file-only).
* **fastmaps is file-only by design (D-4, CLI-SPEC §1 rule 0 exception).** Do NOT add `-`/stdin support; do not re-propose it. The v1 wrapper goes to a PATH, never stdout.
* **Do NOT touch the spiller** (8 bytes/read temp files, preview payloads in RAM — PLAN-QC §1 problems 3-5). That rewrite belongs to fastcheck. This plan only notes it.
* `--dry-run` for fastmaps: **not in this plan** (CLI-SPEC §4 gap 2 is fastprep's). Add a one-line note to the spec doc instead, if anything.

## Commit conventions

Small, imperative, lower-case scope prefix:
`fastmaps: <what>` (e.g. `fastmaps: add --json <PATH> and v1 report wrapper`). Never `git add -A`; stage only the files each task names.

---

### Task 0: Commit the pending -j/--threads change

**Objective:** Land the already-written, verified worker-thread flag so later tasks build on a clean tree.

**Files:**
- Modify (staged only): `tools/fastmaps/src/cli.rs`, `tools/fastmaps/src/main.rs`

**Step 1: Confirm the working tree contains exactly these two fastmaps files**

Run: `cd ~/hermes_scratch/FastSuite && git status --short`
Expected:
```
 M tools/fasthelper/src/compare.rs
 M tools/fastmaps/src/cli.rs
 M tools/fastmaps/src/main.rs
```
(`compare.rs` must be left unstaged and unmodified.)

**Step 2: Confirm the diff is only the threads flag**

Run: `git diff tools/fastmaps/src/cli.rs tools/fastmaps/src/main.rs | grep -E '^[+-]' | grep -v '^[+-][+-]'`
Expected (additions only, nothing removed): `pub threads: usize` arg block (cli.rs) and the `if args.threads > 0 { ... } else { available_parallelism ... }` worker count (main.rs). If anything else shows up, stop and report — this plan assumes the tree is otherwise clean.

**Step 3: Full test run before committing**

Run: `cd tools/fastmaps && cargo test --offline 2>&1 | tail -3`
Expected: `test result: ok. 27 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out`

**Step 4: Commit**

```bash
cd ~/hermes_scratch/FastSuite
git add tools/fastmaps/src/cli.rs tools/fastmaps/src/main.rs
git commit -m "fastmaps: add -j/--threads (0 = auto, bounded by file count)"
git status --short
```
Expected: only ` M tools/fasthelper/src/compare.rs` remains.

**Step 5: Smoke-check the flag**

Run: `cd tools/fastmaps && cargo build --release --offline 2>&1 | tail -1 && ./target/release/fastmaps --help | grep -A1 threads`
Expected: build `Finished`; help shows `-j, --threads <THREADS>`.

---

### Task 1: Corner metrics pure function (`corner_metrics`)

**Objective:** Move the fasthelper-side corner convention into fastmaps as a pure, testable function.

**Files:**
- Create: `tools/fastmaps/src/report.rs`
- Modify: `tools/fastmaps/src/main.rs:8` (add `mod report;`)

**Step 1: Write the failing test**

Append to the new file `tools/fastmaps/src/report.rs`:

```rust
//! The `fastsuite.report.v1` wrapper and the corner metrics it carries.

use std::collections::BTreeMap;

/// Fraction of reads in each heatmap corner, per mate window.
///
/// Reproduces the convention fasthelper historically derived from
/// `heatmap.json` (tools/fasthelper/src/report.rs): a (T+C, T+G) grid of
/// side `size`, corner window `c = max(size/50, 2)`:
///   low_x/low_y  → polyA corner, high_x/high_y → polyT,
///   low_x/high_y → polyG,      high_x/low_y  → polyC,
/// and off-diagonal mass where `|x - y| > size/10`.
///
/// `grids` is `[("all", &counts), ("r1", &counts_r1), ("r2", &counts_r2)]`
/// with flat row-major index `y*size + x`. Keys use the "" / "r1_" / "r2_"
/// prefixes, exactly as fasthelper's adapter did, so consumers see no
/// difference in metric names.
pub fn corner_metrics(size: u32, grids: &[(&str, &[u64])]) -> BTreeMap<String, f64> {
    let c = (size / 50).max(2) as usize;
    let diag_lim = size / 10;
    let mut out = BTreeMap::new();
    for &(label, counts) in grids {
        let (mut poly_a, mut poly_t, mut poly_g, mut poly_c, mut off_diag) =
            (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let mut total = 0.0f64;
        for y in 0..size as usize {
            for x in 0..size as usize {
                let n = counts[y * size as usize + x] as f64;
                if n == 0.0 {
                    continue;
                }
                total += n;
                let (low_x, high_x) = (x <= c, x + c >= size as usize);
                let (low_y, high_y) = (y <= c, y + c >= size as usize);
                if low_x && low_y { poly_a += n; }
                if high_x && high_y { poly_t += n; }
                if low_x && high_y { poly_g += n; }
                if high_x && low_y { poly_c += n; }
                if (x as i64 - y as i64).abs() as u32 > diag_lim { off_diag += n; }
            }
        }
        if total == 0.0 {
            continue;
        }
        let p = |n: f64| n / total;
        let pre = if label.is_empty() { String::new() } else { format!("{label}_") };
        out.insert(format!("{pre}polyA_frac"), p(poly_a));
        out.insert(format!("{pre}polyT_frac"), p(poly_t));
        out.insert(format!("{pre}polyG_frac"), p(poly_g));
        out.insert(format!("{pre}polyC_frac"), p(poly_c));
        out.insert(format!("{pre}off_diagonal_frac"), p(off_diag));
        if label.is_empty() {
            out.insert("total_reads".into(), total);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// size=11, c=2 fixture: 4 reads of length 10.
    /// R1: TTTTTTTTTT (x=10,y=10) and GGGGGGGGGG (x=0,y=10);
    /// R2: GGGGGGGGGG (x=0,y=10) and AAAAAAAAAA (x=0,y=0).
    fn grids() -> (Vec<u64>, Vec<u64>, Vec<u64>) {
        let size = 11;
        let mut all = vec![0u64; size * size];
        let mut r1 = vec![0u64; size * size];
        let mut r2 = vec![0u64; size * size];
        all[10 * size + 10] += 2; // T at (10,10) x2 (one per file)
        all[10 * size + 0] += 2;  // G at (0,10)  x2
        r1[10 * size + 10] += 1;
        r1[10 * size + 0] += 1;
        r2[10 * size + 10] += 1;
        r2[10 * size + 0] += 1;
        all[0 * size + 0] += 1;   // A at (0,0)
        r2[0 * size + 0] += 1;
        (all, r1, r2)
    }

    #[test]
    fn corner_metrics_matches_fasthelper_fixture() {
        let (all, r1, r2) = grids();
        let m = corner_metrics(11, &[("", &all), ("r1", &r1), ("r2", &r2)]);
        let near = |k: &str, want: f64| {
            assert!((m[k] - want).abs() < 1e-9, "{k}: got {}, want {want}", m[k]);
        };
        near("polyA_frac", 0.25);
        near("polyT_frac", 0.25);
        near("polyG_frac", 0.5);
        near("polyC_frac", 0.0);
        near("off_diagonal_frac", 0.5);
        near("total_reads", 4.0);
        near("r1_polyT_frac", 0.5);
        near("r1_polyG_frac", 0.5);
        near("r1_polyA_frac", 0.0);
        near("r2_polyA_frac", 0.5);
        near("r2_polyT_frac", 0.0);
        near("r2_polyG_frac", 0.5);
    }

    #[test]
    fn corner_metrics_empty_grid_is_omitted() {
        let all = vec![0u64; 121];
        let r1 = vec![0u64; 121];
        let r2 = vec![1u64; 121]; // 121 reads uniform: every corner equal, no off-diag
        let m = corner_metrics(11, &[("", &all), ("r1", &r1), ("r2", &r2)]);
        assert!(!m.contains_key("polyT_frac"), "empty window must be skipped");
        assert!(!m.contains_key("r1_polyT_frac"));
        assert!(m.contains_key("r2_polyT_frac"));
        // A uniform grid: no corner is enriched and no pixel is off-diagonal.
        assert!((m["r2_off_diagonal_frac"] - 0.0).abs() < 1e-9);
    }
}
```

**Step 2: Run the test to verify failure**

Run: `cd tools/fastmaps && cargo test --offline report:: 2>&1 | tail -6`
Expected: FAIL — `error[E0433]: failed to resolve: use of undeclared crate or module 'report'` (the module is not yet declared in main.rs).

**Step 3: Register the module**

Edit `tools/fastmaps/src/main.rs` — add `mod report;` in alphabetical position:

```rust
mod cli;
mod encode;
mod file_kind;
mod format;
mod output;
mod report;   // ← new
mod spiller;
mod varint;
mod worker;
```

**Step 4: Run the tests to verify pass**

Run: `cd tools/fastmaps && cargo test --offline report:: 2>&1 | tail -4`
Expected: `test result: ok. 2 passed; 0 failed` (`corner_metrics_matches_fasthelper_fixture`, `corner_metrics_empty_grid_is_omitted`).

**Step 5: Full suite stays green**

Run: `cargo test --offline 2>&1 | tail -3`
Expected: `test result: ok. 29 passed; 0 failed`.

**Step 6: Commit**

```bash
cd ~/hermes_scratch/FastSuite
git add tools/fastmaps/src/report.rs tools/fastmaps/src/main.rs
git commit -m "fastmaps: add corner_metrics (polyA/T/G/C + off-diagonal per mate)"
```

---

### Task 2: `--json <PATH>` flag (parsing + validation only)

**Objective:** Accept `--json <PATH>` and validate its parent directory, without changing any output yet.

**Files:**
- Modify: `tools/fastmaps/src/cli.rs` (new arg)
- Modify: `tools/fastmaps/src/main.rs` (validation block)

**Step 1: Write the failing test**

Append to `tools/fastmaps/src/cli.rs` inside `mod tests`:

```rust
    #[test]
    fn test_json_arg_default_none() {
        let a = Args::parse_from(["fastmaps", "a.fq"]);
        assert_eq!(a.json, None);
    }

    #[test]
    fn test_json_arg_parsed() {
        let a = Args::parse_from(["fastmaps", "a.fq", "--json", "/tmp/r.json"]);
        assert_eq!(a.json.as_deref(), Some("/tmp/r.json"));
    }
```

**Step 2: Run to verify failure**

Run: `cargo test --offline cli:: 2>&1 | tail -6`
Expected: FAIL — compile error: no field `json` on type `Args`.

**Step 3: Add the flag**

In `tools/fastmaps/src/cli.rs`, add before the `threads` field:

```rust
    /// Write the machine-readable FastSuite report (fastsuite.report.v1) to
    /// PATH, alongside the grid JSON. The grid JSON keeps its default
    /// location `<output-dir>/heatmap.json`; this flag only redirects it
    /// and additionally writes the v1 wrapper. `-` is not accepted:
    /// fastmaps is file-only (D-4).
    #[arg(long, value_name = "PATH")]
    pub json: Option<std::path::PathBuf>,
```

Update the `threads` doc comment to record the env default (Task 5 makes it true):

```rust
    /// Number of worker threads to use. Defaults to the number of available
    /// CPU cores (bounded by the file count), or $FASTSUITE_THREADS when set.
    /// Pass 0 for auto.
```

In `tools/fastmaps/src/main.rs`, after the `--preview-reads` validation (line 31-34), add:

```rust
    // ── Validate --json parent directory ───────────────────────────────────
    if let Some(p) = &args.json {
        if p.as_os_str() == "-" {
            eprintln!("Error: fastmaps is file-only (D-4); --json needs a PATH, not '-'");
            std::process::exit(1);
        }
        let parent = p.parent().filter(|d| !d.as_os_str().is_empty());
        if let Some(dir) = parent {
            if !dir.exists() {
                eprintln!(
                    "Error: --json directory {} does not exist",
                    dir.display()
                );
                std::process::exit(1);
            }
        }
    }
```

**Step 4: Run to verify pass**

Run: `cargo test --offline cli:: 2>&1 | tail -4`
Expected: `test result: ok. 6 passed; 0 failed` (4 existing + 2 new).

**Step 5: Full suite**

Run: `cargo test --offline 2>&1 | tail -3`
Expected: `test result: ok. 29 passed; 0 failed`.

**Step 6: Commit**

```bash
cd ~/hermes_scratch/FastSuite
git add tools/fastmaps/src/cli.rs tools/fastmaps/src/main.rs
git commit -m "fastmaps: accept --json <PATH> (validated, not yet written)"
```

---

### Task 3: Emit the v1 wrapper next to `heatmap.json`

**Objective:** When `--json` is given, write the `fastsuite.report.v1` wrapper (schema/tool/tool_version/command/mode/inputs/elapsed_s/metrics + `grid` block pointer) and put the grid JSON at the `--json` path too.

**Files:**
- Modify: `tools/fastmaps/src/report.rs` (add `write_report_v1`)
- Modify: `tools/fastmaps/src/main.rs` (time the run, pass the path to `write_format`)
- Modify: `tools/fastmaps/src/format.rs:135,264` (`write_format` signature + json path)

**Step 1: Write the failing test**

Append to `tools/fastmaps/src/report.rs` tests:

```rust
    #[test]
    fn report_v1_written_and_parseable() {
        let dir = std::env::temp_dir().join(format!("maps_v1_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("r.json");
        let mut metrics = std::collections::BTreeMap::new();
        metrics.insert("total_reads".into(), 4.0);
        metrics.insert("polyT_frac".into(), 0.25);
        metrics.insert("r2_polyG_frac".into(), 0.5);
        super::write_report_v1(
            &path,
            "fastmaps a_R1.fq.gz a_R2.fastq.gz -j 2",
            0.42,
            &["a_R1.fq.gz".into(), "a_R2.fastq.gz".into()],
            &metrics,
        )
        .unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.contains("\"schema\": \"fastsuite.report.v1\""), "{text}");
        assert!(text.contains("\"tool\": \"fastmaps\""));
        assert!(text.contains("\"mode\": \"composition\""));
        assert!(text.contains("\"elapsed_s\": 0.42"));
        assert!(text.contains("\"polyT_frac\": 0.25"));
        assert!(text.contains("\"r2_polyG_frac\": 0.5"));
        assert!(text.contains("\"a_R1.fq.gz\""));
        let _ = std::fs::remove_dir_all(&dir);
    }
```

**Step 2: Run to verify failure**

Run: `cargo test --offline report::report_v1 2>&1 | tail -6`
Expected: FAIL — no function `write_report_v1` in module `report`.

**Step 3: Implement `write_report_v1`**

Add to `tools/fastmaps/src/report.rs` (above `#[cfg(test)]`), using `crate::encode::json_escape` for all string values:

```rust
/// Write the `fastsuite.report.v1` wrapper (docs/REPORT-SCHEMA.md) to `path`.
///
/// fastmaps' grid JSON (`heatmap.json`) keeps its own shape; this wrapper is
/// the machine-readable layer, carrying the corner metrics the tool now owns
/// (they used to be derived by fasthelper from the grid). The grid file's
/// location is recorded under `grid` so a consumer can find both.
pub fn write_report_v1(
    path: &std::path::Path,
    command: &str,
    elapsed_s: f64,
    files: &[String],
    metrics: &BTreeMap<String, f64>,
) -> std::io::Result<()> {
    use crate::encode::json_escape;

    let metrics_json = metrics
        .iter()
        .map(|(k, v)| format!("    {}: {}", json_escape(k), fmt_f64(*v)))
        .collect::<Vec<_>>()
        .join(",\n");
    let files_json = files
        .iter()
        .map(|f| json_escape(f))
        .collect::<Vec<_>>()
        .join(", ");

    let doc = format!(
        concat!(
            "{{\n",
            "  \"schema\": \"fastsuite.report.v1\",\n",
            "  \"tool\": \"fastmaps\",\n",
            "  \"tool_version\": \"{}\",\n",
            "  \"command\": {},\n",
            "  \"mode\": \"composition\",\n",
            "  \"inputs\": {{\"files\": [{files}]}},\n",
            "  \"elapsed_s\": {},\n",
            "  \"metrics\": {{\n{metrics}\n  }},\n",
            "  \"grid\": {{\"format\": \"heatmap.json\", \"see\": \"docs/PLAN-QC.md\"}}\n",
            "}}\n"
        ),
        env!("CARGO_PKG_VERSION"),
        json_escape(command),
        fmt_f64(elapsed_s),
        files = files_json,
        metrics = metrics_json,
    );

    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    use std::io::Write;
    f.write_all(doc.as_bytes())?;
    f.flush()
}

/// Compact float rendering: integers print without a `.0`, everything else
/// with 6 significant decimals. Keeps the hand-written JSON small.
fn fmt_f64(x: f64) -> String {
    if x == x.trunc() && x.abs() < 1e15 {
        format!("{}", x as i64)
    } else {
        let s = format!("{x:.6}");
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    }
}
```

**Step 4: Run to verify pass**

Run: `cargo test --offline report:: 2>&1 | tail -4`
Expected: `test result: ok. 3 passed; 0 failed`.

**Step 5: Wire it into `write_format` + `main`**

In `tools/fastmaps/src/format.rs`, change the signature and the json write:

```rust
pub fn write_format(
    mut state: SpillerState,
    files: &[String],
    file_kinds: &[FileKind],
    preview_n: usize,
    output_dir: &Path,
    json_path: &Path, // ← new: grid JSON destination (default output_dir/heatmap.json)
) -> std::io::Result<()> {
```

and at the former `heatmap.json` write (line 264):

```rust
    let mut jf = BufWriter::new(File::create(json_path)?);
```

In `tools/fastmaps/src/main.rs`:

1. Time the run — add near the top of `main`:
   ```rust
   use std::time::Instant;
   ...
   let started = Instant::now();
   ```
2. Compute the json path after parsing:
   ```rust
   let json_path = args.json.clone().unwrap_or_else(|| output_dir.join("heatmap.json"));
   ```
   (place it after `output_dir` is created, before the spiller is started)
3. Pass it to `write_format` (Task 5 later adds a `quiet` parameter — update this call site then):
   ```rust
   if let Err(e) = write_format(
       state, &files, &file_kinds, args.preview_reads, &output_dir, &json_path,
       started.elapsed().as_secs_f64(),
   ) {
   ```
4. After `write_format` succeeds (and before `write_png`), emit the wrapper when `--json` was given. The corner grids exist inside `write_format` (lines 198-209), so the wrapper is written there too. Final signature (add `json_path` and `elapsed_s`; `quiet` comes in Task 5):

   ```rust
   pub fn write_format(
       mut state: SpillerState,
       files: &[String],
       file_kinds: &[FileKind],
       preview_n: usize,
       output_dir: &Path,
       json_path: &Path,
       elapsed_s: f64,
   ) -> std::io::Result<()> {
   ```

   Right after the three `counts_arr` grids are built (lines 198-209), compute the metrics:

   ```rust
   // Corner metrics (formerly derived by fasthelper from this very file).
   let metrics = crate::report::corner_metrics(
       size as u32,
       &[("", &counts_arr), ("r1", &counts_r1_arr), ("r2", &counts_r2_arr)],
   );
   ```

   Replace the `heatmap.json` write (line 264) with:

   ```rust
   let mut jf = BufWriter::new(File::create(json_path)?);
   ```

   And at the end of `write_format`, before the `println!`:

   ```rust
   // The v1 wrapper lands next to the grid JSON as <stem>.report.json.
   let stem = json_path
       .file_stem()
       .map(|s| s.to_string_lossy().to_string())
       .unwrap_or_else(|| "heatmap".into());
   let mut wrapper_suffix = String::new();
   if let Some(parent) = json_path.parent() {
       let wrapper = parent.join(format!("{stem}.report.json"));
       crate::report::write_report_v1(
           &wrapper,
           &std::env::args().collect::<Vec<_>>().join(" "),
           elapsed_s,
           files,
           &metrics,
       )?;
       wrapper_suffix = format!(" + {stem}.report.json");
   }
   ```

   And the final summary print becomes:

   ```rust
   println!(
       "Saved heatmap.data ({data_written} B), \
        heatmap.preview ({preview_written} B), \
        heatmap.json (size={size}, {} reads){}",
       state.total_reads,
       wrapper_suffix,
   );
   ```

   In `tools/fastmaps/src/main.rs`:

   1. Time the run — add `use std::time::Instant;` to the imports and `let started = Instant::now();` as the first line of `fn main`.
   2. After `output_dir` is created, compute the grid path:
      ```rust
      let json_path = args.json.clone().unwrap_or_else(|| output_dir.join("heatmap.json"));
      ```
   3. Pass it and the elapsed time to `write_format`:
      ```rust
      if let Err(e) = write_format(
          state, &files, &file_kinds, args.preview_reads, &output_dir, &json_path,
          started.elapsed().as_secs_f64(),
      ) {
      ```

**Step 6: Full suite + smoke run**

Run: `cargo test --offline 2>&1 | tail -3`
Expected: `test result: ok. 30 passed; 0 failed`.

Smoke run with the verified fixture:

```bash
cargo build --release --offline
printf '@r1\nTTTTTTTTTT\n+\nIIIIIIIIII\n@r2\nGGGGGGGGGG\n+\nIIIIIIIIII\n' | gzip > /tmp/mapsfix/sampleA_R1.fastq.gz
printf '@r1\nGGGGGGGGGG\n+\nIIIIIIIIII\n@r2\nAAAAAAAAAA\n+\nIIIIIIIIII\n' | gzip > /tmp/mapsfix/sampleA_R2.fastq.gz
rm -rf /tmp/mapsfix/out /tmp/m.json /tmp/m.report.json
./target/release/fastmaps -j2 -o /tmp/mapsfix/out /tmp/mapsfix/sampleA_R1.fastq.gz /tmp/mapsfix/sampleA_R2.fastq.gz --json /tmp/m.json
ls -la /tmp/m.json /tmp/m.report.json
python3 -c "import json; d=json.load(open('/tmp/m.report.json')); m=d['metrics']; print(d['schema'], d['tool'], d['mode']); print({k:round(v,4) for k,v in m.items() if 'frac' in k or k=='total_reads'})"
```
Expected: both files exist; the wrapper prints
`fastsuite.report.v1 fastmaps composition`
and metrics `{'polyA_frac': 0.25, 'polyT_frac': 0.25, 'polyG_frac': 0.5, 'polyC_frac': 0.0, 'off_diagonal_frac': 0.5, 'r1_polyA_frac': 0.0, 'r1_polyT_frac': 0.5, 'r1_polyG_frac': 0.5, 'r1_polyC_frac': 0.0, 'r1_off_diagonal_frac': 0.5, 'r2_polyA_frac': 0.5, 'r2_polyT_frac': 0.0, 'r2_polyG_frac': 0.5, 'r2_polyC_frac': 0.0, 'r2_off_diagonal_frac': 0.5, 'total_reads': 4}` (order per BTreeMap; values must match). Also verify `python3 -c "import json; json.load(open('/tmp/m.json'))['total_reads']"` → `4`, and that **no** `--json` still writes `heatmap.json` into `-o`:
```bash
rm -rf /tmp/mapsfix/out2
./target/release/fastmaps -j2 -o /tmp/mapsfix/out2 /tmp/mapsfix/sampleA_R1.fastq.gz /tmp/mapsfix/sampleA_R2.fastq.gz
ls /tmp/mapsfix/out2/heatmap.json && test ! -f /tmp/mapsfix/out2/heatmap.report.json
```
Expected: `heatmap.json` exists, `heatmap.report.json` does not (no wrapper without `--json`).

**Step 7: Commit**

```bash
cd ~/hermes_scratch/FastSuite
git add tools/fastmaps/src/report.rs tools/fastmaps/src/format.rs tools/fastmaps/src/main.rs
git commit -m "fastmaps: --json <PATH> emits heatmap.json + fastsuite.report.v1 wrapper"
```

---

### Task 4: `--force` (refuse to clobber the output dir)

**Objective:** Without `--force`, exit 1 when the output directory already contains fastmaps outputs; with it, overwrite as today.

**Files:**
- Modify: `tools/fastmaps/src/cli.rs` (flag)
- Modify: `tools/fastmaps/src/main.rs` (guard)

**Step 1: Write the failing test**

Append to `tools/fastmaps/src/cli.rs` tests:

```rust
    #[test]
    fn test_force_default_false() {
        let a = Args::parse_from(["fastmaps", "a.fq"]);
        assert!(!a.force);
    }

    #[test]
    fn test_force_short_and_long() {
        assert!(Args::parse_from(["fastmaps", "a.fq", "--force"]).force);
        assert!(Args::parse_from(["fastmaps", "a.fq", "-f"]).force);
    }
```

**Step 2: Run to verify failure**

Run: `cargo test --offline cli::test_force 2>&1 | tail -6`
Expected: FAIL — no field `force` on `Args`.

**Step 3: Implement**

In `tools/fastmaps/src/cli.rs`:

```rust
    /// Refuse to overwrite an output directory that already contains
    /// fastmaps outputs; without this flag the run exits 1 naming the file.
    #[arg(long, short = 'f')]
    pub force: bool,
```

In `tools/fastmaps/src/main.rs`, replace the output-dir block (lines 36-41) with:

```rust
    // ── Output directory ──────────────────────────────────────────────────
    let output_dir = args.output_dir;
    if !args.force && output_dir.exists() {
        let existing = ["heatmap.json", "heatmap.png", "heatmap.data", "heatmap_3d.html"]
            .iter()
            .map(|f| output_dir.join(f))
            .find(|p| p.exists());
        if let Some(p) = existing {
            eprintln!(
                "Error: {} already exists in {}. Pass --force to overwrite.",
                p.file_name().unwrap().to_string_lossy(),
                output_dir.display()
            );
            std::process::exit(1);
        }
    }
    std::fs::create_dir_all(&output_dir).unwrap_or_else(|e| {
        eprintln!("Cannot create output dir {}: {e}", output_dir.display());
        std::process::exit(1);
    });
```

**Step 4: Run to verify pass**

Run: `cargo test --offline cli:: 2>&1 | tail -4`
Expected: `test result: ok. 8 passed; 0 failed`.

**Step 5: End-to-end verify (both branches)**

```bash
cd tools/fastmaps && cargo build --release --offline
rm -rf /tmp/mapsfix/out3
./target/release/fastmaps -j2 -o /tmp/mapsfix/out3 /tmp/mapsfix/sampleA_R1.fastq.gz /tmp/mapsfix/sampleA_R2.fastq.gz >/dev/null
# second run must refuse:
./target/release/fastmaps -j2 -o /tmp/mapsfix/out3 /tmp/mapsfix/sampleA_R1.fastq.gz /tmp/mapsfix/sampleA_R2.fastq.gz; echo "exit=$?"
./target/release/fastmaps -j2 -o /tmp/mapsfix/out3 /tmp/mapsfix/sampleA_R1.fastq.gz /tmp/mapsfix/sampleA_R2.fastq.gz --force >/dev/null; echo "exit=$?"
```
Expected: first refusal prints `Error: heatmap.json already exists in /tmp/mapsfix/out3. Pass --force to overwrite.` and `exit=1`; the `--force` run prints `exit=0`.

**Step 6: Full suite + commit**

Run: `cargo test --offline 2>&1 | tail -3` → `test result: ok. 29 passed; 0 failed`.

```bash
cd ~/hermes_scratch/FastSuite
git add tools/fastmaps/src/cli.rs tools/fastmaps/src/main.rs
git commit -m "fastmaps: --force/-f refuses to clobber an existing output dir"
```

---

### Task 5: `-q/--quiet`

**Objective:** Suppress the summary `println!` lines (CLI-SPEC §1 rule 5: the single `<tool>: ` summary line). Progress on stderr stays.

**Files:**
- Modify: `tools/fastmaps/src/cli.rs` (flag)
- Modify: `tools/fastmaps/src/main.rs` (gate the three stdout summary prints: lines 200, 233)
- Modify: `tools/fastmaps/src/format.rs:268` (heatmap.json line — takes a `quiet` param)
- Modify: `tools/fastmaps/src/output.rs:147,174,246` (png/html/svg lines)

**Step 1: Write the failing test**

Append to `tools/fastmaps/src/cli.rs` tests:

```rust
    #[test]
    fn test_quiet_default_false() {
        assert!(!Args::parse_from(["fastmaps", "a.fq"]).quiet);
    }

    #[test]
    fn test_quiet_short_and_long() {
        assert!(Args::parse_from(["fastmaps", "a.fq", "-q"]).quiet);
        assert!(Args::parse_from(["fastmaps", "a.fq", "--quiet"]).quiet);
    }
```

**Step 2: Run to verify failure**

Run: `cargo test --offline cli::test_quiet 2>&1 | tail -6`
Expected: FAIL — no field `quiet`.

**Step 3: Implement**

`tools/fastmaps/src/cli.rs`:

```rust
    /// Suppress the end-of-run summary lines on stdout.
    /// Progress and warnings still go to stderr.
    #[arg(long, short = 'q')]
    pub quiet: bool,
```

`tools/fastmaps/src/main.rs` — gate the two stdout prints:

```rust
    if !args.quiet {
        println!(
            "{} unique coordinates, {} total reads",
            state.counts.len(),
            state.total_reads
        );
    }
```

```rust
        if !args.quiet {
            println!("Saved {fname}  ({n} reads extracted)");
        }
```

Pass quiet down to the writers — extend each signature minimally:

```rust
// format.rs
pub fn write_format(
    ...
    quiet: bool,
) -> std::io::Result<()> {
```
```rust
    if !quiet {
        println!(
            "Saved heatmap.data ({data_written} B), \
             heatmap.preview ({preview_written} B), \
             heatmap.json (size={size}, {} reads){}",
            state.total_reads,
            wrapper_suffix,
        );
    }
```
```rust
// output.rs
pub fn write_png(counts: &HashMap<(u32, u32), (u64, f64)>, output_dir: &Path, quiet: bool)
pub fn write_html(counts: &HashMap<(u32, u32), (u64, f64)>, output_dir: &Path, quiet: bool)
pub fn write_svg(counts: &HashMap<(u32, u32), (u64, f64)>, output_dir: &Path, quiet: bool)
```
and in each, wrap the final `println!` in `if !quiet { ... }`.

Update the four call sites in `main.rs` to pass `args.quiet`.

**Step 4: Run to verify pass**

Run: `cargo test --offline cli:: 2>&1 | tail -4`
Expected: `test result: ok. 10 passed; 0 failed`.

**Step 5: End-to-end verify**

```bash
cd tools/fastmaps && cargo build --release --offline
rm -rf /tmp/mapsfix/out4
./target/release/fastmaps -j2 -o /tmp/mapsfix/out4 /tmp/mapsfix/sampleA_R1.fastq.gz /tmp/mapsfix/sampleA_R2.fastq.gz --json /tmp/m4.json -q > /tmp/mapsfix/stdout.txt
wc -l < /tmp/mapsfix/stdout.txt   # stdout must be empty
ls /tmp/m4.json /tmp/m4.report.json  # outputs still written
```
Expected: `0`, both json files exist.

**Step 6: Full suite + commit**

Run: `cargo test --offline 2>&1 | tail -3` → `test result: ok. 29 passed; 0 failed`.

```bash
cd ~/hermes_scratch/FastSuite
git add tools/fastmaps/src/cli.rs tools/fastmaps/src/main.rs tools/fastmaps/src/format.rs tools/fastmaps/src/output.rs
git commit -m "fastmaps: -q/--quiet suppresses the summary lines"
```

---

### Task 6: `FASTSUITE_THREADS` env default + `--threads 0` semantics

**Objective:** `-j` absent (or 0) → `$FASTSUITE_THREADS` if set and valid, else all cores; still bounded by file count. Mirrors fastdedup's `resolve_thread_count`.

**Files:**
- Modify: `tools/fastmaps/src/cli.rs` (add `resolved_threads`)
- Modify: `tools/fastmaps/src/main.rs:112-122` (use it)

**Step 1: Write the failing test**

Append to `tools/fastmaps/src/cli.rs` tests:

```rust
    #[test]
    fn test_resolved_threads_prefers_flag() {
        let a = Args::parse_from(["fastmaps", "a.fq", "b.fq", "-j", "3"]);
        assert_eq!(a.resolved_threads(2), 3);
    }

    #[test]
    fn test_resolved_threads_flag_bounded_by_files() {
        let a = Args::parse_from(["fastmaps", "a.fq", "b.fq", "-j", "8"]);
        assert_eq!(a.resolved_threads(2), 2); // caller clamps by file count
    }

    #[test]
    fn test_resolved_threads_zero_means_auto() {
        let a = Args::parse_from(["fastmaps", "a.fq", "-j", "0"]);
        // No env in tests: falls back to available cores, clamped by file count.
        let n = a.resolved_threads(1);
        assert!(n >= 1 && n <= 16, "got {n}");
    }
```

**Step 2: Run to verify failure**

Run: `cargo test --offline cli::test_resolved 2>&1 | tail -6`
Expected: FAIL — no method `resolved_threads`.

**Step 3: Implement**

In `tools/fastmaps/src/cli.rs`, after the `Args` struct:

```rust
impl Args {
    /// Effective worker count: `-j N` wins; `-j 0`/absent reads
    /// `$FASTSUITE_THREADS` (ignored when unset/unparseable/0), else all
    /// available cores. `max` clamps the result (fastmaps bounds workers by
    /// the input file count). Mirrors fastdedup's `resolve_thread_count`.
    pub fn resolved_threads(&self, max: usize) -> usize {
        let base = if self.threads > 0 {
            self.threads
        } else {
            std::env::var("FASTSUITE_THREADS")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .filter(|&n| n > 0)
                .unwrap_or_else(|| {
                    std::thread::available_parallelism()
                        .map(|n| n.get())
                        .unwrap_or(4)
                })
        };
        base.min(max.max(1))
    }
}
```

In `tools/fastmaps/src/main.rs`, replace lines 112-122:

```rust
    // ── Worker configuration ──────────────────────────────────────────────
    // -j/--threads overrides everything; 0 (or absent) reads
    // $FASTSUITE_THREADS, else all cores. Either way the pool never
    // exceeds the file count (no point in workers with no files).
    let num_workers = args.resolved_threads(total_files);
```

**Step 4: Run to verify pass**

Run: `cargo test --offline cli:: 2>&1 | tail -4`
Expected: `test result: ok. 13 passed; 0 failed`.

**Step 5: End-to-end verify**

```bash
cd tools/fastmaps && cargo build --release --offline
rm -rf /tmp/mapsfix/out5
FASTSUITE_THREADS=3 ./target/release/fastmaps -o /tmp/mapsfix/out5 /tmp/mapsfix/sampleA_R1.fastq.gz /tmp/mapsfix/sampleA_R2.fastq.gz -v 2>&1 | head -2
FASTSUITE_THREADS=bogus ./target/release/fastmaps -o /tmp/mapsfix/out5 --force /tmp/mapsfix/sampleA_R1.fastq.gz /tmp/mapsfix/sampleA_R2.fastq.gz 2>&1 | head -2
```
Expected: both runs succeed (exit 0); the bogus value is silently ignored (no crash, no warning needed — spec says "read when the flag is absent", and an unparseable value is treated as absent).

Direct check that the env actually drives the count (instrument via a debug print is out of scope; instead assert via the resolved_threads unit test, which covers flag/env precedence — the `env::var` path is exercised by the `-j 0` test only when no env is set, which is the test environment).

**Step 6: Full suite + commit**

Run: `cargo test --offline 2>&1 | tail -3` → `test result: ok. 32 passed; 0 failed`.

```bash
cd ~/hermes_scratch/FastSuite
git add tools/fastmaps/src/cli.rs tools/fastmaps/src/main.rs
git commit -m "fastmaps: FASTSUITE_THREADS default, --threads 0 = all cores"
```

---

### Task 7: Documentation — CLI-SPEC fastmaps table + D-4 + spiller note

**Objective:** Make the spec true and record the two decisions that keep getting re-proposed.

**Files:**
- Modify: `docs/CLI-SPEC.md` (§2 fastmaps table, §4 gap list)
- Modify: `docs/DECISIONS.md` (D-4 file-only reinforcement — only if the fastmaps line is missing; it currently exists at D-4, so this file needs **no change** — verify, don't edit)

**Step 1: Update the fastmaps conformance table (§2, lines 108-117)**

Replace the `--json`, `-j --threads`, `--force / -q` rows with:

```markdown
| `--json <PATH>` | ✓ | grid JSON + `fastsuite.report.v1` wrapper with corner metrics |
| `-j --threads` | ✓ | `0` = all cores (bounded by file count); `$FASTSUITE_THREADS` default |
| `--force` / `-f` | ✓ | refuses to clobber an existing output dir |
| `-q --quiet` | ✓ | suppresses the summary lines on stdout |
```

**Step 2: Update the gap list (§4)**

Change items 3, 4, 5, 6 to reflect completion and renumber nothing (the numbering is referenced elsewhere):

```markdown
3. ~~`--json <PATH>` and `-j/--threads` in **fastmaps**.~~ ✅ done (2026-10-12).
4. ~~`--threads 0 = all cores` in **fastprep** and **fastmaps**.~~ fastmaps ✅ done
   (2026-10-12); fastprep still treats 0 as 1.
5. ~~`--quiet` in **fastprep**/**fastmaps**.~~ fastmaps ✅ done (2026-10-12);
   fastprep still open.
6. ~~`FASTSUITE_THREADS` honoured as a default by all four binaries.~~ fastmaps ✅
   done (2026-10-12); fastdedup/fastprep/fastsim still open.
```

**Step 3: Add the two standing notes**

At the end of §2 (after the fastcheck block, before the file's end), append:

```markdown
### fastmaps: standing constraints

* **File-only (D-4).** fastmaps takes positional `FILE...` and classifies
  R1/R2 by filename. Do not add `-`/stdin support: a stream has no name, so a
  pipe would remove information rather than save I/O. `--json` likewise needs a
  PATH.
* **Spiller is deliberately not converged.** fastmaps spills 8 bytes per read
  to per-row temp files and keeps preview read payloads in RAM (PLAN-QC §1,
  problems 3-5). That is fastcheck's dense-grid QC architecture to replace —
  do not rewrite it under a "CLI convergence" plan.
```

**Step 4: Verify DECISIONS.md needs no edit**

Run: `grep -n "fastmaps.*file-only" docs/DECISIONS.md`
Expected: one hit at D-4 (line ~72). If present, leave the file untouched.

**Step 5: Commit**

```bash
cd ~/hermes_scratch/FastSuite
git add docs/CLI-SPEC.md
git commit -m "docs: fastmaps CLI conformance — json/force/quiet/threads/env all done"
```

---

### Task 8: Final verification (full suite + release smoke run)

**Objective:** Prove the whole plan works end-to-end on one command block.

**Step 1: Full offline test suite**

Run:
```bash
cd ~/hermes_scratch/FastSuite/tools/fastmaps && cargo test --offline 2>&1 | tail -3
```
Expected: `test result: ok. 32 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out`

**Step 2: Release build + canonical smoke run**

```bash
cargo build --release --offline
rm -rf /tmp/maps /tmp/maps.json /tmp/maps.report.json
printf '@r1\nTTTTTTTTTT\n+\nIIIIIIIIII\n@r2\nGGGGGGGGGG\n+\nIIIIIIIIII\n' | gzip > /tmp/sampleA_R1.fastq.gz
printf '@r1\nGGGGGGGGGG\n+\nIIIIIIIIII\n@r2\nAAAAAAAAAA\n+\nIIIIIIIIII\n' | gzip > /tmp/sampleA_R2.fastq.gz
./target/release/fastmaps -j2 -o /tmp/maps /tmp/sampleA_R1.fastq.gz /tmp/sampleA_R2.fastq.gz --json /tmp/maps.json
test -f /tmp/maps/heatmap.json && echo GRID-OK
test -f /tmp/maps.json && echo GRID-PATH-OK
test -f /tmp/maps.report.json && echo V1-OK
python3 -c "
import json
d = json.load(open('/tmp/maps.report.json'))
m = d['metrics']
assert d['schema'] == 'fastsuite.report.v1' and d['tool'] == 'fastmaps'
assert abs(m['polyT_frac'] - 0.25) < 1e-6
assert abs(m['polyG_frac'] - 0.5)  < 1e-6
assert abs(m['r2_polyG_frac'] - 0.5) < 1e-6
assert abs(m['off_diagonal_frac'] - 0.5) < 1e-6
assert m['total_reads'] == 4
print('V1-METRICS-OK', {k: round(v,4) for k,v in m.items() if 'frac' in k or k=='total_reads'})
"
```
Expected, in order: `GRID-OK`, `GRID-PATH-OK`, `V1-OK`, and a single
`V1-METRICS-OK {...}` line with the values from the facts table (all:
polyA .25, polyT .25, polyG .5, polyC 0, off_diag .5; r1: polyT .5, polyG .5;
r2: polyA .5, polyG .5; total_reads 4).

**Step 3: Regression guard — fasthelper still reads both shapes**

```bash
cd ~/hermes_scratch/FastSuite/tools/fasthelper && cargo build --release --offline 2>&1 | tail -1
./target/release/fasthelper concatenate /tmp/maps.json -o /tmp/fh_out --formats tsv 2>&1 | tail -2
grep -c polyT_frac /tmp/fh_out/report.tsv
```
Expected: build ok; the tsv contains the `polyT_frac` column (fasthelper reads the
wrapper's flat metrics; its heatmap.json fallback still works unchanged for old
grid files).

**Step 4: Clean up and report**

```bash
rm -rf /tmp/maps /tmp/maps.json /tmp/maps.report.json /tmp/fh_out /tmp/sampleA_R1.fastq.gz /tmp/sampleA_R2.fastq.gz /tmp/mapsfix
cd ~/hermes_scratch/FastSuite && git log --oneline -8 && git status --short
```
Expected: 8 commits (Task 0 through 7); working tree still shows only
` M tools/fasthelper/src/compare.rs` (pre-existing, untouched by this plan).
