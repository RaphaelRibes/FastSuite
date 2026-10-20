# fasthelper — optimization dossier

Tool `tools/fasthelper` 0.1.0. Read-only review: all 11 src files (3,420 lines: runner 1034, report 817,
adapter_guard 567, positions 277, config 199, cli 187, compare 104, check 101, install 77, main 42, lib
15). No src change, no build; numbers are from the pre-built `tools/fasthelper/target/release/fasthelper`
run against scratch trees in `/home/raphael/.hermes/cache/scratch/fh/` (re-run commands at the end).

## Verdict first

At no scale this tool will meet is aggregation a *time* problem. It is a **memory** problem, and two
findings (FH-1, FH-2) are **wrong-output** problems no timing cell would have caught.

Measured, 10k samples x 2 v1 reports (79 MB, 20,000 files), wall / peak RSS / output:

| formats | wall | peak RSS | output |
|---|---|---|---|
| json,tsv,md,html | 1.55 s | 313.9 MB | 49 MB |
| json | 0.91 s | 314.4 MB | 26.5 MB |
| tsv | 0.63 s | 140.0 MB | 4.4 MB |
| md | 0.78 s | 78.6 MB | 5.0 MB |
| html | 0.73 s | 84.8 MB | 12.8 MB |

**Is the 800-sample cell superlinear? No.** Peak RSS at 1k/2k/5k/10k = 33.9/65.2/159.4/314.4 MB (json)
and 17.3/30.8/71.4/140.0 MB (tsv); wall 0.124/0.192/0.524/0.910 s — a straight line in both counters. No
stage is O(n^2): each is O(total report bytes) or O(n_samples x n_metrics), with a large constant made of
copies (FH-3, FH-8). The published 800-sample cell (0.041 s / 3.89 MB, `docs/BENCHMARKS.md:404`) is a
different tree — mine carries 2 reports/sample with 49 and 12 metrics, 3.9 KB each — so neither
extrapolates to the other; use the 10k cell above.

**What breaks first at 10k:** peak RSS, ~32 KB per sample for the default 4-format write — fine on a
workstation, marginal in a 1-2 GB cgroup, 4x the input bytes. Time binds only in the fastmaps grid path
(FH-5, 17.6 s at 10k by linear fit over three measured cells) and the adapter guard (FH-4, 6.3 s per
370k-read mate measured — larger than the whole 10k aggregation).

**Hygiene only:** FH-9 (25 s of `--version` forks at 10k, `run` path only) and FH-10 (thread division,
unmeasured) — nearly free, but they move no percentile. What matters at 10k+ is FH-3..FH-7.

## What is already good

- Self-inclusion is impossible: `report.rs:144-148` rejects any JSON with `generator: fasthelper` or a
  `reports` array. Measured: `concatenate examples/runs/demo`, which holds its own `report.json`, reports
  "4 report(s) from 1 sample(s)".
- Output order is walk-order independent: `BTreeSet` columns (`report.rs:424`), `files.sort()`
  (`report.rs:402`), and `-P` entries tagged `(sample_idx, step_idx)` then sorted (`runner.rs:854-856`),
  so a parallel manifest is byte-identical to a sequential one.
- Guard scratch reports land in `logs/`, which `walk_json` skips (`report.rs:374-377`): the guard cannot
  pollute the report it produces.
- `--positions` refuses to read its own output back (`positions.rs:122-127`) — the one place a matrix
  would silently double every column. One adapter covers three report shapes, so `concatenate` and
  `compare` cannot disagree on a key; `evaluate_thresholds` (`report.rs:745`) names the offending sample,
  so one clean sample cannot mask a dirty one in a 10k gate.

## Ranked propositions — wrong output, then peak memory at 10k, then wall time, then hygiene

### FH-1 — A fastmaps v1 wrapper beside its grid becomes a phantom sample and a duplicate row
`class: correctness-enabling | effort: S | risk: medium | byte-identical: NO (by design)`

`fastmaps --json out/maps/heatmap.json` writes the grid *and* `heatmap.report.json` (`main.rs:84`,
`format.rs:364`); fasthelper normalises both, and the wrapper gets the wrong name because only an exact
stem is tested — `heatmap.report` is not in `GENERIC` (`report.rs:90-94`) — so the walk-up never runs:

```rust
    if GENERIC.contains(&stem.to_ascii_lowercase().as_str()) {   // report.rs:100
```

Measured, 2 samples x {grid, wrapper}: "concatenated 4 report(s) from 3 sample(s)" with rows `alpha`,
`heatmap.report`, `beta`, `heatmap.report` — both libraries collapse to one name, and the two rows per
library carry identical values (`grid_size 151`, `off_diagonal_frac 0.06940`), so sums double-count.

```rust
let stem = stem.strip_suffix(".report").unwrap_or(stem);
reps.retain(|r| r.from_v1 || !v1_keys.contains(&(r.sample.clone(), r.tool.clone())));
```

Effect: a phantom sample and a duplicated metric set disappear, 140 KB of grid parsing per sample stops
(mechanism: one library stops being two rows from two files). Measure: `concatenate probe2 -o out
--format tsv; cut -f1,2 out/report.tsv | sort | uniq -c` — 4 rows / 3 samples before, 2 rows / 2 after.
Risk: a library really named `heatmap.report`; gate the strip on a `GENERIC` prefix.

### FH-2 — Duplicate `[[samples]]` names share one output directory and one manifest key
`class: correctness-enabling | effort: S | risk: low | byte-identical: yes for valid configs`

`config.rs:19` takes the list verbatim; the work dir and the resume key derive from the name alone, so with
`-P 2` both workers `create_dir_all` it, unlink the same `stale` list (`runner.rs:325-329`) and write the
same `fastprep.json`:

```rust
    let work = outdir.join(&sample.name);          // runner.rs:141
```

Measured: the validator meant to catch this class accepts it — `check -c dup.toml` prints
"OK — 2 sample(s), 1 step(s), tools resolved", exit 0, for two entries both named `lib1`.

```rust
let mut seen = BTreeSet::new();
for s in &cfg.samples { if !seen.insert(&s.name) { bail!("duplicate sample name '{}'", s.name); } }
```

Effect: an interleaved-output race becomes a minute-zero error (mechanism: two writers, one path set).
Measure: `check` exits non-zero on `dup.toml`, still 0 on `examples/fastsuite.toml`. Risk: configs that
worked by luck now fail loudly — intended.

### FH-3 — `report.json` is built as a Value tree plus one 26 MB String before any byte is written
`class: memory | effort: M | risk: medium | byte-identical: yes if field order is preserved`

Measured at 10k: json-only peak 314.4 MB for a 26.5 MB document, against 78.6 MB for md-only (same parse
stage, 5.0 MB output). The ~236 MB delta is this function — the Value tree (20k objects x ~30 `Number`s,
each with its own key `String`), the `Vec<Value>`, and the pretty String plus its doubling:

```rust
    let samples: Vec<serde_json::Value> = reports.iter().map(|r| serde_json::json!({ ... })).collect();  // :615-627
    let doc = serde_json::json!({ "schema": ..., "reports": samples });                // :628
    serde_json::to_string_pretty(&doc).unwrap_or_else(|_| "{}".into()) + "\n"          // :635
```

Change: serialize a borrowed view into a buffered writer. Field order must stay alphabetical — `json!`
builds a `BTreeMap` (committed `examples/runs/demo/report.json` opens `generator`, `generator_version`,
`reports`; inside a report `metrics`, `notes`, `sample`, `source`, `tool`, `tool_version`):

```rust
struct Row<'a> { metrics: &'a BTreeMap<String, f64>, notes: &'a [String], sample: &'a str,
                 source: Cow<'a, str>, tool: &'a str, tool_version: Option<&'a str> }
serde_json::to_writer_pretty(&mut BufWriter::new(File::create(p)?), &Doc { reports: rows })?;
```

Effect: peak 314 -> ~80-85 MB at 10k (tree and String gone, 8 KB buffer), ~0.1-0.2 s off the wall from
~590k fewer key allocations. Measure: 10k json cell, VmHWM; gate `sha256sum report.json` against today's
binary. Risk: a field declared out of order silently reorders the JSON — the sha256 gate is the mitigation.

### FH-4 — The adapter guard re-reads the whole input every run, on one thread
`class: io-threads | effort: M | risk: medium | byte-identical: yes (report.*); manifest gains a flag`

`guard_step` measures the untrimmed input, then the trimmed output when adapter was found, and the call
site pins it to one thread while fastprep just ran with `-j threads`:

```rust
                let g = adapter_guard::guard_step(check, &inputs, &outputs, logs, &sample.name,
                                                 guard, &step.args,
                    1,                                            // runner.rs:381
                )?;
```

Measured, one `fastcheck -q --adapter-content` pass, median of 3: 37k reads (1.7 MB gz) 708 ms at `-j 1`
vs 750 ms at `-j 8`; 370k reads (17 MB gz) 6.34 s vs 5.50 s. So one mate of one sample costs 4x the entire
10k aggregation (1.55 s), threading is a 13% lever on real data and zero on small files, and "do not scan
again" is the 100% lever. The key already exists: `input_signatures` at `runner.rs:245`.

```rust
let key = (sample.name.clone(), input_signatures(&cur1, cur2.as_deref()));
if let Some(prev) = cached_guard.get(&key) { return Ok(GuardOutcome::reused(prev)); }
```

Effect: 13% measured from threading; estimated one full pass removed per sample per re-run from the cache
(the verdict is a pure function of the input bytes, keyed on size+mtime). The 10k extrapolation from
6.34 s per 370k reads is an estimate, not a measurement. Measure: a fixed 8-sample config run twice
(second with `--from fastprep`), compare per-step `elapsed_s` in `manifest.json`; `strace -c -f -e
trace=read`. Risk: the cache inherits resume's trust model — an in-place edit preserving size+mtime is
undetectable there already.

### FH-5 — The fastmaps grid is parsed into a Value tree, all 3 x 22,801 cells
`class: algorithmic | effort: M | risk: medium | byte-identical: yes if accumulation order is kept`

A grid is 140,765 B at 151 bp (`examples/runs/demo/sim_E/maps/heatmap.json`) against 3,897 B for a
fastprep report: 36x the bytes and 1.76 ms per file, exactly linear (both `--format tsv`) — 500/1000/2000
files = 0.889/1.745/3.527 s. 2,000 grids (280 MB of JSON) cost 3.53 s; 20,000 ordinary reports (79 MB)
cost 0.63 s, 56x cheaper per file. Every sample gets a grid by default, so a 10k run's aggregation is
~19 s by linear fit, of which 17.6 s is this path — 92% of the stage. The cost is the enum match per cell:

```rust
                    let n = counts[(y * size + x) as usize].as_f64().unwrap_or(0.0);   // :217
```

Change: stream the grid with a visitor that fills `Vec<u64>` and accumulates the five corner sums plus the
total in the same walk:

```rust
struct Grid { size: usize, c: usize, diag: usize, sums: [f64; 5], total: f64 }
impl<'de> Visitor<'de> for Grid {              // visit_seq over counts / counts_r1 / counts_r2
    fn visit_u64<E>(self, n: u64) -> Result<(), E> { self.push(n as f64); Ok(()) } }
```

Effect: estimated 3.53 -> ~1.0-1.4 s at 2k (~10 s at 10k) from dropping 68k `Value` constructions per file;
peak stays ~16 MB either way. Measure: `grid2k` cell, wall + VmHWM, `sha256sum report.tsv`. Risk: must keep
the `n == 0.0 { continue }` shortcut, the `total == 0.0` skip, and x-major accumulation for bit-identical
f64 sums.

### FH-6 — `--all-metrics` builds the whole tidy table in one String
`class: algorithmic | effort: S | risk: low | byte-identical: yes`

Measured, 2k samples x 151 positions x 13 metrics: 1.53 s, 138.2 MB peak, 87.8 MB output. The String's
next power-of-two capacity (~134 MB) *is* the peak. This shape needs no transpose — the loop order is
already input order (`positions.rs:193-202`, `out.push_str` per cell):

```rust
let mut w = BufWriter::new(File::create(&path)?);
for s in &samples { for (pos, values) in &s.rows { for (n, v) in s.metrics.iter().zip(values) {
    writeln!(w, "{}\t{}\t{}\t{}", s.sample, pos, n, fmt(*v))?; } } }
```

Effect: estimated 138 -> ~10 MB at 2k (~0.8 GB -> ~10 MB at 10k by the same slope), and likely slightly
faster with the doubling gone. Measure: `pos2k --all-metrics`, VmHWM + `sha256sum positions_matrix.tsv`.

### FH-7 — The wide matrix keeps all 13 metrics per position, in a `BTreeMap`, to write one column
`class: algorithmic | effort: M | risk: low-medium | byte-identical: yes`

Measured, 2k samples: 0.271 s, 57.5 MB peak, for a 2.1 MB matrix (27x amplification). Twelve of the
thirteen f64 and every BTreeMap node exist only to be discarded, and the wanted column is looked up
inside the inner loop:

```rust
    pub rows: BTreeMap<u32, Vec<f64>>,                                    // positions.rs:34
    if let (Some(idx), Some(values)) = (s.metrics.iter().position(|m| m == metric), s.rows.get(&pos)) {  // :237-240
```

Change: `struct Wide { sample: String, at: usize, max_pos: u32, vals: Vec<f64> }` (NaN = absent) — select
at read time, store densely indexed by `pos - 1`, NaN for a missing cell (the sentinel at
`positions.rs:86`, which `fmt` already renders empty), first-match column semantics by scanning the header
once. Effect: estimated 57.5 -> ~3 MB at 2k (~288 -> 12 MB at 10k). Measure: `pos2k` wide cell, VmHWM +
sha256; also feed a ragged 150/151 bp mix to prove the union and empty cells are unchanged.

### FH-8 — `build_matrix` clones every row's metric map, and `--stdout` builds the matrix twice
`class: allocation | effort: S | risk: low | byte-identical: yes`

Measured at 10k: tsv-only 140.0 MB vs 78.6 MB for md-only (same parse stage, comparable output) — ~61 MB
for the row clone plus the String:

```rust
        rows.push((r.sample.clone(), r.tool.clone(), r.metrics.clone()));      // report.rs:428
            "tsv" => ("report.tsv", matrix_tsv(&build_matrix(reports))),      // :651
        print!("{}", matrix_tsv(&build_matrix(&reports)));                    // :800
```

```rust
pub struct Matrix<'a> { columns: Vec<&'a str>, rows: Vec<(&'a str, &'a str, &'a BTreeMap<String, f64>)> }
write!(s, "{}", fmt_val(*v))?;      // was s.push_str(&fmt_val(*v));  build_matrix once, in cmd_concatenate
```

Same functions, same pass: `columns.extend(r.metrics.keys().cloned())` (`:427`) clones ~500k keys and
drops the duplicates (`BTreeSet<&str>`, own once); `group_by_tool` is rebuilt twice per `write_reports`
(`:471`, `:526`); `Config::step()` (`config.rs:168`) clones a `Vec<String>` per step per sample. Effect:
estimated -61 MB at 10k for tsv, -4.4 MB with `--stdout`. Measure: n10000 tsv cell + sha256; `--stdout |
sha256sum` against `report.tsv`.

### FH-9 — `--version` is forked once per sample per step
`class: io-threads | effort: S | risk: low | byte-identical: yes`

Measured: 300 spawns of `fastcheck --version` = 0.249 s, 0.83 ms each, so 10k x 3 steps = 30,000 forks
~ 25 s by multiplication — paid on `--dry-run` too, since the resume branch precedes the `dry_run` guard:

```rust
                    && prev.version == tool_version(&bin)                    // runner.rs:283
                "tool_version": tool_version(&bin),                       // :348
```

Change: resolve once where the toolset is resolved (`BTreeMap<&str, Option<String>>`), pass it down.
Effect: estimated 25 s at 10k (a version string cannot change inside one run; it is already recorded per
step). Measure: `--dry-run` over a 10k config, wall + `strace -f -c -e trace=execve`.

### FH-10 — `-P` drops the thread remainder and spawns workers with no work
`class: io-threads (scheduling) | effort: S | risk: low | byte-identical: yes`

`-j 8 -P 3` gives 2 threads per sample: 6 of 8 cores, permanently. `-j 4 -P 16` runs 16 one-thread
processes against a stated budget of 4; `-P 64` on 3 samples spawns 64 threads:

```rust
    let per_sample_threads = (cfg.run.threads / args.jobs.max(1)).max(1);   // runner.rs:810
        for _ in 0..args.jobs.max(1) {                                       // :820
```

Change: `let pool = args.jobs.min(samples.len()).max(1); let (base, extra) = (threads / pool, threads %
pool);` with worker `i` taking `base + (i < extra)`. Effect: estimated +33% for `-j 8 -P 3` (the tools are
CPU-bound at `-j threads`; two idle cores is a 25% budget loss). **Unmeasured** — it needs a run that
exercises the tools. Measure: 6-sample config at `-j 8 -P 1/2/3`, then `sum(step elapsed_s)/wall/jobs`
from `manifest.json`; ceiling 1.0.

## Dead ends / explicitly rejected

1. **Parallelise the file parse.** The whole parse stage is 0.63 s at 10k for 79 MB (measured, tsv cell) —
   the small half of the problem, and it needs an ordered merge to keep `files.sort()` determinism.
2. **Swap `BTreeMap<String, f64>` for `HashMap`.** Sorted key order is load-bearing: it makes
   `report.json` field order and TSV column order walk-order independent (`report.rs:424`, `:448`).
3. **Pipeline samples so fastprep(i) overlaps fastdedup(i+1).** Every tool is already multi-threaded at
   `-j`; this trades a bubble that does not exist (overhead 0.667 vs 0.647 s, `docs/BENCHMARKS.md:386`) for
   oversubscription and worse peak RSS.
4. **Hold tool reports in memory instead of re-walking `outdir`.** They are written by child processes;
   holding them parses the same bytes earlier, not fewer, and breaks `concatenate`-parity.

## Re-running these measurements

```
BIN=tools/fasthelper/target/release/fasthelper
$BIN concatenate scale10k -o out --format json     # 10k dirs x 2 v1 reports, 79 MB; peak via VmHWM
$BIN concatenate grid2k    -o out --format tsv     # 2000 hardlinked heatmap.json -> 1.76 ms/file
$BIN concatenate pos2k     -o out --positions [--all-metrics]
$BIN concatenate probe2    -o out --format tsv ; cut -f1,2 out/report.tsv | sort | uniq -c
$BIN check -c dup.toml --tools-dir tools/          # duplicate sample names: exits 0 today
```

Scratch trees and outputs live in `/home/raphael/.hermes/cache/scratch/fh/` (pruned after 24 h idle;
regenerate from `examples/runs/demo`). Peak RSS is the kernel high-water mark, `VmHWM`.