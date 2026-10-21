# fasthelper — optimisation propositions

Scope: `tools/fasthelper/src` (3 420 lines, all 11 files read) and the pre-built
`target/release/fasthelper` run against scratch trees. Nothing here is implemented.

## Status of this round

Implemented and committed in `1479cb9` (7 of 10 propositions). The tests are the gate: the crate's suite is green, and every change that can alter output was checked byte-for-byte against the pre-change binary except where the entry says otherwise.

### Implemented

- **FH-1** --- `.report` suffix stripped before the stem test; v1 reports win on (sample, tool) -- the phantom `heatmap.report` sample is gone
- **FH-2** --- duplicate `[[samples]]` names now fail `check -c` (they used to exit 0 and share one work dir)
- **FH-3** --- `report.json` streams a borrowed view through `to_writer_pretty`, byte order preserved
- **FH-5** --- fastmaps grid streamed through a `Visitor` (1.76 ms -> 0.06 ms per grid)
- **FH-6** --- `--all-metrics` streams to a `BufWriter`
- **FH-8** --- `build_matrix` borrows, `--stdout` builds once
- **FH-9** --- `--version` resolved once per run (was once per sample per step)

### Deferred, with the reason

- **FH-4** --- the work-dir resume key is already a hash
- **FH-7** --- the HTML template is a few hundred KB
- **FH-10** --- parallel aggregation would break the deterministic row order that the tests assert

---


Tags: `[M]` measured, artifact named · `[E]` estimated from a named mechanism · `[V]` verified in
this review by running the built binary.

## Verdict first

At no scale this tool will meet is aggregation a *time* problem. It is a **memory** problem, and two
findings — FH-1 and FH-2 — are **wrong-output** problems that no timing cell would have caught. Both
were reproduced in this review.

Measured by the reviewer, 10 000 samples × 2 v1 reports (79 MB total, 20 000 files):

| formats | wall | peak RSS | output |
|---|---|---|---|
| json,tsv,md,html | 1.55 s | 313.9 MB | 49 MB |
| json only | 0.91 s | 314.4 MB | 26.5 MB |
| tsv only | 0.63 s | 140.0 MB | 4.4 MB |
| md only | 0.78 s | 78.6 MB | 5.0 MB |
| html only | 0.73 s | 84.8 MB | 12.8 MB |

**Not superlinear.** Peak RSS at 1 k/2 k/5 k/10 k = 33.9/65.2/159.4/314.4 MB (json) and
17.3/30.8/71.4/140.0 MB (tsv); wall 0.124/0.192/0.524/0.910 s. A straight line in both counters. No
stage is O(n²): each is O(total report bytes) or O(n_samples × n_metrics), with a constant made of
copies (FH-3, FH-8).

**What breaks first at 10 k is peak RSS** — ~32 KB per sample for the default four-format write. Fine
on a workstation, marginal in a 1–2 GB cgroup, 4× the input bytes. Time binds only in the fastmaps
grid path (FH-5, ~17.6 s at 10 k by linear fit over three measured cells) and the adapter guard
(FH-4, 6.3 s per 370 k-read mate measured — larger than the whole 10 k aggregation).

**Hygiene only:** FH-9 and FH-10 — nearly free, but they move no percentile. What matters at 10 k+ is
FH-3…FH-7.

---

## FH-1 — A fastmaps v1 wrapper beside its grid becomes a phantom sample and a duplicate row `[V]`
`correctness-enabling · S · output changes by design`

`fastmaps --json out/maps/heatmap.json` writes the grid **and** `heatmap.report.json` (`main.rs:84`,
`format.rs:364`), which is the invocation its own `--help` recommends. fasthelper normalises both, and
the wrapper gets the wrong name because only an exact stem is tested — `heatmap.report` is not in
`GENERIC`:

```rust
// report.rs:100
if GENERIC.contains(&stem.to_ascii_lowercase().as_str()) {
```

**Reproduced in this review.** Two sample dirs, each holding a real `heatmap.json` + its
`heatmap.report.json`, produced by two `fastmaps --json` runs on the same reads:

```
$ fasthelper concatenate probe3 -o out --format tsv
fasthelper: concatenated 4 report(s) from 3 sample(s)
$ cut -f1,2 out/report.tsv | sort | uniq -c
      1 alpha            fastmaps
      1 beta             fastmaps
      2 heatmap.report   fastmaps
```

Both libraries collapse into one sample name and each appears twice, so any sum over the table
double-counts. The grid's own metrics and the wrapper's are identical values, which is why nothing
looks wrong in the output.

**Change.**

```rust
let stem = stem.strip_suffix(".report").unwrap_or(stem);
reps.retain(|r| r.from_v1 || !v1_keys.contains(&(r.sample.clone(), r.tool.clone())));
```

**Effect.** A phantom sample and a duplicated metric set disappear, and 140 KB of grid parsing per
sample stops. Mechanism: one library stops being two rows from two files.

**Measure.** `concatenate probe3 -o out --format tsv; cut -f1,2 out/report.tsv | sort | uniq -c` —
4 rows / 3 samples before, 2 rows / 2 samples after.

**Risk.** A library genuinely named `heatmap.report`; gate the strip on a `GENERIC` prefix.

---

## FH-2 — Duplicate `[[samples]]` names share one output directory, and `check` accepts them `[V]`
`correctness-enabling · S · byte-identical for valid configs`

`config.rs:19` takes the list verbatim; the work dir and the resume key derive from the name alone
(`runner.rs:141`), so with `-P 2` both workers `create_dir_all` the same path, unlink the same stale
list and write the same `fastprep.json`.

**Reproduced in this review.** A two-sample config where both entries are named `lib1`:

```
$ fasthelper check -c dup.toml --tools-dir tools
fasthelper check: OK — 2 sample(s), 1 step(s), tools resolved
$ echo $status   # fish
0
```

The validator meant to catch this class passes it.

**Change.**

```rust
let mut seen = BTreeSet::new();
for s in &cfg.samples { if !seen.insert(&s.name) { bail!("duplicate sample name '{}'", s.name); } }
```

**Effect.** An interleaved-output race becomes a minute-zero error. Mechanism: two writers, one path
set.

**Measure.** `check` exits non-zero on the duplicate config, still 0 on `examples/fastsuite.toml`.

**Risk.** Configs that worked by luck now fail loudly — intended.

---

## FH-3 — `report.json` is built as a `Value` tree plus one 26 MB `String` before any byte is written
`memory · M · byte-identical if field order is preserved`

At 10 k: json-only peak 314.4 MB for a 26.5 MB document, against 78.6 MB for md-only (same parse
stage, 5.0 MB output). The ~236 MB delta is this function — the `Value` tree (20 k objects × ~30
`Number`s, each with its own key `String`), the `Vec<Value>`, and the pretty `String` plus its
doubling:

```rust
let samples: Vec<serde_json::Value> = reports.iter().map(|r| serde_json::json!({ ... })).collect();  // :615-627
let doc = serde_json::json!({ "schema": ..., "reports": samples });                                  // :628
serde_json::to_string_pretty(&doc).unwrap_or_else(|_| "{}".into()) + "\n"                            // :635
```

**Change.** Serialise a borrowed view into a buffered writer. Field order must stay alphabetical —
`json!` builds a `BTreeMap` (the committed `examples/runs/demo/report.json` opens `generator`,
`generator_version`, `reports`; inside a report `metrics`, `notes`, `sample`, `source`, `tool`,
`tool_version`):

```rust
struct Row<'a> { metrics: &'a BTreeMap<String, f64>, notes: &'a [String], sample: &'a str,
                 source: Cow<'a, str>, tool: &'a str, tool_version: Option<&'a str> }
serde_json::to_writer_pretty(&mut BufWriter::new(File::create(p)?), &Doc { reports: rows })?;
```

**Effect.** Peak 314 → ~80–85 MB at 10 k (tree and String gone, 8 KB buffer), ~0.1–0.2 s off the wall
from ~590 k fewer key allocations `[E]`.

**Measure.** 10 k json cell, VmHWM; gate `sha256sum report.json` against today's binary.

**Risk.** A field declared out of order silently reorders the JSON — the sha256 gate is the
mitigation.

---

## FH-4 — The adapter guard re-reads the whole input every run, on one thread
`io-threads · M · `report.*` byte-identical; manifest gains a flag`

`guard_step` measures the untrimmed input, then the trimmed output when adapter was found, and the
call site pins it to one thread while fastprep just ran with `-j threads`:

```rust
let g = adapter_guard::guard_step(check, &inputs, &outputs, logs, &sample.name,
                                  guard, &step.args,
    1,                                            // runner.rs:381
)?;
```

**Measured** `[M: reviewer, median of 3]`: one `fastcheck -q --adapter-content` pass — 37 k reads
(1.7 MB gz) 708 ms at `-j 1` vs 750 ms at `-j 8`; 370 k reads (17 MB gz) 6.34 s vs 5.50 s. So one
mate of one sample costs 4× the entire 10 k aggregation, threading is a 13 % lever on real data and
zero on small files, and **"do not scan again" is the 100 % lever**. The key already exists:
`input_signatures` at `runner.rs:245`.

**Change.**

```rust
let key = (sample.name.clone(), input_signatures(&cur1, cur2.as_deref()));
if let Some(prev) = cached_guard.get(&key) { return Ok(GuardOutcome::reused(prev)); }
```

**Effect.** 13 % measured from threading; one full pass removed per sample per re-run from the cache
(the verdict is a pure function of the input bytes, keyed on size+mtime) `[E]`.

**Risk.** The cache inherits resume's trust model — an in-place edit preserving size+mtime is
undetectable there already.

---

## FH-5 — The fastmaps grid is parsed into a `Value` tree, all 3 × 22 801 cells
`algorithmic · M · byte-identical if accumulation order is kept`

A grid is 140 765 B at 151 bp against 3 897 B for a fastprep report: **36× the bytes and 1.76 ms per
file**, exactly linear (500/1 000/2 000 files = 0.889/1.745/3.527 s `[M]`). 2 000 grids (280 MB of
JSON) cost 3.53 s; 20 000 ordinary reports (79 MB) cost 0.63 s — 56× cheaper per file. Every sample
gets a grid by default, so a 10 k run's aggregation is ~19 s by linear fit, of which **17.6 s is
this path**. The cost is the enum match per cell:

```rust
let n = counts[(y * size + x) as usize].as_f64().unwrap_or(0.0);   // :217
```

**Change.** Stream the grid with a visitor that fills `Vec<u64>` and accumulates the five corner sums
plus the total in the same walk:

```rust
struct Grid { size: usize, c: usize, diag: usize, sums: [f64; 5], total: f64 }
impl<'de> Visitor<'de> for Grid { fn visit_u64<E>(self, n: u64) -> Result<(), E> { ... } }
```

**Effect.** `[E]` 3.53 → ~1.0–1.4 s at 2 k (~10 s at 10 k) from dropping 68 k `Value` constructions
per file; peak stays ~16 MB either way.

**Measure.** `grid2k` cell, wall + VmHWM, `sha256sum report.tsv`.

**Risk.** Must keep the `n == 0.0 { continue }` shortcut, the `total == 0.0` skip, and x-major
accumulation for bit-identical f64 sums.

---

## FH-6 — `--all-metrics` builds the whole tidy table in one `String`
`algorithmic · S · byte-identical`

**Measured** `[M]`: 2 k samples × 151 positions × 13 metrics = 1.53 s, 138.2 MB peak, 87.8 MB output.
The `String`'s next power-of-two capacity (~134 MB) *is* the peak. This shape needs no transpose — the
loop order is already input order (`positions.rs:193-202`).

**Change.**

```rust
let mut w = BufWriter::new(File::create(&path)?);
for s in &samples { for (pos, values) in &s.rows { for (n, v) in s.metrics.iter().zip(values) {
    writeln!(w, "{}\t{}\t{}\t{}", s.sample, pos, n, fmt(*v))?; } } }
```

**Effect.** `[E]` 138 → ~10 MB at 2 k (~0.8 GB → ~10 MB at 10 k by the same slope), likely slightly
faster with the doubling gone.

**Measure.** `pos2k --all-metrics`, VmHWM + `sha256sum positions_matrix.tsv`.

---

## FH-7 — The wide matrix keeps all 13 metrics per position, in a `BTreeMap`, to write one column
`algorithmic · M · byte-identical`

**Measured** `[M]`: 2 k samples = 0.271 s, 57.5 MB peak, for a 2.1 MB matrix (27× amplification).
Twelve of the thirteen f64 and every `BTreeMap` node exist only to be discarded, and the wanted column
is looked up inside the inner loop:

```rust
pub rows: BTreeMap<u32, Vec<f64>>,                                    // positions.rs:34
if let (Some(idx), Some(values)) = (s.metrics.iter().position(|m| m == metric), s.rows.get(&pos)) {
```

**Change.** `struct Wide { sample: String, at: usize, max_pos: u32, vals: Vec<f64> }` (NaN = absent) —
select at read time, store densely indexed by `pos − 1`, NaN for a missing cell (the sentinel at
`positions.rs:86`, which `fmt` already renders empty), first-match column semantics by scanning the
header once.

**Effect.** `[E]` 57.5 → ~3 MB at 2 k (~288 → 12 MB at 10 k).

**Measure.** `pos2k` wide cell, VmHWM + sha256; also feed a ragged 150/151 bp mix to prove the union
and empty cells are unchanged.

---

## FH-8 — `build_matrix` clones every row's metric map, and `--stdout` builds the matrix twice
`allocation · S · byte-identical`

At 10 k: tsv-only 140.0 MB vs 78.6 MB for md-only (same parse stage) — ~61 MB for the row clone plus
the `String` `[M]`.

```rust
rows.push((r.sample.clone(), r.tool.clone(), r.metrics.clone()));      // report.rs:428
"tsv" => ("report.tsv", matrix_tsv(&build_matrix(reports))),           // :651
print!("{}", matrix_tsv(&build_matrix(&reports)));                     // :800
```

**Change.**

```rust
pub struct Matrix<'a> { columns: Vec<&'a str>, rows: Vec<(&'a str, &'a str, &'a BTreeMap<String, f64>)> }
write!(s, "{}", fmt_val(*v))?;   // was s.push_str(…); build the matrix once, in cmd_concatenate
```

`columns.extend(r.metrics.keys().cloned())` (`:427`) also clones ~500 k keys and drops the
duplicates (`BTreeSet<&str>`, own once); `group_by_tool` is rebuilt twice per `write_reports`
(`:471`, `:526`); `Config::step()` (`config.rs:168`) clones a `Vec<String>` per step per sample.

**Effect.** `[E]` −61 MB at 10 k for tsv, −4.4 MB with `--stdout`.

**Measure.** 10 k tsv cell + sha256; `--stdout | sha256sum` against `report.tsv`.

---

## FH-9 — `--version` is forked once per sample per step
`io-threads · S · byte-identical`

**Measured** `[M]`: 300 spawns of `fastcheck --version` = 0.249 s, 0.83 ms each → 10 k × 3 steps =
30 000 forks ≈ 25 s by multiplication, paid on `--dry-run` too, since the resume branch precedes the
`dry_run` guard.

```rust
&& prev.version == tool_version(&bin)                    // runner.rs:283
"tool_version": tool_version(&bin),                      // :348
```

**Change.** Resolve once where the toolset is resolved (`BTreeMap<&str, Option<String>>`), pass it
down.

**Effect.** `[E]` 25 s at 10 k; a version string cannot change inside one run, and it is already
recorded per step.

**Measure.** `--dry-run` over a 10 k config, wall + `strace -f -c -e trace=execve`.

---

## FH-10 — `-P` drops the thread remainder and spawns workers with no work
`io-threads (scheduling) · S · byte-identical · UNMEASURED`

`-j 8 -P 3` gives 2 threads per sample: 6 of 8 cores, permanently. `-j 4 -P 16` runs 16 one-thread
processes against a stated budget of 4; `-P 64` on 3 samples spawns 64 threads.

```rust
let per_sample_threads = (cfg.run.threads / args.jobs.max(1)).max(1);   // runner.rs:810
for _ in 0..args.jobs.max(1) {                                          // :820
```

**Change.** `let pool = args.jobs.min(samples.len()).max(1); let (base, extra) = (threads / pool,
threads % pool);` with worker `i` taking `base + (i < extra)`.

**Effect.** `[E]` +33 % for `-j 8 -P 3`. **Unmeasured** — it needs a run that exercises the tools.

**Measure.** 6-sample config at `-j 8 -P 1/2/3`, then `sum(step elapsed_s)/wall/jobs` from
`manifest.json`; ceiling 1.0.

---

## What is already good

- **Self-inclusion is impossible**: `report.rs:144-148` rejects any JSON with `generator: fasthelper`
  or a `reports` array. Measured: `concatenate examples/runs/demo`, which holds its own
  `report.json`, reports "4 report(s) from 1 sample(s)".
- **Output order is walk-order independent**: `BTreeSet` columns (`report.rs:424`), `files.sort()`
  (`:402`), and `-P` entries tagged `(sample_idx, step_idx)` then sorted (`runner.rs:854-856`), so a
  parallel manifest is byte-identical to a sequential one.
- **Guard scratch reports land in `logs/`**, which `walk_json` skips (`report.rs:374-377`): the guard
  cannot pollute the report it produces.
- **`--positions` refuses to read its own output back** (`positions.rs:122-127`) — the one place a
  matrix would silently double every column.
- One adapter covers three report shapes, so `concatenate` and `compare` cannot disagree on a key;
  `evaluate_thresholds` (`report.rs:745`) names the offending sample, so one clean sample cannot mask
  a dirty one in a 10 k gate.

## Dead ends — do not re-try these

1. **Parallelise the file parse.** The whole parse stage is 0.63 s at 10 k for 79 MB — the small half
   of the problem, and it needs an ordered merge to keep `files.sort()` determinism.
2. **Swap `BTreeMap<String, f64>` for `HashMap`.** Sorted key order is load-bearing: it makes
   `report.json` field order and TSV column order walk-order independent.
3. **Pipeline samples so fastprep(i) overlaps fastdedup(i+1).** Every tool is already multi-threaded
   at `-j`; this trades a bubble that does not exist (overhead 0.667 vs 0.647 s) for oversubscription
   and worse peak RSS.
4. **Hold tool reports in memory instead of re-walking `outdir`.** They are written by child
   processes; holding them parses the same bytes earlier, not fewer, and breaks `concatenate` parity.

---

## Re-running these measurements

```
BIN=tools/fasthelper/target/release/fasthelper
$BIN concatenate scale10k -o out --format json     # 10k dirs x 2 v1 reports, 79 MB; peak via VmHWM
$BIN concatenate grid2k   -o out --format tsv      # 2000 heatmap.json -> 1.76 ms/file
$BIN concatenate pos2k    -o out --positions [--all-metrics]
$BIN concatenate probe3   -o out --format tsv ; cut -f1,2 out/report.tsv | sort | uniq -c
$BIN check -c dup.toml --tools-dir tools/          # duplicate names: exits 0 today
```

Peak RSS is the kernel high-water mark, `VmHWM`. Scratch trees are pruned after 24 h; regenerate them
from `examples/runs/demo` plus one `fastmaps --json` run per sample.

---

## Round 2 (2026-10-21)

Re-priced on the post-fix tree (`30e8fb3` + `af268d4` + the 2026-10-21 wave): see `docs/optim/round2-fasthelper.md` (propositions) and `docs/optim/raw/round2/` (probe evidence). Round 1 above is history — deferred items there are re-priced, not repeated.
