# fastcheck 0.1.0 — ranked optimization dossier

All of `tools/fastcheck/src/*.rs` (8321 lines) was read; nothing was built, changed or committed.
Numbers are tagged **[measured]** (artifact named), **[counted]** (counted/derived from code constants,
arithmetic shown) or **[estimated]** (mechanism stated, no artifact). Cells: `C1M`/`C5M`/`C20M`
(single-end 151 bp), threads pinned with `taskset`, `-j 1|2|4` only.

## Baseline, and where the time is

- **[measured]** `docs/BENCHMARKS.md:277-282`: 20 M reads = 55.68 s, 138.8 MB peak — **2784 ns/read,
  54.2 Mbp/s** (derived), 3.35x slower than Falco's 16.63 s.
- **[counted]** the default path walks the sequence **four times**: composition (`accum.rs:218`), heatmap
  (`accum.rs:324`), xxh3 (`accum.rs:331`), 12-mer windows (`accum.rs:352`) — the last touches each base
  **11.1x** (140 windows x 12 B = 1680 decodes for 151 bp). Plus one quality pass (`accum.rs:258`).
- **[counted]** 20 M x 140 = **2.8e9** 12-mer inserts, each a `HashMap<u64,u32>` probe with the **default
  SipHash** hasher (`repeated.rs:15,25`), into a table that saturates at **26.0 MiB** (`memory::project(1e6,
  12)` = 2^21 buckets x 13 B). 2.8e9 probes is the same order as the whole measured runtime.
- **[counted]** `qual_pos` = 151 rows x 2048 B = **302 KiB** (`accum.rs:76`); a read touches 151 distinct rows,
  one line each, so ~151 L2 accesses per read by construction. `base_pos` = 151 x 40 B.

## FC-1 — Exact 12-mer counts in 32 MiB: a dense counter replaces the bounded hash map
`class: algorithmic · effort: M · output: k-mer block changes (approximate -> exact); all other blocks identical`

```rust
// accum.rs:352-357
if seq.len() >= 12 {
    for w in seq.windows(12) {
        if let Some(pk) = repeated::pack_kmer(w) {
            self.kmer_counts.insert(pk);
```

The counter is `std::collections::HashMap<u64,u32>` (`repeated.rs:25`) capped at 1e6 keys with amortized
singleton eviction (`repeated.rs:56-63`); the report labels the result approximate, with its failure mode
(`rich.rs:146-152`). The bound for exactness comes from the code's own constants: k is fixed at 12
(`accum.rs:353`, `rich.rs:137`), so the key space is 4^12 = **16,777,216** [counted]. Dense `Vec<u16>` =
**32 MiB** [counted], input-size independent, exact below 65,536; hot 12-mers above that (adapters, polyA)
overflow into a small `FxHashMap<u32,u32>` that keeps counting, so the count stays exact with only the
*width* bounded.

```rust
struct KmerTable { c: Vec<u16>, hot: FxHashMap<u32, u32> } // 32 MiB + O(hot)
fn add(&mut self, k: u32) {
    let v = &mut self.c[k as usize];
    if *v == u16::MAX { *self.hot.entry(k).or_insert(u16::MAX as u32) += 1 } else { *v += 1 }
}
```

- Effect: removes SipHash and the probe (one indexed RMW) and removes eviction, so the spectrum is exact and
  `distinct` (today unknowable from a capped map) becomes countable from non-zero cells. **[counted]**
  memory is *below* the map it replaces at saturation (32 MiB vs 26.0 MiB plus 1.25e6 entries of slack);
  **[estimated]** net RSS +6 to +10 MiB, still flat. `Vec<u8>` = 16 MiB with the same rule is also exact, but
  the mean count per 12-mer on `C20M` is ~167 [counted: 2.8e9/16.78 M], so a u8 table would route millions
  of keys into the overflow map. u16 is the width that holds.
- Measure: `C5M --only kmer -j 1` wall before/after; VmHWM on `C1M` and `C20M` (must stay flat); `perf stat
  -e cache-misses,LLC-load-misses`; correctness by extending `tests/reference_qc.py` with a dense 4^12 count
  on a fixture, and flipping `approximate` only after that agrees.
- Risk: low-medium. Bonus: today the eviction pattern depends on read order, so `kmer_top_count` can differ
  between `-j 1` and `-j 4` on diverse libraries; dense counting makes the block thread-count-invariant.

## FC-2 — One pass per read: fold heatmap, the rolling 12-mer, polyx and frag-GC into the base loop
`class: algorithmic · effort: M · output: identical`

```rust
// accum.rs:324-325, then accum.rs:331, then accum.rs:352
let (hx, hy) = heatmap::encode(seq);
self.heatmap[heatmap::index(hx, hy)] += 1;
```

`heatmap::encode` re-decodes every base (`heatmap.rs:14-28`) to get x = #(T+C), y = #(T+G) — both derivable
from the `code` the base loop already has (x = [code==1]+[code==3], y = [code==2]+[code==3]).
`pack_kmer` re-decodes 12 bytes per window (`repeated.rs:95-108`), and `measure_polyx` walks the read twice
more (`polyx.rs:23-58`). One loop can yield `code` once and update everything, with a rolling 24-bit k-mer:

```rust
let mut win: u32 = 0; let mut valid = 0u8;         // 12 bases = 24 bits
for (i, &raw) in seq.iter().enumerate() {
    let code = CODE[raw as usize];                  // 256-byte LUT, as adapter.rs:87 already does
    win = ((win << 2) | code as u32) & 0x00FF_FFFF;
    valid = if code == 4 { 0 } else { (valid + 1).min(12) };
    if valid == 12 { kmers.add(win); }              // same 12-mers, one decode per base
    // composition cell, heatmap x/y, N-run state, polyx state ...
}
```

- Effect per base **[counted]**: 11 byte-decodes + 11 match arms leave the k-mer path (1680 -> 151 decodes
  per read), one whole sequence pass leaves (`heatmap::encode`), two `measure_polyx` walks leave.
  **[estimated]** 15-25% of the sequence-side work, ~300-600 ns of the measured 2784 ns/read — and the k-mer
  *probe* survives, which is why FC-1 ranks above this.
- **Answer to "can polyx and the per-position arrays be fused into one pass, and what does it save per
  base":** they fuse trivially but polyx contributes almost nothing — its loops break at the first
  transition, so it costs 2-3 byte compares *per read*, not per base. The per-base saving is the heatmap
  pass (1 decode + 1 match per base) plus the 11x k-mer amplification. The quality plane cannot join this
  loop: it is a second array with a different stride, and interleaving writes two streams (see Dead ends).
- Measure: `C1M`/`C5M -j 1`, `--only bases,heatmap,kmer` before/after; `tests/differential.py` for the full
  per-position table plus both histograms.
- Risk: medium — `run_q_head/run_q_tail` and `w.kmer_windows` must hold their exact values, and the N-run
  state machine (`accum.rs:226-248`) is the delicate part. `gated_all_modules_is_exactly_observe` pins the
  gated twin against `observe`.

## FC-3 — Batch reads so the position rows are walked in blocks, not one line at a time
`class: algorithmic/memory · effort: L · output: identical`

```rust
// accum.rs:242-243 and accum.rs:261-262 — one row touch per base, 2 KiB apart
if i < per_pos {
    self.base_pos[i][code] += 1;
```

Not a retry of the plane-major layout that lost by 25% (`accum.rs:8-15`, measured): rows stay position-major,
only the iteration order across reads changes. Stage B=64 reads into an SoA buffer (codes, quals), then walk
positions in blocks of 16 and update each row for all 64 staged reads before moving on.

- Mechanism **[counted]**: a `qual_pos` row is 2048 B and a read touches each of 151 rows once — line
  fetched, one u64 written, line dead: 151 L2/L3 accesses per read for one increment each. With 64 reads
  staged, the same line takes 64 increments before eviction, so L2 traffic per read falls from ~151 lines to
  ~2.4; `base_pos` (2 rows/line) gains too. **[estimated]** 8-20% of wall clock, bounded by the position
  writes' share of the 2784 ns/read. Sequence FC behind FC-1/FC-2 and re-decide from the `--only bases`
  and `--only quality` deltas.
- Measure: `C5M -j 1`, `perf stat -e L1-dcache-load-misses,LLC-load-misses,cache-misses` plus wall; the
  counter that must move is misses per read, not only seconds.
- Risk: high-ish (`observe` becomes `observe_batch`; per-read tables — fingerprint, k-mer, adapter_pos — stay
  per read). Guard with the gated-equality test and `tests/differential.py`.

## FC-4 — A single plain FASTQ uses exactly one core, whatever `-j` says
`class: io-threads · effort: L · output: identical`

```rust
// main.rs:240-243
let reports: Vec<FileReport> = if files.len() == 1 || threads <= 1 {
    for path in &files {
        out.push(scan_one(path, threads.max(1), ...
```

`gz_threads` only reaches `reader::open`'s BGZF branch (`reader.rs:126-135`), so for a plain `.fastq` the
parameter is dead: `-j 4` on one file is one thread. For BGZF the inflate is parallel but the accumulator is
a single serial consumer — the module comment says so (`main.rs:238-239`).

- Change: split a plain file by byte ranges snapped to record boundaries (seek, then find 4 newlines), scan
  each range into its own `Accum`, merge. `Accum::merge` is exact and associative for every field including
  the fingerprint table (`accum.rs:362-432`), so this is the cheapest parallelism in the codebase to justify.
  For BGZF, chunk the inflated stream the same way to decouple accumulation from the inflate pool.
- Effect **[estimated]**: plain cells 1.0x -> ~3.5x at `-j 4`, bounded by the merge (2 MB of counter adds per
  split at 151 positions — negligible); gz cells gain only where the accumulator is the bottleneck, which
  `-j 1` gz-vs-plain at equal decompressed size already measures.
- Measure: `C5M` plain and gz, `taskset -c 0-3`, `-j 1|2|4`, medians of 3; report the speedup curve, not one
  point; diff full JSON against `-j 1`.
- Risk: medium. Two ordering hazards: the k-mer block is order-dependent today (FC-1 fixes it), and the
  budgets (`--reads`, `--bytes`) are per file, so chunking must not make a budgeted run stop at a different
  record or report a different `partial` reason.

## FC-5 — The adapter matcher recomputes constants per position, and its seeded path is dead code
`class: algorithmic · effort: M · output: identical (semantics pinned by tests)`

```rust
// adapter.rs:166-168 — inside a per-read, per-position loop
for off in 0..=a.saturating_sub(6) {
    if seed_code(&adapter.seq[off..off + 6]) != code {
```

`seed_code` runs for every adapter offset at **every read position** — 28 offsets x 6 decodes = 168 ops per
position per adapter [counted] — for values constant for the whole run. That is `match_at`, which nothing in
`src/` calls (grep: only `adapter_pos.rs:40` uses `all_matches`). What `--adapter-content` really runs is the
plain scan (`adapter.rs:233-243`): ~136 positions x ~8 compares per adapter for a clean 151 bp read.
Precompute the seed index once per run and share it: `all_matches` walks the read's 6-mers once and tests
only the offsets the index reports; `match_at` reuses the same structure.

- Effect **[estimated]**: ~4 adapters x 136 x k byte-compares per read -> ~150 table lookups; 10-20x on that
  module, which is opt-in (`cli.rs:118-119`).
- Measure: `C1M --adapter-content -j 1`; correctness gate `adapter::seeded_matcher_agrees_with_the_exhaustive_
  reference` (3000 random cases) must stay green for both matchers.
- Risk: medium — the seed-soundness argument is load-bearing (`adapter.rs:14-17`); keep the exhaustive
  fallback for the short-overlap/high-error region unchanged.

## FC-6 — `--profile` readthrough search: ~3600 branchy byte compares per read, replaceable with memchr
`class: algorithmic (SIMD) · effort: S · output: identical`

```rust
// accum.rs:280-281, --profile path
for (i, (_, pat)) in ADAPTER_PREFIXES.iter().enumerate() {
    if contains_ci(seq, pat.as_bytes()) {
```

`contains_ci` (`accum.rs:127-134`) is `windows(13).any(...)` with `to_ascii_uppercase` per byte, twice:
2 x 138 windows x 13 B = **~3600** operations per read [counted], plus a second full pass for the fragment GC
(`accum.rs:288-294`). `--profile` is the fastsim fitter path — the coupling the tool exists for.

- Change: two-byte SIMD prefilter (`memchr2`-style; a `u64` word loop with a 0x0101.. mask suffices, no new
  dependency), then verify the 13-byte pattern case-insensitively only at candidates. Fold the frag-GC loop
  into the base pass: `start` (first non-G position) is exactly what that loop already walks past, so
  accumulate `frag_gc/frag_bases` for `i >= start` there.
- Effect **[counted]**: 3600 -> ~150-300 ops per read, one sequence pass removed; **[estimated]** 5-10% of a
  `--profile` run, zero without the flag.
- Measure: `C1M --profile -j 1` wall, plus byte-identity of the produced profile JSON (the fit must not move).
- Risk: low. Case-folding must stay exact (`accum.rs:132`) or `adapter_hits` -> `readthrough` shifts.

## FC-7 — SWAR the count reductions in the base and quality loops
`class: code-level (SIMD) · effort: S · output: identical`

```rust
// accum.rs:257-261 — per-read sum plus byte census over 151 bytes
for (i, &q) in qual.iter().enumerate() {
    self.qual_byte[q as usize] += 1;
    qsum += q as u64;
```

The scatter halves (`base_pos[i][code]`, `qual_pos[i][q]`) cannot vectorize — data-dependent scatter, and a
gather/scatter emulation loses to the scalar RMW. The reductions can: `acgt`, `gc`, `qsum` are sums over 151
bytes, expressible as 8-byte words with a byte-equality mask and `count_ones`; the byte census stays scalar
(histogram) but `qsum` falls out of the same word loop.

- Effect **[estimated]**: 2-3x on the reduction half of the base loop and on `qsum`, i.e. 3-8% of the base
  loop and less of the run. Small; a day's work with no semantic risk.
- Measure: `--only bases` and `--only quality` at `-j 1`, `C5M`; `perf stat -e instructions,cycles` (the
  ratio must move, not just the clock).
- Risk: low (integer sums, no FP reordering).

## FC-8 — `qual_pos` rows are 2048 B for a 94-value alphabet
`class: memory/cache layout · effort: S · output: identical for valid FASTQ`

```rust
// accum.rs:26, accum.rs:77
pub const QBINS: usize = 256;
pub qual_pos: Vec<[u64; QBINS]>,
```

Legal Phred bytes are 33..=126 (`accum.rs:490-494`), so 128 bins cover every valid byte: 151 rows x 128 x 8 B
= **151 KiB** instead of 302 KiB [counted]. Bytes >= 128 (already flagged `malformed`, `accum.rs:497`) go to
a small out-of-range census instead of into the rows.

- Effect **[estimated]**: halves the hottest array; the touched-line count per read is unchanged (151), but
  the whole footprint now fits a 151 KiB region that survives across reads, so steady-state L2 hits improve.
  Low single digits alone; compounds with FC-3.
- Measure: `C5M -j 1`, `perf stat -e L1-dcache-load-misses,LLC-load-misses`, RSS on `C1M`.
- Risk: low. Identity holds unless a file carries bytes >= 128, where per-position stats would exclude them
  — say so in the note the report already prints for `malformed`.

## FC-9 — `--tiles` allocates a `Vec` per read and re-checks the row length per base
`class: allocation · effort: S · output: identical`

```rust
// tile.rs:24, tile.rs:46-50
let fields: Vec<&[u8]> = header[..end].split(|&b| b == b':').collect();
for (i, &q) in qual.iter().enumerate() {
    if rows.len() <= i { rows.resize(i + 1, [0; QBINS]); }
```

One heap allocation per read (20 M on `C20M`) [counted], plus a length test and possible `resize` inside the
per-base loop, for a header field that is `fields[4]`.

- Change: walk the 4th/5th colon with an index loop (no `Vec`), `rows.resize(qual.len(), ..)` once before the
  loop, branch-free index inside.
- Effect **[counted]** 1 allocation/read -> 0, 1 branch/base -> 0; **[estimated]** 2-4% of a `--tiles` run,
  none otherwise. Opt-in flag, so ranked last among the scan-path items.
- Measure: `C1M --tiles -j 1`, wall + `perf stat -e instructions`; contract pinned by `tile.rs:143-153`.
- Risk: low.

## FC-10 — `--overrep-mode table` costs ~1.7x its published projection
`class: memory/allocation · effort: M · output: identical`

```rust
// overrep.rs:53 and overrep.rs:311-312 — the structure and the projection it is admitted against
pub struct Packed { len: u16, words: Vec<u64> }
let per = (read_len.div_ceil(4) + 8) as u64;
```

A 151 bp read packs to 3 `u64` = 24 B payload, but each `Packed` also carries a 24-byte `Vec` header, its own
malloc (16-32 B allocator overhead) and a 13 B hashbrown bucket: ~80-90 B resident per unique read against
47 B charged [counted]. The mode's selling point is refusing rather than being OOM-killed
(`main.rs:168-187`), so a 1.7x understatement is a defect in the refusal, not just in the memory.

- Change: arena of packed payloads (one bump `Vec<u8>`, `u32` offsets) plus open-addressed
  `FxHashMap<u64, (u32 off, u32 len, u32 count)>` keyed by xxh3 of the payload: 24 B payload + 13 B bucket,
  which is what the projection already promises.
- Effect **[counted]**: ~1.7x lower resident cost, an honest projection, and SipHash off that path.
- Measure: `C1M --overrep-mode table`, VmHWM against the projection printed on stderr (`main.rs:163-167`) —
  today print and RSS disagree, and `table_projection_matches_the_plan_quote` only pins the formula.
- Risk: medium (a second hash set to keep consistent with `Packed`'s `Hash`); do not change the refusal
  threshold and the layout in the same commit.

## FC-11 — Insert pass: pack the anchor; stop sampling only the first 20 000 pairs
`class: algorithmic/correctness-enabling · effort: S · output: search identical, sampling changes (behind a flag)`

```rust
// insert.rs:167-169
let mut p = 0usize;
while p + ANCHOR <= bio1 {
    if anchor_matches(&s1[p..p + ANCHOR], anchor) {
```

`anchor_matches` uppercases per byte (`insert.rs:185-190`), so a censored pair (the majority for long inserts)
pays ~240 positions x 12 B = **2900** compares [counted]. `revcomp` allocates a `Vec` per pair
(`insert.rs:131-133`, called at 155), and `observe_pair` recomputes `trailing_run` for both mates
(`insert.rs:259-260`) that `overlap` already computed.

- Change: precompute the anchor as a 24-bit code, compare against a rolling code of `s1` (one shift-or and one
  compare per position, ~240 ops); reuse two scratch buffers for the reverse complements; return `t1/t2` from
  `overlap`. Sampling: `scan_pair` reads pairs `0..max_pairs` in file order (`insert.rs:570-586`) — the first
  tiles of the flowcell. A strided or reservoir sample makes the fitted distribution representative; that
  changes `insert_hist`, so it must be `--insert-sample first|strided`, default `first`.
- Effect **[estimated]**: 3-10x on the per-pair search, which at the default 20 000 pairs (`cli.rs:90`) is tens
  of milliseconds of a run — not where the wall clock is. The sampling change is about the fit, not speed.
- Measure: `--insert-hist --insert-pairs 200000` on `C1M -j 1`: wall, and the `measured` /
  `censored_no_overlap` counts must be identical for the search changes.
- Risk: low for the search; the sampling change is report-visible.

## FC-12 — `--overrep-mode second-pass` re-reads the file with one thread, whatever `-j` says
`class: io-threads · effort: S · output: identical`

```rust
// report.rs:369 — thread count is the literal 1
"overrep": crate::overrep::report(overrep, total, &files_paths(files), 1, overrep_min_count)?,
```

`overrep::report` forwards `threads` to `second_pass` -> `reader::open(path, gz_threads)`
(`overrep.rs:280-286`), so the BGZF pool is built with 1 thread: the one mode that reads the file twice
decompresses the second pass serially.

- Change: thread the resolved `cli.thread_count()` through `write_report`. One argument.
- Effect **[estimated]**: second pass on `C5M` gz goes from ~1x to ~2-3x inflate; the absolute saving is one
  pass over the file, in the only mode that pays one.
- Measure: `C1M --overrep-mode second-pass` gz at `-j 1|2|4` against `--overrep-mode fingerprint` on the same
  file (which does not re-read) to isolate the pass.
- Risk: none beyond plumbing; report bytes unchanged (`passes: 2` either way).

## FC-13 — Per-record parsing goes through `&mut dyn BufRead` four times and copies each line
`class: io/code-level · effort: M · output: identical`

```rust
// reader.rs:160-162
pub(crate) fn read_line(rdr: &mut dyn BufRead, out: &mut Vec<u8>) -> io::Result<usize> {
    out.clear();
    let n = rdr.read_until(b'\n', out)?;
```

Four virtual calls and four copies per record (80 M of each on `C20M`) [counted]. The `BufReader` is 1 MiB
(`reader.rs:26`), so there is no syscall problem — the cost is the non-inlinable `read_until` and the copy.

- Change: a `Records<R: BufRead>` iterator over the 1 MiB buffer that finds the 4 newlines with a memchr and
  yields `&[u8]` slices into it (one straddling record per buffer handled by refill). No copy in the common
  case.
- Effect **[estimated]**: 5-10% of a gz cell once FC-1/FC-2 have thinned the accumulator, less before; also
  removes `Vec::clear` + capacity churn from the hot path.
- Measure: `C5M -j 1`, wall + `perf stat -e instructions,branch-misses`; identity via
  `tests/differential.py`, framing via the wrapped/short-record refusals (`reader.rs:396-402`).
- Risk: medium — a buffer-boundary bug is silent, so the differential run must include a file whose records
  straddle 1 MiB.

## Protocol (before coding any of this)

1. Cells: `C1M`, `C5M`, `C20M`, plain and gz, `taskset -c 0`/`-c 0,1`/`-c 0-3` with `-j 1|2|4`, median of 3, one
   tool at a time on an idle machine. Nothing else is quotable.
2. Attribution without guessing: `--only <module>` runs only that module's loops and prints
   `fastcheck: measurement work: ...` (`gated.rs:113-126`); a module's cost is the delta between `--only`
   cells, and a zero counter proves a loop did not run. Use it to settle FC-1 vs FC-2 order before writing
   either.
3. Counters: `perf stat -e cycles,instructions,cache-misses,L1-dcache-load-misses,LLC-load-misses` for the
   cache claims (FC-3, FC-8); `/usr/bin/time -v` VmHWM for memory claims (FC-1, FC-10) — VmHWM must stay flat
   from `C1M` to `C20M`, which is the property the tool is sold on.
4. Identity gate: every prop marked `output: identical` must produce a byte-identical JSON on a fixture and a
   clean `tests/differential.py` run (it compares every metric, both histograms and the whole per-position
   table against `tests/reference_qc.py`).

## What is already good

- Memory is a function of read length, not of input: **[measured]** 133-140 MB from 100 k to 20 M pairs.
- The fingerprint table is projected by an identity over the real entry width and admission-checked, with the
  projection pinned against a real hashbrown allocation (`memory.rs:392-416`). Refusing before allocating is
  rare and right. Counters are integers everywhere and `merge` is exact for every field, so file-level
  parallelism costs nothing in correctness (`accum.rs:362-432`).
- The Phred offset is decided from an exact byte census, and per-position rows are indexed by the raw byte, so
  a wrong +33/+64 cannot silently propagate (`accum.rs:20-23`, `488-545`).
- Gating is enforced in the accumulator and the `Work` counters make "did it skip the work" answerable from
  outside the process (`gated.rs:80-127`) — the benchmarking infrastructure this dossier leans on.
- The main path allocates once per file, not per record: `read_line` reuses four `Vec`s, there are no temp
  files and no second pass (`reader.rs:209-214`).
- The k-mer approximation is labelled with its exact failure mode and bound (`rich.rs:146-152`), its worst
  pathology (a full `retain` per insert at saturation) is already fixed (`repeated.rs:44-64`), and a losing
  layout is recorded with its measurement so it cannot be re-tried (`accum.rs:8-15`).

## Dead ends / explicitly rejected

1. **Preallocate the fingerprint table to its projection.** `main.rs:112-117` computes a capacity from file
   size / 350 B while the table starts at 1024 and grows (`accum.rs:163`), which looks like a miss. It is not:
   at 20 M distinct reads the projection is **544 MiB** [counted, `project(2e7, 16)`] against a **[measured]**
   138.8 MB peak, so real distinctness is a small fraction of it and preallocating would quadruple RSS to save
   ~14 rehashes over 20 M inserts.
2. **Plane-major `qual_pos[byte][pos]`.** Already measured 25% slower (5.14 s vs 4.10 s, 2 M x 150 bp) and
   recorded at `accum.rs:8-15`. Do not re-derive it; position-major won.
3. **Emit only the non-zero heatmap cells.** 65,536 numbers per JSON (`report.rs:371`, `rich.rs:175`) for at
   most a few thousand non-zero cells, while every other block skips zeros (`report.rs:79-87`). It breaks
   fasthelper, which requires a dense array and indexes it positionally (`if counts.len() < (size * size)`
   then `counts[(y * size + x)]`, `fasthelper/src/report.rs:209-217`), and D-1 forbids shipping the schema
   change alongside. The block is dense by contract.
4. **Replace the per-position histograms with moments or a quantile sketch.** `pos_stats` needs exact
   q10/q50/q90 per position (`report.rs:96-116`) and the accumulator's claim is that nothing is estimated; a
   t-digest or log-binned histogram would halve the array and break that claim. The same reasoning rejects
   merging the base and quality loops into one body: two arrays, two strides, and the merge buys loop
   bookkeeping while costing each chain its locality.
