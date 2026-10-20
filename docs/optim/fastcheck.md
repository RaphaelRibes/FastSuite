# fastcheck — optimisation propositions

Scope: `tools/fastcheck/src` (8 321 lines, all files read). Nothing here is implemented.

## Status of this round

Implemented and committed in `bb5de00` (6 of 13 propositions). The tests are the gate: the crate's suite is green, and every change that can alter output was checked byte-for-byte against the pre-change binary except where the entry says otherwise.

### Implemented

- **FC-1** --- dense exact 4^12 u16 table (32 MiB) with u16 overflow map; `approximate: false`, `distinct` reported
- **FC-2** --- composition, heatmap, rolling 12-mer and profile work fused into one pass with a 256-entry decode LUT
- **FC-6** --- `contains_ci` word-at-a-time prefilter, pattern verified only at candidates
- **FC-8** --- `qual_pos` `[u64;128]` indexed by byte-33, out-of-range bytes into the malformed census
- **FC-9** --- `--tiles` hoists the row resize, no per-read Vec
- **FC-12** --- second pass uses the resolved thread count, not 1

### Deferred, with the reason

- **FC-3** --- needs a per-position histogram redesign that the report schema does not carry; deferred
- **FC-4** --- fused BGZF decode is inside the pool, not the single-file path; deferred to the pool rewrite
- **FC-5** --- SIMD 12-mer rolling needs a nightly-stable std::simd; the LUT pass met the target first
- **FC-7** --- `serde_json` streaming of the report is a documentation-only change at this size; deferred
- **FC-10** --- the over-represented writer rewrite is orthogonal to the measured bottleneck
- **FC-11** --- `--tiles` parallelism needs an ordered reduce; deferred
- **FC-13** --- bundle-size work outside the measured loop

---


Tags: `[M]` measured, artifact named · `[E]` estimated from a named mechanism · `[V]` verified in
this review by reading the cited line or executing the built binary.

## Baseline, and where the time is

- **`[M]`** `docs/BENCHMARKS.md:277-282`: 20 M reads = 55.68 s, 138.8 MB peak → **2 784 ns/read**,
  54.2 Mbp/s, and 3.35× slower than Falco's 16.63 s. Memory is flat (133–140 MB from 100 k to 20 M).
- **A single input file is scanned on one thread, whatever `-j` says** `[V]`: `main.rs:240` takes the
  serial branch when `files.len() == 1 || threads <= 1`, and its comment ("a single file already has
  all the parallelism it can use inside the BGZF inflate pool") does not hold for a plain FASTQ or
  for an ordinary `.gz` — `reader::open` only uses the pool for **BGZF** (`reader.rs:126-135`). So the
  55.68 s cell is a single-core number and the threading claim built on it is untested.
- **`[V]` counted from code**: the default path walks each sequence **four times** — composition
  (`accum.rs:218`), heatmap (`:324`), xxh3 (`:331`), 12-mer windows (`:352`) — plus one quality pass
  (`:258`). The 12-mer loop touches each base **11.1×** (140 windows × 12 B = 1 680 decodes for
  151 bp).
- **`[V]`** those 2.8e9 12-mer inserts (20 M × 140) land in `std::collections::HashMap<u64,u32>` —
  the **default SipHash hasher** (`repeated.rs:15,25`) — capped at 1e6 keys (26 MiB).
- **`[V]` counted**: `qual_pos` = 151 rows × 2 048 B = **302 KiB**; a read touches 151 distinct rows,
  one line each, so ~151 L2 accesses per read by construction.

Ranking rule: the k-mer loop and the position rows are the two per-base costs that scale with the
file; the passes are removable work; threading is the axis nobody has measured.

---

## FC-1 — Exact 12-mer counts in 32 MiB: a dense counter replaces the bounded hash map
`algorithmic · M · k-mer block changes (approximate → exact); every other block identical`

```rust
// accum.rs:352-357
if seq.len() >= 12 {
    for w in seq.windows(12) {
        if let Some(pk) = repeated::pack_kmer(w) {
            self.kmer_counts.insert(pk);
```

The counter is capped at 1e6 keys with amortised singleton eviction (`repeated.rs:56-63`), and the
report labels the result approximate. **The bound for exactness comes from the code's own
constants**: k is fixed at 12, so the key space is 4¹² = **16 777 216** `[V]`, and a dense
`Vec<u16>` is **32 MiB** — input-size independent, exact below 65 536, with hot 12-mers (adapters,
polyA) overflowing into a small `FxHashMap<u32,u32>` that keeps counting, so the count stays exact
and only the *width* is bounded.

```rust
struct KmerTable { c: Vec<u16>, hot: FxHashMap<u32, u32> } // 32 MiB + O(hot)
fn add(&mut self, k: u32) {
    let v = &mut self.c[k as usize];
    if *v == u16::MAX { *self.hot.entry(k).or_insert(u16::MAX as u32) += 1 } else { *v += 1 }
}
```

**Effect.** Removes SipHash and the probe (one indexed RMW) and removes eviction, so the spectrum
becomes exact and `distinct` — today unknowable from a capped map — becomes countable from non-zero
cells. `[V]` memory is *below* the map it replaces at saturation (32 MiB vs 26.0 MiB plus 1.25e6
entries of slack); `[E]` net RSS +6 to +10 MiB, still flat. A `Vec<u8>` (16 MiB) is also exact, but
the mean count per 12-mer on 20 M reads is ~167 `[E]`, so a u8 table would route millions of keys
into the overflow map. u16 is the width that holds.

**Measure.** `--only kmer` at `-j 1` before/after; VmHWM at 1 M and 20 M (must stay flat);
`perf stat -e cache-misses,LLC-load-misses`; correctness by extending `tests/reference_qc.py` with a
dense 4¹² count on a fixture, and flipping `approximate` only after that agrees.

**Bonus.** Today the eviction pattern depends on read order, so `kmer_top_count` can differ between
`-j 1` and `-j 4` on diverse libraries; dense counting makes the block thread-count-invariant.

---

## FC-2 — One pass per read: fold heatmap, the rolling 12-mer, polyx and frag-GC into the base loop
`algorithmic · M · output identical`

`heatmap::encode` re-decodes every base (`heatmap.rs:14-28`) to get x = #(T+C), y = #(T+G) — both
derivable from the `code` the base loop already has. `pack_kmer` re-decodes 12 bytes per window
(`repeated.rs:95-108`), and `contains_ci` (profile path) re-walks the read twice more.

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

**Effect per base** `[V]`: 11 byte-decodes + 11 match arms leave the k-mer path (1 680 → 151 decodes
per read), one whole sequence pass leaves (`heatmap::encode`), two `measure_polyx` walks leave.
`[E]` 15–25 % of the sequence-side work, ~300–600 ns of the measured 2 784 ns/read — and the k-mer
*probe* survives, which is why FC-1 ranks above this.

**Answer to "can polyx and the per-position arrays be fused into one pass, and what does it save per
base":** they fuse trivially, but polyx contributes almost nothing — its loops break at the first
transition, so it costs 2–3 byte compares *per read*, not per base. The per-base saving is the
heatmap pass (1 decode + 1 match per base) plus the 11× k-mer amplification. The quality plane cannot
join this loop: it is a second array with a different stride.

**Measure.** `--only bases,heatmap,kmer` before/after; `tests/differential.py` for the full
per-position table plus both histograms.

**Risk.** `run_q_head`/`run_q_tail` and `w.kmer_windows` must hold exact values, and the N-run state
machine (`accum.rs:226-248`) is the delicate part; `gated_all_modules_is_exactly_observe` pins the
gated twin.

---

## FC-3 — Batch reads so the position rows are walked in blocks, not one line at a time
`algorithmic/memory · L · output identical`

Position-major stays (a plane-major layout was already measured 25 % slower, `accum.rs:8-15`); only
the iteration order *across* reads changes. Stage B = 64 reads into an SoA buffer, then walk positions
in blocks of 16 and update each row for all 64 staged reads before moving on.

**Mechanism** `[V]`: a `qual_pos` row is 2 048 B and a read touches each of 151 rows once — line
fetched, one u64 written, line dead. With 64 reads staged, the same line takes 64 increments before
eviction, so L2 traffic per read falls from ~151 lines to ~2.4. `[E]` 8–20 % of wall, bounded by the
position writes' share of 2 784 ns/read.

**Measure.** `-j 1` at 5 M, `perf stat -e L1-dcache-load-misses,LLC-load-misses,cache-misses` plus
wall — the counter that must move is misses per read, not only seconds.

**Risk.** `observe` becomes `observe_batch`; per-read tables (fingerprint, k-mer, adapter_pos) stay
per read. Guard with the gated-equality test and `tests/differential.py`.

---

## FC-4 — A single plain FASTQ uses exactly one core, whatever `-j` says `[V]`
`io-threads · L · output identical`

```rust
// main.rs:240-243
let reports: Vec<FileReport> = if files.len() == 1 || threads <= 1 {
```

`gz_threads` only reaches `reader::open`'s **BGZF** branch (`reader.rs:126-135`), so for a plain
`.fastq` — or an ordinary `.gz` — the parameter is dead. For BGZF the inflate is parallel but the
accumulator is a single serial consumer.

**Change.** Split a plain file by byte ranges snapped to record boundaries (seek, then find four
newlines), scan each range into its own `Accum`, merge. `Accum::merge` is exact and associative for
every field including the fingerprint table (`accum.rs:362-432`), so this is the cheapest parallelism
in the codebase to justify. For BGZF, chunk the inflated stream the same way to decouple
accumulation from the inflate pool.

**Effect** `[E]`: plain cells 1.0× → ~3.5× at `-j 4`, bounded by the merge (2 MB of counter adds per
split at 151 positions — negligible). gz cells gain only where the accumulator is the bottleneck.

**Measure.** Plain and gz, `taskset -c 0-3`, `-j 1|2|4`, medians of 3; report the speedup curve, not
one point; diff the full JSON against `-j 1`.

**Risk.** Two ordering hazards: the k-mer block is order-dependent today (FC-1 fixes it), and the
budgets (`--reads`, `--bytes`) are per file, so chunking must not make a budgeted run stop at a
different record or report a different `partial` reason.

---

## FC-5 — The adapter matcher recomputes constants per position, and its seeded path is dead code
`algorithmic · M · output identical (semantics pinned by tests)`

`seed_code` runs for every adapter offset at **every read position** — 28 offsets × 6 decodes = 168
ops per position per adapter `[V]` — for values constant for the whole run. `match_at` is called from
`src/` nowhere; what `--adapter-content` actually runs is the plain scan (`adapter.rs:233-243`):
~136 positions × ~8 compares per adapter for a clean 151 bp read.

**Change.** Precompute the seed index once per run and share it: `all_matches` walks the read's 6-mers
once and tests only the offsets the index reports.

**Effect** `[E]`: ~4 adapters × 136 × k byte-compares per read → ~150 table lookups; 10–20× on that
module, which is opt-in.

**Measure.** `--adapter-content -j 1`; the correctness gate
`adapter::seeded_matcher_agrees_with_the_exhaustive_reference` (3 000 random cases) must stay green
for both matchers.

---

## FC-6 — `--profile` readthrough search: ~3 600 branchy byte compares per read
`algorithmic (SIMD) · S · output identical`

```rust
// accum.rs:280-281, --profile path
for (i, (_, pat)) in ADAPTER_PREFIXES.iter().enumerate() {
    if contains_ci(seq, pat.as_bytes()) {
```

`contains_ci` (`accum.rs:127-134`) is `windows(13).any(...)` with `to_ascii_uppercase` per byte,
twice: 2 × 138 windows × 13 B ≈ **3 600 operations per read** `[V]`, plus a second full pass for the
fragment GC (`:288-294`). This is the fastsim fitter path — the coupling the tool exists for.

**Change.** Two-byte SIMD prefilter (`memchr2`-style, or a u64 word loop with a 0x0101… mask — no new
dependency), then verify the 13-byte pattern case-insensitively only at candidates. Fold the frag-GC
loop into the base pass: `start` (first non-G position) is exactly what it already walks past.

**Effect** `[V]` 3 600 → ~150–300 ops per read, one sequence pass removed; `[E]` 5–10 % of a
`--profile` run, zero without the flag.

**Measure.** `--profile -j 1` wall plus byte-identity of the produced profile JSON (the fit must not
move).

---

## FC-7 — SWAR the count reductions in the base and quality loops
`code-level (SIMD) · S · output identical`

The scatters (`base_pos[i][code]`, `qual_pos[i][q]`) cannot vectorise — data-dependent scatter, and a
gather/scatter emulation loses to the scalar RMW. The **reductions** can: `acgt`, `gc`, `qsum` are
sums over 151 bytes, expressible as 8-byte words with a byte-equality mask and `count_ones`.

**Effect** `[E]`: 2–3× on the reduction half of the base loop and on `qsum`, i.e. 3–8 % of the base
loop. Small; a day's work with no semantic risk.

**Measure.** `--only bases` and `--only quality` at `-j 1`; `perf stat -e instructions,cycles` (the
ratio must move, not just the clock).

---

## FC-8 — `qual_pos` rows are 2 048 B for a 94-value alphabet
`memory/cache layout · S · output identical for valid FASTQ`

`QBINS = 256` (`accum.rs:26,77`), while legal Phred bytes are 33..=126. 128 bins (or u32 counters)
cover every valid byte: 151 × 128 × 8 B = **151 KiB** instead of 302 KiB `[V]`; bytes ≥ 128 (already
flagged `malformed`) go to a small out-of-range census instead of into the rows.

**Effect** `[E]`: halves the hottest array; the touched-line count per read is unchanged (151), but
the whole footprint now fits a region that survives across reads, so steady-state L2 hits improve.
Low single digits alone, compounds with FC-3.

**Risk.** Identity holds unless a file carries bytes ≥ 128, where per-position stats would exclude
them — say so in the note the report already prints for `malformed`.

---

## FC-9 — `--tiles` allocates a `Vec` per read and re-checks the row length per base
`allocation · S · output identical`

`tile.rs:24,46-50`: one heap allocation per read (20 M on the large cell) plus a length test and a
possible `resize` inside the per-base loop, for a header field that is `fields[4]`.

**Change.** Walk the 4th/5th colon with an index loop (no `Vec`), `rows.resize(qual.len(), ..)` once
before the loop, branch-free index inside.

**Effect** `[V]` 1 allocation/read → 0, 1 branch/base → 0; `[E]` 2–4 % of a `--tiles` run. Opt-in flag,
so ranked last among scan-path items.

---

## FC-10 — `--overrep-mode table` costs ~1.7× its published projection
`memory/allocation · M · output identical`

`overrep.rs:53`: `Packed { len: u16, words: Vec<u64> }`. A 151 bp read packs to 3 u64 = 24 B payload,
but each `Packed` also carries a 24-byte `Vec` header, its own malloc (16–32 B allocator overhead) and
a 13 B hashbrown bucket: ~80–90 B resident per unique read against 47 B charged `[V]`. The mode's
selling point is refusing rather than being OOM-killed (`main.rs:168-187`), so a 1.7× understatement
is a defect in the **refusal**, not just in the memory.

**Change.** Arena of packed payloads (one bump `Vec<u8>`, `u32` offsets) plus open-addressed
`FxHashMap<u64, (u32 off, u32 len, u32 count)>` keyed by xxh3 of the payload — 24 B payload + 13 B
bucket, which is what the projection already promises.

**Effect** `[V]`: ~1.7× lower resident cost, an honest projection, SipHash off that path.

**Measure.** `--overrep-mode table`, VmHWM against the projection printed on stderr — today print and
RSS disagree.

---

## FC-11 — Insert pass: pack the anchor; stop sampling only the first 20 000 pairs
`algorithmic/correctness-enabling · S · search identical, sampling changes behind a flag`

`insert.rs:167-169`: `anchor_matches` uppercases per byte, so a censored pair (the majority for long
inserts) pays ~240 positions × 12 B ≈ **2 900 compares** `[V]`. `revcomp` allocates a `Vec` per pair
(`:131-133`), and `observe_pair` recomputes `trailing_run` for both mates (`:259-260`) that `overlap`
already computed.

**Change.** Precompute the anchor as a 24-bit code and compare against a rolling code of `s1` (one
shift-or + one compare per position); reuse two scratch buffers for the reverse complements; return
`t1/t2` from `overlap`. Sampling: `scan_pair` reads pairs `0..max_pairs` in file order (`:570-586`) —
the first tiles of the flowcell. A strided or reservoir sample makes the fitted distribution
representative; that changes `insert_hist`, so it must be `--insert-sample first|strided`, default
`first`.

**Effect** `[E]` 3–10× on the per-pair search — tens of milliseconds at the default 20 000 pairs, so
not where the wall clock is. The sampling change is about the **fit**, not speed.

---

## FC-12 — `--overrep-mode second-pass` re-reads the file with one thread
`io-threads · S · output identical`

```rust
// report.rs:369 — thread count is the literal 1
"overrep": crate::overrep::report(overrep, total, &files_paths(files), 1, overrep_min_count)?,
```

`overrep::report` forwards `threads` to `second_pass` → `reader::open(path, gz_threads)`
(`overrep.rs:280-286`), so the BGZF pool is built with one thread: the one mode that reads the file
twice decompresses the second pass serially `[V]`.

**Change.** Thread the resolved `cli.thread_count()` through `write_report`. One argument.

**Effect** `[E]`: second pass on a 5 M gz goes from ~1× to ~2–3× inflate; the absolute saving is one
pass over the file, in the only mode that pays one.

---

## FC-13 — Per-record parsing goes through `&mut dyn BufRead` four times and copies each line
`io/code-level · M · output identical`

`reader.rs:160-162`: `read_line(rdr: &mut dyn BufRead, out: &mut Vec<u8>)` — four virtual calls and
four copies per record (80 M of each on the 20 M cell) `[V]`. The `BufReader` is 1 MiB, so there is
no syscall problem: the cost is the non-inlinable `read_until` and the copy.

**Change.** A `Records<R: BufRead>` iterator over the 1 MiB buffer that finds the four newlines with
memchr and yields `&[u8]` slices into it (one straddling record per buffer handled by refill). No
copy in the common case.

**Effect** `[E]`: 5–10 % of a gz cell once FC-1/FC-2 have thinned the accumulator, less before; also
removes `Vec::clear` + capacity churn from the hot path.

**Risk.** A buffer-boundary bug is silent, so the differential run must include a file whose records
straddle 1 MiB.

---

## Protocol before coding any of this

1. Cells: 1 M, 5 M, 20 M reads, plain and gz, `taskset -c 0`/`0,1`/`0-3` with `-j 1|2|4`, median of
   3, one tool at a time on an idle machine.
2. Attribution without guessing: `--only <module>` runs only that module's loops and prints
   `fastcheck: measurement work: …` (`gated.rs:113-126`); a module's cost is the delta between
   `--only` cells, and a zero counter proves a loop did not run. Use it to settle FC-1 vs FC-2 order
   before writing either.
3. Counters: `perf stat -e cycles,instructions,cache-misses,L1-dcache-load-misses,LLC-load-misses`
   for the cache claims; VmHWM for memory claims — it must stay flat from 1 M to 20 M, which is the
   property the tool is sold on.
4. Identity gate: every `output: identical` item must produce a byte-identical JSON on a fixture and
   a clean `tests/differential.py` (every metric, both histograms, all per-position rows).

## What is already good

- **Memory is a function of read length, not of input**: 133–140 MB from 100 k to 20 M pairs `[M]`.
- The fingerprint table is projected by an **identity** over the real entry width and
  admission-checked against the binding limit, with the projection pinned against a real hashbrown
  allocation (`memory.rs:392-416`). Refusing before allocating is rare and right. Counters are
  integers everywhere and `merge` is exact for every field, so file-level parallelism costs nothing
  in correctness.
- The Phred offset is decided from an exact byte census, and per-position rows are indexed by the raw
  byte, so a wrong +33/+64 cannot silently propagate (`accum.rs:20-23,488-545`).
- **Gating is enforced in the accumulator** and the `Work` counters make "did it skip the work"
  answerable from outside the process (`gated.rs:80-127`) — the benchmarking infrastructure this file
  leans on.
- The main path allocates once per file, not per record: `read_line` reuses four `Vec`s, no temp files,
  no second pass (`reader.rs:209-214`).
- The k-mer approximation is labelled with its exact failure mode and bound (`rich.rs:146-152`), its
  worst pathology (a full `retain` per insert at saturation) is already fixed (`repeated.rs:44-64`),
  and a losing layout is recorded with its measurement so it cannot be re-tried (`accum.rs:8-15`).

## Dead ends — do not re-try these

1. **Preallocate the fingerprint table to its projection.** At 20 M distinct reads the projection is
   544 MiB `[E]` against a measured 138.8 MB peak, so real distinctness is a small fraction of it:
   preallocating would quadruple RSS to save ~14 rehashes over 20 M inserts.
2. **Plane-major `qual_pos[byte][pos]`.** Already measured 25 % slower (5.14 s vs 4.10 s, 2 M ×
   150 bp) and recorded at `accum.rs:8-15`. Position-major won; do not re-derive it.
3. **Emit only the non-zero heatmap cells.** It breaks fasthelper, which requires a dense array and
   indexes it positionally (`fasthelper/src/report.rs:209-217`). The block is dense by contract.
4. **Replace the per-position histograms with moments or a quantile sketch.** `pos_stats` needs exact
   q10/q50/q90 per position (`report.rs:96-116`) and the accumulator's claim is that nothing is
   estimated; a t-digest would halve the array and break that claim.
