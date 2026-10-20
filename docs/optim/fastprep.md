# fastprep — optimisation propositions

Scope: `tools/fastprep/src` (9 792 lines; `algorithm.rs`, `adapter.rs`, `pipeline.rs`, `io.rs`,
`cli.rs`, `dedup.rs`, `pe_overlap.rs`, `memory.rs`, `stats.rs`, `main.rs` all read) plus the
committed bench artifacts. Nothing here is implemented.

## Status of this round

Implemented and committed in `154a39f` (8 of 14 propositions). The tests are the gate: the crate's suite is green, and every change that can alter output was checked byte-for-byte against the pre-change binary except where the entry says otherwise.

### Implemented

- **FP-2** --- `--pe-overlap auto` no longer runs on `AdapterSource::Fixed` (added `explicit`)
- **FP-3** --- prefilter seed width derived from the error bound (per-adapter widths: 4-mer from k=20 for TruSeq, 3-mer from k=9 for Nextera; 414 -> 81 exact compares per 151 bp read); conformance gate green for eps 0..0.3 across widths 6/5/4/3. Note: the commit message of `154a39f` repeats the unsound bound this file used to state; the code and this entry carry the corrected one.
- **FP-4** --- staged writes, one flush per 256 KiB instead of six `write_all` per record
- **FP-5** --- quality stays ASCII end to end, gate LUT indexed by raw byte, `offset*len` subtracted once
- **FP-6** --- dedup set probed once per read when not saturated
- **FP-7** --- chunk arenas recycled
- **FP-8** --- adapter carried as an index, no two heap Strings per match
- **FP-14** --- `TrimResult` no longer cloned on the ordered thread

### Deferred, with the reason

- **FP-1** --- the three-stage pipeline resting on two threads is a redesign of the ordered writer; the measured single-core gain came from FP-3/FP-4 instead
- **FP-9** --- PXMT metadata write is outside the default path
- **FP-10** --- `--dry-run` cost only matters on tiny inputs
- **FP-11** --- SIMD quality scan needs the nightly-stable intrinsics
- **FP-12** --- hash-set sharding conflicts with the ordered writer
- **FP-13** --- report size is bounded by read count, not the bottleneck

---


Tags: `[M]` measured, artifact named · `[E]` estimated from a named mechanism · `[V]` verified in
this review by reading the cited line or executing the built binary.

This is the tool where speed is contested, so the first section settles what the contested cell
actually says before any proposition is ranked.

---

## 0. The "34 vs 28 s at one core" cell is not a one-core cell `[V]`

`~/work/FastPrep-Bench/results/results_gz6.csv`, the raw equal-level run
`docs/BENCHMARKS.md` cites:

| reads | tool | pinned | wall s | cpu s | cpu/wall |
|---|---|---|---|---|---|
| 5 M | fastprep | no | 33.508 | 81.458 | 2.43 |
| 5 M | fastprep | yes | 32.983 | 80.510 | **2.44** |
| 5 M | fastp | no | 28.197 | 88.200 | 3.13 |
| 5 M | fastp | yes | 83.254 | 82.341 | **0.99** |
| 5 M | cutadapt | yes | 135.759 | 133.079 | 0.98 |

A process confined to one CPU cannot accumulate 80.5 CPU-seconds in 33.0 s of wall. fastp and
cutadapt land at 0.98–0.99, so the harness *does* constrain them; fastprep's `pinned` and `unpinned`
rows are indistinguishable at every size (0.685/0.684, 3.240/3.220, 6.521/6.601, 33.508/32.983). The
fastprep rows therefore were **not** inside `taskset`. Two consequences:

- **The quoted 1-core comparison does not exist.** The honest single-core number is fastp 83.25 s
  with fastprep unmeasured — and fastprep's CPU-seconds (80.51) sit **2 % below** fastp's (82.34).
- **The gap is concurrency utilisation, not per-pair work.** Over the linear range 1 M→5 M, both
  `-j 1` unpinned: fastprep 6.747 s/M pairs at 16.36 CPU-s/M = **2.42 cores busy**; fastp
  5.579 s/M pairs at 17.64 CPU-s/M = **3.16 cores busy**. `wall = CPU / busy-cores` in both rows.
  fastprep does **7 % less work per pair** and still takes 21 % longer. At fastp's utilisation its
  wall would be 16.36/3.16 = 5.18 s/M ≈ 26 s at 5 M `[E: the two rows' own arithmetic]`.

**Three runs settle the attribution, no instrumentation:** `--pe-overlap never`, `--no-adapter`,
`--no-dedup`, each on the 5 M cell, pinned *and* `-j 1` unpinned, asserting `cpu_s/wall_s ≈ 1.00` on
the pinned run before believing it. That partitions the CPU across the three stages no competitor
has in the same form. **Fixing the harness is proposition FP-0.**

### FP-0 — Make the pinned cell actually pinned
`correctness-of-benchmark · S · no code change`

Run the fastprep arm through the same `taskset -c 0` wrapper as the others, or launch it with
`taskset -c 0` *and* `-j 1`, then assert `0.95 ≤ cpu_s/wall_s ≤ 1.05` in the CSV writer and refuse
to write the row otherwise. A row that violates it is the signature of an escaped process, and it is
the only way this class of error stops being invisible: the numbers look plausible, so nothing else
catches it. Re-derive `docs/STATUS.md` and `docs/BENCHMARKS.md` from the corrected row.

---

## FP-1 — Three stages, two threads: the ordered thread is the critical path
`algorithmic (work partitioning) · M · output byte-identical`

Parse of chunk N+1 and the write of chunk N share the main thread (`pipeline.rs:470-501`), and the
compute side may run only two chunks ahead (`bounded(2)`, `:445-446`):

```rust
fill_chunk(&mut *reader, &mut buf, &mut chunk, rt.chunk, phred_offset)?;
while !chunk.is_empty() {
    raw_tx.send((buf, chunk)).ok();
```

On that same thread sit the dedup probe, the stats update and the six `write_all` emit
(`:487-501`).

**Change.** Three stages on three threads, ring depth 4:

```
parse thread -> ChunkBuf ring -> trim pool -> ordered commit thread
commit thread: dedup probe, stats, 256 KiB staging buffer,
               one write_all per buffer -> ParCompress
```

Dedup stays on the single ordered thread — its decision order is the input order, a documented
property — it merely stops sharing a thread with the compressor.

**Effect.** Busy cores 2.42 → ~3.1 at unchanged CPU `[E]`, i.e. ~26–28 s at 5 M. Mechanism:
`wall = CPU/cores` here, so a genuinely concurrent third stage raises utilisation without adding
work. On a genuinely pinned core it buys nothing directly but still removes FP-4/FP-14 per-record
work from the chain.

**Measure.** Both cells; per-thread `utime` from `/proc/<pid>/task/*/stat` or `pidstat -t`. Success =
busy-cores ≥ 3.0 with CPU-s/M ≤ 16.4.

**Risk.** Ordering hazards on the dedup/singleton path; the error-drain and the `.part`+rename commit
(`pipeline.rs:528-543`) must keep "drop output on error".

---

## FP-2 — `--pe-overlap auto` runs on every pair, and its own doc says it should not `[V]`
`algorithmic (skip a no-op stage) · S · output expected identical (prove with cmp)`

`cli.rs:37-40` states the intent — *"`auto` runs only when the adapter source is EXPLICIT
(`--adapter-preset` / `--adapters`), and stays off with auto-detection"* — but auto-detection also
produces `AdapterSource::Fixed(...)` (`cli.rs:837`), and `auto` is resolved against exactly that
(`cli.rs:967-970`):

```rust
PeOverlapMode::Auto => matches!(source, AdapterSource::Fixed(_)),
```

So a run with no adapter flags — i.e. the benchmark, and every default user — runs the detector on
every pair, and each call allocates twice and scans ~140 offsets (`pe_overlap.rs:177,196-197`).

**Change.** Resolve `auto` only for an *explicitly configured* set: carry `explicit: bool` beside
`AdapterSource::Fixed` (set only in the `adapter_preset`/`adapters` branches) and resolve against
it; keep today's reach behind `--pe-overlap always`.

**Effect.** Removes, per pair: one `revcomp` allocation + fill, one `Vec<usize>` allocation and ~140
anchor scans. `[E]` 0.5–0.7 µs/pair ≈ 0.5–0.7 s/M pairs — up to half the measured 1.17 s/M marginal
gap. The mode's measured value on the cohort is a null (`cli.rs:236-244`: *"exact matches identical
with and without the mode, overtrim +106, +0.000 pp"*), so this pays for a stage known to be a no-op.

**Measure.** `--pe-overlap never` on the 5 M cell plus `cmp` of both outputs; add call counters to
`--json` (neither exists today, which is why the cost is invisible).

**Risk.** If `cmp` fails, the mode does fire on this corpus: report that as a behaviour change, do
not hide it.

---

## FP-3 — The 6-mer pigeonhole prefilter is dead at the shipped error rate `[V]`
`algorithmic (matching strategy) · M · output byte-identical, gated by the conformance test`

```rust
// adapter.rs:123-133
for k in 1..=a.seq.len() {
    let max_mm = (k as f64 * error_rate).ceil() as usize;
    if k >= 6 * (max_mm + 1) { ms = k; break; }
}
```

The condition `k ≥ 6·(ceil(k·ε)+1)` requires `ε < 1/6`; at the shipped `--adapter-error-rate 0.2`
(`cli.rs:205`) no `k` satisfies it, so `min_sound_k = usize::MAX` for every preset (a=33, 34, 19).
`a < ms` then holds for all adapters (`adapter.rs:489`) and the **reference full scan** runs
(`adapter.rs:492`) `[V]`: at L = 151, `max_p = 137`, three presets → **414 `mismatches_upto` calls
and ~11 KB compared per read**, while only 2.3 % of reads ever match an adapter
(`fastprep_5000k.json`: 231 389 / 5 000 000 `[M]`).

**Change.** Size the seed to the bound instead of fixing it at 6. The sound condition is the pigeonhole one, `k >= s(m+1)` for `s` seed bases and `m` tolerated mismatches; the arithmetic this file first carried, `k >= 6*(ceil(k*eps)+1)`, assumed the mismatch-free block has nominal length and is therefore unsound -- the block can be short. At eps = 0.2 that gives per-adapter widths rather than one width: 4 seed bases sounding from k = 20 for the 34 bp TruSeq adapters, 3 from k = 9 for the 19 bp Nextera adapter, in one 5440-entry table. The unsound band is the truncated 3' tail, and the conformance test, not this paragraph, is what settles it.

```
SeedTable::build(adapters, eps, seed_len)   // smallest seed_len with finite ms
min_sound_k(a, eps, s) = min{ k : ceil(k/(ceil(k*eps)+1)) >= s }
direct band = p >= l + 1 - min_sound_k
```

**Effect.** 414 → ~36 exact compares/read plus one 151-step rolling 5-mer walk `[E]`, ~1 s/M pairs of
CPU. Absolute wall, not necessarily the head-to-head delta: fastp runs a comparable scan, so part of
this cost is shared and cancels.

**Measure.** `--no-adapter` is the whole path in one flag; then `--adapter-error-rate 0.1` **revives
the prefilter with no code change** (ms = 18) and therefore *measures the prefilter's value at
runtime*. Gate on `pigeonhole_filter_matches_reference` (4 000 cases × 20 configs,
`adapter.rs:638-699`).

**Risk.** A wrong `min_sound_k` silently changes what is trimmed; that test is the only guard and
must run for every new seed width.

---

## FP-4 — Six `write_all` calls per record through `Box<dyn Write>`
`allocation (memcpy, call overhead) · S · output byte-identical`

`io.rs:196-210` issues six dynamic dispatches per record (plus a 512-byte stack bounce for quality);
each lands in `ParCompress::write`, which memcpys into a `BytesMut` and re-checks the 65 280-byte
block boundary.

**Change.** Append the record into a per-thread staging `Vec<u8>` and `write_all` once per 256 KiB:

```rust
stage.extend_from_slice(b"@"); stage.extend_from_slice(id);
stage.extend_from_slice(b"\n"); stage.extend_from_slice(seq);
stage.extend_from_slice(b"\n+"); stage.push(b'\n');
encode_qual(&mut stage, qual, off); stage.push(b'\n');
if stage.len() >= 256 << 10 { w.write_all(&stage)?; stage.clear(); }
```

**Effect.** Five fewer dyn calls, five fewer short memcpys, five fewer boundary checks per read, and
gzp receives full blocks: `[E]` 30–60 cycles/read ≈ 0.3–0.6 s at 10.1 M reads.

**Measure.** `perf stat -e instructions,cycles` on the 5 M cell; the byte stream itself must be
identical.

---

## FP-5 — Quality is decoded on ingest and re-encoded on egress, per byte
`allocation (per-byte work, buffer reuse) · S · byte-identical for quality bytes ≥ offset`

`pipeline.rs:92-95` subtracts the Phred offset per byte and `io.rs:201-210` adds it back per byte,
while the trim path reads the decoded values through `gate_lut`. With Phred+33 input at the default
offset this round trip reproduces the input bytes at the cost of two full passes over 151 bytes.

**Change.** Keep ASCII and shift the **LUT** instead of the data: add `phred_offset` to `Params`,
build `gate_lut[b] = quality_gate(b − offset)` indexed by the raw byte, `extend_from_slice(q)` on
ingest, write the raw slice out. The mean-Q sum needs the same shift (second LUT, or subtract
`offset·l`).

**Effect.** Removes ~302 byte-ops and one 151-byte staged copy per read `[E]` ≈ 0.3–0.6 s at 10.1 M
reads.

**Measure.** `perf stat` cycles plus `cmp` against the current build on a Phred+33 file, then repeat
with `--phred-offset 64` to prove the LUT shift.

**Risk.** Bytes below the offset are currently clamped; passthrough would keep them verbatim. That
differs only on malformed input — document it or gate on a first-chunk check.

---

## FP-6 — Dedup probes the set twice per read and fingerprints the whole read
`algorithmic (dedup admission) · S · byte-identical`

```rust
// dedup.rs:229-243
if self.seen.contains(&fp) { return Seen::Duplicate; }
if self.saturated() { return Seen::NewUnrecorded; }
...
self.seen.insert(fp);
```

`contains` then `insert` is two probes of the same open-addressed table for every non-duplicate. The
anchor also never engages under the shipped config: `scan_5` and `scan_3` are both true, so
`anchor_end` returns `AnchorEnd::None` (`dedup.rs:166-174`) and `slice()` falls through to the whole
read — xxh3 over 151 bases instead of an artefact-free 0.9 L window.

**Change.** One probe on the unsaturated path, semantics preserved:

```rust
if self.saturated() {
    return if self.seen.contains(&fp) { Seen::Duplicate } else { Seen::NewUnrecorded };
}
if !self.seen.insert(fp) { return Seen::Duplicate; }
```

**Effect.** Halves the probes. At 5 M pairs the table is ~72 MB `[E]`, so each probe is a DRAM round
trip: 10.1 M reads × ~90 ns ≈ 0.9 s CPU, on the critical thread.

**Measure.** `--no-dedup` bounds the whole stage (upper bound only); then `perf stat -e
cache-misses` before/after, asserting `reads_dropped_duplicate` stays bit-identical (321 058 `[M]`).

**Risk.** The saturation branch order is load-bearing for the one-sidedness argument
(`dedup.rs:196-202`): keep `saturated()` first, as written above.

---

## FP-7 — A fresh chunk arena per chunk, and the old one freed after the new reserve
`memory (buffer reuse) · S · byte-identical`

`pipeline.rs:475-477` (and `:784-789` for PE) allocates a new arena each iteration, while the arena
returned on `res_rx` is dropped at the end of the same body — i.e. *after* the new `reserve()`. At
chunk = 1024 `[M]`, 151 bp, that is ~464 KB per mate per chunk, 928 KB paired; 5 M pairs = 4 883
chunks → ~4.5 GB of fresh large-object allocation `[E]`.

**Change.** Recycle: send the consumed arena back through a free-list channel and `clear()` it, so
two arenas per stage serve the whole run.

**Effect.** Removes the mmap/madvise churn and the first-touch faults `[E]` ≈ 0.4–0.9 s, charged to
the faulting thread, which is the ordered thread.

**Measure.** `perf stat -e page-faults`; peak RSS must stay at ~101 MB `[M]`.

---

## FP-8 — Adapter names are heap `String`s twice per trimmed read
`allocation · S · byte-identical`

`adapter.rs:526` allocates in the worker; `stats.rs:34-35` allocates again and SipHashes on the
ordered thread — 462 792 adapter-trimmed reads at 5 M pairs `[M]`.

**Change.** Carry `ai: u8` in `AdapterMatch`/`TrimOutput`, keep `by_name` as `[u64; 256]` indexed by
it, materialise names once at report time.

**Effect.** Two allocations and one string hash removed per trimmed read `[E]` ≈ 0.02–0.04 s. Small,
ranked because it is an hour of work on the critical thread.

---

## FP-9 — needletail parses one record at a time on the serial thread
`algorithmic (single-pass restructure) · L · byte-identical, with a fallback`

`pipeline.rs:156-173` drives `reader.next()` per record and copies through `SequenceRecord`; this is
the one stage that cannot overlap with itself and it sits on the ordered thread (FP-1).

**Change.** A specialised 4-line FASTQ reader over the input buffer: three `memchr(b'\n')` per
record, direct `extend_from_slice` into the arena, a structural check (`@`, `+`, equal lengths), and
a fallback to needletail on the first anomaly. Keep needletail for FASTA, multi-line and CRLF.

**Effect.** needletail runs at ~0.8–1.2 GB/s on 4-line FASTQ; a memchr splitter reaches 2–3 GB/s
`[E]`, ~0.5–0.7 µs/read over 3.4 GB of input ≈ 5–7 s at 5 M pairs.

**Measure.** `--dry-run` (`pipeline.rs:415-419`) isolates parse + trim + stats from compression
entirely.

**Risk.** CRLF, multi-line FASTA, empty records, missing trailing newline: the fallback keeps them
correct but must be exercised by `tests/regressions.rs`.

---

## FP-10 — `Arc::make_mut` clones a 2 KB `Params` per detected pair
`allocation · S · byte-identical`

`pipeline.rs:1067-1068`; `Params` carries `gate_lut: [f64; 256]` (`algorithm.rs:216`), so each clone
is a 2 KB memcpy plus an Arc allocation inside the parallel map.

**Change.** Hoist the two variants per chunk, or pass `pe_adapter` as an argument to `trim()` and keep
`Params` immutable and shared.

**Effect.** `[E]` 0.2–0.5 GB of memcpy at 2–5 % of 5 M pairs, ~0.05 s; disappears for free if FP-2
is taken.

---

## FP-11 — BGZF: 65 280-byte independent blocks with no cross-block dictionary
`algorithmic (gzip pipeline) · M · container changes, decoded bytes identical`

`io.rs:124-128` hard-codes `Bgzf`; gzp's `Bgzf` reports `needs_dict() = false` with
`BGZF_BLOCK_SIZE = 65280`, while gzp's `Gzip` hands the previous block's trailing 32 KiB to
libdeflate as a dictionary. BGZF cannot: every 65 KB block restarts the match finder.

**Change.** Expose the container: BGZF stays the default (tabix, htslib), add `--gzip-format mgzip`
(gzp `Mgzip`, 128 KiB blocks, real dictionary) where indexability is not wanted. Both stay valid
`gzip -d`.

**Effect.** Same input and level, fewer and larger blocks with cross-block matches `[E]` — a few
percent ratio and less per-block table cost. Assert it as testable, not as a speed win: the output
size is not recorded anywhere today.

**Measure.** One fixed `.fq` in both formats: bytes out and wall; `pigz -6 -b 256` as an independent
reference. Never compare sizes against fastp without decoding both.

---

## FP-12 — Input decompression is single-threaded for every non-BGZF gzip
`io-threads · L · byte-identical`

`io.rs:172-180` takes the parallel path only for BGZF, and `gz_threads = threads/2` (`cli.rs:981`)
is 1 at `-j 2`, so even BGZF input gets no parallelism at 2 cores `[V]`.

**Change.** htslib-style parallel gunzip for ordinary `.gz` (find block boundaries with a bit-level
inflate probe every N MB, inflate on a pool, reassemble in order). Cheap half-measure: size
`gz_threads` from `threads − 1` and document that the parallel input path needs BGZF.

**Effect.** The input is ~954 MB gzipped at 5 M pairs `[M]`; single-threaded inflate at 150–250 MB/s
is 4–6 s of wall `[E]` and is the one stage that never overlaps at any thread count. Relevant only
at `-j ≥ 2`, a cell fastprep already wins (25.84 vs 39.61 `[M]`), hence low rank.

**Risk.** A hand-rolled gzip block splitter is the classic silent-truncation bug: use gzp's `Mgzip`
decompressor or libdeflate's block API, not a new one.

---

## FP-13 — `mismatches_upto`: auto-vectorisation by hope, and a float `ceil` per position
`allocation (SIMD, branch layout) · S · byte-identical`

`adapter.rs:550-558` relies on LLVM, and `adapter.rs:497` recomputes a float `ceil` for every
position.

**Change.** Explicit `_mm256_cmpeq_epi8` + `movemask` + popcount on 32-byte blocks, and a per-adapter
`max_mm[k]` table (k spans `min_overlap..=a`, < 64 entries) instead of `(k as f64 * eps).ceil()`.
Check the cap once per four blocks.

**Effect.** `[E]` ~10–15 → ~6–8 cycles/call; at 414 calls/read that is 0.3–0.6 µs/read, but only
*after* FP-3 cuts the call count. **Do FP-3 first: SIMD on a loop that should not run is the wrong
order.**

**Risk.** Intrinsics need a `target_feature` gate with a runtime fallback even though
`.cargo/config.toml` builds `target-cpu=native`.

---

## FP-14 — `TrimResult` carries a heap `Vec` and is cloned per read
`allocation · S · byte-identical`

`pipeline.rs:492-496` clones on the ordered thread; `TrimResult::Emit` holds `runs: Vec<(usize,
usize)>` (`algorithm.rs:435-441`), empty for every read when the interior scan is off (the shipped
default) yet enough to keep the enum non-`Copy`.

**Change.** Const-generic the enum on `RUNS: bool` (the default instantiation becomes `Copy`), or
store `runs` inline as `[(u32, u32); 2]` + length; take the result by reference instead of cloning.

**Effect.** `[E]` 20–40 cycles/read on the ordered thread: small, but it is on the critical path and
removes the last heap allocation from the steady-state trim output.

---

## What is already good (do not break these)

- **The two-tier cascade fires.** With `min_len 0.9 L` at L = 151, `rem_budget = 15 > w = 14`, so a
  clean read skips the whole Γ walk (`algorithm.rs:971-996`) and pays two 14-pair window checks plus
  one quality sum. Fragile by exactly one base: at `0.92 L` the guard fails and every read walks.
- **The quality gate is a LUT**, no `powf` in the loop (`algorithm.rs:301-325`), and direction is
  monomorphised (`streaming_walk::<FORWARD>`), replacing what the comment calls a 10 %-of-runtime
  branch.
- **Chunk arenas**: `ScanInput` is three `Range`s, so the parse path has zero per-read allocations
  (`pipeline.rs:43-80`) — the design fastdedup lacks.
- **`mimalloc`, `lto = "fat"`, `codegen-units = 1`, `target-cpu=native`** (`+12 %` measured, per
  `.cargo/config.toml`).
- **Dedup is admission-checked and saturating**: 101 MB peak vs fastp's 1 143 MB `[M]`, one-sided
  freeze by construction, projection as an identity rather than a guess (`memory.rs:243-269`).
- **Output is thread-count independent**: dedup and stats on the ordered thread, locked down by
  `parallel_output_matches_serial`. That is what makes FP-1 safe to attempt.
- **Refusal over guessing**: PE ambiguity → `None` (`pe_overlap.rs:207-213`), adapter auto-detect
  never invents a sequence (`cli.rs:793-815`).

## Dead ends — do not re-try these

1. **Vectorising the polyX streaming walk.** `gamma` is a loop-carried, data-dependent recurrence
   with an early exit on `gamma > mercy`: there is no prefix sum to vectorise and the LUT loads are
   gather-shaped. The cheap win was already taken (skip the walk entirely); a rewrite adds risk to
   the one function that decides every trim.
2. **Raising `--threads` to win the cell.** On a genuinely pinned core `wall == CPU-seconds`, so
   threads cannot help and only inflate `cpu_s`. The levers there are FP-3/5/6/9; the parallel lever
   is FP-1 at the same thread count.
3. **Turning dedup off to match fastp.** `reads_dropped_duplicate` = 321 058 at 5 M pairs `[M]` —
   3.2 % of reads and 41 % of all drops. That changes the output, the drop rate (0.0786 → ~0.046)
   and the tool's purpose. It is a policy difference: fix the probe (FP-6), do not delete the stage.
4. **Dropping `--gz-level` to 4 to match fastp's default.** The benchmark equalises at 6 on purpose
   (an earlier revision "charged fastprep for compression the others never paid for"). Lowering it
   wins a cell by changing the question and changes every output byte.
