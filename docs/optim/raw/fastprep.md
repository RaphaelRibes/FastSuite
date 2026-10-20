# fastprep 0.3.0 — ranked optimization dossier

Read-only static review of `tools/fastprep/src` (9,792 lines) plus the committed benchmark artifacts. No build, no run, no source edit. Every number is tagged `[M]` measured (artifact named) or `[E]` estimated (mechanism named). Cycle counts assume ~3.5 GHz. Files read in full: `algorithm.rs`, `adapter.rs`, `pipeline.rs`, `io.rs`, `cli.rs`, `dedup.rs`, `pe_overlap.rs`, `memory.rs`, `stats.rs`, `main.rs`; also `docs/BENCHMARKS.md`, `~/work/FastPrep-Bench/scripts/run_bench.py`, `~/work/FastPrep-Bench/results/{results_gz6.csv,fastprep_5000k.json}`, and gzp-2.0.2 in the cargo registry.

## 0. The 34-vs-28 s cell: the artifacts do not say what the table says

`results_gz6.csv` (the raw equal-level run `docs/BENCHMARKS.md` cites) `[M]`:

| reads | tool | pinned | wall s | cpu s | cpu/wall |
|---|---|---|---|---|---|
| 5,000,000 | fastprep | no | 33.508 | 81.458 | 2.43 |
| 5,000,000 | fastprep | yes | 32.983 | 80.510 | 2.44 |
| 5,000,000 | fastp | no | 28.197 | 88.200 | 3.13 |
| 5,000,000 | fastp | yes | 83.254 | 82.341 | 0.99 |

**(a) The "pinned 1 core" fastprep rows are not pinned.** A process confined to one CPU cannot accumulate 80.5 CPU-seconds in 33.0 s of wall. fastp's pinned rows sit at 0.99-1.00; fastprep's pinned and unpinned rows are indistinguishable at every size (0.685/0.684, 3.240/3.220, 6.521/6.601, 33.508/32.983). So the fastprep rows came from outside `taskset` — `run_bench.py:106-150` wraps both tools identically, which makes this stale rows, not a code path. The 34.37-vs-28.50 pair is therefore `-j 1` on 16 cores, and the genuinely pinned 1-core cell is fastp 83.25 s with fastprep unmeasured and its CPU-seconds within 2% of fastp's (80.5 vs 82.3). A true 1-core comparison is a near-tie.

**(b) The gap is concurrency utilisation, not per-pair work.** Over the linear range 1M->5M pairs, both `-j 1` unpinned `[M]`: fastprep 6.747 s/M pairs at 16.36 CPU-s/M pairs = 2.42 cores busy, fastp 5.579 s/M pairs at 17.64 CPU-s/M pairs = 3.16 cores busy. `wall = CPU / busy-cores` in both rows (16.36/2.42 = 6.76; 17.64/3.16 = 5.58). fastprep does **7% less work per pair** and still takes 21% longer. At the same CPU total and fastp's utilisation its wall would be 16.36/3.16 = **5.18 s/M pairs, ~26 s at 5M** `[E: the two rows' own arithmetic]`. Nothing in the per-base arithmetic explains the gap; 0.74 core of parallelism does, and FP-1 shows where it is left idle.

Corollary for ranking: on the cell as the brief states it the lever is CPU-seconds per read (FP-3, FP-5, FP-6, FP-9, FP-12); on the cell as shipped it is FP-1. Both are covered.

**Three runs settle it, no instrumentation:** `--pe-overlap never`, `--no-adapter`, `--no-dedup`, each on the 5M cell, pinned *and* `-j 1` unpinned, asserting `cpu_s/wall_s = 1.00` on the pinned run before believing it. Those three flags partition the trim-path CPU across the three stages no competitor has in the same form.

## 1. Ranked propositions

### FP-1 — Three stages, two threads: the ordered thread is the critical path
**class:** algorithmic (work partitioning) · **effort:** M · **byte-identical:** yes

**Evidence** — `pipeline.rs:470-486`: parse of chunk N+1 and the write of chunk N share the main thread, in lockstep:
```rust
fill_chunk(&mut *reader, &mut buf, &mut chunk, rt.chunk, phred_offset)?;
while !chunk.is_empty() {
    raw_tx.send((buf, chunk)).ok();
    let mut next_buf = ChunkBuf::default();
    ...
    let (res_buf, results) = res_rx.recv()...
```
and `pipeline.rs:487-501` puts the dedup probe, the stats update and the six-`write_all` emit on that same thread, while the compute thread may run only 2 chunks ahead (`bounded(2)`, lines 445-446):
```rust
let dup = dedup.as_mut().is_some_and(|d| d.is_duplicate(seq, None));
stats.record(seq.len(), &result, p.output.adapter_name.as_deref());
if !dup { if let Some(w) = writer.as_mut() { p.emit_to(&res_buf, w, phred_offset)?; } }
```

**Change** — three stages on three threads, ring depth 4:
```
parse thread  -> ChunkBuf ring -> trim pool -> ordered commit thread
commit thread: dedup probe, stats, batch into a 256 KB staging buffer,
               one write_all per buffer -> ParCompress
```
Dedup stays on the single ordered thread (its decision order is the input order, a documented property); it merely stops sharing a thread with the compressor.

**Effect** — busy cores 2.42 -> ~3.1 at unchanged CPU `[E]`, i.e. ~5.2-5.5 s/M pairs, ~26-28 s at 5M. Mechanism: wall = CPU/cores here, so a real third concurrent stage raises utilisation without adding work. On a genuinely pinned core it buys nothing directly (wall = CPU) but still removes FP-4's per-record work from the chain.

**Measure** — both cells; `pidstat -t -p <pid> 1` or sample `/proc/<pid>/task/*/stat` for per-thread `utime`; report busy-cores. Success = busy-cores >= 3.0 with CPU-s/M pairs <= 16.4.

**Risk** — ordering hazards on the dedup/singleton path; the error-drain and `.part`+rename commit (`pipeline.rs:528-543`) must keep "drop output on error".

### FP-2 — `--pe-overlap auto` runs on every pair; its own docs say it should not
**class:** algorithmic (skip a no-op stage), also correctness-enabling · **effort:** S · **byte-identical:** expected yes, must be proven with `cmp`

**Evidence** — `cli.rs:37-40` states the intent:
```rust
/// `auto` ... runs only when the adapter source is EXPLICIT
/// (`--adapter-preset` / `--adapters`), and stays off with auto-detection
```
but detection also yields `AdapterSource::Fixed` (`cli.rs:837`), `auto` is resolved against `Fixed` (`cli.rs:967-970`), and `runs()` is `!matches!(self, Never)` (`cli.rs:72-74`). So the benchmark, which passes no adapter flags and auto-detects, runs the detector on every pair (`pipeline.rs:1060-1061`), and each call allocates twice and scans 140 offsets (`pe_overlap.rs:177,196-197`):
```rust
let rc2 = revcomp(r2);
let mut admissible: Vec<usize> = Vec::new();
for o in (anchor..=max_insert).rev() {
```

**Change** — resolve `auto` only for an explicitly configured set: carry `explicit: bool` alongside `AdapterSource::Fixed` (set only in the `adapter_preset` / `adapters` branches) and resolve `auto` against it; keep today's reach behind `--pe-overlap always`.

**Effect** — removes, per pair: one 151-byte `revcomp` allocation + fill, one `Vec<usize>` allocation, 140 `anchor_matches` calls (~3 compares each on random sequence). Estimated 0.5-0.7 us/pair `[E: 140 x ~12 cycles + 2 malloc/free + revcomp ~450 cycles]` = 2.5-3.5 s CPU at 5M = 0.5-0.7 s/M pairs, i.e. up to half the 1.17 s/M pairs marginal gap `[M]`. The mode's measured value on the cohort is a null (`cli.rs:236-244`: "exact matches identical with and without the mode, overtrim +106, +0.000 pp"), so this is paying for a stage known to be a no-op.

**Measure** — `--pe-overlap never` on the 5M cell plus `cmp` of both outputs; add counters for `detect_pe_overlap_trim` calls and ambiguity refusals to `--json` (neither exists today, which is why the cost is invisible).

**Risk** — if `cmp` fails the mode does fire on this corpus and the difference must be reported as a behaviour change, not hidden.

### FP-3 — The 6-mer pigeonhole prefilter is dead at the shipped error rate
**class:** algorithmic (matching strategy) · **effort:** M · **byte-identical:** yes, gated by the existing conformance test

**Evidence** — `adapter.rs:122-134`, soundness needs a 6-base exact run, which no preset has once `eps = 0.2` (the shipped default, `cli.rs:205`):
```rust
// Smallest k with k >= 6 * (max_mm + 1), max_mm = ceil(k·e).
for k in 1..=a.seq.len() {
    let max_mm = (k as f64 * error_rate).ceil() as usize;
    if k >= 6 * (max_mm + 1) { ms = k; break; }
}
```
For a in {33,34,19} at eps 0.2 there is no such k, so `ms = usize::MAX` (by inspection: eps 0.1 -> 18,18,18; eps 0.2 -> none). `ms = MAX` makes `a < ms` true for every adapter (`adapter.rs:489`) and takes the reference-style full scan (`adapter.rs:492`):
```rust
if a < ms {
    for p in 0..=max_p {
        let k = a.min(l - p);
        if k < min_overlap { break; }
        let max_mm = (k as f64 * set.error_rate).ceil() as usize;
```
At L = 151, `max_p = 137`, three presets: **414 `mismatches_upto` calls and ~11 KB compared per read, on every read**, while only 2.3% of reads ever match (`fastprep_5000k.json`: `r1_reads_adapter_trimmed` 231,389 / 5,000,000 `[M]`). The seed pass at `adapter.rs:441-477` is never entered.

**Change** — size the seed to the bound instead of fixing it at 6. A 5-mer table needs an exact run of >= 5, which holds from k = 25 at eps 0.2, so the unsound band shrinks to the 3' tail (k < 25) — 11 positions per adapter, not 138:
```
SeedTable::build(adapters, eps, seed_len)   // smallest seed_len with finite ms
min_sound_k(a, eps, s) = min{ k : ceil(k/(ceil(k*eps)+1)) >= s }
direct band            = p >= l + 1 - min_sound_k
```

**Effect** — 414 -> ~36 exact compares/read plus one 151-step rolling 5-mer walk `[E: arithmetic above]`; saves 0.9-1.2 us/read = 9-12 s CPU at 10.1M reads, ~1.0-1.2 s/M pairs. Absolute wall, not necessarily the head-to-head delta: fastp runs a comparable scan, so part of this cost is shared and cancels.

**Measure** — `--no-adapter` is the whole path in one flag; then `--adapter-error-rate 0.1`, which revives the prefilter with no code change (ms = 18) and so *measures the prefilter's value* at runtime. Gate on `pigeonhole_filter_matches_reference` (4,000 cases x 20 configs, `adapter.rs:638-699`).

**Risk** — a wrong `min_sound_k` silently changes what is trimmed; that test is the only guard and must be run for each new seed width.

### FP-4 — Six `write_all` calls per record through `Box<dyn Write>`
**class:** allocation (memcpy, call overhead) · **effort:** S · **byte-identical:** yes

**Evidence** — `io.rs:196-200` (plus a 512-byte stack bounce for quality at 201-210); every call is a dynamic dispatch into `ParCompress::write`, which memcpy's into a `BytesMut` and re-checks the 65,280-byte block boundary (`gzp/src/par/compress.rs:413-415`):
```rust
w.write_all(b"@")?;
w.write_all(id_trimmed)?;
w.write_all(b"\n")?;
w.write_all(seq)?;
w.write_all(b"\n+\n")?;
```

**Change** — append the record into a per-thread staging `Vec<u8>`, one `write_all` per 256 KB:
```rust
stage.extend_from_slice(b"@"); stage.extend_from_slice(id);
stage.extend_from_slice(b"\n"); stage.extend_from_slice(seq);
stage.extend_from_slice(b"\n+\n"); encode_qual(&mut stage, qual, off);
stage.push(b'\n');
if stage.len() >= 256 << 10 { w.write_all(&stage)?; stage.clear(); }
```

**Effect** — ~5 fewer dyn calls, 5 fewer short memcpys and 5 fewer boundary checks per read; one contiguous copy per block `[E: 30-60 cycles/read = 0.3-0.6 s at 10.1M reads]`. Mechanism: fewer call boundaries, and gzp receives full blocks.

**Measure** — `perf stat -e instructions,cycles` on the 5M cell; assert the output size is identical to the byte (the byte stream itself is unchanged).

**Risk** — staging memory must be flushed on the error path and before commit.

### FP-5 — Quality is decoded on ingest and re-encoded on egress, per byte
**class:** allocation (per-byte work, buffer reuse) · **effort:** S · **byte-identical:** for quality bytes >= offset (all well-formed FASTQ; see risk)

**Evidence** — `pipeline.rs:92-95` and `io.rs:201-210`:
```rust
buf.qual.extend(q.iter().map(|&b| b.saturating_sub(phred_offset)));
...
for (dst, &q) in buf[..n].iter_mut().zip(rest.iter()) {
    *dst = q.saturating_add(phred_offset);
}
```
The trim path then reads decoded values through `gate_lut` (`algorithm.rs:313-316`). With Phred+33 input at the default offset 33 this round trip reproduces the input bytes at the cost of two full passes over 151 bytes plus a 512-byte bounce, per read.

**Change** — keep ASCII and shift the LUT instead of the data: add `phred_offset` to `Params`, build `gate_lut[b] = quality_gate(b - offset)` indexed by the raw byte, `extend_from_slice(q)` on ingest, write `raw_qual` straight out in `write_fastq`. The mean-q sum needs the same shift (second LUT, or subtract `offset * l`).

**Effect** — removes ~302 byte-ops and one 151-byte staged copy per read `[E: 100-200 cycles/read = 0.3-0.6 s at 10.1M reads]`. Mechanism: two whole passes over the quality string collapse to one memcpy.

**Measure** — `perf stat` cycles plus `cmp` against the current build on a Phred+33 file, then repeat with `--phred-offset 64` to prove the LUT shift.

**Risk** — bytes below the offset are currently clamped (`sub` then `add` yields exactly `offset`); passthrough would keep them verbatim. Differs only on malformed input: document it, or gate on a minimum check over the first chunk.

### FP-6 — Dedup probes the set twice per read and fingerprints the whole read
**class:** algorithmic (dedup admission) · **effort:** S · **byte-identical:** yes

**Evidence** — `dedup.rs:229-243`: `contains` then `insert` is two probes of the same open-addressed table for every non-duplicate:
```rust
if self.seen.contains(&fp) {
    return Seen::Duplicate;
}
if self.saturated() { return Seen::NewUnrecorded; }
...
self.seen.insert(fp);
```
and the anchor never engages under the shipped config: `scan_5` and `scan_3` are both true, so `anchor_end` returns `AnchorEnd::None` (`dedup.rs:166-174`) and `slice()` falls through to the whole read (`dedup.rs:156`) — xxh3 over 151 bases instead of an artefact-free 0.9L window.

**Change** — one probe on the unsaturated path, semantics preserved:
```rust
if self.saturated() {
    return if self.seen.contains(&fp) { Seen::Duplicate } else { Seen::NewUnrecorded };
}
if !self.seen.insert(fp) { return Seen::Duplicate; }
```
Optional (the code argues against it at `pipeline.rs:238-242`, and that argument holds): fingerprint in the compute stage and carry `fp` in `Processed`.

**Effect** — halves the probes. At 5M pairs the table is ~72 MB (8M buckets x 9 B `[E]`, consistent with the measured 101.4 MB peak `[M]`), so each probe is a DRAM round trip: 10.1M reads x ~90 ns = ~0.9 s CPU `[E]`, on the critical thread.

**Measure** — `--no-dedup` bounds the whole stage (it changes drop counts, so upper bound only); then `perf stat -e cache-misses` before/after and assert `reads_dropped_duplicate` is bit-identical (321,058 `[M]`).

**Risk** — the saturation branch order is load-bearing for the one-sidedness argument (`dedup.rs:196-202`); keep `saturated()` first, as written above.

### FP-7 — A fresh chunk arena per chunk, and the old one freed after the new reserve
**class:** memory (buffer reuse) · **effort:** S · **byte-identical:** yes

**Evidence** — `pipeline.rs:475-477` (and 784-789 for PE) allocates a new arena each iteration, while the arena returned on `res_rx` is dropped at the end of the same body, i.e. after the new `reserve`:
```rust
let mut next_buf = ChunkBuf::default();
next_buf.reserve(rt.chunk, 151);
let mut next = Vec::with_capacity(rt.chunk);
```
At chunk = 1024 `[M, cli.rs:390]`, 151 bp, that is ~464 KB per mate per chunk, 928 KB paired; 5M pairs = 4,883 chunks -> ~4.5 GB of fresh large-object allocation `[E]`.

**Change** — recycle: send the consumed arena back through a free-list channel and `clear()` it, so two arenas per stage serve the whole run.

**Effect** — removes ~4.5 GB of mmap/madvise churn and the first-touch faults `[E: ~1.1M faults x 0.4-0.8 us = 0.4-0.9 s]`. Mechanism: faults are charged to the faulting thread, which is the ordered thread.

**Measure** — `perf stat -e page-faults` before/after; peak RSS must stay at ~101 MB `[M]`.

**Risk** — a returned arena must be provably unborrowed; `ScanInput` ranges make that checkable, but the drain path needs care.

### FP-8 — Adapter names are heap Strings twice per trimmed read
**class:** allocation · **effort:** S · **byte-identical:** yes

**Evidence** — `adapter.rs:526` allocates in the worker; `stats.rs:34-35` allocates again and SipHashes on the ordered thread:
```rust
name: set.adapters[ai].name.clone(),
...
*self.by_name.entry(name.to_string()).or_insert(0) += 1;
```
462,792 adapter-trimmed reads at 5M pairs `[M, fastprep_5000k.json]`.

**Change** — carry `ai: u8` in `AdapterMatch`/`TrimOutput`, keep `by_name` as a `[u64; 256]` indexed by it, materialise names once at report time.

**Effect** — two allocations and one string hash removed per trimmed read `[E: 40-80 ns each = 0.02-0.04 s]`. Small; ranked here because it is an hour of work on the critical thread.

**Measure** — `perf stat -e instructions`; assert `--json`'s adapter block is identical.

**Risk** — the JSON must keep name keys: render the array into the same map.

### FP-9 — needletail parses one record at a time on the serial thread
**class:** algorithmic (single-pass restructure) · **effort:** L · **byte-identical:** yes, with a FASTQ-only fast path and the existing fallback

**Evidence** — `pipeline.rs:156-173` drives `reader.next()` per record and copies through `SequenceRecord`; `io.rs:92-95` then maps the quality slice again. This is the one stage that cannot overlap with itself, and it sits on the ordered thread (FP-1).

**Change** — a specialized 4-line FASTQ reader over the input buffer: three `memchr(b'\n')` per record, direct `extend_from_slice` into the arena, and a structural check (`@`, `+`, equal lengths) that falls back to needletail on the first anomaly; keep needletail for FASTA, multi-line and `\r\n` input.

**Effect** — needletail runs at roughly 0.8-1.2 GB/s on 4-line FASTQ; a memchr splitter reaches 2-3 GB/s `[E]`, i.e. ~0.5-0.7 us/read saved over 3.4 GB of input = 5-7 s at 5M pairs `[E]`. Mechanism: three scans instead of a state machine plus `Cow` and validation layers.

**Measure** — `--dry-run` (`pipeline.rs:415-419`) isolates parse + trim + stats from compression entirely; compare it to the full run to size everything that is not the writer.

**Risk** — CRLF, multi-line FASTA, empty records, missing trailing newline: the fallback keeps them correct but must be exercised by `tests/regressions.rs`.

### FP-10 — `Arc::make_mut` clones a 2 KB `Params` per detected pair
**class:** allocation · **effort:** S · **byte-identical:** yes

**Evidence** — `pipeline.rs:1067-1068`; `Params` carries `gate_lut: [f64; 256]` (`algorithm.rs:216`), so each clone is a 2 KB memcpy plus an Arc allocation, inside the parallel map:
```rust
Arc::make_mut(&mut a).pe_adapter = Some(t1);
Arc::make_mut(&mut b).pe_adapter = Some(t2);
```

**Change** — hoist the two variants per chunk, or pass `pe_adapter` as an argument to `trim()` and keep `Params` immutable and shared.

**Effect** — `[E: 2 KB memcpy + alloc per firing pair; at 2-5% of 5M pairs = 0.2-0.5 GB of memcpy, ~0.05 s]`. Low, and it disappears for free if FP-2 is taken.

**Measure** — `--pe-overlap always` on a library with real readthrough (the Cell B cohort named at `cli.rs:236-244`); the bench corpus keeps it silent.

**Risk** — none; `pe_adapter` semantics unchanged.

### FP-11 — BGZF: 65,280-byte independent blocks with no cross-block dictionary
**class:** algorithmic (gzip pipeline), io-threads · **effort:** M · **byte-identical:** no (container changes; decoded bytes identical)

**Evidence** — `io.rs:124-128` hardcodes `Bgzf`; gzp's `Bgzf` reports `needs_dict() = false` (`gzp/src/deflate.rs:608`) with `BGZF_BLOCK_SIZE = 65280` (`gzp/src/bgzf.rs:20`), while gzp's `Gzip` returns `cfg!(feature = "any_zlib")` — it hands the previous block's trailing 32 KiB to libdeflate as a dictionary. BGZF cannot: every 65 KB block restarts the match finder.
```rust
let enc = ParCompressBuilder::<Bgzf>::new()
    .num_threads(levels.gz_threads.max(1))
    .compression_level(GzCompression::new(levels.gz))
    .from_writer(buf);
```

**Change** — expose the container: BGZF stays the default (tabix, htslib, the tool's own readers), add `--gzip-format mgzip` (gzp `Mgzip`, 128 KiB blocks, real dictionary) for when indexability is not wanted. Both stay valid `gzip -d`.

**Effect** — same input and level, fewer and larger blocks with cross-block matches `[E: a few percent ratio and less per-block table cost; io.rs:118-120 already records "a few percent larger" as BGZF's price]`. Asserted as testable, not as a speed win: the output size is not recorded anywhere today (`file_mb` in the CSV is the input size, identical across all three tools `[M]`).

**Measure** — one fixed `.fq` through the tool in both formats: bytes out and wall; `pigz -6 -b 256` on the same file as an independent reference. Never compare sizes against fastp without decoding both.

**Risk** — a second container is a support surface; keep BGZF the default.

### FP-12 — Input decompression is single-threaded for every non-BGZF gzip
**class:** io-threads · **effort:** L · **byte-identical:** yes (the byte stream into the parser is identical)

**Evidence** — `io.rs:172-180` takes the parallel path only for BGZF, and `gz_threads = threads/2` (`cli.rs:981`) is 1 at `-j 2`, so even BGZF input gets no parallelism at 2 cores:
```rust
if looks_like_bgzf(head) && gz_threads > 1 {
    let dec = ParDecompressBuilder::<Bgzf>::new().num_threads(gz_threads)...
} else {
    Ok(Box::new(BufReader::with_capacity(INPUT_BUFFER, file)))
}
```

**Change** — htslib-style parallel gunzip for ordinary `.gz` (find member/block boundaries with a bit-level inflate probe every N MB, inflate on a pool, reassemble in order). Cheap half-measure: size `gz_threads` from `threads - 1` and document that the parallel input path needs BGZF.

**Effect** — the input is ~954 MB gzipped at 5M pairs `[M]`; single-threaded inflate at 150-250 MB/s is 4-6 s of wall `[E]` and is the one stage that never overlaps with anything at any thread count. Relevant only at `-j 2`+, a cell fastprep already wins (25.84 vs 39.61 `[M, docs/BENCHMARKS.md]`), hence low rank.

**Measure** — 2-core cell only: `perf stat -e cycles` on the reader thread, or BGZF vs ordinary gzip input of identical content at `-j 2`.

**Risk** — a hand-rolled gzip block splitter is the classic silent-truncation bug; use gzp's `Mgzip` decompressor or libdeflate's block API, not a new one.

### FP-13 — `mismatches_upto`: auto-vectorisation by hope, and a float `ceil` per position
**class:** allocation (SIMD, branch layout) · **effort:** S · **byte-identical:** yes

**Evidence** — `adapter.rs:550-558` relies on LLVM, and `adapter.rs:497` recomputes a float `ceil` for every position:
```rust
for (ca, cb) in a.chunks(32).zip(b.chunks(32)) {
    for i in 0..ca.len() {
        mm += (ca[i] != cb[i]) as usize;
    }
    if mm > cap { return mm; }
}
```

**Change** — explicit `_mm256_cmpeq_epi8` + `movemask` + popcount on 32-byte blocks, and a per-adapter `max_mm[k]` table (k spans `min_overlap..=a`, < 64 entries) instead of `(k as f64 * eps).ceil()`. Check the cap once per 4 blocks.

**Effect** — `[E: ~10-15 cycles/call -> ~6-8]`; at 414 calls/read that is 0.3-0.6 us/read, but only *after* FP-3 cuts the call count. Do FP-3 first: SIMD on a loop that should not run is the wrong order.

**Measure** — `perf stat -e instructions` per read on a `--dry-run` cell; keep the 4,000-case reference test green.

**Risk** — intrinsics need a `target_feature` gate with a runtime-detect fallback, even though `.cargo/config.toml` already builds `target-cpu=native`.

### FP-14 — `TrimResult` carries a heap `Vec` and is cloned per read
**class:** allocation · **effort:** S · **byte-identical:** yes

**Evidence** — `pipeline.rs:492-496` clones on the ordered thread, and `TrimResult::Emit` holds `runs: Vec<(usize, usize)>` (`algorithm.rs:435-441`), empty for every read when the interior scan is off (the shipped default, `cli.rs:320-323`) yet enough to keep the enum non-`Copy` and to drag a `map().collect()` path at `algorithm.rs:1169-1175`:
```rust
let result = if dup {
    duplicate_result()
} else {
    p.output.result.clone()
};
```

**Change** — const-generic the enum on `RUNS: bool` (the default instantiation then becomes `Copy`), or store `runs` inline as `[(u32, u32); 2]` + length; take the result by reference instead of cloning.

**Effect** — `[E: 20-40 cycles/read on the ordered thread]`; small, but it is on the critical path and removes the last heap allocation from the steady-state trim output.

**Measure** — `perf stat -e instructions`, plus `--json` `bases_interior_mass` unchanged.

**Risk** — `--interior` must keep working; both instantiations need tests.

## 2. What is already good (and worth not breaking)

* **The two-tier cascade fires.** With `min_len 0.9L` at L = 151,
`rem_budget = 15 > w = 14`, so a clean read skips the whole Γ walk (`algorithm.rs:971-996`) and pays two 14-pair window checks plus one quality sum. Fragile by exactly one base: at `0.92L` the guard fails and every read walks. Anyone touching `min_len` should know that cliff exists.
* **The quality gate is a LUT**, no `powf` in the loop (`algorithm.rs:301-325`),
and direction is monomorphised (`streaming_walk::<FORWARD>`), replacing what the comment calls a 10%-of-runtime branch.
* **Chunk arenas**: `ScanInput` is three ranges, so the parse path has zero
per-read allocations (`pipeline.rs:43-80`).
* **`mimalloc`, `lto = "fat"`, `codegen-units = 1`, `target-cpu=native`**
(measured +12% in the comment in `.cargo/config.toml`).
* **Dedup is admission-checked and saturating**: 101 MB peak vs fastp's 1,143 MB
`[M]`, one-sided freeze by construction (`dedup.rs:196-245`), and a projection that is an identity rather than a guess (`memory.rs:243-269`).
* **Output is thread-count independent**: dedup and stats on the ordered thread,
locked down by `parallel_output_matches_serial`. That is what makes FP-1 safe to attempt.
* **Refusal over guessing**: PE ambiguity -> `None` (`pe_overlap.rs:207-213`),
adapter auto-detect never invents a sequence (`cli.rs:793-815`).

## 3. Dead ends / explicitly rejected

1. **Vectorising the polyX streaming walk.** `gamma` is a loop-carried,
data-dependent recurrence with an early exit on `gamma > mercy`: there is no prefix sum to vectorise and the LUT loads are gather-shaped. The cheap win was already taken (skip the walk entirely, tier 1), and a rewrite would add risk to the one function that decides every trim.
2. **Raising `--threads` (or the shipped default) to win the cell.** The cell
fixes `-j 1`; on a genuinely pinned core wall == CPU-seconds, so threads cannot help and would only inflate `cpu_s`. The lever there is FP-3/5/6/9; the parallel lever is FP-1 at the same thread count.
3. **Turning dedup off to match fastp.** `reads_dropped_duplicate` = 321,058 at
5M pairs `[M]` — 3.2% of reads and 41% of all drops. Removing it changes the output, the drop rate (0.0786 -> ~0.046) and the tool's purpose. It is a policy difference; the honest move is FP-6's single probe, not deletion.
4. **Dropping `--gz-level` to 4 to match fastp's default.** The benchmark
equalises at 6 on purpose (`docs/BENCHMARKS.md`: an earlier revision "charged fastprep for compression the others never paid for"). Lowering it wins a cell by changing the question and changes every output byte.
