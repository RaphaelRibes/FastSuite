# fastdedup — optimisation propositions

Scope: `tools/fastdedup/src` (3 376 lines, all files read) + the three integration suites.
The propositions below are Round 1's analysis; what landed is listed under Status of this round.

## Status of this round

Implemented and committed in `2dfaff2` (5 of 11 propositions). The tests are the gate: the crate's suite is green, and every change that can alter output was checked byte-for-byte against the pre-change binary except where the entry says otherwise.

### Implemented

- **FD-2** --- parallel gzip (gzippy) on the `-o -` stream path
- **FD-3** --- `worker_threads = num_workers.saturating_sub(1).max(1)` budgets the pool (was 1 + j + j + 1 threads)
- **FD-7** --- `BufWriter::with_capacity(1 MiB)` before the gzip encoder on the file path
- **FD-8** --- `ChanWriter` moves whole chunks instead of two copies
- **FD-9** --- `RawRecord.crlf` + `ending()`: serial and parallel output now agree on CRLF input

### Deferred, with the reason

- **FD-1** --- stop offloading the hash (it is < 1 % of a pair) and parallelise the serialisation instead: the ordered writer is the design constraint, a restructure not a patch
- **FD-4** --- SIMD hash is pointless until FD-1 moves the work
- **FD-5** --- lock-free set only matters once several threads insert
- **FD-6** --- read-side zero-copy is small against the 21 MB/s ceiling
- **FD-10** --- `--json` streaming is a documentation-size change
- **FD-11** --- mmap output conflicts with the ordered-writer design

---


Tags used throughout:
`[M]` measured, artifact named · `[E]` estimated from a named mechanism · `[V]` verified in this
review by executing the built binary or reading the cited line.

No fastdedup throughput artifact exists in this repo (`docs/BENCHMARKS.md` carries zero fastdedup
rows), so the only timed anchor is the suite-wide figure: **~15–17 µs per pair per core** `[M]`.
At 150 bp that is ~330 B per pair, i.e. ~21 MB/s — gzip-class. **Consequence: removing work that
touches bytes outranks removing instructions**, and the ranked list follows that rule.

## Where the time can possibly be

The threaded path is a three-stage channel pipeline (`parallel.rs:388-496`): one reader thread
copies each record into three `Vec`s, N workers hash it, the calling thread inserts into the table
and writes. Everything that touches bytes — parse, serialise, deflate — is on one thread or split
across a thread boundary that copies.

---

## FD-1 — The pipeline parallelises the one stage that cannot pay for itself
`algorithmic · L · output byte-identical`

Workers do nothing but hash; the table insert *and* the write are on the writer thread.

```rust
// parallel.rs:447-452 — the entire worker body
while let Ok(chunk) = chunk_rx.recv() {
    for record in chunk {
        let hash = T::hash_sequence(&record.seq);
```

xxh3 over 330 B is tens of nanoseconds against 15–17 µs — **under 1 %** `[E: 330 B at
memory-bandwidth rates]` — yet the pipeline pays 3 allocations per record, 2 channel hops and two
cache-cold handoffs to offload it, while parse, insert, serialisation and the compression feed stay
serial. `ShardedHashVerifier` (`parallel.rs:25-53`) is **only ever called from that writer loop**
(`parallel.rs:498`), so the sharding is dead weight plus a second hash per record
(`parallel.rs:51`) `[V]`.

**Change.** Cut the pipeline at the ordered-commit boundary: the commit thread runs parse → hash →
ordered insert → serialise survivors into a chunk arena (one `Vec<u8>` + `Vec<Span>`); workers do
the order-free per-byte work on the previous chunk. Verdicts never leave the commit thread, so
input order and the surviving representative are unchanged by construction.

**Effect.** Removes 6 allocations and 2 channel sends per pair `[E]`; removes the *unbounded*
`hash_rx` (`parallel.rs:391`), which can buffer a whole file when the compressor lags; makes
per-byte work parallel. New floor is the commit thread: parse + one table probe.

**Measure.** `taskset -c 0-3`, `-j {1,2,4}`, wall + VmHWM + `perf stat -e cycles,cache-misses,
context-switches`; a debug counter on `hash_tx.send` must fall from 2/pair to 0.

**Risk.** 135 tests pin this path; the ordered-commit form keeps `tests/cli.rs:858`
(`contenu_ser == contenu_par`) green by construction.

---

## FD-2 — Parallel gzip never runs on the stream path
`algorithmic · S · gzip bytes differ, decoded bytes identical`

```rust
// parallel.rs:192-195
if crate::stream::is_stream(path) {
    // D-4: stdout cannot be seeked or reopened, so there is no gzippy thread
    let writer = crate::writer::prepare_stream_writer(out_compression, compression_level);
```

The stated reason does not apply: `gzippy::compress_to_writer_with_threads` writes to any
`W: Write + Send` and never seeks, and the file branch already proves the channel pattern
(`parallel.rs:208-224`) `[V: signature read]`. So `-1 - -o - --out-compression gz -j 8` — the
invocation that replaces `fastp … | gzip` — compresses on **one** flate2 thread at any `-j`.

**Change.** Route streams through the same gzippy pipe, moving `stdout()` into the compression
thread, keeping the `level >= 6` guard (below 6 gzippy emits a non-standard multi-block format,
`parallel.rs:205`).

**Effect.** Stream compression goes from 1 deflate worker to `j − 1` over disjoint blocks. If
compression is the ~21 MB/s stage the µs/B ratio suggests `[E]`, this is the largest single-thread
cost in the tool, on a code path it already ships for files.

**Measure.** `-o - --out-compression gz | cat > /dev/null` at `-j 1/2/4` pinned; report the
`-j4/-j1` ratio and `perf stat -e task-clock` (>1 core busy); `gunzip -t` plus a byte compare of
the decoded output against a `-j 1` run.

---

## FD-3 — The thread budget oversubscribes the cores the benchmarks pin
`io-threads · S · output unchanged`

```rust
// parallel.rs:465
let gzip_threads = if num_workers > 1 { num_workers } else { 0 };
```

At `-j 4` single-end that is 1 reader (`:404`) + 3 workers (`:443-455`) + 4 gzippy workers + 1
gzippy coordinator (`:213`) = **8–9 runnable threads against a 4-core budget** `[V]`. Paired-end
adds `(num_workers / 2).max(1)` compression threads per output (`:698-702`).

**Change.** Budget the commit thread: `gzip_threads = num_workers.saturating_sub(1).max(1)`, or an
explicit `--compress-threads` with that default; under FD-1 the budget is `1 commit + (j−1) compress`.

**Effect.** Never more runnable threads than CPUs, so the thread holding a cache-hot output buffer
is not preempted, and deflate stops competing with the table insert for a core.

**Measure.** `perf stat -e context-switches,cpu-migrations` and `ctxt` in `/proc/self/status` around
a pinned `-j 4` gz-in/gz-out run; expect fewer switches and no wall regression.

---

## FD-4 — `--abundance` touches the data three times, and its "fast" path is the slower one
`algorithmic · L · output byte-identical`

Abundance needs the final count before the first header, so *some* second traversal is inherent —
but both current options are the worst ones: without `--temp-dir` the input gzip is re-opened and
re-decompressed (`processor.rs:556,580`); with it, uniques are re-serialised uncompressed behind a
256 KiB buffer (`processor.rs:443`), then re-parsed, re-allocated (3 `Vec`s per record,
`processor.rs:753-755`) and re-hashed `[V]`.

**Change.** Tier 1: pass 1 appends each unique record to an in-RAM arena (`Vec<u8>` + `(span, hash)`
spans); pass 2 walks spans — header formatted once, no parse, no alloc, no re-hash. Tier 2: at the
arena cap, flush to the temp file and keep appending; pass 2 copies spans with `read_at()` into the
same formatter. `--abundance-ram` (default a fraction of the discovered limit) bounds tier 1.

**Effect.** Tier 1 deletes one full gzip decompression and one full FASTX parse of the input; tier 2
deletes the second parse and second hash from the temp path.

**Measure.** `--abundance -j 1` on a 1 M-pair corpus, wall + VmHWM, both paths; `strace -c -f -e
trace=read,write` (tier 1 collapses the syscall count); byte-compare against the current binary.

---

## FD-5 — The classic abundance path keeps a second full-size hash set
`memory · M · output unchanged`

```rust
// processor.rs:581,595
let mut written = rustc_hash::FxHashSet::default();
...
if !written.insert(hash) {
```

`written` only remembers which uniques were emitted, after `tracker` already counted them:
`|unique|` entries grown by doubling over the same key set, duplicated in the paired classic path
(`processor.rs:890,915`) `[V]`.

**Change.** Use the tracker value as the emitted flag: read the count, emit, set the entry to 0; a
later occurrence reading 0 is a duplicate already emitted, since a tracked key starts at count ≥ 1
and the pass is a single forward traversal.

**Effect.** Removes one hash table of `|unique|` entries — `[E]` ~9 B/entry at 7/8 load, ≈180 MB at
20 M uniques before power-of-two rounding.

**Measure.** `/usr/bin/time -v` max RSS on a many-unique `--abundance` run without `--temp-dir`;
`heaptrack`/`massif` to attribute the drop.

**Risk.** `get_count` also returns 0 for absent keys, so the sentinel needs a distinct accessor and a
test that an absent key is never silently skipped.

---

## FD-6 — Table sizing comes from a file-size estimate; on stdin it is zero
`memory · M · `hash_bits`/`hash_capacity_estimate` in `--json` change`

```rust
// utils.rs:60-63
let path_ref = path.as_ref();
if !path_ref.exists() { return Ok(0); }
```

`-1 -` is not a path, so the estimate is 0, `get_hash_method(0, _)` returns **64-bit**
(`utils.rs:15-21` ← `main.rs:295-298`) and `try_reserve(0)` reserves nothing: a streamed run
silently gets the weakest fingerprint and grows by doubling `[V]`. Every run also reserves the full
estimate in one `try_reserve` before the first record (cap 500 M entries, `utils.rs:8`), so peak RSS
is paid at t=0.

**Change, three separable pieces.** (1) streams default to 128-bit and say so in `--json`;
(2) reserve in fixed chunks (e.g. 4 M entries) and top up at 7/8 load, so peak RSS tracks the real
key count; (3) print the projection and refuse above the discovered limit, copying fastprep's model
(`tools/fastprep/src/memory.rs:60`, `cli.rs:894-906`).

**Effect.** (1) removes a silent collision class on the pipe path; (2) removes first-touch faults for
the gap between estimate and reality (on a 4× repetitive input the divisor at `utils.rs:48-55` is 4×
too large `[E]`); (3) turns a late OOM kill into an upfront refusal with a number.

---

## FD-7 — `.gz` output: flate2 writes into a bare `File`
`io · S · output unchanged`

```rust
// writer.rs:101-104
let writer: Box<dyn Write> = if is_gz {
    Box::new(GzEncoder::new(file, Compression::new(compression_level)))
```

The plain branch is wrapped by `BufWriter` in the caller (`processor.rs:144`); the gz branch is not,
and `GzEncoder` dumps each finished block into the inner writer.

**Change.** `GzEncoder::new(BufWriter::with_capacity(1 << 20, file), …)`, and the same for gzippy's
output file (`parallel.rs:219`).

**Effect.** `[E]` deflate emits in ~32 KiB units, so a 1 MiB buffer cuts write syscalls and
short-write memmoves ~32× on compressed output.

**Measure.** `strace -c -f -e trace=write` (writes per MiB) + wall; `gunzip -t` unchanged.

---

## FD-8 — The gzippy channel copies every output byte twice and mallocs per chunk
`allocation · S · output unchanged`

```rust
// parallel.rs:94-95
fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
    self.tx.send(buf.to_vec()).map_err(|_| {
```

`ChanReader` then copies again into gzip's read buffer (`parallel.rs:123-126`) `[V]`. A 175 KiB
`BufWriter` hands over a slice, so every output byte is copied twice plus one 175 KiB malloc/free
per chunk, when ownership could simply move.

**Change.** Send `struct Chunk { buf: Vec<u8>, len: usize }`, `mem::take` it reader-side and slice
without copying, and return emptied `Vec`s through a small free-list channel.

**Effect.** Removes 2 memcpys of every output byte and one large allocation per buffer `[E]` — a
per-byte tax on the compression path, hence ranked below FD-2/FD-3.

---

## FD-9 — Serial and threaded output disagree on line endings `[V]`
`correctness-enabling · S · output changes for CRLF input only`

The serial path writes through needletail (`processor.rs:175-177`), which re-emits the record's own
line ending; the threaded writer hard-codes `\n` (`parallel.rs:506-510`).

Verified on a CRLF FASTQ with the built binary: `-j 1` emits `@r1\r\n…`, `-j 4` emits `@r1\n…`;
`cmp` fails at byte 4. `tests/cli.rs:858` asserts serial == parallel, but its fixtures are LF-only,
so the gap is untested.

**Change.** Carry the record's `LineEnding` (one byte in span metadata) into the threaded writer, or
normalise both paths to LF and document it. The first keeps every existing output byte-identical.

**Effect.** Removes a data-dependent `-j 1` vs `-j >1` difference, and keeps future
parallel-formatting work (FD-1) under an honest parity claim.

**Measure.** Add a CRLF fixture to `tests/cli.rs` and run the existing serial-vs-parallel comparison
on it — it fails today.

---

## FD-10 — `--dry-run` silently runs the serial engine, so `-s -j 4` measures one thread
`io-threads · S · output unchanged`

```rust
// parallel.rs:375-378
if cfg.dry_run {
    return crate::processor::execute_deduplication::<T>(...);
```

The comment defends correctness, which is right, but the consequence is that the cheapest way to
measure the engine (no writer) is the one path that ignores `-j`, and any benchmark that subtracts
`-s` from a full run compares 1 thread against N.

**Change.** Run the same pipeline with a sink at the commit stage, so `-s` exercises the real
reader/commit path at any `-j`.

**Effect.** Turns `-s` into a true engine measurement and gives the protocol below a clean
separator for the write+compress stage.

---

## FD-11 — One cache-cold table probe per record, exposed serially
`algorithmic · M · output byte-identical · GATED on the protocol below`

One `FxHashSet` insert per record: hashbrown re-hashes an already-uniform fingerprint and touches a
bucket that misses cache — and often the TLB — once the table exceeds L3. The surrounding loads and
stores overlap; the miss does not.

**Change.** A purpose-built open-addressed table (`Vec<u64>`/`Vec<u128>`, power-of-two, linear
probing, load 0.75), index taken directly from the fingerprint's low bits, inserts batched 16 at a
time with a software prefetch of the next batch's buckets. Fingerprints stay xxh3, so verdicts do
not move.

**Effect.** 16 serial DRAM latencies become 16 overlapped ones per batch `[E]`; the price is
density — ~10.7 B/entry at 64-bit against ~9 B/entry for hashbrown `[E]`.

**Measure.** Only once the protocol shows the insert dominates the `-s` run:
`perf stat -e cache-misses,LLC-load-misses,cycles` per record.

**Risk.** A hand-rolled table is where a dedup tool gets silent wrong answers: property-test it
against `FxHashSet` on random keys, keep the 135 tests, and report the memory regression honestly.

---

## Protocol — run this before coding any of it

No throughput fixture exists in `bench/fixtures` (largest is 2.5 MB). Build one: 1 M and 10 M pairs,
150 bp, ~30 % duplicates, from `bench/fixtures/make_fixtures.py` or `fastsim`. Same machine,
`taskset -c 0-{0,1,3}`, three repetitions, median.

| run | isolates |
|---|---|
| `-s -j N` | parse + hash + insert (no writer) |
| `-j N` plain out | + serialisation + write |
| `-j N` gz out | + compression |
| `-j N` gz out, `-o -` | + the stream path (FD-2) |

Record wall, VmHWM, minor faults, `perf stat -e cycles,instructions,cache-misses,context-switches,
page-faults`, `strace -c -e trace=read,write`, and `--json` (`reads_in`, `reads_out`,
`duplicate_rate`, `hash_bits`, `collision_estimate`) to prove the verdict did not change. A
proposition that cannot show its mechanism in this table does not ship.

---

## What is already good

- The **serial** hot loop is allocation-free per record: `seq()` is a `Cow` into needletail's buffer,
  then hash, insert, write (`processor.rs:156-184`). Leave it alone `[V]`.
- `try_reserve` with an explicit OOM message instead of a panic (`hasher.rs:21-27`), plus a
  documented pre-allocation cap (`utils.rs:6-8`).
- Collision risk is a dial with a printed probability, not a hidden assumption (`utils.rs:10-21`,
  `--json.collision_estimate`).
- Resume measures its truncation offset with a no-op writer (`writer.rs:199-212`) instead of
  re-serialising to disk.
- Refusals are specific, ordered and individually tested (`main.rs:194-228`,
  `tests/streaming_tests.rs:345-369`).
- Multi-thread gzip is gated at level ≥ 6 to avoid gzippy's non-standard format (`parallel.rs:205`):
  a compatibility landmine correctly avoided.

## Dead ends — do not re-try these

1. **SIMD or 2-bit-pack the hasher.** The hash is <1 % of a pair and xxh3 already runs at
   memory-bandwidth rates; it optimises the cheapest stage and adds an encoding that must stay in
   sync with every hash change.
2. **Bloom filter or sampled table to cut table memory.** A filter answers "possibly seen", not
   "seen": every false positive drops a unique read, silently, with a load-dependent loss. Output
   stops being byte-identical *or* reproducible. Legitimate only under `--dry-run`.
3. **Parallel-decompress the input with `gzippy::decompress_with_threads`.** Its API takes the whole
   archive as one `&[u8]`, so a 2.7 GB input must be resident before decoding starts, and stdin
   cannot be handed to it at all — the constant-memory promise traded for a stage that is not the
   bottleneck.
4. **Make `ShardedHashVerifier` real (N workers, one shard each) without FD-1.** Insert order across
   shards becomes nondeterministic, and duplicates differ in quality string
   (`tests/streaming_tests.rs:13`), so the surviving representative varies: output is neither
   byte-identical to `-j 1` nor stable across runs. Today the code is accidentally safe because only
   the writer calls `verify`.

---

## Round 2 (2026-10-21)

Re-priced on the post-fix tree (`30e8fb3` + `af268d4` + the 2026-10-21 wave): see `docs/optim/round2-fastdedup.md` (propositions) and `docs/optim/raw/round2/` (probe evidence). Round 1 above is history — deferred items there are re-priced, not repeated.
