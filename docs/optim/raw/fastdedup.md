# fastdedup 1.3.0 — ranked optimization dossier

All of `tools/fastdedup/src/*.rs` (3376 lines) and the three integration suites were read; nothing was built,
changed or committed. Figures are tagged **[measured]**, **[counted]** (counted in code, not timed) or
**[estimated]** with the mechanism that produces them. No fastdedup throughput artifact exists in this repo
(`docs/BENCHMARKS.md` has zero fastdedup rows), so the only timed anchor is the suite figure below and every
proposition carries its own measurement recipe (§Protocol).

## Baseline, and what it implies

- **[measured]** ~15–17 µs per pair per core across the suite (anchor supplied with this task).
- **[counted]** ~330 B per pair at 150 bp (`utils.rs:39`: `2 * read_length + 42`). 15.5 µs over 330 B is ~21
  MB/s — gzip-class — so removing *bytes* of work outranks removing *instructions*; §Protocol decides which
  branch applies before anything is coded.
- **[counted]** the threaded path allocates 3 `Vec`s per record (`parallel.rs:419-421`): 6 per pair, 8 for
  paired-end (`parallel.rs:642-649` plus two `base_read_id().to_vec()` at 630-631), two channel sends each.
- **[counted]** serial abundance re-reads and re-decompresses the input twice (`processor.rs:556`, `580`);
  `--temp-dir` instead writes a full uncompressed copy and re-parses it.
- **[estimated]** the verdict cannot be reordered: the survivor is the first in input order, and duplicates
  differ in quality (`tests/streaming_tests.rs:13`: same `ACGTACGT`, `IIIIIIII` vs `!!!!!!!!`), so "any member
  of the duplicate class" is not byte-identical. Hence D4.

## FD-1 — The pipeline parallelizes the one step that cannot pay for itself
`class: algorithmic · effort: L · output: byte-identical`

Workers only hash; the table insert and the write are both on the writer thread (`parallel.rs:496-498`).

```rust
// parallel.rs:447-452 — the whole worker body
while let Ok(chunk) = chunk_rx.recv() {
    for record in chunk {
        let hash = T::hash_sequence(&record.seq);
        if hash_tx.send(HashedRecord { hash, record }).is_err() {
```

Hashing ~330 B at memory-bandwidth rates is tens of nanoseconds against **[measured]** 15–17 µs, under 1% of
a pair — yet the pipeline spends 6 allocations, 2 channel sends and two cache-cold handoffs to offload it,
while parse, insert, serialization and the compression feed stay serial. `ShardedHashVerifier`
(`parallel.rs:25-53`) is never touched concurrently (`verify` is called only in that writer loop), so the
sharding is dead weight plus a second hash per record (`parallel.rs:51`).

- Change: cut the pipeline at the ordered boundary. Commit thread: parse → hash → ordered insert → serialize
  survivors into a chunk arena (one `Vec<u8>` + `Vec<Span>` per chunk). Workers: speculative hash of chunk
  k+1, then the order-free per-byte work on chunk k (header/qual formatting, hand-off to the compressor).
  Verdict never leaves the commit thread, so order and representative are unchanged by construction. Cheaper
  intermediate: keep parse+hash+insert on one thread and give the whole `-j` budget to compression (FD-2/3).
- Effect: removes 6 allocations and 2 channel sends per pair **[counted]**; removes the unbounded `hash_rx`
  (`parallel.rs:391`), which can buffer a whole file when the compressor lags; makes per-byte work parallel.
  New floor = commit thread = parse + one table probe.
- Measure: `taskset -c 0-3`, `-j {1,2,4}`, wall clock, max RSS, `perf stat -e cycles,cache-misses,
  context-switches`, and a debug counter of `hash_tx.send` (must go 2/pair → 0).
- Risk: medium — 104 tests pin this path; the ordered-commit form keeps `tests/cli.rs:858`
  (`contenu_ser == contenu_par`) green by construction.

## FD-2 — Parallel gzip never runs on the stream path, which is the fastp-replacement path
`class: algorithmic · effort: S · output: gzip bytes differ (same decoded bytes)`

```rust
// parallel.rs:192-195
if crate::stream::is_stream(path) {
    // D-4: stdout cannot be seeked or reopened, so there is no gzippy
    // thread and no resume — serial flate2 gzip when asked for.
    let writer = crate::writer::prepare_stream_writer(out_compression, compression_level);
```

The stated reason does not apply: `gzippy::compress_to_writer_with_threads` writes to any `W: Write + Send`
and never seeks (`gzippy-0.8.0/src/lib.rs:142`), and the file branch already proves the channel pattern
(`parallel.rs:208-224`). So `-1 - -o - --out-compression gz -j 8`, the invocation that replaces
`fastp … | gzip`, compresses on one flate2 thread at any `-j`.

- Change: route streams through the same gzippy pipe, moving `std::io::stdout()` into the compression thread,
  keeping the level >= 6 guard (`parallel.rs:205`, rationale at 168-172: below 6 with threads, gzippy emits a
  non-standard multi-block format).
- Effect: stream compression goes from 1 deflate worker to `j - 1` over disjoint blocks. If compression
  dominates at ~21 MB/s **[estimated from the µs/B ratio above]**, this is the largest single-thread cost in
  the tool that already has a vendored parallel implementation.
- Measure: `-o - --out-compression gz | cat > /dev/null` at `-j 1/2/4` pinned; report the -j4/-j1 ratio and
  `perf stat -e task-clock` (>1 core busy); `gunzip -t` plus byte-compare decoded output with `-j 1`.
- Risk: low — same crate, same guard, a thread that only holds stdout; `ChanWriter`'s 16-slot bounded channel
  is the intended backpressure when the consumer is slow.

## FD-3 — Thread budget oversubscribes the cores the benchmarks pin
`class: io-threads · effort: S · output: unchanged`

```rust
// parallel.rs:465
let gzip_threads = if num_workers > 1 { num_workers } else { 0 };
```

At `-j 4` single-end that is 1 reader (`parallel.rs:404`) + 3 workers (`parallel.rs:443-456`) + 4 gzippy
workers + 1 gzippy coordinator (`parallel.rs:213`): 8–9 runnable threads against a budget the suite pins to
1/2/4 cores with `taskset`. Paired-end adds `(num_workers / 2).max(1)` per output (`parallel.rs:698-702`).

- Change: budget the reader/commit thread — `gzip_threads = num_workers.saturating_sub(1).max(1)`, or an
  explicit `--compress-threads` with that default; under FD-1 the budget is `1 commit + (j-1) compress`.
- Effect: never more runnable threads than CPUs, so no involuntary preemption of a thread holding a cache-hot
  output buffer, and the deflate workers stop competing with the table insert for a core.
- Measure: `perf stat -e context-switches,cpu-migrations` and `ctxt` in `/proc/self/status` around a pinned
  `-j 4` gz-in/gz-out run; expect fewer switches, no wall-clock regression.
- Risk: low; neutral on unpinned runs with spare cores.

## FD-4 — `--abundance` touches the data three times, and its "fast" path is the slower one
`class: algorithmic · effort: L · output: byte-identical`

The temp path re-serializes every unique record in pass 1, then re-parses, re-allocates and re-hashes it in
pass 2 (`processor.rs:778-779`).

```rust
// processor.rs:753-755
id: rec.id().to_vec(),
seq: rec.seq().to_vec(),
qual: rec.qual().map(|q| q.to_vec()),
```

Abundance needs the final count before the first header, so *some* second traversal is inherent — but both
current options are the worst ones: without `--temp-dir` the input gzip is re-opened and re-decompressed
(`processor.rs:556`, `580`); with it, uniques are written uncompressed behind a 256 KiB buffer
(`processor.rs:443`) then parsed again, hashed again, allocated again (6 `Vec`s per pair).

- Change, tier 1: pass 1 appends each unique record to an in-RAM arena (`Vec<u8>`) with `(span, hash)` spans;
  if the arena stays under budget, pass 2 walks spans — header formatted once, no parse, no alloc, no re-hash.
  Tier 2: at the arena cap, flush to the temp file and keep appending; pass 2 copies spans with `read_at()`
  into the same formatter. Both emit today's bytes (same records, same first-occurrence order, same `;size=`).
- Effect: tier 1 deletes one full gzip decompression and one full FASTX parse of the input; tier 2 deletes the
  second parse and second hash from the temp path (6 allocations + 1 xxh3 per pair). A new `--abundance-ram`
  (default a fraction of the discovered limit, FD-6) bounds tier 1, so memory is not an unbounded regression.
- Measure: `--abundance -j 1` on the 1 M-pair corpus, wall clock + max RSS, both paths; `strace -c -f -e
  trace=read,write` (tier 1 collapses the syscall count); byte-compare against the current binary on the
  `tests/abundance_tests.rs` corpus.
- Risk: medium — single-threaded by design (`main.rs:327-331`) so no concurrency risk, but the spill boundary
  is a new failure mode (disk full mid-arena).

## FD-5 — The classic abundance path keeps a second full-size hash set
`class: memory · effort: M · output: unchanged`

```rust
// processor.rs:581, 595
let mut written = rustc_hash::FxHashSet::default();
...
if !written.insert(hash) {
```

`written` only remembers which uniques were emitted, after `tracker` already counted them: `|unique|` entries
grown by doubling over the same key set, duplicated in the paired classic path (`processor.rs:890`, `915`).

- Change: use the tracker value as the emitted flag — read the count, emit, set the entry to 0; a later
  occurrence reading 0 is a duplicate already emitted, since a tracked key starts at count >= 1 and the pass
  is a single forward traversal: `let count = tracker.take_count(hash); if count == 0 { continue; }`
- Effect: removes one hash table of `|unique|` entries. **[estimated]** at ~9 B/entry (hashbrown keeps one
  control byte per slot at 7/8 load — the same `bits / 8 + 1` fastprep projects, `tools/fastprep/src/cli.rs:906`),
  20 M uniques ≈ 180 MB, ~360 MB once the table rounds up to a power of two.
- Measure: `/usr/bin/time -v` max RSS on a many-unique `--abundance` run without `--temp-dir`;
  `heaptrack`/`massif` to attribute the drop to the removed set.
- Risk: low but subtle — `get_count` also returns 0 for absent keys, so the sentinel needs a distinct accessor
  and a test that an absent key is never silently skipped.

## FD-6 — Table sizing comes from a file-size estimate; for stdin it is zero, forcing 64-bit
`class: memory · effort: M · output: hash_bits / hash_capacity_estimate in --json change`

```rust
// utils.rs:60-63
let path_ref = path.as_ref();
if !path_ref.exists() {
    return Ok(0);
}
```

`-1 -` is not a path, so the estimate is 0, `get_hash_method(0, _)` returns 64-bit (`utils.rs:15-21`, called
from `main.rs:295-298`) and `try_reserve(0)` (`hasher.rs:22-24`) reserves nothing: a streamed run silently
gets the weakest fingerprint and grows by doubling. Every run also reserves the *full* estimate in one
`try_reserve` before the first record (cap 500 M entries, `utils.rs:8`), so peak RSS is paid at t=0.

- Change, three separable pieces: (1) streams default to 128-bit and say so in `--json`; (2) reserve in fixed
  chunks (e.g. 4 M entries) and top up at 7/8 load, so peak RSS tracks the real key count; (3) print the
  projection and refuse above the discovered limit, copying fastprep's model (`tools/fastprep/src/memory.rs:60`,
  `tools/fastprep/src/cli.rs:894-906`) per D-1.
- Effect: (1) removes a silent collision class on the pipe path; (2) removes first-touch faults for the gap
  between estimate and reality — on a 4x-repetitive 20 GB input the divisor at `utils.rs:48-55` is 4x too large
  **[estimated]**; (3) turns a late OOM kill into an upfront refusal with a number.
- Measure: max RSS and `perf stat -e page-faults` at three duplication levels on the same corpus; the `--json`
  diff is the only output change.
- Risk: medium — chunked growth costs a few rehashes, and the refusal is a new user-visible error that must be
  tested as one (cf. `tools/fastcheck/tests/cli.rs:383`).

## FD-7 — `.gz` output: flate2 writes into a bare `File`, one syscall and one `drain` per block
`class: io · effort: S · output: unchanged`

```rust
// writer.rs:101-104
let writer: Box<dyn Write> = if is_gz {
    Box::new(GzEncoder::new(file, Compression::new(compression_level)))
} else {
```

The plain branch is wrapped by `BufWriter::with_capacity(buf_size, …)` in the caller (`processor.rs:144`); the
gz branch is not. flate2's `GzEncoder` dumps its output buffer straight into the inner writer and then
`self.buf.drain(..n)` (`flate2-1.1.9/src/zio.rs:238-252`) — with an unbuffered `File` every block is a syscall
and a short write leaves a memmove of the remainder.

- Change: `GzEncoder::new(BufWriter::with_capacity(1 << 20, file), …)` in `prepare_writer`, and the same for the
  gzippy output file (`parallel.rs:219`) if it is fed small blocks.
- Effect: **[estimated]** deflate emits in ~32 KiB units, so a 1 MiB `BufWriter` cuts write syscalls and
  short-write memmoves ~32x on compressed output and removes the partial-write retry path.
- Measure: `strace -c -f -e trace=write` (writes per MiB) plus wall clock; `gunzip -t` unchanged.
- Risk: very low; 1 MiB RSS per open output (2 MiB paired-end).

## FD-8 — The gzippy channel copies every output byte twice and mallocs per chunk
`class: allocation · effort: S · output: unchanged`

```rust
// parallel.rs:94-95
fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
    self.tx.send(buf.to_vec()).map_err(|_| {
```

`ChanReader` then copies again into gzip's read buffer (`parallel.rs:123-126`, `buf[..n].copy_from_slice(…)`).
`BufWriter` (175 KiB at 150 bp, `utils.rs:38-42`) hands over a slice, so the output stream is copied twice plus
one 175 KiB malloc/free per chunk, when ownership could simply move.

- Change: send `struct Chunk { buf: Vec<u8>, len: usize }`, `mem::take` it reader-side and slice without
  copying, and return emptied `Vec`s through a small free-list channel so the allocator reaches steady state.
- Effect: removes 2 memcpys of every output byte and one 175 KiB allocation per 175 KiB of output **[counted]** —
  a per-byte tax on the compression path, hence ranked below FD-2/FD-3.
- Measure: `perf stat -e cycles,instructions` and wall clock on gz output at `-j 2/4`; a debug counter of `Chunk`
  allocations must reach zero once the free list warms.
- Risk: low — `Read` semantics (partial fills, EOF) must be exact; `tests/streaming_tests.rs:325` covers stdout.

## FD-9 — Serial and threaded output disagree on line endings
`class: correctness-enabling · effort: S · output: yes, for CRLF input only`

```rust
// parallel.rs:506-507
buffered_writer.write_all(b"\n")?;
buffered_writer.write_all(&hr.record.seq)?;
```

The serial path calls `seq_record.write(&mut buffered_writer, None)` (`processor.rs:175-177`), and needletail
re-emits the record's own `line_ending` (`record.rs:202-220`). `tests/cli.rs:858` asserts serial == parallel
byte-for-byte, but its fixtures are LF-only, so the gap is untested: CRLF input round-trips as CRLF under `-j 1`
and as LF under `-j 4`. Id, seq, `+` and qual already match.

- Change: carry the record's `LineEnding` (one byte in span metadata) into the threaded writer; or normalise both
  paths to LF and document it. The first keeps every existing output byte-identical.
- Effect: removes a data-dependent `-j 1` vs `-j >1` output difference for one byte per record, and keeps future
  parallel-formatting work (FD-1) under an honest parity claim.
- Measure: add a CRLF fixture to `tests/cli.rs` and run the existing serial-vs-parallel comparison on it (it fails
  today); `cmp` the two outputs.
- Risk: low; LF normalisation instead would be an observable output change needing a release note.

## FD-10 — `--dry-run` silently runs the serial engine, so `-s -j 4` measures one thread
`class: io-threads · effort: S · output: unchanged`

```rust
// parallel.rs:375-378
// --dry-run writes nothing: the serial counting pass reaches the identical
// verdict (same refusals, same duplicate count) without creating a file.
if cfg.dry_run {
```

The comment defends correctness, which is right, but the consequence is that the cheapest way to measure the
engine (no writer) is the one path that ignores `-j`: users scanning a big file with `-s` get one thread, and any
benchmark that subtracts `-s` from a full run compares 1 thread against N.

- Change: run the same pipeline with a sink at the commit stage (drop winners and losers after the insert), so
  `-s` exercises the real reader/commit path at any `-j`.
- Effect: turns `-s` into a true engine measurement and gives §Protocol a clean separator for the write+compress
  stage (full run minus `-s` run). No output change.
- Measure: `-s -j {1,4}` must scale, and `--json` metrics must stay identical to the serial run — that identity is
  the correctness gate.
- Risk: low; stream refusals (`--pre-count` on stdin, `main.rs:209-213`) are untouched.

## FD-11 — One cache-cold table probe per record, exposed serially (gated on §Protocol)
`class: algorithmic · effort: M · output: byte-identical`

```rust
// hasher.rs:31-33
pub fn verify(&mut self, hash: T) -> bool {
    self.memory.insert(hash)
}
```

One `FxHashSet` insert per record: hashbrown re-hashes an already-uniform fingerprint (the same double hash as
`parallel.rs:51`) and touches a bucket that misses cache and often the TLB once the table exceeds L3. The
surrounding loads and stores overlap; the miss does not.

- Change: a purpose-built open-addressed table — `Vec<u64>`/`Vec<u128>`, power-of-two, linear probing, load 0.75
  (0.875 at 128-bit), index taken directly from the fingerprint's low bits, inserts batched 16 at a time with a
  software prefetch of the next batch's buckets. Fingerprints stay xxh3, so verdicts do not move.
- Effect: 16 serial DRAM latencies become 16 overlapped ones per batch **[estimated]**; the price is density —
  10.7 B/entry at 64-bit vs ~9 B/entry for hashbrown **[estimated]**, ~15% more memory at 64-bit, 1.8x at 128-bit
  if the load factor were left at 0.75.
- Measure: only once §Protocol shows the insert dominates the `-s` run — `perf stat -e cache-misses,
  LLC-load-misses,cycles` per record.
- Risk: medium — a hand-rolled table is where a dedup tool gets silent wrong answers; needs a property test against
  `FxHashSet` on random keys, the existing 104 tests, and an honest memory regression number.

## Protocol (run before coding any of this)

No throughput fixture exists in `bench/fixtures` (largest is 2.5 MB), so build one: 1 M and 10 M pairs, 150 bp,
~30% duplicates, from `bench/fixtures/make_fixtures.py` patterns or `fastsim`. Same machine, `taskset -c 0-{0,1,3}`,
three repetitions, median.

| run | isolates |
|---|---|
| `-s -j N` | parse + hash + insert (no writer) |
| `-j N` plain out | + serialization + write |
| `-j N` gz out | + compression |
| `-j N` gz out, `-o -` | + the stream path (FD-2) |

Per run record wall clock, `/usr/bin/time -v` max RSS and minor faults, `perf stat -e
cycles,instructions,cache-misses,context-switches,page-faults`, `strace -c -e trace=read,write`, and `--json`
(`reads_in`, `reads_out`, `duplicate_rate`, `hash_bits`, `collision_estimate`) to prove the verdict did not change.
A proposition that cannot show its mechanism in this table does not ship.

## What is already good

- The serial hot loop is allocation-free per record: `seq()` is a `Cow` into needletail's 64 KiB buffer, then hash,
  insert, write (`processor.rs:156-184`). Leave it alone.
- `try_reserve` with an explicit OOM message instead of a panic (`hasher.rs:21-27`), plus a documented
  pre-allocation cap (`utils.rs:6-8`).
- Collision risk is a dial with a printed probability, not a hidden assumption (`utils.rs:10-21`, `main.rs:400-403`,
  `--json.collision_estimate`).
- Resume measures its truncation offset with a no-op writer (`writer.rs:199-212`) instead of re-serializing to disk:
  cheap and correct.
- Refusals are specific, ordered and individually tested so a generic error cannot shadow a real one
  (`main.rs:194-228`, `tests/streaming_tests.rs:345-369`).
- Multi-thread gzip is gated at level >= 6 to avoid gzippy's non-standard multi-block format (`parallel.rs:205`):
  a compatibility landmine correctly avoided.
- Write-buffer sizing follows `--read-length` with sane floor/ceiling, which matters for ONT (`utils.rs:38-42`).

## Dead ends / explicitly rejected

**D1. SIMD or 2-bit-pack the hasher.** The hash is <1% of a pair (~330 B against **[measured]** 15–17 µs) and xxh3
already runs at memory-bandwidth rates; it would optimize the cheapest stage and add an encoding that must stay in
sync with every hash change.

**D2. Bloom filter or sampled table to cut table memory.** A filter answers "possibly seen", not "seen": every false
positive drops a unique read, silently, with a loss that depends on load — output stops being byte-identical or even
reproducible. Legitimate only inside `--dry-run`, where the estimate is the product.

**D3. Parallel-decompress the input with `gzippy::decompress_with_threads`.** The API takes the whole archive as one
`&[u8]` (`gzippy-0.8.0/src/lib.rs:190-215`), so a 2.7 GB gz input must be resident before decoding starts, and stdin
cannot be handed to it at all — the streaming, constant-memory promise traded for a stage that is not the bottleneck.

**D4. Use `ShardedHashVerifier` as intended: N workers, one shard each.** Insert order across shards is then
nondeterministic, so which member of a duplicate class survives varies — and duplicates differ in quality string
(`tests/streaming_tests.rs:13`), so output is neither byte-identical to `-j 1` nor stable across runs. Today the code
is accidentally safe because only the writer calls `verify` (`parallel.rs:498`); making the shards real without FD-1's
ordered commit converts that safety into a bug.
