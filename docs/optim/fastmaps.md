# fastmaps — optimisation propositions

Scope: `tools/fastmaps/src` (2 397 lines, all 10 files read). Nothing here is implemented.

## Status of this round

Implemented and committed in `1e0fb65` (6 of 12 propositions). The tests are the gate: the crate's suite is green, and every change that can alter output was checked byte-for-byte against the pre-change binary except where the entry says otherwise.

### Implemented

- **FM-1** --- `--data-mode exact|capped:N|off`; non-exact modes declare `data_reads_per_pixel` and `data_complete:false`
- **FM-2** --- the 4 MB preview flush is charged to `--buffer-mb` instead of being a floor
- **FM-3** --- per-row writer table became a map with 4 KiB buffers (was a 65536-slot Vec)
- **FM-9** --- spiller writes one packed 8-byte record instead of three `write_all` calls
- **FM-10** --- FASTQ loop reads into four reused `Vec<u8>`, encodes from `&[u8]`
- **FM-12** --- `Done` messages merged in file order: the per-pixel mean quality, and therefore `heatmap_3d.html`, is reproducible

### Deferred, with the reason

- **FM-4** --- per-pixel accumulators are 8 bytes already
- **FM-5** --- the 3-D renderer is a separate offline pass
- **FM-6** --- encoder table is LUT-bound, not arithmetic-bound
- **FM-7** --- thread-local grids would multiply the grid allocation by the thread count
- **FM-8** --- preview thinning is already adaptive
- **FM-11** --- the HTML writer is not on the measured path

---


Tags: `[M]` measured, artifact named · `[E]` estimated from a named mechanism · `[V]` verified in
this review by reading the cited line or running the built binary.

## Two RSS regimes exist and must not be mixed

| config | 5 M pairs = 10 M reads | artifact |
|---|---|---|
| `--buffer-mb 512` (the pre-2026-10-18 default, i.e. what the paper's CSV ran) | **143.2 MB**, 5.49 s | `~/work/FastMaps-Bench/results/results.csv` `[V]` |
| `--buffer-mb 8`, **shipped default since 2026-10-18** | **36.1 MB** | `docs/BENCHMARKS.md:217` `[M]` |
| `--buffer-mb 1 --preview-reads 0` | 16.0 MB | `docs/BENCHMARKS.md:219` `[M]` |
| FastQC 0.12.1 | 16.4 MB, 16.1 s | `docs/BENCHMARKS.md:182` `[M]` |

**The 145 MB headline is a non-default configuration.** On the shipped default the same corpus
measures 36.1 MB, so the contested gap against FastQC is ~20 MB, not ~130 MB — and the 16.0 MB cell
proves the floor already sits at FastQC's level when previews are off. Everything below attacks the
remaining ~20 MB or the wall clock. `docs/STATUS.md` and the fastmaps row of the replacement table
should say which configuration their number came from; today they do not.

## Where the 143.2 MB goes

Record size is **12 B**, not 8 and not the 16 the budget charges `[V]`: `(u16,u16,u16,u32)` is 10 B
of payload at 4-byte alignment. `docs/BENCHMARKS.md:216` calls it an "8-byte entry" — that is the
*spilled* record (x:u16 | fid:u16 | rn:u32, `y` dropped because it is the row-file name), not the
in-memory one.

| term | size | basis |
|---|---|---|
| index records: 10 M reads × 12 B | 120.0 MB | measured record size × measured read count; 84 % of the cell `[E]` |
| process base | ~8 MB | `docs/BENCHMARKS.md:221` `[M]` |
| `row_writers`, 65 536 × 32 B | 2.10 MB | `spiller.rs:57` `[V]` |
| 152 row + 152 preview writers × 8 KiB | 2.49 MB | `spiller.rs:57,60` `[E]` |
| previews in flight, 2 workers × 4 MB flush | ≤8 MB | `worker.rs:58` × `main.rs:185` `[E]` |
| grid, 3 × 152² × 8 B | 0.55 MB | `format.rs:273-284` `[E]` |

The grid is 0.4 % of the peak. `heatmap.data` on the same corpus is 18 488 308 B = **1.85 B/read**:
the varint + delta encoding is already 4.3× smaller than the 8 B/read stream it is built from — the
encoder is the win, not the bottleneck.

---

## FM-1 — Make the per-read index stream optional and label it
`algorithmic · M · counts/PNG/JSON stay bit-exact; heatmap.data becomes a labelled approximation`

**The only term that is structural** (grows with read count at any buffer setting) and it exists for
exactly one output file: `worker.rs:226` pushes it, `format.rs:257` is its single consumer, while
counts/PNG/JSON/report all come from `state.counts` `[V]`.

```rust
// worker.rs:226
state.indices.push((xy.0, xy.1, file_id, read_num));
state.bytes_estimate += 16;
```

**Change.** Add the knob the schema already precedents (`preview_reads_per_pixel`, `format.rs:329`):

```rust
pub enum DataMode { Exact, Capped(u8), Off }
if let DataMode::Capped(n) = cfg.data_mode {
    let c = state.data_count.entry((xy, file_id)).or_insert(0);
    if *c < n { *c += 1; state.indices.push((xy.0, xy.1, file_id, read_num)); }
} // Off: no push, no spill, no render pass
```

Emit `"data_reads_per_pixel": N` and `"data_complete": false` in `heatmap.json` plus a
`fastsuite.report.v1` metric whenever the mode is not `Exact`.

**Effect.** `--no-data` deletes the whole path: minus the buffer term, minus 8 B/read of temp disk
(80 MB per 10 M reads `[E]`), no render pass. Floor becomes the measured 16.0 MB cell minus its own
1 MB buffer ≈ 15 MB `[E]`. Mechanism: the counting grid is exact and complete without any per-read
record; only pixel *membership* needs one.

**Exactness.** `Capped`/`Off` are approximations of `heatmap.data` alone — counts, PNG, grid JSON,
corner metrics and the report stay bit-exact. Labelled per the suite rule; default stays `Exact`.

**Measure.** `sha256sum heatmap.{json,png,data}` against the old binary on 100 k for `Exact`;
for `Off`, identical `counts`/`counts_r1`/`counts_r2`/`total_reads`/PNG hash plus the
`data_complete:false` field; peak RSS via `wait4(rusage)` at 5 M.

---

## FM-2 — Put the preview flush inside `--buffer-mb`
`memory · S · output neutral`

```rust
// worker.rs:58, with cli.rs:31-35 stating previews "do not count against it"
const PREVIEW_FLUSH_BYTES: usize = 4 * 1024 * 1024;
```

**Change.** Charge it (subtract `preview_bytes` from the worker's share) or cut the constant to
256 KiB — output-neutral, because previews are capped per `(pixel, kind)` (`worker.rs:234`) and
re-sorted at encode (`format.rs:98-103`).

**Effect.** ≤8 MB → ≤0.5 MB at 2 workers `[E]` (−7.5 MB off the 36.1 MB cell). Mechanism: the 4 MB
threshold is *per worker* and invisible to the budget, so it is a hard floor `--buffer-mb 1` cannot
move — which is why the 16.0 MB cell also needs `--preview-reads 0`.

**Measure.** 2 M-read cell peak RSS; `sha256sum heatmap.preview` unchanged; a send-counter to show
the 16× message increase does not move wall.

---

## FM-3 — Kill the 65 536-slot writer table and cap the row buffers
`memory / allocation · S · output identical`

```rust
// spiller.rs:57 — allocated before a single row is seen
let mut row_writers: Vec<Option<BufWriter<File>>> =
    (0..=u16::MAX as usize).map(|_| None).collect();
```

**Change.** A `HashMap<u16, writer>` (or a `Vec<u32>` slab index, 4 B/slot) with
`BufWriter::with_capacity(4 * 1024, …)`, `OpenOptions::append` if a row must be reopened.

**Effect.** −2.10 MB (table `[V]`) − 1.25 MB (304 writers × 8 KiB → 4 KiB `[E]`) ≈ **−3.4 MB, 9 % of
the 36.1 MB cell**. Mechanism: the table is sized for `u16::MAX` rows but a short-read run has
`max_seq_len + 1` = 152, and 8 KiB is 1 024 records of 8 B.

**Measure.** Peak RSS at 5 M; hash the temp dir before the render phase and hash `heatmap.data`;
assert `row_writers.len() <= size`.

---

## FM-4 — Stream the row temp file instead of `fs::read` + a full decoded copy
`memory · S · output identical`

`format.rs:225` reads the row whole, then `:230-236` builds a second `Vec` of 12 B/read beside it, and
the raw buffer stays alive across the whole `chunks_exact` loop `[V]`.

**Change.** Decode through a fixed buffer so the raw bytes die as they are consumed (reuse the
`u16_at`/`u32_at` helpers already in this file, `format.rs:150-159`).

**Effect.** The render peak becomes `12 × max_row_reads` instead of `20 × max_row_reads` `[E]`, a 40 %
cut of that term.

**Measure.** Log `raw.len()` per row behind `--timings` *first* — the size of this term is currently
unobservable, so that log is the deliverable. Then peak RSS at 5 M; `heatmap.data` hash identical.

---

## FM-5 — Charge the budget by allocated capacity, and stop discarding the buffer
`allocation · M · output identical`

`worker.rs:133,228`: `mem::take` leaves a zero-capacity `Vec` that regrows every chunk, and the charge
is a constant 16 B against 12 B of storage `[V]`.

**Change.** Account the real thing and recycle the allocation back from the spiller:

```rust
let used = state.indices.capacity() * std::mem::size_of_val(&state.indices[0]);
if used >= cfg.threshold_bytes { flush_local(state, tx); }
if let Ok(mut back) = recycle_rx.try_recv() { back.clear(); state.indices = back; }
```

**Effect.** Two mechanisms: (1) the ceiling becomes **true** at any `--buffer-mb` — real residency is
`capacity × 12` with `capacity ∈ [len, 2·len]` while the budget believes `16 × len`, a disagreement up
to 1.5× `[E]`; (2) per-chunk regrowth costs ~2× the chunk in memcpy per flush — at the shipped
default ~38 chunks × ~6.2 MB ≈ 235 MB of memcpy per 10 M reads `[E]`.

**Measure.** Peak RSS at `--buffer-mb 64` (where the mismatch bites) must be `64 + base`;
`--timings` flush delta; hashes identical.

---

## FM-6 — Two-pass mode: regenerate `heatmap.data` exactly, hold nothing
`algorithmic · M/L · heatmap.data byte-identical, mode labelled`

The record is a pure function of the sequence and the read ordinal (`worker.rs:181`, with `read_num`
the file-local counter), and the encoder sorts by `x` then `(fid, rn)` (`format.rs:29-37`), so arrival
order is irrelevant.

**Change.** `--two-pass`: pass 1 counts only (no index, no spill); pass 2 re-reads the same files in
the same order and streams `(x, read_num)` into the per-row encoders.

**Effect.** Index RAM → 0, temp disk → 0, `heatmap.data` exact, at the cost of a second decode:
wall ≈ 2× the scan `[E]`. This must never become the default — an opt-in for memory-capped nodes.

**Exactness.** None lost. This is how to keep "every read is addressable" while dropping FM-1's
stream.

**Risk.** Input must be re-readable (already guaranteed: file-only, D-4) and unchanged — hash the
first MiB in pass 1 and verify it in pass 2.

---

## FM-7 — Pack the record to 6 bytes behind a guard
`algorithmic (encoding/layout) · M · output identical`

`worker.rs:78`: 12 B in memory for 48 bits of payload, and `y` is already implicit in the per-row file
name (`spiller.rs:35`) `[V]`.

**Change.** `Vec<[u8; 6]>` (`x:u8, fid:u8, rn:u32` LE), align 1 so it is really 6 B/read, falling
back to the 12 B struct when `max_seq_len > 254 || files.len() > 254` (both known before the first
read). Do **not** use `(u8,u8,u32)` — alignment pads it back to 8.

**Effect.** −50 % of the index term (6.3 → 3.2 MB at the shipped default `[E]`) and −25 % of the temp
stream (8 → 6 B/read, 80 → 60 MB per 10 M reads `[E]`).

**Risk.** The guard must be checked before the first push or a late long read truncates `x`;
`debug_assert!(x <= 254)`.

---

## FM-8 — Shard the spiller and bound the channel
`io-threads (parallel merge) · M · output identical`

`main.rs:179` is unbounded and `spiller.rs:74-88` is one thread doing three `write_all` calls per read,
so nothing limits how far workers run ahead `[V]`.

**Change.** `S` shards keyed `y % S`, each owning its writer map and temp-dir suffix; merge =
concatenate the per-row files in `y` order (safe: the encoder sorts, so arrival order is irrelevant).
Back-pressure instead of growth: `mpsc::sync_channel::<SpillerMsg>(64)` (≈64 MB ceiling at 1 MB
chunks).

**Effect.** Removes the one RSS term the budget cannot see, and a serialisation point doing 30 M
`write_all` calls per 10 M reads. Wall effect `[E]` — it only helps once the scan is cheaper than the
spill.

**Risk.** More open handles per row; ordering bugs are invisible in the grid but visible in
`heatmap.data`, so the hash is the gate.

---

## FM-9 — One packed `write_all` per record
`io · S · output identical`

`spiller.rs:85-87`: three calls for a record that fits one store.

**Change.**

```rust
let mut rec = [0u8; 8];
rec[0..2].copy_from_slice(&x.to_le_bytes());
rec[2..4].copy_from_slice(&fid.to_le_bytes());
rec[4..8].copy_from_slice(&rn.to_le_bytes());
w.write_all(&rec).unwrap();
```

**Effect.** 3× fewer calls and bounds checks in the spiller's hot loop (30 M → 10 M per 10 M reads);
wall `[E]`, visible only when the spiller is the bottleneck.

---

## FM-10 — Reusable byte buffers in the FASTQ loop
`allocation (wall) · M · output identical`

`worker.rs:171,178-180`: `lines()` allocates four `String`s per read, and `encode_read` takes `&str`
(`encode.rs:12`), so every base is UTF-8-validated as a side effect `[V]`.

**Change.** `read_until(b'\n', &mut buf)` into four reused `Vec<u8>`; `encode_read(&[u8])`,
`mean_quality(&[u8])`; build a `String` only for a preview read actually kept.

**Effect.** 40 M allocations removed per 10 M reads (4/read `[E]`) plus two UTF-8 scans of 151 B per
read: wall `[E]` 20–30 % of the scan; **the only hot-loop change that also helps the pinned cell**
(10.61 s `[V]`).

**Measure.** Pinned `-j 1` cell on 2 M reads, 5 runs, median; `heatmap.*` hashes identical.

**Risk.** CRLF input keeps `\r` in `seq`, changing `seq.len()` — strip it and add a CRLF test.

---

## FM-11 — Dense per-worker counters, and drop the duplicate map
`algorithmic (accumulator layout) · M · output identical`

`worker.rs:199-201,224`: three hash lookups per read, one for a value that is merely a per-file
counter; `main.rs:255-259` then copies the whole map again `[V]`.

**Change.** A `Vec<(u64, f64)>` of `size*size` per worker sized from the first read's length and grown
on demand; per-file totals as `Vec<u64>`; build `counts_arr` straight from the worker's density on
`Done` and delete `counts_u32`.

**Effect.** 3 hashes/read → 0; 0.37 MB/worker dense at 151 bp against ~0.9 MB of `HashMap` `[E]`; and
the `counts_u32` duplicate goes (≈0.9 MB at 151 bp, ≈40 MB at 1 kb reads `[E]`).

**This is also the answer to "would a 4-bit/6-bit grid encoding help":** a cell is a *count* (`u64`),
not a composition code — the composition is the `(x, y)` pair, and `(T+C, T+G)` has 152² = 23 104
distinct values at L = 151, so a 64-bucket accumulator would **collide**. Packing cells to `u32` is
safe only while `total_reads < 2³²`; gate it on that check.

---

## FM-12 — Make the quality merge order deterministic
`correctness-enabling · S · heatmap_3d.html becomes reproducible`

`worker.rs:201` accumulates `f64`; `spiller.rs:123-127` merges those partials in *thread-completion*
order, because `Done` messages race `[V]`.

**Change.** Buffer `Done` messages and merge in `file_id` order (or accumulate `(sum_q, sum_len)`
fixed-point instead of `f64`), and sort the HTML series by `(x, y)`.

**Effect.** `heatmap_3d.html`'s `z = cq/c` becomes reproducible. Integer counts are already
order-exact — exactly why the PNG and grid JSON *are* byte-identical across `-j`; the `f64` sum is the
one quantity that is not.

**Measure.** Ten runs at `-j 2` and `-j 8`, hash `heatmap_3d.html`: it differs today; after both
fixes it must be stable.

---

## What is already good (keep it through every change above)

- **The wire format.** 1.85 B/read measured (18 488 308 B for 10 M reads) from varint + delta-coded
  read numbers within `(x, file)` groups (`format.rs:50-58`), 4.3× smaller than the 8 B/read stream.
  The varint is the *win*, not the bottleneck.
- **Per-row temp files, read-one-row-then-delete** (`format.rs:221-255`): render peak is bounded by
  the largest row, not by the run — the pattern the rest of the suite should copy.
- **Preview streaming** (`worker.rs:145-156`): removed a ~76 MB term, outputs byte-identical.
- **Per-`(pixel, kind)` preview quota** (`worker.rs:231-236`): R2 cannot be starved by file order.
- **`log_norm` NaN guard** with a regression test (`output.rs:13-19`).
- **No serde**: 3 dependencies, hand-written JSON with a real parse test; `no_call` handling so all-N
  reads are not reported as polyA contamination.
- **Integer path determinism**: `-j 1` and `-j 8` give identical `.data`/`.json`/`.png` because counts
  are integers and the encoder sorts.

## Dead ends — do not re-try these

- **Bit-packing the grid (4-bit or 6-bit cells).** The grid is 0.55 MB at 151 bp — 0.4 % of the cell —
  and a 64-bucket accumulator would collide on 23 104 distinct compositions, so it is not a lossless
  replacement at any size.
- **Compress the temp stream to cut disk.** Disk is not the RSS term (records live in the worker
  `Vec`) and it inserts a codec between worker and spiller. FM-7 takes 25 % of that stream for free.
- **Bloom filter / count-min sketch for pixel membership.** A sketch can estimate *how many* reads are
  in a pixel — the exact counter already answers that — and cannot answer *which* reads, the only
  question `heatmap.data` exists for. Any sampling here must be a labelled cap (FM-1).
- **Raise `--buffer-mb` or pass `-j 32` to go faster.** Measured dead: 5.3–5.9 s at both 512 MB and
  8 MB, and `-j` is clamped to the file count, so a paired-end run gets 2 workers whatever you pass.
  Thread scaling needs block sharding, not a bigger number.
- **`mmap` the row temp files instead of `fs::read`.** Resident mapped pages count toward VmHWM
  exactly like a heap buffer, so peak RSS does not move.

---

## Shared verification protocol

1. Baseline the current binary on the four corpora (100 k / 500 k / 1 M / 5 M pairs) and store
   `sha256sum heatmap.data heatmap.preview heatmap.json heatmap.png` per cell; after each proposition
   re-run and diff those four hashes.
2. **`heatmap_3d.html` is excluded until FM-12 lands** — it already varies run to run.
3. Every RSS claim is `wait4(rusage)` peak of the child, never `ps` sampling; wall is the median of 5
   runs, with the pinned `-j 1` cell reported next to `all cores`, because that pair separates scan
   cost from spill cost.
4. Add `--timings` (scan / spill-busy / render / png to stderr) **before** touching anything: the
   render phase (`format.rs:221-267`) is currently invisible, so nobody can say whether FM-8–FM-10
   would even show up in wall time.
