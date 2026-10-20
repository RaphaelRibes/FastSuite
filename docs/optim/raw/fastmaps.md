# fastmaps 0.1.0 — ranked optimization dossier

Read-only review of `tools/fastmaps/src/` (2 397 lines, all 10 files read). Nothing under `src/` was modified, nothing built, nothing committed. Every number is tagged **[measured]** (artifact named) or **[est.]** (mechanism stated). Two measured RSS regimes exist and must not be mixed:

| config | 5 M pairs = 10 M reads | artifact |
|---|---|---|
| `--buffer-mb 512` (the bench/paper binary) | 143.2 MB, 5.49 s (pinned: 142.7 MB, 10.61 s) | `~/work/FastMaps-Bench/results/results.csv` **[measured]** |
| `--buffer-mb 8`, **shipped default since 2026-10-18** | 36.1 MB | `docs/BENCHMARKS.md:217` **[measured]** |
| `--buffer-mb 1 --preview-reads 0` | 16.0 MB | `docs/BENCHMARKS.md:219` **[measured]** |
| FastQC 0.12.1 | 16.4 MB, 16.1 s | `docs/BENCHMARKS.md:182` **[measured]** |

The 145 MB headline is the 512 MB-buffer configuration; the shipped default measures 36.1 MB on the same corpus, so the contested gap is ~20 MB, not ~130 MB — and the 16.0 MB cell proves the floor already sits at FastQC's level. Everything below attacks the remaining ~20 MB or the wall clock.

## Where the 143.2 MB goes

Record size is 12 B, not 8 and not the 16 the budget charges; `Option<BufWriter<File>>` = 32 B, `PreviewRead` = 88 B **[measured]**, scratch probe `~/.hermes/cache/scratch/fmprobe/probe.rs`.

| term | size | basis |
|---|---|---|
| index records: 10 M reads × 12 B | 120.0 MB | measured record size × measured 10 M reads (2 × 20 M lines); 84 % of the cell **[est.]** |
| process base | ~8 MB | `docs/BENCHMARKS.md:221` **[measured]** |
| `row_writers`, 65 536 × 32 B | 2.10 MB | `spiller.rs:57` **[measured]** |
| 152 row + 152 preview writers × 8 KiB | 2.49 MB | `spiller.rs:57,60` **[est.]** |
| previews in flight, 2 workers × 4 MB flush | ≤8 MB | `worker.rs:58` × `main.rs:185` **[est.]** |
| grid, 3 × 152² × 8 B | 0.55 MB | `format.rs:273-284`, 152 = 151 bp + 1 **[measured]** |

The grid is 0.4 % of the peak. `heatmap.data` on the same corpus is 18 488 308 B = **1.85 B/read** **[measured]** (`.bench_tmp/5000k/fastmaps/default/heatmap.data`): the varint+delta encoding is already 4.3× smaller than the 8 B/read stream it is built from.

---

## FM-1 — Make the per-read index stream optional and label it
**class:** algorithmic — the only term that is structural (grows with read count at any buffer setting) and it exists for one output file.

**evidence** `src/worker.rs:226`; `format.rs:257` is its single consumer, while counts/PNG/JSON/report all come from `state.counts`:
```rust
// Index entry for heatmap.data.
state.indices.push((xy.0, xy.1, file_id, read_num));
state.bytes_estimate += 16;
```
**change** add the knob the schema already precedents (`preview_reads_per_pixel`, `format.rs:329`):
```rust
pub enum DataMode { Exact, Capped(u8), Off }
if let DataMode::Capped(n) = cfg.data_mode {
    let c = state.data_count.entry((xy, file_id)).or_insert(0);
    if *c < n { *c += 1; state.indices.push((xy.0, xy.1, file_id, read_num)); }
} // Off: no push, no spill, no render pass
```
Emit `"data_reads_per_pixel": N` and `"data_complete": false` in `heatmap.json`, plus a `fastsuite.report.v1` metric, whenever the mode is not `Exact`.

**effect** `--no-data` deletes the whole path: −(buffer term), −8 B/read of temp disk (80 MB per 10 M reads **[est.]**), no render pass. Floor becomes the measured 16.0 MB cell minus its own 1 MB buffer ≈ 15 MB **[est.]**. Mechanism: the counting grid is exact and complete without any per-read record; only pixel *membership* needs one.

**exactness** `Capped`/`Off` are approximations of `heatmap.data` alone — counts, PNG, grid JSON, corner metrics and the report stay bit-exact. Labelled per the suite rule; default stays `Exact`.

**measure** `sha256sum heatmap.{json,png,data}` vs the old binary on 100 k for `Exact` (identical); for `Off`, identical `counts`/`counts_r1`/`counts_r2`/`total_reads`/PNG hash plus the `data_complete:false` field; peak RSS via `wait4(rusage)` at 5 M.
**risk** notebook consumers must read the new field; default unchanged, so low. **effort:** M

## FM-2 — Put the preview flush inside `--buffer-mb`
**class:** memory
**evidence** `src/worker.rs:58`, with `src/cli.rs:31-35` stating previews "do not count against it":
```rust
const PREVIEW_FLUSH_BYTES: usize = 4 * 1024 * 1024;
```
**change** charge it (subtract `preview_bytes` from the worker's share) or cut the constant to 256 KiB — output-neutral, because previews are capped per `(pixel, kind)` (`worker.rs:234`) and re-sorted at encode (`format.rs:98-103`).
**effect** ≤8 MB → ≤0.5 MB at 2 workers **[est.]** (−7.5 MB off the 36.1 MB cell). Mechanism: the 4 MB threshold is *per worker* and invisible to the budget, so it is a hard floor `--buffer-mb 1` cannot move — which is why the 16.0 MB cell also needs `--preview-reads 0`.
**measure** 2 M-read cell peak RSS; `sha256sum heatmap.preview` unchanged; a send-counter under `--timings` to show the 16× message increase does not move wall.
**risk** more channel traffic on the single spiller — verify wall. **effort:** S

## FM-3 — Kill the 65 536-slot writer table and cap the row buffers
**class:** memory / allocation
**evidence** `src/spiller.rs:57` — allocated before a single row is seen:
```rust
let mut row_writers: Vec<Option<BufWriter<File>>> =
    (0..=u16::MAX as usize).map(|_| None).collect();
```
**change** a `HashMap<u16, writer>` (or a `Vec<u32>` slab index, 4 B/slot) with `BufWriter::with_capacity(4 * 1024, …)`, `OpenOptions::append` if a row must be reopened:
```rust
let w = row_writers.entry(y).or_insert_with(|| BufWriter::with_capacity(
    4 * 1024, File::create(row_temp_path(&temp_dir, y)).unwrap()));
```
**effect** −2.10 MB (table **[measured]**) − 1.25 MB (304 writers × 8 KiB → 4 KiB; 8 KiB is the `BufWriter` default and is oversized for 8 B records **[est.]**) ≈ −3.4 MB, 9 % of the 36.1 MB cell. Mechanism: the table is sized for `u16::MAX` rows but a short-read run has `max_seq_len + 1` = 152, and 8 KiB is 1 024 records of 8 B.
**measure** peak RSS at 5 M; hash the temp dir before the render phase and hash `heatmap.data`; assert `row_writers.len() <= size`.
**risk** none to output; 2× more flushes per row — which is what FM-9 pays for. **effort:** S

## FM-4 — Stream the row temp file instead of `fs::read` + a full decoded copy
**class:** memory
**evidence** `src/format.rs:225` then `230-236`: the row is read whole, then a second `Vec` of 12 B/read is built beside it.
```rust
let raw = std::fs::read(&temp)?;
let mut v = Vec::with_capacity(raw.len() / 8);
for chunk in raw.chunks_exact(8) { ... v.push((x, fid, rn)); }
```
**change** decode through a fixed buffer so the raw bytes die as they are consumed (reuse the `u16_at`/`u32_at` helpers already in this file, `format.rs:150-159`):
```rust
let mut f = BufReader::with_capacity(256 * 1024, File::open(&temp)?);
let mut rec = [0u8; 8];
while f.read_exact(&mut rec).is_ok() {
    v.push((le16(&rec[0..2]), le16(&rec[2..4]), le32(&rec[4..8])));
}
```
**effect** the render peak becomes `12 × max_row_reads` instead of `20 × max_row_reads` **[est.]**, a 40 % cut of that term. Mechanism: `fs::read` materialises the row twice and the raw buffer stays alive across the whole `chunks_exact` loop.
**measure** log `raw.len()` per row behind `--timings` *first* — the size of this term is currently unobservable, so that log is the deliverable. Then peak RSS at 5 M; `heatmap.data` hash identical.
**risk** low; keep an explicit parity/length check so a truncated row still fails loudly. **effort:** S

## FM-5 — Charge the budget by allocated capacity, and stop discarding the buffer
**class:** allocation
**evidence** `src/worker.rs:133,228` — `mem::take` leaves a zero-capacity `Vec` that regrows every chunk, and the charge is a constant 16 B against 12 B of storage:
```rust
indices: std::mem::take(&mut state.indices),
...
state.bytes_estimate += 16;
```
**change** account the real thing and recycle the allocation back from the spiller:
```rust
let used = state.indices.capacity() * std::mem::size_of_val(&state.indices[0]);
if used >= cfg.threshold_bytes { flush_local(state, tx); }
if let Ok(mut back) = recycle_rx.try_recv() { back.clear(); state.indices = back; }
```
**effect** two mechanisms: (1) the ceiling becomes true at any `--buffer-mb` — real residency is `capacity × 12` with `capacity ∈ [len, 2·len]` (probe: 5 M elements → capacity 8 388 608) while the budget believes `16 × len`, a disagreement up to 1.5×; (2) per-chunk regrowth costs ~2× the chunk in `memcpy` per flush — at the shipped default that is ~38 chunks × ~6.2 MB ≈ 235 MB of `memcpy` per 10 M reads **[est.]**.
**measure** peak RSS at `--buffer-mb 64` (where the mismatch bites) must be `64 + base`; `--timings` flush delta; hashes identical.
**risk** the recycle channel adds a sync point; `reserve_exact` before each chunk is the zero-risk fallback. **effort:** M

## FM-6 — Two-pass mode: regenerate `heatmap.data` exactly, hold nothing
**class:** algorithmic
**evidence** the record is a pure function of the sequence and the read ordinal — `src/worker.rs:181`, with `read_num` the file-local counter:
```rust
let (cx, cy, n_acgt) = encode_read(&seq);
let mq = mean_quality(&qual);
```
**change** `--two-pass`: pass 1 counts only (no index, no spill); pass 2 re-reads the same files in the same order and streams `(x, read_num)` into the per-row encoders. Row order stays sorted by `x` then `(fid, rn)` (`format.rs:29-37`), so the bytes are unchanged.
```rust
let (x, y, _) = encode_read(&seq);
rows[y].entry(x).or_default().push((fid, read_num));
```
**effect** index RAM → 0, temp disk → 0, `heatmap.data` exact, at the cost of a second decode: wall ≈ 2× the scan **[est.]** (only anchor: 10.61 s pinned vs 5.49 s at 2 threads). This must never become the default — an opt-in for memory-capped nodes.
**exactness** none lost. This is how to keep "every read is addressable" while dropping FM-1's stream.
**measure** `sha256sum heatmap.data` equal to the single-pass `Exact` output on all four corpus sizes; `mode: "two-pass"` in the wrapper.
**risk** input must be re-readable (already guaranteed: file-only, D-4) and unchanged — hash the first MiB in pass 1 and verify it in pass 2. **effort:** M/L

## FM-7 — Pack the record to 6 bytes behind a guard
**class:** algorithmic (encoding/layout), bit-packing
**evidence** `src/worker.rs:78` — 12 B in memory for 48 bits of payload; `y` is already implicit in the per-row file name (`spiller.rs:35`):
```rust
pub indices: Vec<(u16, u16, u16, u32)>,
```
**change** `Vec<[u8; 6]>` (`x:u8, fid:u8, rn:u32 LE`), align 1 so it is really 6 B/read, falling back to the 12 B struct when `max_seq_len > 254 || files.len() > 254` (both known before the first read). Do **not** use `(u8,u8,u32)` — alignment pads it back to 8.
```rust
let mut b = [0u8; 6]; b[0] = x; b[1] = fid; b[2..].copy_from_slice(&rn.to_le_bytes());
```
**effect** −50 % of the index term (6.3 → 3.2 MB at the shipped default **[est.]**) and −25 % of the temp stream (8 → 6 B/read, 80 → 60 MB per 10 M reads **[est.]**). Mechanism: `y` redundant, `x ≤ 150`, `fid < 4` for short-read runs.
**measure** peak RSS + temp-dir bytes at 5 M; `heatmap.data` hash identical.
**risk** the guard must be checked before the first push or a late long read truncates `x` — `debug_assert!(x <= 254)`. **effort:** M

## FM-8 — Shard the spiller and bound the channel
**class:** io-threads (parallel merge)
**evidence** `src/main.rs:179` is unbounded and `src/spiller.rs:74-88` is one thread doing three `write_all` calls per read, so nothing limits how far workers run ahead:
```rust
let (spill_tx, spill_rx) = mpsc::channel::<SpillerMsg>();
w.write_all(&x.to_le_bytes()).unwrap();
w.write_all(&fid.to_le_bytes()).unwrap();
w.write_all(&rn.to_le_bytes()).unwrap();
```
**change** `S` shards keyed `y % S`, each owning its writer map and temp-dir suffix; merge = concatenate the per-row files in `y` order (safe: the encoder sorts by `x` then `(fid, rn)`, `format.rs:29-37`, so arrival order is irrelevant). Back-pressure instead of growth:
```rust
let (tx, rx) = mpsc::sync_channel::<SpillerMsg>(64);   // ~64 MB ceiling at 1 MB chunks
```
**effect** removes the one RSS term the budget cannot see, and a serialization point doing 30 M `write_all` calls per 10 M reads. Wall effect **[est.]** — it only helps once the scan is cheaper than the spill (FM-10).
**measure** max channel depth (counter under `--timings`) ≤ cap; peak RSS at `-j 8`; `heatmap.data` hash identical at `-j 1/2/8` and `S=1/4`.
**risk** more open handles per row; ordering bugs are invisible in the grid but visible in `heatmap.data`, so the hash is the gate. **effort:** M

## FM-9 — One packed `write_all` per record
**class:** io
**evidence** `src/spiller.rs:85-87` (quoted in FM-8) — three calls for a record that fits one store.
**change**
```rust
let mut rec = [0u8; 8];
rec[0..2].copy_from_slice(&x.to_le_bytes());
rec[2..4].copy_from_slice(&fid.to_le_bytes());
rec[4..8].copy_from_slice(&rn.to_le_bytes());
w.write_all(&rec).unwrap();
```
**effect** 3× fewer calls and bounds checks in the spiller's hot loop (30 M → 10 M per 10 M reads); wall **[est.]**, visible only when the spiller is the bottleneck.
**measure** `--timings` spiller-busy at `-j 1` (where the scan cannot hide it); byte-identical temp files. **risk** none. **effort:** S

## FM-10 — Reusable byte buffers in the FASTQ loop
**class:** allocation (wall)
**evidence** `src/worker.rs:171,178-180` — `lines()` allocates four `String`s per read, and `encode_read` takes `&str` (`encode.rs:12`), so every base is UTF-8-validated as a side effect:
```rust
let mut lines = reader.lines();
while let (Some(Ok(header)), Some(Ok(seq)), Some(Ok(plus)), Some(Ok(qual))) =
    (lines.next(), lines.next(), lines.next(), lines.next())
```
**change** `read_until(b'\n', &mut buf)` into four reused `Vec<u8>`; `encode_read(&[u8])`, `mean_quality(&[u8])`; build a `String` only for a preview read actually kept:
```rust
if reader.read_until(b'\n', &mut hdr)? == 0 { break }
reader.read_until(b'\n', &mut seq)?; /* then plus, qual */
let (x, y, n) = encode_read(&seq);
```
**effect** 40 M allocations removed per 10 M reads (4/read **[est.]**) plus two UTF-8 scans of 151 B per read. Wall **[est.]** 20–30 % of the scan; the only hot-loop change that also helps the pinned cell (17.91 s).
**measure** pinned `-j 1` cell on 2 M reads, 5 runs, median; `heatmap.*` hashes identical.
**risk** CRLF input keeps `\r` in `seq`, changing `seq.len()` — strip it and add a CRLF test. **effort:** M

## FM-11 — Dense per-worker counters, and drop the duplicate map
**class:** algorithmic (accumulator layout)
**evidence** `src/worker.rs:199-201,224` — 3 hash lookups per read, one for a value that is merely a per-file counter; `src/main.rs:255-259` then copies the whole map again:
```rust
let e = state.counts.entry(xy).or_insert((0, 0.0));
e.0 += 1; e.1 += mq;
*state.per_file_counts.entry(file_id).or_insert(0) += 1;
```
**change** a `Vec<(u64, f64)>` of `size*size` per worker sized from the first read's length and grown on demand; per-file totals as `Vec<u64>` indexed by `file_id`; build `counts_arr` straight from the worker's density on `Done` and delete `counts_u32`:
```rust
let idx = cy as usize * size + cx as usize;
state.counts[idx].0 += 1;
state.counts[idx].1 += mq;
```
**effect** 3 hashes/read → 0; 0.37 MB/worker dense at 151 bp against ~0.9 MB of `HashMap` **[est.]**; and the `counts_u32` duplicate (≈0.9 MB at 151 bp, ≈40 MB at 1 kb reads **[est.]**) goes. This is also the answer to the encoding question: a cell is a *count* (`u64`), not a 4-bit composition code — the composition is the `(x, y)` pair, needing `2 × ceil(log2(L+1))` = 16 bits per read at L = 151. A 6-bit (64-bucket) hash accumulator would **collide**: `(T+C, T+G)` has 152² = 23 104 distinct values. Packing cells to `u32` is safe only while `total_reads < 2^32`; gate it on that check.
**measure** wall on the pinned 2 M cell; the three count arrays in `heatmap.json` byte-identical; peak RSS.
**risk** the dense array must grow when a longer read appears mid-file, and at 1 kb reads it is 16 MB/worker — it needs FM-5's admission control. **effort:** M

## FM-12 — Make the quality merge order deterministic
**class:** correctness-enabling
**evidence** `src/worker.rs:201` accumulates `f64`; `src/spiller.rs:123-127` merges those partials in *thread-completion* order, because `Done` messages race:
```rust
for (k, (cnt, q)) in c {
    let e = counts.entry(k).or_insert((0, 0.0));
    e.0 += cnt; e.1 += q;
}
```
**change** buffer `Done` messages and merge in `file_id` order (or accumulate `(sum_q, sum_len)` fixed-point instead of `f64`), and sort the HTML series by `(x, y)`:
```rust
dones.sort_unstable_by_key(|d| d.first_file_id);
```
**effect** `heatmap_3d.html`'s `z = cq/c` becomes reproducible. Integer counts are already order-exact — exactly why the PNG and grid JSON *are* byte-identical across `-j`; the `f64` sum is the one quantity that is not. No RSS effect.
**measure** ten runs at `-j 2` and `-j 8`, hash `heatmap_3d.html`: it differs today (doc blames `HashMap` order, `docs/BENCHMARKS.md:238`); after both fixes it must be stable.
**risk** none to other outputs. **effort:** S

---

## What is already good (keep it through every change above)

- **The wire format.** 1.85 B/read measured (18 488 308 B for 10 M reads) from varint + delta-coded read numbers within `(x, file)` groups (`format.rs:50-58`), 4.3× smaller than the 8 B/read stream. The varint is the *win*, not the bottleneck: ~4 pushes per read in a pass that is already I/O bound, and it is what makes `heatmap.data` 18 MB instead of 80 MB.
- **Per-row temp files, read-one-row-then-delete** (`format.rs:221-255`): render peak is bounded by the largest row, not by the run — the pattern the rest of the suite should copy.
- **Preview streaming** (`worker.rs:145-156`, `spiller.rs:89-103`): removed a ~76 MB term **[measured]** (`docs/BENCHMARKS.md:232-239`), outputs byte-identical.
- **Per-`(pixel, kind)` preview quota** (`worker.rs:231-236`): R2 cannot be starved by file order.
- **`log_norm` NaN guard** with a regression test (`output.rs:13-19`) — the one division by `ln(1)` that could have blacked out the PNG.
- **No serde**: 3 dependencies, hand-written JSON with a real parse test (`report.rs:243`); **`no_call` handling** (`report.rs:51-56`, `tests/nocall_poly_a_test.rs`) so all-N reads are not reported as polyA contamination.
- **Integer path determinism**: `-j 1` and `-j 8` give identical `.data`/`.json`/`.png` because counts are integers and the encoder sorts (`format.rs:29-37`).

## Dead ends / explicitly rejected

- **Bit-packing the grid (4-bit or 6-bit cells).** The grid is 0.55 MB at 151 bp — 0.4 % of the 143.2 MB cell — and a 64-bucket accumulator would *collide* on 23 104 distinct compositions, so it is not a lossless replacement at any size. At L = 1 000 the three `u64` arrays are 24 MB and FM-11's dense-`u32` option halves that; nothing here justifies a lossy encoding.
- **Compress the temp stream (flate2/zstd) to cut disk.** Disk is not the RSS term (records live in the worker `Vec`, not on disk) and it inserts a codec between worker and spiller. The stream is 80 MB per 10 M reads **[est.]**; FM-7 takes 25 % of that for free.
- **Bloom filter / count-min sketch for pixel membership.** A sketch can estimate *how many* reads are in a pixel — the exact counter already answers that — and cannot answer *which* reads, the only question `heatmap.data` exists for. Any sampling here must be a labelled cap (FM-1), never a sketch presented as exact.
- **Raise `--buffer-mb` or pass `-j 32` to go faster.** Measured dead: 5.3–5.9 s at both 512 MB and 8 MB (`docs/BENCHMARKS.md:242`), and `-j` is clamped to the file count (`cli.rs:120`), so a paired-end run gets 2 workers whatever you pass. Thread scaling needs block sharding, not a bigger number.
- **`mmap` the row temp files instead of `fs::read`.** Resident mapped pages count toward VmHWM exactly like a heap buffer, so peak RSS does not move; it only avoids the copy FM-4 already removes.

## How to prove nothing broke (shared protocol)

1. Baseline the current binary on the four bench corpora (100 k / 500 k / 1 M / 5 M pairs) and store `sha256sum heatmap.data heatmap.preview heatmap.json heatmap.png` per cell. After each proposition, re-run and diff those four hashes.
2. **`heatmap_3d.html` is excluded until FM-12 lands** — it serialises a `HashMap` and already varies run to run (`BENCHMARKS.md:238`).
3. Every RSS claim is `wait4(rusage)` peak of the child, never `ps` sampling; wall is the median of 5 runs, with the pinned `-j 1` cell reported next to `all cores` because that pair separates scan cost from spill cost.
4. Add `--timings` (scan / spill-busy / render / png to stderr) **before** touching anything: the render phase (`format.rs:221-267`) is currently invisible, so no one can say whether FM-8–FM-10 would even show up in wall time.
