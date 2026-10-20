# fastsim — ranked optimization dossier

Tool: `tools/fastsim` 0.1.0 (ex-novaseq-sim), the simulator that closes the QC -> sim loop
(`fastcheck --profile` fits, `fastsim --profile` consumes). Read-only review: no source
modified, nothing built or committed. Every file under `src/` was read (3681 lines). Probes used
`target/release/fastsim`, mtime 2026-10-20 09:22 — newer than the newest source (`src/main.rs`,
09:09) — so the numbers describe this tree without a rebuild. Every figure is `[measured]`
(artifact named) or `[estimated]` (mechanism named): artifacts are
`docs/optim/raw/fastsim_probe.txt` (raw readings + exact commands) and `docs/PLAN-ML.md:195`.
These are single-host probes on a contended 15 W part (Ryzen 7 7730U, load avg 1.9), not the
suite benchmark, which is 1/2/4 threads pinned with `taskset`.

## Cost model — what everything below is ranked against

| cell, 200 000 pairs, mandrill_E358, `taskset -c 0`, best of 2 | wall | pairs/s | peak RSS |
|---|---|---|---|
| default (gzip, level 6) | 14.34 s `[measured]` | 13 950 | 139 MB |
| `--no-gzip` | 3.49 s `[measured]` | 57 285 | 139 MB |
| gzip, `--dup-pool 4096` | 14.36 s `[measured]` | 13 932 | 9 MB |

Deflate is **75.7 % of wall** (10.85 s of 14.34 s = 54.2 us/pair `[measured]`); generation is
17.5 us/pair `[measured]`; the duplicate pool is 130 MB of RSS and 0 s of time `[measured]`.
`docs/PLAN-ML.md:195` records 19.9 s for the same unit on the workstation (10 050 pairs/s
`[measured]`) — same shape, slower clock. Ranking rule: **a 2x faster generator is worth 1.22x
until the compression is fixed; the compression change is worth 3x on its own** — so per-base
work is ranked *under* the two compression items.

## FS-1 — Deflate level 6 is three quarters of the runtime and buys 17 % of file size

**class** io-threads · **effort** S · **stream** fastq/GT identical (deflate is lossless), `.gz` bytes differ.
**evidence** `src/fastq.rs:34-40` — the level is a constant:
```rust
let gz = GzEncoder::new(Box::new(f) as Box<dyn Write>, Compression::new(6));
```
`[measured]` sweep over 14.12 MB of this tool's own output (artifact C): level 1 = 108.5 MB/s,
4.51 MB; level 4 = 59.8 MB/s, 4.09 MB; level 6 = 9.3 MB/s, 3.84 MB; level 9 = 2.9 MB/s, 3.74 MB.

**change** parameterise it, default 1, plus `--gz-level` (document 4 as the balanced choice):
```rust
pub fn open_gz(path: &Path, level: u32) -> Result<Self> {
    let gz = GzEncoder::new(Box::new(f) as Box<dyn Write>, Compression::new(level));
```

**effect** `[estimated]` 3.1-3.3x on the default cell: deflate 10.85 s -> ~0.93 s at the
measured level-1 rate, generation unchanged, 14.34 s -> ~4.4 s. **Mechanism**: level 6 runs
the hash-chain search with lazy matching, level 1 the greedy path with a 4-byte lookahead;
the cost is per-byte search and the saving is match quality — and FASTQ quality strings are
near-incompressible noise, so the search mostly fails. The proxy is the C `gzip` CLI;
in-process miniz_oxide is *faster* than it at level 6 (13 vs 9.3 MB/s `[measured]`,
artifact E), so hold the estimate at ±25 % and confirm against artifact A's second row.

**measure** new cell `fastsim_perf` under `rule fastsim_sweep`
(`bench/snakefile/Snakefile:230`): 200k pairs, `taskset -c 0`, best of 3, recording wall,
`VmHWM`, `wc -c` of both `.gz`; the `--json` `elapsed_s` (`main.rs:311`) is the same
quantity. **risk** +17.4 % `.gz` size `[measured]` — ~3.4 MB per 100k-pair unit (70.7 MB plain ->
19.2 -> 22.5 MB gz `[estimated]` from that ratio); nothing downstream parses size.

## FS-2 — Parallel emission: per-pair seeding, block decomposition, per-block gzip members

**class** algorithmic (parallel emission + determinism) · **effort** L · **stream** identical across `-j`; **changed** against 0.1.0.
**evidence** one loop, one stream: `main.rs:130` `ChaCha8Rng::seed_from_u64(cfg.seed)`,
`main.rs:142` `for pair_idx in 0..cfg.pairs`; `-j` is parsed then discarded —
`src/cli.rs:251-257` "simulation is single-threaded; --threads N accepted ... and ignored".
Per-pair RNG consumption is **not constant**, which kills any offset-based split: a
duplicate draws one coin and skips generation (`main.rs:163-176`, then `continue`), the N
overlay draws only when the rate fires (`quality.rs:157`), the boundary rewrite draws only
when a biology base is G (`artefact.rs:196-213`).

**change**
```rust
let plan = plan_pass(&profile, &cfg, &mut ChaCha8Rng::seed_from_u64(cfg.seed));
// content pass, per block: pair-local stream, dup regenerates its source
let mut rng = ChaCha8Rng::seed_from_u64(mix(cfg.seed, pair_idx as u64));
```
1. *Plan pass, serial, cheap*: per pair draw only `is_dup`/`is_drop` and, for a duplicate, the
   pool draw. The pool's identity mapping depends on the classification only, never on bytes, so
   the whole cross-pair dependency fits in 2-3 f64 per pair (<1 % of the work):
   `enum Plan { Fresh, Drop, Dup { src: u32 } }`.
2. *Content pass, parallel over blocks*: pair `i` gets its own stream,
   `ChaCha8Rng::seed_from_u64(mix(cfg.seed, i as u64))`, and `Dup { src }` **regenerates** pair
   `src` from its seed — a pure function of `(seed, src)` — instead of reaching into another
   worker's buffers (1.4 % extra pairs at the measured `dup_rate 0.0141`): the last cross-thread
   dependency gone.
3. *Emission*: `BLOCK = 8192` pairs, a function of `(pairs, BLOCK)` only. Each worker owns one
   block and writes one gzip member per output file at a fixed level, concatenated in block
   order (RFC 1952 multi-member streams are legal; `gzip -d`, `zcat`, `samtools` concatenate
   them). **`-j1` uses the same decomposition**, so threading never moves a block boundary.

**determinism** `-j1`/`-j2`/`-j4` byte-identical in all three artefacts *by construction*: pair
bytes are a pure function of `(seed, pair_idx)`, the plan pass is serial and thread-count
independent, block boundaries ignore the thread count, members concatenate in block-index order.
Verification: `sha256sum` of `_R1`, `_R2`, `.gt.tsv` at `-j1/-j2/-j4` under `taskset -c 0-3` plus
`cmp` after `zcat`. **This is a stream change against 0.1.0** — today's draws are interleaved from
one stream, so per-pair streams cannot reproduce them: keep `--legacy-stream` (current path, `-j1`
only) and a golden-hash test; do not claim the old bytes survive.

**effect** `[estimated]` after FS-1: generation 17.5 + deflate 6.7 us/pair (measured 108.5 MB/s
over 141 MB) = 24.2 us/pair serial; at 4 workers with deflate included, 17.5/4 + 6.7/4 = 6.1
us/pair -> 1.2 s per 200k, i.e. **3.2-3.6x on top of the post-FS-1 cell**. **If FS-1 is skipped
this item is worth 1.22x** — Amdahl from the measured split (54.2 of 71.7 us/pair is serial
deflate). **What limits it** is the deflate, and only after FS-1 the generator: measured
4-process deflate scaling here is 1.8x at level 6 and 2.0x at level 1 (artifact D, a contended
15 W part — re-measure on a pinned server before quoting).

**measure** `fastsim_perf` at `-j1/-j2/-j4`, `taskset -c 0-3`, best of 3: wall, pairs/s,
core-seconds (`getrusage(RUSAGE_CHILDREN)`; no `/usr/bin/time` on this host), `VmHWM`, plus
the three sha256s as the same cell's correctness assertion. **risk** documented stream
change; multi-member gzip less portable for a few consumers; `--legacy-stream` doubles the
emit path for one release.

## FS-3 — The duplicate pool caches bytes it never reads, and those bytes block parallelism

**class** algorithmic (duplicate realisation carried across copies) / memory · **effort** M · **stream** identical.
**evidence** the entry carries `insert` (`src/artefact.rs:50-57`) but the replay path never
reads it — `main.rs:163-176` uses `src.artefact`, `src.r1`, `src.r2`, `src.insert_len` — while
the push copies it every pair:
```rust
// main.rs:255-261
hist.push(PooledPair {
    insert: insert.clone(),
```
`[measured]` RSS 139 MB at 200k / `dup_pool 131072` vs 9 MB at `dup_pool 4096`, wall unchanged
(14.34 vs 14.36 s): the pool is 130 MB of the 139 MB and none of the time.

**change** (1) delete the dead `PooledPair::insert`; (2) replace the byte cache with FS-2's
plan — store `(source_pair_idx)` and regenerate on replay, so `History` becomes a `Vec<u32>`
of fresh pair indices plus a ring cursor and the `r1`/`r2` copies disappear.
**effect** `[measured]` RSS 139 -> ~9 MB (the `dup_pool 4096` cell is the proxy for a
byte-free pool), and it removes the cross-thread dependency that makes FS-2 hard. **Time
unchanged `[measured]` — do not sell it as a speedup.**
**measure** `/proc/<pid>/VmHWM`; correctness by `pooled_pair_replays_identical_bytes`
(`artefact.rs:569`) plus the closed-loop duplicate rate (the 1.0-of-requested property in
`artefact.rs:38-48`). **risk** regenerating `src` costs a full pair of work (1.4 % of pairs)
and reproduces the original exactly only with FS-2's seeding — land them together.

## FS-4 — Per-base quality calls `exp`/`erf`/`powf`/`sqrt` for numbers the profile already fixes

**class** algorithmic (profile-derived tables — fastsim's analogue of "fit once") · **effort** M · **stream** **byte-identical** (same floats, same draw order): the one item here with no stream risk.
**evidence** `src/quality.rs:117-123`:
```rust
let sd = p.qual_noise;
if sd > 0.0 {
    let z = sample_gauss(0.0, sd, -20.0, 20.0, rng);
    v += cap_loss(q, sd) + z;
```
`cap_loss` costs one `exp` plus one `erf` (`quality.rs:38`, `:52`), `apply_errors` one `powf`
(`:193`), `Normal::new` one `sqrt` — reconstructed per base, not per run (`randutil.rs:44`) — and
`read_meanq_mean` is recomputed per mate per read (`profile.rs:594`).

**change** build the tables once at profile load:
```rust
struct QTables { p_err: [f64; 42], cap_loss: [f64; 42],              // sd constant per profile
                 pre_noise: HashMap<(u8, u8, u32, u32), Box<[u8]>> } // (mate, off, h, t)
```
`pre_noise` is the integer quality vector *before* the noise draw, and it is small: the offset is
quantised by the measured per-read histogram — mandrill_E358 has **6 `read_meanq` bins per mate**
`[measured]`, so `off` takes <=6 values — and `(h, t)` come from 2 head bins and 9-11 tail bins.

**effect** `[estimated]` removes 2 `exp`, 1 `sqrt`, 1 `powf` per base while keeping the one real
draw: 40-60 % of `qualify`, ~15-25 % of the plain 17.5 us/pair, ~2-3 % of today's wall.
**Mechanism**: every term except `z` is deterministic given `(profile, mate, offset, h, t)`.
**measure** `fastsim_perf --no-gzip` plus the assertions that move if the model moved:
`substitution_rate_matches_q30` (`quality.rs:501`), per-position mean Q, byte bisect against the
parent binary. **risk** the key must cover everything the vector depends on — a missing `t`
silently clips tails; guard with a cached-vs-recomputed assertion over 1000 pairs.

## FS-5 — One uniform draw per base to almost never make an error

**class** algorithmic · **effort** S · **stream** not identical (draw counts change).
**evidence** `src/quality.rs:182-198`:
```rust
let p = 10f64.powf(-(q as f64) / 10.0);
let u: f64 = rng.gen();
if u < p && seq[i] != b'N' {
```
At the measured mean quality (~Q32) `p` is ~6e-4, so 99.94 % of these draws and all of the
`powf` calls decide *not* to act: ~300 draws + 300 `powf` per pair `[counted]`.

**change** geometric gaps — exactly equivalent for a Bernoulli process, since `P(no error
over k bases) = (1-p)^k`; `p_err[q]` and `ln_q[q] = ln(1-p_err[q])` precomputed, Q0 keeps
today's always-substitute branch:
```rust
let mut i = 0;
loop { i += (rng.gen::<f64>().ln() / ln_q).floor() as usize;   // 1 draw per ERROR
       if i >= seq.len() { break } out[i] = substitute(out[i], rng); i += 1; }
```
**effect** `[estimated]` 300 draws + 300 `powf` per pair -> ~1 draw: most of `apply_errors`,
~10-20 % of the plain 17.5 us/pair, ~2 % of today's wall. **Mechanism**: draws collapse from bases
to errors, and the error *count* distribution is unchanged by construction.
**measure** `substitution_rate_matches_q30` plus a `substitutions` counter in the `--json` metrics
(`main.rs:312-328`) — the count, not just the rate. **risk** a sign slip turns 1e-4 into 1.0; the
existing 1e6-base test catches it and `ln(1-p)` stays finite at Q41.

## FS-6 — 300 one-byte `write_all` calls per pair through a `dyn Write`

**class** allocation / io · **effort** S · **stream** identical.
**evidence** `src/fastq.rs:80-84`:
```rust
for &q in &mr.qual {
    self.out.write_all(&[q + PHIRED_OFFSET])?;
}
```
`out` is `BufWriter<Box<dyn Write>>` (`fastq.rs:24`), so each byte costs a vtable hop, a capacity
check and a 1-byte memcpy: 300 per pair, 60 M per 200k-pair unit `[counted]`.

**change** encode into a reused buffer, emit once:
```rust
const ASCII: [u8; 42] = /* q + 33 */;
self.qbuf.clear(); self.qbuf.extend(mr.qual.iter().map(|&q| ASCII[q as usize]));
self.out.write_all(&self.qbuf)?;  self.out.write_all(b"\n")?;  // 1 call, not 150
```
**effect** `[estimated]` 300 dispatches -> 1 write per read: 2-5 % of the plain path, 3-6 %
of the gzip path (fragmented input also costs capacity checks and deflate restarts at buffer
edges). **Mechanism**: per-call fixed cost, not per-byte bandwidth.
**measure** `fastsim_perf` gz and plain cells; acceptance is `sha256sum` equality with the
parent binary at the same seed — this change *must* be byte-identical. **risk** none beyond
the table bounds (`q <= 41` is enforced at `quality.rs:123`).

## FS-7 — `Box<dyn Write>` in the writer type erases every inlining opportunity

**class** allocation / codegen · **effort** S · **stream** identical.
**evidence** `src/fastq.rs:24` and `:37-40`:
```rust
pub struct FastqWriter { out: BufWriter<Box<dyn Write>>, .. }
let gz = GzEncoder::new(Box::new(f) as Box<dyn Write>, Compression::new(6));
```
**change** monomorphise or use inlined enum dispatch: `pub enum Sink { Plain(BufWriter<File>),
Gz(BufWriter<GzEncoder<File>>) }` with `#[inline] fn put(&mut self, b: &[u8])`.
**effect** `[estimated]` 1-3 % alone, but it is what lets FS-6's single `write_all` inline into the
copy loop — land them together. **measure** `fastsim_perf` plain cell. **risk** low, crate-internal.

## FS-8 — Insert generation: one f64 draw and a four-way branch per base

**class** algorithmic (insert sampling) · **effort** S · **stream** not identical.
**evidence** `src/biology.rs:81-91`:
```rust
let u: f64 = rng.gen();
out.push(if u < gc / 2.0 { b'C' } else if u < gc { b'G' }
         else if u < gc + (1.0 - gc) / 2.0 { b'A' } else { b'T' });
```
Insert *lengths* are already one CDF search (146 bins in mandrill_E358 `[measured]`); the cost is
per-base content, ~300 bases/pair `[counted]`.

**change** one `u64` per 8 bases through a 256-entry table built from `gc` at load
(byte-uniform -> exactly 64 slots per base, so the marginal is exact, not approximate); for
biased positions pre-convert each `[f32;4]` row to `[u32;4]` cumulative and compare one `u32`
draw instead of `gen::<f32>()`: `for (c, w) in out.chunks_mut(8).zip(rng.gen::<u64>().to_le_bytes())`.
**effect** `[estimated]` 8x fewer draws on the unbiased share (1/base -> 1/8 bases) and no `f64`
compare chain the compiler cannot vectorise: ~5-10 % of the plain path, ~1 % of today's wall. **measure** `fastsim_perf --no-gzip`, `random_dna_gc_matches_request`
(`biology.rs:245`), closed-loop `g_frac_max_diff`. **risk** a LUT built from a default rather
than `profile.gc` shifts GC by a constant and fails every base-composition gate.

## FS-9 — ~10 heap allocations per pair for buffers whose shape never changes

**class** allocation · **effort** M · **stream** identical.
**evidence** per pair, `[counted]`: 2 read seqs (`artefact.rs:161`), 2 qual vecs
(`quality.rs:65`), 2 error seqs (`quality.rs:183`, `seq.to_vec()`), the insert
(`biology.rs:55`), the insert clone (`main.rs:256`), 2 header `String`s (`fastq.rs:68`).
**change** one `Scratch` owned by the loop, passed as `&mut`:
`struct Scratch { insert: Vec<u8>, seq: [Vec<u8>;2], qual: [Vec<u8>;2], line: String }` with
`qualify_into(&p, mate, &read, &mut qual, &mut seq, rng)` returning `()`.
**effect** `[estimated]` 2 M allocations per 200k-pair unit -> 0; 1-3 % of wall (malloc/free
~20-40 ns against 17.5 us/pair) — hygiene, worth doing after FS-4/FS-5 remove the real work.
**measure** `fastsim_perf --no-gzip`, `perf stat -e page-faults` if available.
**risk** borrow-checker churn in `main.rs`; the duplicate path hands `MateRead` to the pool, so
FS-3 must land first or the pool copies.

## FS-10 — Two formatting sites put `fmt` dispatch in the per-pair path

**class** formatting / allocation · **effort** S · **stream** identical.
**evidence** `src/fastq.rs:64-72` returns a fresh `String` per record (the caller borrows the
temporary at `main.rs:243`, `w1.write_record(&w1.header(pair_idx, 1), &m1)`), and
`src/truth.rs:103-112` formats a GT row through 12 `fmt::Display` arguments into a `BufWriter`
(`truth.rs:41`). `[measured]` the GT is 6.1 MB per 200 000 pairs, so this is instruction count,
not bytes.
**change** append integers into two reused 96-byte buffers — header (6 of its 7 fields constant per
run) and GT row — one `write_all` each, no `fmt`.
**effect** `[estimated]` 4 allocations + 14 `fmt` dispatches per pair -> 0; <2.5 % of wall, and it
removes the last two `format!`/`writeln!` sites from the hot loop, which is what makes FS-9 (an
allocation-free loop) verifiable. **measure** `header_format_matches_novaseq_style`
(`fastq.rs:114`), `row_count_matches_fastq_pairs` (`truth.rs:178`), byte equality of all three
artefacts. **risk** none.

## FS-11 — No counter exists to prove a parallel or table-driven run is equivalent

**class** correctness-enabling (instrument) · **effort** S · **stream** identical.
**evidence** the counters at `main.rs:136-141` and the `--json` block (`main.rs:312-328`) carry
rates, not the primitives these items move: no RNG-draw, substitution, mean-Q or gzip-byte count.
**change** add `rng_draws` (instrumented RNG wrapper behind `--stats`, `#[inline]` no-op
otherwise), `substitutions`, `mean_q`, `gz_bytes`, and `tests/golden_stream.rs` pinning `sha256` of
the three artefacts for `--seed 42 --pairs 1000` per committed profile.
**effect** `[estimated]` 0 % on the default path; it is what makes FS-2's `-j1 == -j4` and FS-5's
"same error count" checkable rather than asserted. **measure** itself — the golden test is the cell.
**risk** the wrapper must be unreachable from the default path or it taxes every draw.

---

## What is already good (do not regress)

- **The insert ceiling is geometry, not a model gap**: lengths above `2L-12` collapse into one
  saturated bin (`profile.rs:56-64`, `cli.rs:58-61`); no proposal here touches it.
- **Duplicate replay is byte-exact and measured.** `PooledPair` stores finished mates and
  `main.rs:163-176` re-emits them without calling `qualify`; `artefact.rs:38-48` records the
  measured failure of the older design (0.35 of requested, independent of pool size: 1024 ->
  0.322, 1e6 -> 0.346, which ruled the ring buffer out). The brief's "a duplicate re-draws
  substitution errors per copy" is the *pre-fix* design, not this HEAD — if a scorer still sees
  <1.0 of requested, look at the drop-class exclusion (`main.rs:254`) or the fingerprint's anchor
  requirement, not at `qualify`.
- **RSS is flat in `--pairs`**: saturates with the pool (139 MB at 200k vs 9 MB at `dup_pool 4096`
  `[measured]`); nothing else grows with run length.
- **Buffering is layered in the right order**: `BufWriter` above the gzip encoder
  (`fastq.rs:38-40`), `BufWriter` on the GT (`truth.rs:41`), and R1/R2 are two independent
  streams — 2-way deflate overlap is free the moment FS-2 exists.
- **No per-base indirection in the profile**: `pos_q()` returns a slice (`profile.rs:564`),
  histograms are flattened at load, and the only `HashMap` (base planes, `profile.rs:649`) is
  read once at load — no per-base hashing to remove.
- **The hot loop has no `String` formatting, no logging, no `Instant::now()`**; the single
  `t0` is at `main.rs:48`.
- **Determinism is already real**: one seed, one stream, no time/thread/address input — FS-2 starts
  from a reproducible baseline.

## Dead ends / explicitly rejected

1. **Swap flate2's backend (zlib / zlib-ng) for speed.** `Cargo.lock:149-155` pins `miniz_oxide`,
   whose measured in-process rate (~13 MB/s at level 6, artifact A) is *faster* than system C zlib
   (9.3 MB/s, artifact C): at level 6 the cost is the algorithm, not the implementation, and after
   FS-1 both are noise.
2. **Split the RNG by pair-index arithmetic** (one stream, per-block offsets). Impossible: a
   duplicate draws one coin and skips generation (`main.rs:163-176`), the N overlay is
   conditional (`quality.rs:157`), the boundary rewrite is conditional
   (`artefact.rs:196-213`). Any offset table is wrong after the first duplicate.
3. **Optimise the pool's sampling path** (`gen_range`, ring index, the *time* cost of
   `insert.clone()`). Measured: `--dup-rate 0` vs default 1.47 vs 1.48 s at 20k, `dup_pool` 131072
   vs 4096 14.34 vs 14.36 s at 200k — the pool is a memory item (FS-3), never a time one.
4. **Parallelise only the gzip writer, keep the generator serial.** At level 6 that is 1.0x —
   deflate *is* the serial part. At level 1 it caps at 2.96x (3.49 s generation + 1.34 s deflate
   `[estimated]`), after which the generator is the ceiling: the same wall FS-2 hits from the other
   side, without the determinism work done.
5. **Memoise whole reads by `(insert_len, h, t)` to skip generation.** The run-length marginals
   are the point of the model; a memo makes reads identical by construction, which is exactly what
   the duplicate fingerprint measures — it would inflate the duplicate rate and break the loop.

## Verification protocol (every item)

1. Cell `fastsim_perf`: 200 000 pairs, mandrill_E358, `taskset -c 0..3`, one warm-up then best of
   3; gz and `--no-gzip` sub-cells; `-j1/-j2/-j4` once FS-2 lands. Record wall, pairs/s,
   core-seconds, `VmHWM`, `wc -c` of both `.gz`, and `--json` `elapsed_s` (the same quantity).
2. Gates that must not move: closed-loop QC parity, duplicate rate vs requested,
   `g_frac_max_diff`, substitution rate at Q30, `sha256` of the three artefacts across `-j`.
3. Anything derived from the level-sweep proxy is `[estimated]` with its mechanism, never mixed in.
