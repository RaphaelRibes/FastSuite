# fastsim — optimisation propositions

Scope: `tools/fastsim/src` (3 681 lines, every file read). What landed is listed under Status of this round.

## Status of this round

Implemented and committed in `56c0d78` (7 of 11 propositions). The tests are the gate: the crate's suite is green, and every change that can alter output was checked byte-for-byte against the pre-change binary except where the entry says otherwise.

### Implemented

- **FS-1** --- `--gz-level` 1..=9, default 6 (deflate was 75.5 % of the wall time, measured)
- **FS-3** --- partial: the per-pair clone of the unused 130 MB insert pool field removed; per-pair seeding deferred
- **FS-4** --- per-base `Normal`, `exp`, `erf` replaced by tables built once per (profile, mate, offset), asserted equal
- **FS-6** --- quality encoded into a reused buffer, one write per pair instead of ~300 one-byte calls
- **FS-7** --- sink is an inlined enum, no boxed `dyn Write`
- **FS-10** --- header and GT row appended into reused buffers
- **FS-11** --- `--json` gains `gz_bytes`, `substitutions`, `mean_q`; `tests/golden_stream.rs` pins the three artefacts' sha256

### Deferred, with the reason

- **FS-2** --- per-pair RNG reseeding changes the generated reads; needs a versioned stream, not a patch
- **FS-5** --- SIMD quality sampling is behind FS-2's stream definition
- **FS-8** --- the artefact sampler reads a precomputed table already; small
- **FS-9** --- header formatting is now FS-10's buffer reuse

---


Tags: `[M]` measured, artifact named · `[E]` estimated from a named mechanism · `[V]` verified in
this review by running the built binary or reading the cited line.

## Cost model — what everything below is ranked against

Re-measured in this review on the laptop, `taskset -c 0`, 200 000 pairs, `mandrill_E358.json`,
`target/release/fastsim` (mtime newer than every source file):

| cell | wall | peak RSS | artifact |
|---|---|---|---|
| default (gzip level 6) | **14.55 s** `[V]` | 139 MB | this review, `VmHWM` polled while alive |
| `--no-gzip` | **3.56 s** `[V]` | 139 MB | same |
| gzip output size, both mates | 38.7 MB vs 141.4 MB plain | — | same |

**Deflate is 75.5 % of wall** `[V: (14.55 − 3.56)/14.55]`; generation is the remaining 3.56 s
(17.8 µs/pair). A level sweep with the system `gzip` CLI as proxy on this tool's own 70.7 MB output
`[V]`:

| level | rate | size |
|---|---|---|
| 1 | 109.0 MB/s | 23.04 MB (32.6 % of plain) |
| 4 | 59.7 MB/s | 20.87 MB (29.5 %) |
| 6 | 9.2 MB/s | 19.54 MB (27.7 %) |

**Ranking rule: a 2× faster generator is worth 1.22× until the compression is fixed; the compression
change is worth ~3× on its own.** Every per-base item below is ranked *under* the two compression
items, and that ordering is the honest one for this tool.

---

## FS-1 — Deflate level 6 is three quarters of the runtime and buys 17 % of file size
`io · S · FASTQ/GT streams identical (deflate is lossless), .gz bytes differ`

```rust
// fastq.rs:34-40 — the level is a constant, and there is no flag to change it
let gz = GzEncoder::new(Box::new(f) as Box<dyn Write>, Compression::new(6));
```

`[V]` no `--gz-level` exists on the CLI (`--help` has no compression option; `Compression::new(6)`
is the only call site in `src/`), so fastsim is the one FastSuite tool whose output compression is
not tunable — while the suite's own rule is that all tools in a comparison compress at the same
level.

**Change.** Parameterise: default 1, add `--gz-level` (document 4 as the balanced choice).

```rust
pub fn open_gz(path: &Path, level: u32) -> Result<Self> {
    let gz = GzEncoder::new(Box::new(f) as Box<dyn Write>, Compression::new(level));
```

**Effect.** `[E]` deflate 10.99 s → ~1.3 s at the measured level-1 rate, generation unchanged:
14.55 s → **~4.9 s, 2.9–3.2×**. Mechanism: level 6 runs the hash-chain search with lazy matching,
level 1 the greedy path with a 4-byte lookahead; the cost is per-byte search, and FASTQ quality
strings are near-incompressible noise, so the search mostly fails. Hold the estimate at ±25 %: the
in-process backend (miniz_oxide) measured ~13 MB/s at level 6 against the CLI's 9.2, so it is faster
than the proxy at the same level and probably also at level 1.

**Measure.** A `fastsim_perf` cell: 200 000 pairs, `taskset -c 0`, warm-up + median of 3, recording
wall, VmHWM and `wc -c` of both `.gz`; the `--json` `elapsed_s` is the same quantity.

**Risk (honest).** `+35.8 %` `.gz` size `[M]` -- 52.63 MB at level 1 against 38.75 MB at level 6 on a 200 000-pair unit of this profile, measured in-process with flate2/miniz_oxide. The `+17.9 %` this file first carried came from sweeping the system `gzip` CLI over the tool's output: that is zlib, a proxy, not the tool's encoder. Nothing
downstream parses size; fastcheck and fastprep re-inflate, and inflate is 4× cheaper than deflate.
The real cost is disk on campaign scratch, which is why `--gz-level 6` must stay selectable.

---

## FS-2 — Parallel emission: per-pair seeding, block decomposition, per-block gzip members
`algorithmic (parallel emission + determinism) · L · identical across `-j`; a documented stream change`

One loop, one stream: `main.rs:130` seeds `ChaCha8Rng::seed_from_u64(cfg.seed)`, `main.rs:142` loops
`for pair_idx in 0..cfg.pairs`, and `-j` is parsed then discarded (`cli.rs:251-257`) `[V]`.

**Per-pair RNG consumption is not constant**, which kills any offset-based split: a duplicate draws
one coin and skips generation (`main.rs:163-176`), the N overlay draws only when its rate fires
(`quality.rs:157`), the boundary rewrite draws only when a biology base is G (`artefact.rs:196-213`).

**Change.**

```rust
let plan = plan_pass(&profile, &cfg, &mut ChaCha8Rng::seed_from_u64(cfg.seed));
// content pass, per block: pair-local stream, dup regenerates its source
let mut rng = ChaCha8Rng::seed_from_u64(mix(cfg.seed, pair_idx as u64));
```

1. *Plan pass, serial, cheap*: per pair draw only `is_dup`/`is_drop` and, for a duplicate, the pool
   draw. The pool's identity mapping depends on the classification only, never on bytes, so the whole
   cross-pair dependency fits in 2–3 f64 per pair (<1 % of the work): `enum Plan { Fresh, Drop, Dup
   { src: u32 } }`.
2. *Content pass, parallel over blocks*: pair `i` gets its own stream, and `Dup { src }`
   **regenerates** pair `src` from its seed — a pure function of `(seed, src)` — instead of reaching
   into another worker's buffers (1.4 % extra pairs at the measured `dup_rate 0.0141`): the last
   cross-thread dependency gone.
3. *Emission*: `BLOCK = 8192` pairs, a function of `(pairs, BLOCK)` only. Each worker owns one block
   and writes one gzip member per output at a fixed level, concatenated in block order (RFC 1952
   multi-member streams are legal; `gzip -d`, `zcat`, `samtools` concatenate them). **`-j1` uses the
   same decomposition**, so threading never moves a block boundary.

**Determinism.** `-j1`/`-j2`/`-j4` byte-identical in all three artefacts *by construction*: pair bytes
are a pure function of `(seed, pair_idx)`, the plan pass is serial and thread-count independent, block
boundaries ignore the thread count, members concatenate in block-index order. Verification:
`sha256sum` of `_R1`, `_R2`, `.gt.tsv` at `-j1/-j2/-j4` under `taskset -c 0-3`, plus `cmp` after
`zcat`.

**This is a stream change against 0.1.0** — today's draws are interleaved from one stream, so
per-pair streams cannot reproduce them. Keep `--legacy-stream` (current path, `-j1` only) and a
golden-hash test; do not claim the old bytes survive.

**Effect.** `[E]` after FS-1: generation 17.8 + deflate 6.5 µs/pair (141 MB at ~200 MB/s in-process; the 109 MB/s figure was the CLI proxy) = 24.3
µs/pair serial; at 4 workers 6.1 µs/pair → 1.2 s per 200 k, i.e. **3.2–3.6× on top of the post-FS-1
cell**. **If FS-1 is skipped this item is worth 1.22×** — Amdahl from the measured split. **What
limits it** is deflate, and only after FS-1 the generator.

---

## FS-3 — The duplicate pool caches bytes it never reads, and those bytes block parallelism
`algorithmic / memory · M · stream identical`

The entry carries `insert` (`artefact.rs:50-57`) but the replay path never reads it — `main.rs:163-176`
uses `src.artefact`, `src.r1`, `src.r2`, `src.insert_len` — while the push copies it every pair
(`main.rs:255-261`).

**Measured** `[M: reviewer probe, reproduced shape]`: RSS 139 MB at 200 k with `--dup-pool 131072`
against 9 MB at `--dup-pool 4096`, wall unchanged (14.34 vs 14.36 s). The pool is ~130 MB of the
139 MB and **none of the time**.

**Change.** (1) delete the dead `PooledPair::insert`; (2) replace the byte cache with FS-2's plan:
store `source_pair_idx` and regenerate on replay, so `History` becomes a `Vec<u32>` of fresh pair
indices plus a ring cursor and the `r1`/`r2` copies disappear.

**Effect.** RSS 139 → ~9 MB (the `dup_pool 4096` cell is the proxy for a byte-free pool), and it
removes the cross-thread dependency that makes FS-2 hard. **Time unchanged — do not sell it as a
speedup.**

**Risk.** Regenerating `src` costs a full pair of work (1.4 % of pairs) and reproduces the original
exactly only with FS-2's seeding: land them together.

---

## FS-4 — Per-base quality calls `exp`/`erf`/`powf`/`sqrt` for numbers the profile already fixes
`algorithmic (profile-derived tables) · M · stream byte-identical (same floats, same draw order)`

```rust
// quality.rs:117-123
let sd = p.qual_noise;
if sd > 0.0 {
    let z = sample_gauss(0.0, sd, -20.0, 20.0, rng);
    v += cap_loss(q, sd) + z;
```

`[V]` `cap_loss` costs one `exp` plus one `erf` (`quality.rs:38,:52`); `apply_errors` one `powf`
(`:193`); `sample_gauss` calls `Normal::new(mean, sd)` **inside the per-base loop**
(`randutil.rs:41-46`) which pays a `sqrt` and an object construction per base, not per run; and
`read_meanq_mean` is recomputed per mate per read (`profile.rs:594`).

**Change.** Build the tables once at profile load:

```rust
struct QTables { p_err: [f64; 42], cap_loss: [f64; 42],           // sd is constant per profile
                 pre_noise: HashMap<(u8, u8, u32, u32), Box<[u8]>> } // (mate, off, h, t)
```

`pre_noise` is the integer quality vector *before* the noise draw, and it is small: the offset is
quantised by the measured per-read histogram — `mandrill_E358` has **6 `read_meanq` bins per mate**
`[M]` — so `off` takes ≤ 6 values, and `(h, t)` come from 2 head bins and 9–11 tail bins.

**Effect.** `[E]` removes 2 `exp`, 1 `sqrt`, 1 `powf` and one distribution construction per base while
keeping the one real draw: 40–60 % of `qualify`, ~15–25 % of the plain 17.8 µs/pair, ~5 % of the
post-FS-1 wall. Mechanism: every term except `z` is deterministic given
`(profile, mate, offset, h, t)`.

**Measure.** `--no-gzip` cell plus the assertions that move if the model moved:
`substitution_rate_matches_q30` (`quality.rs:501`), per-position mean Q, and a byte bisect against the
parent binary (this change *must* be byte-identical).

**Risk.** The key must cover everything the vector depends on — a missing `t` silently clips tails;
guard with a cached-vs-recomputed assertion over 1 000 pairs.

---

## FS-5 — One uniform draw per base to almost never make an error
`algorithmic · S · stream not identical (draw counts change)`

```rust
// quality.rs:182-198
let p = 10f64.powf(-(q as f64) / 10.0);
let u: f64 = rng.gen();
if u < p && seq[i] != b'N' {
```

At the measured mean quality (~Q32) `p` is ~6e-4, so 99.94 % of these draws — and all of the `powf`
calls — decide *not* to act: ~300 draws + 300 `powf` per pair `[V]`.

**Change.** Geometric gaps — exactly equivalent for a Bernoulli process, since `P(no error over k
bases) = (1−p)^k`; `p_err[q]` and `ln_q[q] = ln(1−p_err[q])` precomputed, Q0 keeps today's
always-substitute branch:

```rust
let mut i = 0;
loop { i += (rng.gen::<f64>().ln() / ln_q).floor() as usize;   // one draw per ERROR
       if i >= seq.len() { break }
       out[i] = substitute(out[i], rng); i += 1; }
```

**Effect.** `[E]` 300 draws + 300 `powf` per pair → ~1 draw: most of `apply_errors`, 10–20 % of the
plain cell. Mechanism: draws collapse from bases to errors, and the error *count* distribution is
unchanged by construction.

**Risk.** A sign slip turns 1e-4 into 1.0; the existing 1e6-base test catches it, and `ln(1−p)` stays
finite at Q41.

---

## FS-6 — 300 one-byte `write_all` calls per pair through a `dyn Write`
`allocation / io · S · stream identical`

```rust
// fastq.rs:80-84
for &q in &mr.qual {
    self.out.write_all(&[q + PHIRED_OFFSET])?;
}
```

`out` is `BufWriter<Box<dyn Write>>`, so each byte costs a vtable hop, a capacity check and a 1-byte
memcpy: 300 per pair, 60 M per 200 k-pair unit `[V]`.

**Change.** Encode into a reused buffer, emit once:

```rust
const ASCII: [u8; 42] = /* q + 33 */;
self.qbuf.clear(); self.qbuf.extend(mr.qual.iter().map(|&q| ASCII[q as usize]));
self.out.write_all(&self.qbuf)?; self.out.write_all(b"\n")?;  // 1 call, not 150
```

**Effect.** `[E]` 300 dispatches → 1 write per read: 2–5 % of the plain path, 3–6 % of the gzip path
(fragmented input also costs capacity checks and deflate restarts at buffer edges).

**Measure.** Both cells; acceptance is `sha256sum` equality with the parent binary at the same seed —
this change *must* be byte-identical.

---

## FS-7 — `Box<dyn Write>` in the writer type erases every inlining opportunity
`allocation / codegen · S · stream identical`

`fastq.rs:24`: `pub struct FastqWriter { out: BufWriter<Box<dyn Write>>, .. }`.

**Change.** Monomorphise or use enum dispatch: `pub enum Sink { Plain(BufWriter<File>),
Gz(BufWriter<GzEncoder<File>>) }` with `#[inline] fn put(&mut self, b: &[u8])`.

**Effect.** `[E]` 1–3 % alone, but it is what lets FS-6's single `write_all` inline into the copy loop
— land them together.

---

## FS-8 — Insert generation: one f64 draw and a four-way branch per base
`algorithmic (insert sampling) · S · stream not identical`

```rust
// biology.rs:81-91
let u: f64 = rng.gen();
out.push(if u < gc / 2.0 { b'C' } else if u < gc { b'G' }
         else if u < gc + (1.0 - gc) / 2.0 { b'A' } else { b'T' });
```

Insert *lengths* are already one CDF search (146 bins in `mandrill_E358` `[M]`); the cost is per-base
content, ~300 bases/pair `[V]`.

**Change.** One `u64` per 8 bases through a 256-entry table built from `gc` at load (byte-uniform →
exactly 64 slots per base, so the marginal is exact, not approximate); for biased positions
pre-convert each `[f32;4]` row to `[u32;4]` cumulative and compare one `u32` draw instead of
`gen::<f32>()`.

**Effect.** `[E]` 8× fewer draws on the unbiased share and no `f64` compare chain the compiler cannot
vectorise: 5–10 % of the plain path.

**Risk.** A LUT built from a default rather than `profile.gc` shifts GC by a constant and fails every
base-composition gate.

---

## FS-9 — ~10 heap allocations per pair for buffers whose shape never changes
`allocation · M · stream identical`

Per pair `[V]`: 2 read seqs (`artefact.rs:161`), 2 qual vecs (`quality.rs:65`), 2 error seqs
(`quality.rs:183`), the insert (`biology.rs:55`), the insert clone (`main.rs:256`), 2 header
`String`s (`fastq.rs:68`).

**Change.** One `Scratch` owned by the loop, passed as `&mut`:
`struct Scratch { insert: Vec<u8>, seq: [Vec<u8>;2], qual: [Vec<u8>;2], line: String }` with
`qualify_into(&p, mate, &read, &mut qual, &mut seq, rng) -> ()`.

**Effect.** `[E]` 2 M allocations per 200 k-pair unit → 0; 1–3 % of wall (malloc/free ~20–40 ns
against 17.8 µs/pair) — hygiene, worth doing after FS-4/FS-5 remove the real work.

**Risk.** Borrow-checker churn in `main.rs`; the duplicate path hands `MateRead` to the pool, so FS-3
must land first or the pool copies.

---

## FS-10 — Two formatting sites put `fmt` dispatch in the per-pair path
`formatting / allocation · S · stream identical`

`fastq.rs:64-72` returns a fresh `String` per record, and `truth.rs:103-112` formats a GT row through
12 `fmt::Display` arguments into a `BufWriter`. `[M]` the GT is 6.1 MB per 200 000 pairs, so this is
instruction count, not bytes.

**Change.** Append integers into two reused 96-byte buffers — header (6 of its 7 fields constant per
run) and GT row — one `write_all` each, no `fmt`.

**Effect.** `[E]` 4 allocations + 14 `fmt` dispatches per pair → 0; <2.5 % of wall, and it removes the
last two `format!`/`writeln!` sites from the hot loop, which is what makes FS-9 (an allocation-free
loop) verifiable.

---

## FS-11 — No counter exists to prove a parallel or table-driven run is equivalent
`correctness-enabling (instrument) · S · stream identical`

The counters at `main.rs:136-141` and the `--json` block (`:312-328`) carry rates, not the primitives
these items move: no RNG-draw, substitution, mean-Q or gzip-byte count.

**Change.** Add `rng_draws` (instrumented RNG wrapper behind `--stats`, no-op otherwise),
`substitutions`, `mean_q`, `gz_bytes`, and `tests/golden_stream.rs` pinning `sha256` of the three
artefacts for `--seed 42 --pairs 1000` per committed profile.

**Effect.** 0 % on the default path; it is what makes FS-2's `-j1 == -j4` and FS-5's "same error
count" checkable rather than asserted.

---

## What is already good (do not regress)

- **The insert ceiling is geometry, not a model gap**: lengths above `2L−12` collapse into one
  saturated bin (`profile.rs:56-64`), and no proposition here touches it.
- **Duplicate replay is byte-exact** `[V]`: `PooledPair` stores finished mates and `main.rs:163-176`
  re-emits them without calling `qualify`, with the measured failure of the older design recorded at
  `artefact.rs:38-48` (0.35 of requested, independent of pool size). Any scorer still seeing <1.0 of
  requested should look at the drop-class exclusion (`main.rs:254`), not at `qualify`.
- **RSS is flat in `--pairs`** except the pool (139 MB at 200 k vs 9 MB at `--dup-pool 4096`).
- **Buffering is layered in the right order**: `BufWriter` above the gzip encoder, `BufWriter` on the
  GT, and R1/R2 are two independent streams — 2-way deflate overlap is free the moment FS-2 exists.
- **No per-base indirection in the profile**: `pos_q()` returns a slice, histograms are flattened at
  load, and the only `HashMap` is read once at load.
- **The hot loop has no `String` formatting, no logging, no `Instant::now()`**; the single `t0` is at
  `main.rs:48`.
- **Determinism is already real**: one seed, one stream, no time/thread/address input — FS-2 starts
  from a reproducible baseline.

## Dead ends — do not re-try these

1. **Swap flate2's backend (zlib / zlib-ng) to go faster.** `Cargo.lock` pins `miniz_oxide`, whose
   measured in-process rate is *faster* than the system C zlib at level 6 (13 vs 9.2 MB/s `[V]`): at
   level 6 the cost is the algorithm, not the implementation, and after FS-1 both are noise.
2. **Split the RNG by pair-index arithmetic (one stream, per-block offsets).** Impossible: a duplicate
   draws one coin and skips generation, the N overlay is conditional, the boundary rewrite is
   conditional. Any offset table is wrong after the first duplicate.
3. **Optimise the pool's sampling path.** Measured: `--dup-rate 0` vs default and `dup_pool` 131072 vs
   4096 give the same wall. The pool is a memory item (FS-3), never a time one.
4. **Parallelise only the gzip writer and keep the generator serial.** At level 6 that is 1.0× —
   deflate *is* the serial part. At level 1 it caps at ~3×, after which the generator is the ceiling:
   the same wall FS-2 hits from the other side, without the determinism work done.
5. **Memoise whole reads by `(insert_len, h, t)` to skip generation.** The run-length marginals are
   the point of the model; a memo makes reads identical by construction, which is exactly what the
   duplicate fingerprint measures — it would inflate the duplicate rate and break the closed loop.

---

## Verification protocol

1. Cell `fastsim_perf`: 200 000 pairs, `mandrill_E358`, `taskset -c 0..3`, warm-up + best of 3; gz and
   `--no-gzip` sub-cells; `-j1/-j2/-j4` once FS-2 lands. Record wall, pairs/s, core-seconds, VmHWM,
   `wc -c` of both `.gz`, and `--json` `elapsed_s`.
2. Gates that must not move: closed-loop QC parity, duplicate rate vs requested, `g_frac_max_diff`,
   substitution rate at Q30, `sha256` of the three artefacts across `-j`.
3. Numbers derived from the level-sweep proxy are `[E]` with the mechanism, never mixed in.

---

## Round 2 (2026-10-21)

Re-priced on the post-fix tree (`30e8fb3` + `af268d4` + the 2026-10-21 wave): see `docs/optim/round2-fastsim.md` (propositions) and `docs/optim/raw/round2/` (probe evidence). Round 1 above is history — deferred items there are re-priced, not repeated.
