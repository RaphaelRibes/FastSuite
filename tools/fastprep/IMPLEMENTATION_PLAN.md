# fastprep implementation plan

Target: a FASTX preprocessor that deduplicates, adapter-trims, polyX-trims and
quality-filters in a single pass, and is **not slower than fastp** on any
tested workload while holding a flat memory profile outside the
deduplication stage.

Priority order throughout: **speed first, memory second, features third.**
Where the two conflict, the tie-breaker is that memory must stay bounded and
predictable, because an OOM on an HPC node is a total failure while a 10%
slowdown is a regression.

---

## Ground rules

1. **No optimisation lands without a before/after measurement.** The current
   numbers in the paper's optimisation table are the only measured ones, and
   they predate everything below.
2. **Output-identical by default.** Every speed change is either provably
   semantics-preserving or lives behind a flag. The one deliberate exception
   is deduplication, which is new behaviour and ships on by default
   (`--no-dedup` turns it off).
3. **One golden-output test guards all of it.** Freeze the output of the
   current build on a fixed 100k-read input; every phase must reproduce it
   byte for byte unless it changes semantics on purpose.
4. **Surgical changes.** Each phase touches the modules named and nothing
   else. If a phase grows past its listed files, stop and re-scope.

---

## Phase 0: Measurement harness (gate: nothing else starts without this)

Everything below is an estimate derived from reading the code. Phase 0 turns
estimates into numbers, and may reorder the phases.

**Tasks**

- Add `--profile` emitting per-stage nanosecond totals: parse, dedup,
  adapter, scan, write. Accumulate per chunk, merge with `Stats`.
- Benchmark matrix: `{plain, gz, zst} × {SE, PE} × {1, 8, 96} threads` on a
  10M-read slice, driven by `hyperfine`, wired into the existing pixi task.
- Isolation runs, in this order:
  ```
  fastprep -1 in.fq  -o /dev/null --no-adapter --no-5 --no-3 --no-interior
  fastprep -1 in.fq  -o /dev/null --no-adapter
  fastprep -1 in.fq  -o /dev/null
  fastprep -1 in.fq.gz -o out.fq.gz
  ```
- `perf stat -e cycles,instructions,branch-misses,cache-misses,LLC-load-misses`
  on the scan-only configuration.

**Verify:** the four isolation runs' deltas sum to within 10% of total wall
clock. If they do not, the instrumentation is wrong, not the pipeline.

**Decision point:** if run 1 is already ~80% of run 4's time, Phases 3 to 5
are near-worthless and Phase 1 is the entire project.

---

## Phase 1: Pipeline architecture and allocation

**Expected: 3–5× wall clock. This is almost certainly the largest single win.**

*(Landed in part: `run_single` streams recycled `ChunkBuf` arenas from the
reading thread to a compute thread that runs rayon, with dedup and writing on
the ordered thread. The zstd default (1.2) and byte-sized chunks (1.4) were not
adopted: output compression still follows the extension and `--chunk` still
counts records.)*

The pre-phase code did decompress, parse and three heap allocations per read on
the main thread, then compress and write on the same thread. Only `trim` ran on
rayon. If trimming is 20% of the work, adding threads buys almost nothing,
which matches the observed "slower than fastp".

**1.1 Chunk arena (memory + speed)**

Replace the three `Vec<u8>` per read with one arena per chunk:

```rust
struct Chunk {
    bytes: Vec<u8>,                 // one allocation, reused across chunks
    recs: Vec<Rec>,                 // { id: Range<u32>, seq: Range<u32>, qual: Range<u32> }
}
```

At 10⁸ reads this removes ~3×10⁸ allocations. `mimalloc` currently reduces
the cost of each; this removes the count. Arenas are recycled through the
return channel rather than freed, so steady-state allocation is zero.

- Files: `pipeline.rs` (`ChunkBuf`; no separate `chunk.rs`).
- Verify: `valgrind --tool=dhat` or `heaptrack` shows allocation count flat in
  read count; golden output unchanged.

**1.2 Get compression off the critical thread**

gzip level 6 runs at roughly 30–50 MB/s on one core. If the tool writes
`.fq.gz` from one thread, that number *is* the ceiling and no other work in
this document matters.

- Default output codec becomes zstd, using
  `zstd::stream::write::Encoder::multithread(n)`.
- For gzip output, use `gzp` for block-parallel deflate.
- Lower the default gzip level from 6 to 2 and document the ratio cost
  (small on FASTQ).
- Files: `io.rs`, `cli.rs`.
- Verify: output decompresses identically; single-thread wall clock unchanged;
  8-thread wall clock scales.

**1.3 Decompression and parsing off the main thread**

Reader thread produces raw byte blocks split at record boundaries; workers
parse and process; writer owns compression. The main thread coordinates only.

- Files: `pipeline.rs`.
- Verify: PE record ordering preserved (see the ordering test in Phase 2.3);
  CPU utilisation approaches `min(threads, cores)` on uncompressed input.

**1.4 Chunk sizing**

`chunk = 1024` reads is too small for the per-chunk rayon dispatch and channel
round trip. Size chunks **by bytes** (2–4 MB) rather than read count, so
behaviour is stable across read lengths from 150 bp Illumina to 100 kb ONT.

- Verify: throughput curve vs chunk size is flat above 2 MB.

**Memory budget after Phase 1:** `channel_depth × chunk_bytes × 2` for input
and output, plus one arena per in-flight chunk. With depth 4 and 4 MB chunks
that is ~64 MB regardless of input size. **Assert this in a test.**

---

## Phase 2: Deduplication

New capability. Sequenced here because it must be built against Phase 1's
final thread topology, not because it is the next-largest speed win. Phases
3–5 are independent of it and can proceed in parallel.

The value is not the ~10% of trim work skipped on duplicates. It is that
users currently run a separate dedup tool, costing one full
decompress-parse-write-compress round trip. On compressed data that round trip
is more than 90% of the dedup tool's own wall clock. Fusing it removes that
pass entirely.

**2.1 Fingerprinting**

- Depend on `xxhash-rust` and `rustc-hash` directly. Do **not** depend on the
  `fastdedup` crate (binary crate, `hasher.rs` not public) and do not vendor
  `processor.rs`.
- Port the `SequenceHasher` trait and `HashVerifier` only, ~60 lines.
- Anchor window per the specification: `a = T5_hard`,
  `k = min(L_min, L - T5_hard - T3_hard)`, taken from the end that the
  configured scan does **not** touch.
- PE: hash both mates' anchor windows, combined as
  `h1 ^ rot(h2, N/2)`.

**2.2 Threading**

`HashVerifier::verify` takes `&mut self` and cannot go inside
`into_par_iter().map()`. A `Mutex` would mean one lock per read at 10⁸ reads.
Split it:

- **Parallel (rayon, alongside `trim`)**: compute the fingerprint. Pure
  function, embarrassingly parallel, cheap next to the polyX scan.
- **Sequential (the ordered thread from Phase 1.3)**: `verify(hash)`, then
  drop or emit.

No locks, ordering preserved. Cost estimate: at 3×10⁸ reads the table is
~4.6 GB, so essentially every lookup is a DRAM round trip at ~80 ns, giving
~24 s of added serial time against a ~450 s baseline. Acceptable, but it lands
on the serial thread, so re-measure after Phase 1.

**2.3 Semantics and plumbing**

- `DropReason::Duplicate`; PXMT status byte 7; **format version bump to 3**;
  `fastprep-meta` viewer updated to read v2 and v3.
- `Stats::reads_dropped_duplicate` plus a row in the `causes` array.
- Pair-level decision: a duplicate pair drops both mates, produces no
  singleton. Singletons are not deduplicated against a SE fingerprint space.

**Verify (correctness of the stage itself):**
- Round-trip: a file duplicated `n` times yields exactly the unique set.
- Anchor invariance: appending polyG tails of differing length to two copies of
  one read still yields one survivor.
- PE ordering: R1[i] pairs with R2[i] at `-j 1` and `-j 16`.
- PXMT: a v3 file round-trips through `fastprep-meta`, and v2 files still read.

**2.4 Memory projection (this is the only unbounded stage)**

The projection is an identity, not an estimate. `FxHashSet` is hashbrown, and
hashbrown's sizing is deterministic:

```rust
fn table_bytes<T>(capacity: usize) -> u64 {
    let buckets = ((capacity * 8) / 7).next_power_of_two();  // load factor 7/8
    (buckets * (size_of::<T>() + 1)) as u64                  // +1 control byte
}
```

Checked against the FastDedup benchmark table at 300M reads:

| | buckets | predicted | measured |
|---|---|---|---|
| 64-bit | 2²⁹ | 4608 MB | 4616 ± 1 MB |
| 128-bit | 2²⁹ | 8704 MB | 8712 ± 1 MB |

Within 8 MB, which is the rest of the pipeline. Because it is exact, refusing a
run on the strength of it is defensible.

| reads | unique | 64-bit | 128-bit |
|---|---|---|---|
| 10⁷ | 9×10⁶ | ~150 MB | ~290 MB |
| 10⁸ | 9×10⁷ | ~1.5 GB | ~2.9 GB |
| 3×10⁸ | 2.7×10⁸ | ~4.6 GB | ~8.7 GB |

Two consequences of `next_power_of_two`, both easy to get wrong. **Note the
boundaries are at `7 * 2^k / 8`, not at powers of two** (verified numerically):

| buckets | max capacity | 64-bit size |
|---|---|---|
| 2²⁸ | 234,881,024 | 2304 MB |
| 2²⁹ | 469,762,048 | 4608 MB |
| 2³⁰ | 939,524,096 | 9216 MB |

- **It is a step function.** Capacity 469,762,048 costs 4608 MB; 469,762,049
  costs 9216 MB. Report the projected **bucket count**, not only the byte
  total, so a user marginally over a boundary can act on it.
- **`MAX_PREALLOC_ENTRIES` was `500_000_000`, badly placed.** It sat 6% above
  the 2²⁹ boundary, so it rounded to 2³⁰ buckets and the clamp itself projected
  9216 MB where 4608 MB would do. It is now `469_762_048`
  (`memory.rs::capacity_ceiling(29)`). Rounding it to a
  power of two, which an earlier draft of this plan advised, is *worse*: 2²⁹
  as a capacity also rounds up to 2³⁰ buckets.

Free win: `estimate_records` sizes for **all** reads, but the table only holds
**unique** ones. Scale by `(1 - expected_dup_rate)`. That ~10% reduction will
often cross a bucket boundary downwards and halve the allocation. Take the rate
from a flag, or from the sampled prefix already read for adapter auto-detection.

- Verify: `table_bytes` matches measured RSS within 2% at 10⁷, 10⁸ and 3×10⁸.

**2.5 Admission control: detect the enforced limit, not the machine**

Total system RAM is the wrong number on a cluster. A SLURM job with `--mem=8G`
on a 1.5 TB node sees 1.5 TB from `sysinfo`, allocates happily, and is killed by
the cgroup. The FastDedup benchmarks ran under exactly that setup (33 GB
allotted, killed above 32 GB), so this must read the enforcement boundary.

Take the **minimum** over everything discoverable:

```rust
#[cfg(target_os = "linux")]
pub fn memory_limit() -> Option<u64> {
    [cgroup_limit(), slurm_limit(), rlimit_as(), mem_available()]
        .into_iter().flatten().min()
}
```

- **cgroup v2**: parse `/proc/self/cgroup` (format `0::/path`), read
  `/sys/fs/cgroup/<path>/memory.max`. **Walk the whole ancestor chain and take
  the minimum**, because the limit is frequently set on a parent slice rather
  than the leaf. Value is a byte count or the literal `max`.
- **cgroup v1**: `/sys/fs/cgroup/memory/memory.limit_in_bytes`, path from the
  `:memory:` line. The unlimited sentinel is `9223372036854771712`; treat
  anything above physical RAM as absent.
- **SLURM**: `SLURM_MEM_PER_NODE` (MB), or `SLURM_MEM_PER_CPU` ×
  `SLURM_CPUS_ON_NODE`. Read it even under cgroup enforcement, since it is what
  the scheduler kills against, and some sites set the env without cgroups.
- **`getrlimit(RLIMIT_AS)`**: respect if set, but do not treat as an RSS proxy.
  It bounds virtual address space, and mimalloc reserves far more VA than it
  commits, so a naive check refuses jobs that would have fit.
- **`/proc/meminfo`**: use `MemAvailable`, not `MemTotal`. MemAvailable accounts
  for reclaimable page cache and is the honest figure.

No new dependency: ~120 lines of file parsing plus one `libc::getrlimit`.
`sysinfo` pulls in a lot for one number and does not walk nested cgroup slices
correctly. Gate the module behind `cfg(target_os = "linux")` and return `None`
elsewhere, which downgrades the policy to warn-only.

> Housekeeping: the pixi manifest is `linux-64` only, but the paper's §8.1
> fragment claims `linux-64` and `osx-arm64`. Reconcile before this lands, since
> the platform list now has behavioural consequences.

**Policy**, with `reserve = max(512 MB, 10% of limit)`:

```
projected + reserve  <  0.6 × limit   proceed
                     <  0.9 × limit   proceed, warn with the numbers
                     >= 0.9 × limit   refuse before allocating
```

Before refusing, retry the projection at 64-bit if the width was 128. If that
fits, downgrade with a warning stating the resulting collision probability from
the birthday bound. A quantified degradation beats a failure.

The refusal message is the deliverable here, so make it prescriptive and name
the **source** of the limit:

```
error: deduplication needs ~4.9 GB (2^29 buckets x 9 B, 64-bit fingerprints)
       but this job is limited to 4.0 GB (cgroup v2: /slurm/uid_1234/job_98765)

  --dedup-bits 64          already selected
  --no-dedup               run without deduplication
  sbatch --mem=6G          rerun with enough memory
```

When someone is confused about why a 1.5 TB node refuses a 5 GB allocation,
"cgroup v2" is the answer they need.

**Verify:**
- Unit tests against synthetic `/proc` and `/sys` trees for cgroup v1, cgroup
  v2, nested v2 slices where the parent is the binding limit, and the v1
  unlimited sentinel.
- Integration: run under `systemd-run --user --scope -p MemoryMax=1G` and assert
  the tool refuses rather than being OOM-killed.
- Assert the SLURM env path is exercised with cgroups absent.

**2.6 Saturation: the runtime safety net**

The startup projection depends on the file-size estimate, and compression ratio
varies enough that it can be wrong. Poll actual usage every N chunks via cgroup
`memory.current` (what the killer actually looks at), falling back to
`/proc/self/statm` field 2.

On crossing the threshold, **stop inserting but keep querying**. The set
freezes; existing fingerprints still catch duplicates, newly seen reads are no
longer recorded. This is one-sided: it can only miss duplicates, never remove a
read it should have kept. Proved in §6.1.2 of the paper by induction on
`S'_t ⊆ S_t`. Same error structure as the cascade prefilter in §3.11.

The alternative, aborting cleanly, is also safe (the `TmpFile` guard removes the
`.part` files) but throws away the whole run. Ship saturation as the default and
`--dedup-on-oom=abort` as the opt-out.

**Saturation must be loud.** Silent partial deduplication is exactly what this
project criticises in probabilistic filters; the only difference here is that
the error is one-sided and declared, and that difference evaporates if it is not
declared:

- warning on stderr at the moment of freezing;
- `dedup_saturated: true` and `dedup_saturated_at_read: <index>` in the JSON
  stats;
- a line in the report block;
- the reported duplication rate flagged as a lower bound.

**Verify:**
- Force saturation with an artificially low `--dedup-max-mem` and assert the
  survivor set is a **superset** of the unsaturated run's survivors (one-sided
  property, tested directly rather than assumed).
- Assert `dedup_saturated` appears in JSON and the exit code stays 0.
- Assert RSS plateaus after saturation rather than continuing to climb.

**2.7 Ordering experiment (paper material)**

Deduplicate before trimming (anchor window), before trimming (full read), and
after trimming, on the same NovaSeq lane. Report the duplicate-count delta.
This quantifies whether the anchor window buys the sensitivity the
specification claims. Either answer is publishable.

---

## Phase 3: Adapter path

**Expected: large, possibly the second-biggest win.** On reading the code this
component plausibly costs more per read than the entire polyX scan.

**3.1 Revert to suffix-anchored alignments**

`match_adapters` searches every internal start position (~148 per adapter at
L=150), for an event the contamination model excludes: readthrough always
places the adapter as a suffix. Restrict to the ~30 suffix alignments the
specification already describes.

- Files: `adapter.rs`.
- Verify: adapter-trimmed base counts unchanged on a contaminated test set
  (the internal-position generalisation should never have been firing on real
  readthrough anyway); measure the delta if it does change.

**3.2 Seed prefilter**

Build one `memchr::memmem::Finder` per adapter over its first 10–12 bases,
stored in `AdapterSet`. One SIMD substring search locates the candidate
offset; verify a single suffix alignment instead of 30.

- Verify: identical output to 3.1; per-read adapter time drops by ~10×.

**3.3 PE overlap: vectorise the mismatch count**

`detect_pe_overlap` is O(L²) in the **common** case (clean long-insert pairs
admit no overlap, so all ~120 offsets are tried, each accumulating
`⌈Iε⌉+1 ≈ 16` mismatches at ~21 comparisons). Compare 32 bases at a time
against the reverse-complemented mate and accumulate popcounts, with the
saturation test hoisted out of the inner loop.

- Verify: `mm_PE(I)` identical for all I on randomised pairs; ~10× on the
  component.
- Do **not** ship the sampled probe gate as a default. It is not sound: an
  admissible offset may carry up to `⌈Iε⌉` mismatches, and nothing stops them
  landing in the sample. If shipped at all it is a flag, excluded from
  conformance testing, and reported in the stats block.

---

## Phase 4: PolyX composition prefilter (the cascade)

**Expected: 10–50× on the scan component, on clean libraries.**

Specified in §3.11 of the paper. Two propositions make it output-identical:
`mm(K) > ⌊Δ_m⌋ ⟹ Γ > Δ_m`, so counting gives an upper bound on `K*`; and the
rejection test with window `w = m_h + ⌊Δ_m⌋ - 1`.

**Tasks**

- `fn probe_end(bases, anchor, params) -> Option<usize>` returning `l_max`, or
  `None` when the end provably resolves to `T = 0`.
- Anchored case: byte-compare a 32-byte window against a broadcast anchor,
  extract the mask, popcount. Generic case: compare the window against itself
  shifted by one lane.
- Call it before `scan` / `streaming_walk`, using `l_max` to tighten `upper`.
- Round `w` up to 32; widening is sound (it can only raise `c_X`).

**Verify (this is the whole justification, so make it a property test):**
```rust
proptest! {
    #[test]
    fn cascade_matches_plain(read in arb_read(), p in arb_params()) {
        prop_assert_eq!(trim_cascade(&read, &p), trim_plain(&read, &p));
    }
}
```

**Also in this phase: unfuse the mean-quality pass.** `streaming_walk` drags
the scalar per-step machinery across every base of every read purely to obtain
a quality sum that vectorises to ~1 vector add per 32 bytes. The measured 1.4×
for fusion was against *two scalar passes*, which is the wrong baseline. Split
them, vectorise the sum, and add the early exit: if
`partial + 41·(L-i) < θ_q·L`, the read is condemned and the rest of the pass
can be skipped.

Unfusing also frees the scan order, allowing the 3′ end to be scanned first on
two-colour chemistry where the polyG tail is the expected artefact.

---

## Phase 5: Γ inner loop

**Expected: 3–5× on whatever scan work survives Phase 4.** Lower absolute
value than Phase 4, since it accelerates work rather than removing it. Do it
after, and only if Phase 0's numbers still justify it.

- **Branchless per-step gap.** `if matches` is data-dependent and close to
  unpredictable on non-homopolymeric sequence; a mispredict costs more than the
  arithmetic it guards. Both branches are affine in W, so:
  `ADD = [1.0, μ]`, `MUL = [ν, -μ]`, `gap = ADD[m] + MUL[m] * w`.
- **Fixed-point accumulator.** Γ is the sole loop-carried dependency, so the
  loop cannot retire faster than its accumulation latency: 4 cycles for `f64`,
  1 for `u16` fixed point. This is a hard 4× ceiling improvement.
- **Quality binning.** RTA3 (NovaSeq) emits four distinct Phred values, HiSeq X
  eight. A 256-byte binning LUT collapses the per-step table to
  4×4×4×4 = 256 entries, register-resident rather than L2-resident. Detect the
  number of distinct values during the sampling pass already used for adapter
  auto-detection; keep the unbinned path as fallback.
- **Bounds checks and rule dispatch.** Monomorphise on `MatchRule` the way
  direction already is; iterate with zipped slice iterators so the checks
  actually vanish rather than relying on LLVM to hoist them.

**Verify:** golden output unchanged; `perf stat` branch-miss rate on the scan
configuration drops below 1%.

---

## Phase 6: SIMD Γ (deferred)

Gather-based vectorisation of the scan body, as described in §3.6. Deferred
indefinitely: after Phase 4 the scan runs on a small fraction of reads, so the
absolute gain shrinks in proportion. Revisit only if a profile shows scan
compute still dominant, which would be surprising.

---

## Memory budget (target end state)

| Component | Bound | Notes |
|---|---|---|
| Chunk arenas | `depth × chunk_bytes` | ~64 MB at depth 4, 4 MB chunks |
| Compression buffers | `threads × 1 MiB` | zstd multithread internal |
| Gate / score LUTs | < 1 MB | shared read-only, per-CCD replication optional |
| `AdapterSet` + finders | < 100 kB | `Arc`, one per process |
| Fingerprint set | **Θ(unique reads)** | only unbounded stage; on by default (`--no-dedup` opts out), admission-checked against the enforced limit (2.5), saturating rather than OOM-ing (2.6) |

**Invariant to assert in CI:** with `--no-dedup`, peak RSS is independent of
input size. Run 1M and 100M reads, assert RSS ratio < 1.2.

---

## Non-goals

- **Do not port FastDedup's resume machinery.** Hash preloading, corrupted-file
  detection and `set_len` truncation contradict fastprep's crash-safety model,
  which writes to `.part` and deletes on any error path. The two designs answer
  the same question in opposite ways; gluing them together produces a tool that
  sometimes resumes and sometimes discards depending on which stage failed.
- No quality-based sliding-window trimming. Out of scope, and Trimmomatic and
  fastp already do it.
- No indel-tolerant adapter alignment as a default. Illumina indel rate is
  ~10⁻⁶; the banded DP exists for chemistries this tool does not target. Keep
  it available as a third tier behind the seed if sensitivity parity with
  cutadapt is ever needed for a benchmark.
- No multi-file input concatenation, no output splitting, no HTML reports.

---

## Separate task: two bugs in FastDedup

Not blocking any phase above, since only `hasher.rs` is being ported, but both
are live in the released tool.

1. **Desync check inside the dry-run loop.** In
   `execute_paired_deduplication`, the
   `if reader_r1.next().is_some() || reader_r2.next().is_some()` block sits
   before the closing brace of the `while let`, so it runs every iteration,
   consumes the next pair, and bails. `--dry-run` fails on any paired file with
   two or more pairs. The non-dry-run path has the same check correctly placed
   after the loop.
2. **Main write loop does not call `base_read_id`.** The dry-run path uses it
   (strips whitespace, tabs, and the `/1` / `/2` mate suffix); the write loop
   inlines `.split(|&b| b == b' ')` instead. On pre-CASAVA-1.8 headers or
   `fastq-dump --split-files` output, the IDs never compare equal and the tool
   bails on the first pair. `test_base_read_id_strips_mate_suffix` passes
   because the function is tested and then not called.

Also worth checking: when the readers run out unevenly, the `while let` tuple
consumes a record from the longer file before the pattern fails, so the
post-loop check is off by one and a single trailing extra record goes
undetected.

---

## Sequencing summary

```
Phase 0  measurement harness          gate, ~1 day
Phase 1  architecture + allocation    est. 3-5x        ← largest win
Phase 2  deduplication                feature          ← depends on Phase 1
Phase 3  adapter path                 est. 2-10x on that component
Phase 4  polyX cascade                est. 10-50x on the scan
Phase 5  Γ inner loop                 est. 3-5x on residual scan
Phase 6  SIMD Γ                       deferred
```

Phases 3, 4 and 5 are independent of each other and of Phase 2. Phase 1 gates
everything, and Phase 0 gates Phase 1.

**Definition of done:** fastprep is not slower than fastp on any cell of the
Phase 0 benchmark matrix, peak RSS with `--no-dedup` is independent of input
size, and the property test asserting cascade/plain equivalence passes on
10⁵ generated cases.
