# Plan — `fastcheck`: a real QC binary

Replace the two-wheeled approach (fastmaps for composition + a hypothetical
"add the FastQC modules later" phase 4) with one QC binary that answers every
short-read QC question **exactly**, in **one pass**, with a **printed memory
projection** and **no hidden thresholds**.

Non-goals: PASS/WARN/FAIL verdicts from magic constants; sampling presented as
measurement; "enrichment" scores whose formula lives only in the source; a
zipped HTML bundle.

---


## 0. Phase 1–2 status: implemented, measured (`tools/fastcheck`)

**Update 2026-10-12 — what shipped after this plan was written:** the rich
report (polyX run-length distributions per base, bounded repeated-sequence
and 12-mer k-mer counters, the 64×64 composition heatmap, and before/after
side-by-side via `--processed-1/--processed-2`), the `--profile` fitter (a
fastsim profile written directly from the run's measurements — the QC→sim
coupling this plan's §6 anticipated, now in Rust with no Python fit step),
and the MultiQC-style HTML report (inline SVG, no JS, accessible tabs).
fastcheck now has 54 tests plus the Python differential harness below.
**What remains is phases 3–5 of this document** — `--overrep-mode`,
`--max-mem` admission, adapter-per-position, per-tile modules, `--only`
gating and the ground-truth recovery test — now written up as TDD tasks in
`docs/plans/2026-10-12_fastcheck-exactness.md`.

Built: one-pass reader (plain / `.gz` / BGZF block-parallel via `gzp`), exact
integer counters, basic statistics, per-position base + N + quality matrices,
per-read mean-Q and length histograms, exact quality-byte census with an
explicit Phred-offset verdict, `--reads`/`--bytes` budgets with self-declared
partial runs, `fastsuite.report.v1` JSON, per-position TSV, memory projection
printed before the numbers. 23 tests, plus a differential harness against an
independent Python implementation (`tests/reference_qc.py` +
`tests/differential.py`) that compares **every** metric, both histograms and all
150 per-position rows: identical on the 6-read fixture and on 20 000 real reads
from fastsim, in both the serial-gzip and the BGZF-parallel paths.

Measured on 2 M reads x 150 bp (16 cores, /tmp = tmpfs, so RAM-backed I/O):

| run | wall | peak RSS |
|---|---|---|
| `fastcheck -j1`, uncompressed (pure scan, no inflate) | 4.51 s | 15.5 MB |
| `fastcheck -j1`, single-member BGZF | 5.49 s | 15.8 MB |
| `fastcheck -j16`, single-member BGZF | 4.54 s | 15.6 MB |
| `fastmaps` (current), same BGZF file | 4.76 s | 33.2 MB |
| `fastcheck --reads 5000` on a 20 k-read file (budget path) | 10 ms | 17 MB |

What that says honestly:

* **The analysis is the bottleneck, not inflate.** Uncompressed input (zero
  inflate cost) still takes 4.51 s, so ~75 % of the BGZF wall clock is the
  counter scan; parallel inflate accounts for 5.49 s -> 4.54 s.
* **Single-file scan is single-threaded.** `-j` parallelises inflation and
  multiple files, not one file's scan. Closing that gap is the next task:
  dispatch BGZF blocks to workers, each with its own `Accum`, merged at the end
  (the merge is exact and cheap: 5 x R + 256 x R counters plus histograms).
  Expected 3-4x on one big file, to be *measured*, not assumed.
* **Layout experiment, negative result.** Storing the per-position planes
  plane-major (`plane[code * max_pos + pos]`) to improve locality made the scan
  **25 % slower** (5.14 s vs 4.10 s) because the quality plane is
  `256 x max_read_len x 8 B` = 2 MB at the default `--max-read-len 1000`.
  Position-major keeps the working set at `read_len x 2 KB` = 300 KB. Reverted;
  the reasoning is recorded in `src/accum.rs` so it is not retried blindly.
* Memory is flat in dataset size: 15.6 MB peak RSS regardless of file size, and
  the accumulator state is printed (0.3 MB per file at R=150).

Still open from the phase list: duplication/over-representation/adapters/k-mers
/ per-tile modules (phases 3-4), and the `fasthelper` step + HTML wiring
(phase 5).

---

## 1. Why the current fastmaps cannot just grow into this

Evidence, `tools/fastmaps/src/`:

| # | problem | evidence |
|---|---|---|
| 1 | **gzip decode is single-threaded per file** — the actual bottleneck | `encode.rs`: `get_reader()` uses `flate2::read::MultiGzDecoder`. One inflate stream per file, serial, regardless of cores. |
| 2 | **parallelism is per file, not per core** | `main.rs`: `num_workers = available_parallelism().min(total_files)` with a work queue of *files*. R1+R2 → 2 workers on a 32-core box. |
| 3 | **writes 8 bytes per read to temp files** | `worker.rs` `SpillChunk.indices: Vec<(u16,u16,u16,u32)>`, consumed by `spiller.rs` → `heatmap_row_*.tmp`. 3×10⁸ reads ⇒ **2.4 GB written and re-read**, for a grid that is 178 KB dense. |
| 4 | **keeps read payloads in RAM** | `worker.rs` `PreviewRead { header, seq, qual }`, `preview_n = 5` per pixel, and `spiller.rs` holds `preview_rows: HashMap<u16, Vec<PreviewRead>>` until the end. **Fixed 2026-10-18**: previews stream to per-row temp files (spiller) instead of being retained; see `docs/BENCHMARKS.md` §2. |
| 5 | **counters are HashMaps keyed by (x,y)** | `spiller.rs` `counts/counts_r1/counts_r2: HashMap<(u16,u16), …>` — hashing per read, and no dense-array fast path. |
| 6 | **no per-position accumulator at all** | there is no per-position base or quality array anywhere in the crate; that is the core of a FastQC replacement. |
| 7 | **short-read assumption baked in** | grid is `(max_seq_len+1)²`; a 10 kb ONT read ⇒ 10⁸ pixels. |

Points 1–5 are the "I/O is the bottleneck" problem stated directly: today the
tool decodes one stream serially and pays a disk round-trip per read to build a
structure that fits in L2 cache. The rewrite fixes the data-structure choice
first (dense counters), then the reader.

---

## 2. Design principles (non-negotiable)

1. **One pass, ever.** Every module shares one decode. A second pass exists only
   where exactness requires it (see over-representation, §5.9) and is an
   explicit flag that the report records.
2. **Exact integers, no estimates.** Counts are `u64`. Rates are computed from
   exact counts. No Bloom filters, no "up to 100 000 sequences" sampling.
3. **No hidden parameters.** Every threshold/tolerance is a flag, its value is
   echoed into the report, and the default is documented there.
4. **No verdicts.** Measured values + evidence only. Pass/fail comes from the
   user's own `--fail-on metric>value` (fasthelper already supports this), and
   the report says which rule was applied.
5. **Label the models.** Where a quantity is not a measurement (k-mer
   expectation, insert-size bound, GC expectation), the JSON carries
   `"model": {…}` naming the assumption and its inputs.
6. **Memory is projected before the run** and the run refuses to start above the
   ceiling (`--max-mem`, defaults to the process limit) — the fastdedup
   admission model, reused.
7. **Deterministic.** Same input ⇒ byte-identical JSON except timings. One test
   asserts it.
8. **Partial runs are self-declaring.** If the user bounds the run
   (`--reads N`, `--bytes N`, `--first N`), output carries
   `"partial": {"reads_analyzed": …, "bytes_read": …, "reason": "budget"}`.
   Never a silent sample.

---

## 3. The reader (the I/O work)

Data flow: `files → block sources → worker pool → per-worker accumulators → merge → report`.

| input | strategy |
|---|---|
| **BGZF `.gz`** (Illumina FASTQ is BGZF: gzip header FEXTRA + `BC` subfield) | block-parallel inflate. Detect with the existing 4-byte test from `fastprep/src/io.rs::looks_like_bgzf`, then `gzp` parallel inflate. Blocks are independent ⇒ real core scaling on one file. This is the single biggest win and the code already exists next door. |
| plain `.gz` (non-BGZF, e.g. pigz output) | sequential inflate per file, but N files inflate concurrently; report the mode per file so a slow case is visible. |
| uncompressed | `mmap` + chunk boundaries aligned to record starts (`\n@` scan forward from a guessed offset), workers take chunks. No copy of the payload. |
| `.zst/.bz2/.xz` | same fallback ladder (compressed ⇒ serial per file unless the format has independent blocks). |

Scheduling: a work queue of **blocks**, not files, so `-j 32` uses 32 cores on a
single file. R1/R2/lanes are just more sources.

Record handling: one reusable `Record` buffer per worker (header/seq/qual
slices into it), parsed with a hand-rolled FASTQ scanner — no `String`
allocation per read, no `needletail` iterator overhead in the hot loop
(needletail stays for format detection and for the plain-gzip path).

Module gating: each module declares `needs()` — `Bases`, `Quals`, `Length`,
`Header`. The parser only materialises what the union of enabled modules needs:

* `--no-quality` ⇒ quality bytes skipped (still parsed for framing), no quality
  accumulators, no quality decode cost.
* `--only dup` ⇒ bases hashed into fingerprints, nothing else touched.
* `--no-heatmap` ⇒ the 2-D grid is never built (grid is the only structure that
  needs `max_seq_len` up front).

Budgeting (all explicit, all reported):

* `--reads N` — stop after N records. The report states exactly how many were
  read and that the run was partial.
* `--bytes N` — stop after N compressed bytes (useful when the input is a 34 TB
  directory and you want a 1-minute answer; the JSON says so).
* both stop the *sources*, not the analysis.

**No spill in the default path.** The heatmap grid becomes a dense
`u64[(R+2)²]` (178 KB at R=150) plus a small overflow map only if a read
exceeds `--max-read-len`. Corner extraction (the reads inside a region) stays
opt-in (`--preview`, `--extract-corners`) because it is the only feature that
needs payloads retained, and it is what forces today's spill+RAM cost.

Expected outcome: for a 10 M-pair BGZF file, decode across all cores, zero temp
I/O, and a resident set in the low tens of MB even with every QC module on.
Benchmark phase 5 must show this as wall-clock and RSS, not as a claim.

---

## 4. Accumulators (per mate unless stated)

| accumulator | exact definition | size at R=150 |
|---|---|---|
| `n_reads`, `n_bases` | counters | 16 B |
| `length_hist[len]` | exact length distribution | (R+1)×8 B = 1.2 KB |
| `base_pos[pos][A/C/G/T/N]` | per-position base counts | 5(R+2)×8 ≈ 6 KB |
| `qual_pos[pos][q]` | per-position quality histogram, `q ∈ 0..=QMAX`, QMAX from `--qmax` (default 41; 93 for legacy Phred+33 extremes) | 42×R×8 = 50 KB |
| `qual_byte_census[byte]` | 256-bin histogram of raw quality bytes ⇒ **exact Phred-offset verdict** (`+33`/`+64`), no guessing | 2 KB |
| `gc_hist[gc]`, `meanq_hist[q]` | per-read GC (N excluded, denominator = ACGT count, stated) and integer-floor mean Q | 1.2 KB + 336 B |
| `heatmap[(R+2)²]` | the composition plane, now dense | 178 KB |
| `adapter_pos[pos]` | reads whose adapter match **starts** at `pos`, plus overlap-length sum | R×8 = 1.2 KB |
| `kmer_pos[pos][kmer]` | exact counts, `k` = `--kmer-k` (4 default…7 max) | 4^k×R×8: 0.8 MB @k=4, 3.1 MB @k=5, 12.6 MB @k=6, 50 MB @k=7 (per mate) — printed as the cost of the dial |
| `dup_table` | 64/128-bit fingerprints of the sequence (quality ignored — stated), multiplicity exact | opt-in, projected: `unique×8 B / 0.8`; refused above `--max-mem` |
| `overrep` | sequences with multiplicity ≥ `--overrep-min-count` | see §5.9 |
| `tile_pos[tile][pos][q]` | per-tile per-position quality, only if headers carry instrument coords and `--tiles` is on | num_tiles×R×42×8 ≈ 5 MB per 100 tiles |

Total with everything on except dup/overrep/kmer: **≈ 250 KB per mate**, plus
the k-mer dial and the duplication table. Memory is a function of read length
and enabled modules — never of the dataset size (except the two table-based
modules, whose growth is exactly what the projection reports).

Long reads: arrays grow with `R`; `--max-read-len` (default 1000) bounds them,
reads beyond it are counted in `length_hist` and reported as
`reads_over_max_len`, with their bases excluded from per-position arrays (the
report says so). Binned position mode (`--pos-bins 100`) is available when the
user prefers coarse positions over exclusion — again, labelled.

---

## 5. Modules: FastQC/Falco → exact definition

| # | module | exact definition we compute | one pass? |
|---|---|---|---|
| 1 | Basic statistics | file, format, encoding (from §4 census), total sequences, total bases, length min/max/mean/median (exact from `length_hist`), %GC | ✅ |
| 2 | Per-base sequence content | `base_pos[p][b] / (Σ_b)`, N reported separately, never folded into ACGT | ✅ |
| 3 | Per-sequence GC content | `gc_hist` normalised; the theoretical distribution is *not* drawn unless asked (`--gc-theory`), and then it is labelled as a model with the observed GC as its only input | ✅ |
| 4 | Per-base N content | part of §2, reported as its own curve | ✅ |
| 5 | Per-base quality | `qual_pos[p][q]` ⇒ mean, median, and quantiles per position (all from the exact histogram) | ✅ |
| 6 | Per-sequence quality | `meanq_hist` | ✅ |
| 7 | Per-tile quality *(FastQC's tile modules)* | per-tile × position × quality, from parsed Illumina headers; only when headers have the 7-field form | ✅ (optional) |
| 8 | Sequence length distribution | `length_hist` verbatim | ✅ |
| 9 | **Duplication levels** | exact multiplicity histogram from fingerprints (`count==1: n₁, count==2: n₂, …`), reported as counts and as `% of unique` and `% of total`; "total deduplicated percentage" computed from exact numbers | ✅ |
| 10 | **Over-represented sequences** | exact: sequences whose multiplicity ≥ `--overrep-min-count` (default 2 = "appears more than once"). Exactness needs identity, which needs either the sequence table or a second pass — see below | ⚠ flag |
| 11 | Adapter content per position | for each adapter in the catalogue (reuse `fastprep/src/adapter.rs`: TruSeq R1/R2, Nextera, SmallRNA, user FASTA) and each `p`, the count of reads with an adapter match starting at `p`, using the **documented** matcher (`--min-adapter-overlap`, `--adapter-error-rate`); plus the per-adapter totals and the k-mer evidence that selected them | ✅ |
| 12 | K-mer content | exact `kmer_pos[p][kmer]` for the top-N k-mers, plus for each reported k-mer an expected count from **this run's measured per-position composition** under a stated independence model and the binomial tail probability. Numbers and model both printed; no verdict | ✅ |
| 13 | Quality-offset census | the 256-bin quality-byte histogram; if the mode is `< 33` the report says the file is malformed instead of silently shifting | ✅ |

### 5.9 Over-representation — the one honest trade-off

Exact identity requires one of:

| mode | cost | report field |
|---|---|---|
| `--overrep-mode none` (default) | zero | duplication levels still exact (§9) |
| `--overrep-mode fingerprint` | + `unique×8 B/0.8` RAM; reports multiplicity exactly, but cannot name the sequence | `names_available: false` |
| `--overrep-mode second-pass` | one extra read of the file, only to emit the sequence text for fingerprints above the threshold; memory stays at the fingerprint table | `passes: 2` |
| `--overrep-mode table` | RAM for the unique sequences themselves (2-bit packed, `⌈len/4⌉+8` B each; 3×10⁸ unique @150 bp ≈ 11 GB) — refuses to start if it exceeds `--max-mem` | `unique_sequences_retained: true` |

FastQC picks a fourth option: analyse the first 100 000 sequences and hope. We
report which of the four the user chose — and never present a sample as the
whole.

---

## 6. Output

```
<outdir>/
  report.json          fastsuite.report.v1 — every metric, exact, + provenance
  positions.tsv        per-position matrices (base/qual/N), optional --positions-tsv
  positions.bin        same, binary (varint, the existing fastmaps format) for re-plotting
  heatmap.{png,svg,json,data,preview}   from the composition module (fastmaps behaviour kept)
  report.html          rendered by fasthelper from report.json (single file, no JS)
```

`fastmaps` remains as a thin alias/step for the heatmap + corner-extraction
workflow; the QC binary owns the decode and the metrics, so nothing decodes a
file twice. `fasthelper` gains the module as a step name (`fastcheck`), and
`fastsuite.toml` can request `--only` module subsets per step.

---

## 7. Verification (how "actually good" gets proven)

1. **Hand-computed fixtures.** Every module has a tiny FASTQ with counts written
   out by hand in the test (e.g. a 4-read file where per-position quality is
   known to the unit digit).
2. **Independent reference implementation.** A short Python script in
   `tests/` recomputes the same statistics from the same fixture with a
   different code path (and, where available, with `pysam`/`biopython` as a
   third check). Fixture files + expected JSON are committed.
3. **Ground-truth recovery via fastsim — the decisive test.** `fastsim`
   generates reads from a *known* profile: known per-position quality curves,
   known head/tail G-run histograms, known adapter readthrough rate, known
   duplicate rate, known GC. The QC binary must recover them within stated
   tolerances (e.g. per-position mean Q within ±0.05, duplicate rate within
   ±0.2 pp, adapter rate within ±0.3 pp), and the test fails loudly otherwise.
   This is the difference between "my numbers look plausible" and "my numbers
   are right".
4. **Differential against FastQC/Falco** on the same input for the modules they
   share: exact equality where both are exact (total sequences, length
   distribution, per-position counts, per-base quality histograms, GC), and a
   written explanation where they differ (they sample; we do not).
5. **Determinism:** run twice, `diff` the JSON minus timestamps.
6. **Memory ceiling:** a synthetic 20 M-read run asserted to stay under
   `--max-mem`, and a fingerprint-mode run asserted to refuse cleanly above it.
7. **I/O accounting:** the test suite asserts that the default path creates **no
   temp files** (`temp_dir` empty after a run) — the regression guard for the
   2.4 GB-per-300 M-read spill we are deleting.

---

## 8. Phases

| phase | work | files | estimate |
|---|---|---|---|
| **0** | Decisions (§9) + freeze the metric names in `REPORT-SCHEMA.md` | docs | 1 h |
| **1** | **Skeleton + reader.** New crate `tools/fastcheck`: CLI, record scanner, block sources (BGZF parallel via `gzp`, plain-gzip serial, mmap chunked), worker pool over blocks, merge, JSON writer. No modules yet; the "basic statistics" + length/GC/meanQ modules only. Tests: fixtures + determinism + no-temp-files | new crate, ~1200 LOC | 2 d |
| **2** | **Position accumulators.** base/qual/N per position, quality-byte census, quantiles from histograms. Differential vs the Python reference on 3 fixtures | `src/pos.rs`, tests | 1.5 d |
| **3** | **Table modules.** duplication-level histogram (reuse fastdedup's hasher via the shared crate, §9.2), memory projection + admission, overrep modes (§5.9) | `src/dup.rs`, `src/report.rs` | 2 d |
| **4** | **Sequence modules.** adapter-per-position (reuse fastprep's catalogue/matcher), k-mer with the stated expectation model, per-tile modules | `src/adapter.rs`, `src/kmer.rs`, `src/tile.rs` | 2.5 d |
| **5** | **Integration + proof.** heatmap module folded in (fastmaps kept as alias), `fasthelper` step + HTML rendering, fastsim ground-truth recovery test, differential vs FastQC/Falco, throughput/RSS matrix vs the current fastmaps and vs FastQC | wiring + tests + bench | 2 d |

Total: **~10 working days** with agent support, in the order above. Phases 1–2
alone already produce a usable tool (statistics + per-position + quality
census) and are the minimum to retire Falco.

## 9. Decisions needed before phase 1

1. **Name.** `fastcheck` (my proposal), or `fastqc` (honest, but collides with
   the tool being replaced in every search and on every PATH), or keep growing
   `fastmaps`?
2. **Shared core crate.** Extract `fastsuite-core` (report writer, hasher, BGZF
   reader, FASTQ scanner) out of the five tools, or copy the ~600 lines into
   `fastcheck` and keep the tools independent? Proposal: extract — it is the
   only way `fastsuite.report.v1` stays true across tools, and the BGZF reader
   is needed by three of them.
3. **Fingerprint sharing with fastdedup.** Reuse one implementation (library
   call) so `fastdedup --json` and `fastcheck --modules dup` cannot disagree.
4. **Long reads.** Bound at `--max-read-len 1000` with exclusion reporting
   (proposal), or implement binned positions as the default for R > 1000?
5. **Tile modules.** Worth the header parsing + 5 MB/100 tiles? Only useful for
   raw Illumina runs; proposal: implement, default off.
6. **Second pass.** Allowed at all for exact over-representation, or is a second
   read of a 34 TB dataset forbidden by policy? (Proposal: allowed, off by
   default, `passes: 2` recorded.)

## 10. Risks

| risk | impact | mitigation |
|---|---|---|
| BGZF detection fails on a real dataset (non-BGZF gzip in the wild) | single-threaded decode, the thing we set out to fix | per-file mode is reported; benchmark both paths; file-level parallelism still gives N cores for multi-lane runs |
| k-mer expectations get read as verdicts | exactly the "heuristic bullshit" failure mode | the model and its inputs are in the JSON; no PASS/FAIL field exists in the schema |
| duplication table memory on 3×10⁸ unique reads | OOM at hour 3 | projection printed + admission refusal (fastdedup's proven path); `64` vs `128` bit dial with the collision probability stated |
| per-position arrays for very long reads | RAM blow-up | `--max-read-len` + exclusion reporting, or `--pos-bins` |
| differential vs FastQC disagrees on a shared module | credibility | disagreements are written down with the reason (they sample, or use a different offset), not smoothed over |
