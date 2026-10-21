# Phase 2 — the speed levers, measured

Running evidence log for `docs/plans/2026-10-21_implementation-round2-speed-undertrim.md`
Phase 2. Every entry names the artifact it came from. Claim tiers follow the
usual convention: **[M]** measured with a named artifact, **[E]** estimated,
**[V]** verified against current source or binary.

The reference machine for all speed numbers is the **desktop**
(`cachyos-x8664`, 16 threads, 46 GB RAM, NVMe), per `bench/profile/README.md`.
The laptop is used only for functional verification (tests, byte-identity).

---

## Wave A — free wins

### Task 2.14 — fastmaps `lto = "fat"`, `codegen-units = 1`

**[V] The premise held.** fastmaps' `[profile.release]` carried only
`opt-level = 3`; the other five tools all carried both flags. This was the only
tool left without them.

**Byte-identity [M]** — 2 M-pair corpus (`bench/profile/fastmaps_ab.py`), all
five grid files compared by md5:

```
GRID-IDENTICAL
GRID-MD5-IDENTICAL        (5 files: heatmap.data, heatmap.preview,
                           heatmap.json, heatmap.png, heatmap_3d.html)
```

The `fastmaps.json` report differs **only** in `command` (the binary path) and
`elapsed_s` (a timing field); normalised, the two reports are byte-identical:

```
normalised identical: True
data fields: 971 bytes vs 971 bytes
```

**Binary size [M]:** 1 692 288 B → 1 308 368 B, a **23 % reduction**.

**2 M-pair timing [M]**, 5 reps, wall/cpu `min` of the run set:

| build | wall min | wall median | cpu min | peak RSS |
|---|---:|---:|---:|---:|
| baseline (no lto) | 1.859 s | 1.868 s | 3.610 s | 21.4 MB |
| `lto=fat, cgu=1` | 1.847 s | 1.852 s | 3.591 s | **16.4 MB** |

Wall is **not** meaningfully different (0.6 %, inside run-to-run noise);
**peak RSS falls 23 %** (21.4 → 16.4 MB). This is a memory win and a
code-size win, not a wall-time win, and it is reported as such.

**5 M-pair campaign [M]** — three builds from the same tree state, one shared
corpus (883 MB of trimmed FASTQ), 5 reps each, `-j 8`:

| build | wall min | wall median | cpu min | peak RSS | vs baseline |
|---|---:|---:|---:|---:|---:|
| baseline (no lto, std HashMap) | 4.626 s | 4.629 s | 9.044 s | 26.0 MB | — |
| `+ lto` (2.14) | 4.564 s | 4.598 s | 8.880 s | 26.7 MB | −1.3 % wall |
| `+ lto + fxhash + Vec` (2.13) | **4.457 s** | **4.473 s** | **8.702 s** | **21.4 MB** | **−3.7 % wall, −3.8 % cpu, −17.7 % peak** |

Grid byte-identity across all three: `fm-lto: GRID-IDENTICAL`,
`fm-lto-fx: GRID-IDENTICAL`.

**Read the ratio, not the absolute.** The plan's reference row is 3.856 s wall
at 5 M pairs; this campaign's baseline is 4.626 s because the corpus was
produced at `--gz-level 1` rather than 6 and is a different simulation. Only
the three rows above are comparable to each other.

---

### Task 2.13 — fastmaps: cheap hasher, `Vec` for per-file counts

**What changed [V]:**

* the four pixel maps (`counts`, `counts_r1`, `counts_r2`, and the worker's
  local) move from `std::collections::HashMap` (SipHash-1-3) to
  `rustc_hash::FxHashMap` — the suite's existing crate, the same version
  fastcheck pins, so no second hashing dependency was introduced;
* `per_file_counts` moves from `HashMap<u16, u64>` to `Vec<u64>`. The key is a
  dense file id, which is what an array is for.

**Why it is safe [V]:** the keys are `(u16, u16)` pairs derived from read
content, never untrusted network input, so SipHash's DoS resistance buys
nothing here.

**Byte-identity [M]:** `per_file_counts` is serialised by **index order**
(`format.rs` loops `0..files.len()`), not hash order, so the `Vec` is a
drop-in — and an out-of-range id still serialises as `0`, exactly as the old
`.get().unwrap_or(0)` did. Same md5-identical grid as above.

**The merge stays exact [V]:** `spiller.rs` accumulates by index over the
partial's `Vec`, which is exact in any `Done` arrival order (FM-12) and grows
to fit an id beyond the current length rather than dropping it.

**Tests [M]:** fastmaps 96 pass, 0 fail; full `pixi run ci` green at 765.

---

### Task 2.20 — fasthelper: "parse each report once"

**[V] The premise did not hold on current source — no change was made.**

The claim was that `check → aggregate → concatenate → html` re-reads every
report. It does not:

* `cmd_concatenate` calls `collect(&a.inputs, …)` **once** and hands the
  resulting `Vec<ToolReport>` to `write_reports`;
* `write_reports` builds the matrix at most once (`owned_matrix`) and reuses it
  for the TSV, and renders `json`, `md` and `html` from the same slice;
* `runner::run` likewise calls `collect_named` once, then `write_reports`.

Measured rather than assumed, with a parse-dominated corpus
(`bench/scripts/fh_parse_once_probe.py`: 40 reports, 9.3 MB, 40 000 counts
each) — if each format re-parsed, four formats would cost about four times one:

```
1 format  (json)            : 0.048 s
4 formats (json,tsv,md,html): 0.048 s
ratio 4/1                   : 1.00x

VERDICT: parsed ONCE — adding three output formats costs render only.
```

A re-parse would have shown ~4.0x. **Nothing to optimise here**; the probe is
kept as the regression test for the property.

---

### Task 2.8 — fastprep: the PXMT writer

**[V] The premise did not hold, and the real defect was elsewhere.**

The plan said to "give the PXMT writer the same 256 KiB staging the FASTQ
writer already has". It already has it, and more: `open_output` wraps every
writer in a `BufWriter::with_capacity(1 MiB)` — ten times the proposed staging —
and the metadata writer goes through the same `open_output`. A 256 KiB layer in
front of a 1 MiB layer is pure overhead, so it was **not** added.

**The actual defect [M]:** the metadata file was **never compressed**, even when
named `.pxmt.gz`.

`open_output(tp, …)` inferred the codec from the *temp* path, and `part_path`
appends `.part` — so `m.pxmt.gz` became `m.pxmt.gz.part`, whose extension is
`.part`. `encode_output`'s match fell through to the plain branch. The FASTQ
writer never had this bug because it already used `open_output_as(tmp, final,
…)`, which infers from the **final** name; the metadata writer did not.

Evidence, 5 000 reads: the file began with `PXMT`, not `1f 8b` — raw, at
38.8 bytes/read, which matches the plan's own "38.4 bytes per read" figure and
explains it.

**Fix:** `open_output_as(tp, fp, METADATA_COMPRESSION)`. The `.pxmt.gz` name is
now true; a `.pxmt` name still writes plain, unchanged.

| | before | after |
|---|---:|---:|
| on-disk size, 200 000 reads | 8 088 895 B (38.8 B/read) | **513 326 B (2.57 B/read)** |
| cost of `--metadata` (wall) | +0.11 s (plan, pre-fix) | **+0.012 s** |
| cost of `--metadata` (cpu) | +0.20 s (plan, pre-fix) | **+0.084 s** |

A **15.8× smaller file**, and the wall cost of enabling metadata drops by about
9×. Measured with `bench/scripts/meta_cost_probe.py` (200 000 reads, `-j 8`,
5 reps, median).

**Byte-identity of the payload [M]:** compressing the container does not touch
the records. Two runs of the same binary produce decompressed payloads that are
identical byte for byte (193 895 B, md5 `2b9eb83b…`, starting with `PXMT`).

**The second defect this exposed [M]:** once the file really was compressed, the
viewer failed on it —

```
$ fastprep-meta m.pxmt.gz
Error: failed to fill whole buffer
```

The output is **BGZF**, a chain of ~64 KiB gzip *members*, and the viewer used
`flate2::read::GzDecoder`, which decodes only the **first member**. It now uses
`MultiGzDecoder`. The two bugs had been hiding each other: no compression meant
no multi-member stream to mis-read.

**Regression test [M]:** `a_gzipped_metadata_file_is_really_gzip_and_reads_back_whole`
asserts (a) the `.gz` file starts with `1f 8b`, and (b) the viewer reaches
record `bulk5999`, which lives in the final member. Verified to **fail** against
`GzDecoder` and pass against `MultiGzDecoder`. An earlier draft used a fixture
under 64 KiB and passed against the broken decoder — the trap this guard exists
to catch, so the fixture is now 6 000 reads (~197 KiB of metadata, 3 members).

**Flush errors now propagate.** The finalisation was `mw.flush().ok()`, which
renamed a truncated metadata file into place and exited 0 — a D-6 violation that
a real compressor makes more likely. It is now `?`.

---

## Wave B — the gzip lever

### Task 2.4b — fastdedup: flipping the default to level 1

**[M] The plan's own warning does not materialise, but its byte estimate is off
by a factor of five.**

The plan calls this "the single most likely way the gzip lever backfires": the
parallel compressor (gzippy) engages only at `compression_level >= 6`, so a
level-1 default drops the writer onto the single-threaded flate2 path. That
mechanism is real, and `parallel.rs:222` plus the `--help` text both document
it. What it costs, measured on one 5 M-pair corpus (441.5 MB per mate) at
three reps per arm:

| arm | wall | CPU | busy cores | peak RSS | output |
|---|---:|---:|---:|---:|---:|
| `--gz-level 6 -j 8` (today) | 20.26 s | 124.79 s | **6.16** | **3 780 MB** | 390 MB |
| `--gz-level 1 -j 8` | **17.46 s** | 64.51 s | 3.70 | **152 MB** | 729 MB |
| `--gz-level 6 -j 1` | 84.22 s | 84.06 s | 1.00 | 79.6 MB | 390 MB |
| `--gz-level 1 -j 1` | **12.86 s** | 12.83 s | 1.00 | 79.6 MB | 726 MB |

**It does not backfire.** Level 1 at `-j 8` is **1.16× faster** (17.46 s vs
20.26 s) and uses **25× less peak memory** (3 780 MB → 152 MB). Busy cores fall
from 6.16 to 3.70 — the parallel compressor really is disengaged, exactly as the
code says — but the level-1 single-threaded path is fast enough that it still
wins outright.

**The cost is bytes, and the plan understates it by five times.** The plan
predicted "~18 % larger". Measured: **1.87× larger** (390 MB → 729 MB). A
100 M-pair run would go from roughly 7.8 GB to 14.6 GB of output.

Two further facts the plan did not anticipate:

* **Both outputs are standard gzip.** `gzip -t` accepts the level-1 output, so
  the proprietary-format hazard does not apply on this path — it applies only
  to a hypothetical `-j >1` build that forced gzippy below level 6, which the
  code does not do. It falls back to plain flate2 instead.
* **At level 1, `-j 1` beats `-j 8`** (12.86 s vs 17.46 s). Once compression is
  cheap, the dedup hashing is I/O-bound and extra threads only contend. This
  also explains the long-standing `-j 1..4` flatline the plan asked about in
  Task 2.6: it is a *level-6* phenomenon, and at level 1 the curve inverts.

**Recommendation: do not flip this default unilaterally.** The trade is
genuine — 1.16× wall and 25× memory against 1.87× disk — and neither side
dominates. `--gz-level 1` is now the documented answer for a memory- or
time-constrained run; whether it should be the default is a release decision
(C-8), not a measurement outcome.
