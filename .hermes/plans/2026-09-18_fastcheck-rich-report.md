# fastcheck: rich report (polyX, repeats, heatmap) + before/after side-by-side

> **For Hermes:** use subagent-driven-development to implement task-by-task.

**Goal:** make `fastcheck --html` a complete QC report — polyX run-length distributions
(4 base tabs), repeated sequences (full sequences **and** k-mer spectrum), the 2-colour
polyG-head signal, and the fastmaps-style heatmap — and let the report compare **raw vs
processed** reads side by side, both by feeding processed files into fastcheck directly and
via a new `fasthelper compare` verb that merges two separate runs. Every new number flows
through `fasthelper concatenate` (`report` and `--positions`).

**Architecture:** all new measurements are computed in the existing single pass (new `Accum`
fields), written to the JSON report (`metrics` for per-sample numerics, tool-specific blocks
for grids/distributions) and the positions TSV (new per-position columns), and rendered as
new inline-SVG sections in the HTML. Side-by-side is a second scan of the processed files
(`--processed-1/--processed-2`) whose results are stored alongside the raw ones. No new deps.

**Tech stack:** Rust (fastcheck, fasthelper), inline SVG (no JS), serde_json. Reuses
fastmaps' heatmap encoding (`x = #(T+C)`, `y = #(T+G)`, row-major `index = y*size + x`).

---

## Context / assumptions (verified in code)

- **`concatenate report`** (`fasthelper/src/report.rs::adapt_file`): for a `fastsuite.report.v1`
  JSON it iterates `metrics` and keeps **every numeric value**. → New per-sample numerics in
  `metrics` flow through `concatenate report` for free.
- **`concatenate --positions`** (`fasthelper/src/positions.rs::read_positions`): reads the
  positions TSV header at runtime, every column after `pos` is a metric. → New per-position
  columns become matrix metrics for free. A before/after pair is just two samples in the
  matrix (already supported).
- **Nested structures** (heatmap grid, polyX distributions, repeated-sequence list, k-mer
  spectrum) are arrays/objects → they do NOT flow through `concatenate`; they live in the
  JSON tool-specific blocks + HTML. Their numeric summaries go in `metrics`.
- **Current HTML sections** (`html.rs`): Summary, Per-base sequence content, N content by
  position, Distributions, Insert sizes, Observations. New sections slot in after
  "Per-base sequence content".
- **Heatmap semantics** (fastmaps `encode.rs`/`format.rs`): per read `x = #(T+C)`,
  `y = #(T+G)`, `index = y*size + x`, per mate (`counts_r1`, `counts_r2`). fastcheck uses
  **64×64** (~32 KB per mate).
- **Domain:** raw, uniform-length short reads (PE or SE). Pre-trimmed libraries out of scope.

---

## Compatibility contract (the user's hard requirement)

| New data | Where it lives | Flows through concatenate? |
|---|---|---|
| per-sample numerics (polyX max, head-dropout, repeated-seq count, duplication, heatmap peak, k-mer top) | `metrics` | ✅ `concatenate report` (automatic) |
| per-position G fraction, per-position polyX reach | positions TSV new columns | ✅ `concatenate --positions` (automatic) |
| polyX distributions, repeated-sequence list, k-mer spectrum, heatmap grids | JSON tool-specific blocks + HTML | ❌ carried in JSON/HTML only (documented) |
| before/after pair | two samples in the matrix / `fasthelper compare` | ✅ both verbs |

No change to `fasthelper` is required for the numerics to flow. `fasthelper compare` is a new
verb that reuses the existing `collect`/`adapt_file` machinery.

---

## Design

### New `Accum` fields (one pass, per mate)
- `polyx3: [Vec<u64>; 4]` — trailing (3′) run length per base (A,C,G,T), capped at read_len.
- `polyx5: [Vec<u64>; 4]` — leading (5′) run length per base, capped at read_len.
- `heatmap: Vec<u64>` — 64×64 grid, `index = y*64 + x`, `x = #(T+C)`, `y = #(T+G)`.
- `seq_counts: HashMap<Vec<u8>, u32>` — bounded over-represented full-sequence counter.
- `kmer_counts: HashMap<u64, u32>` — bounded 12-mer counter (packed into a u64: 2 bits/base).
- `distinct_seqs: u64` — distinct full sequences (for the duplication level).

### Bounded counters (new `repeated.rs`)
Streaming with singleton eviction: when distinct entries exceed a cap (default 1,000,000),
evict every entry with count == 1 (a singleton can never be over-represented; losing its
first occurrence costs ≤1 of a count that must reach a large threshold to be reported).
Memory stays bounded on 50 M-read files. Two counters: full sequences (report those with
`count >= max(10, 0.1% of reads)`, top 20) and 12-mers (report top 20 with count + fraction).

### Side-by-side (fastcheck `--processed-1/--processed-2`)
- New optional CLI: `--processed-1 <fq>` / `--processed-2 <fq>` (SE: `--processed-1` only).
- fastcheck scans the processed files with the same Accum, stores them as a second report set.
- JSON: the existing top-level report is the **before**; a new `processed` block holds the
  after metrics + blocks (same schema shape).
- HTML: every section gains a before/after view (two columns or overlaid curves).
- positions TSV: written for both (e.g. `positions.tsv` + `positions.processed.tsv`), so
  `concatenate --positions` sees them as two samples.

### `fasthelper compare` verb (new)
`fasthelper compare --before <dir|json> --after <dir|json> -o out/` → a side-by-side report:
- **HTML**: per-metric table with before | after | delta; per-position curves overlaid;
  heatmaps side by side; polyX/repeats side by side.
- **TSV**: `metric, before, after, delta` (and a per-position long form).
- Reuses `collect`/`adapt_file` (already reads fastcheck JSONs + positions TSVs); no new
  parsing.

### New JSON report content
- `metrics` additions (numeric → concatenate): `polyx3_max`, `polyx5_max`,
  `polyx3_frac_gt_5`, `polyx5_frac_gt_1`, `head_g_frac` (position-1 G fraction),
  `head_dropout` (library G baseline − position-1 G fraction; >0 = 2-colour first-cycle
  dropout), `repeated_seq_count`, `repeated_seq_frac`, `distinct_seqs`, `duplication_level`
  (1 − distinct/reads), `kmer_top_count`, `kmer_top_frac`, `heatmap_peak`, `heatmap_peak_x`,
  `heatmap_peak_y`.
- New tool-specific blocks (outside `metrics`): `polyx` (per-mate `{3p: {base: {len: count}},
  5p: {base: {len: count}}}`), `repeated_sequences` (top list `{seq, count, frac}`),
  `kmer_spectrum` (top list `{kmer, count, frac}`), `heatmap` (per-mate `{size, counts_r1,
  counts_r2}` — same shape as fastmaps `heatmap.json`).
- When `--processed-*` given: a `processed` block mirroring the above.

### New positions TSV columns (per position)
`g_frac` (G fraction), `polyx3_frac` (fraction of reads whose trailing run reaches this
position), `polyx5_frac` (fraction whose leading run reaches this position). These flow
through `concatenate --positions`.

### New HTML sections (inline SVG, no JS)
1. **Repeated sequences** — two views: the top over-represented full sequences (table: seq,
   count, %) and the top 12-mer spectrum (table + bar chart); plus duplication level and
   distinct-sequence count.
2. **polyX** — **4 tabs (A, C, G, T)**, each with the 3′ and 5′ run-length distributions
   (bar charts); a callout for the G-specific head-dropout (2-colour chemistry) with the
   position-1 G fraction vs the library baseline.
3. **Heatmap** — **two tabs**: "side by side" (before/after panels) and "full" (the large
   64×64 heatmap), colour-scaled to the peak, matching fastmaps' orientation.
4. When `--processed-*` given, every existing section (Summary, Per-base content, N, etc.)
   gains a before/after view.

---

## Tasks (TDD, bite-sized)

### Task 1: polyX run-length measurement
**Files:** `tools/fastcheck/src/accum.rs`, new `tools/fastcheck/src/polyx.rs`.
- **Step 1 (test):** reads `ACGTGGG` (trailing G run 3), `GGGACGT` (leading G run 3),
  `AAAACGT` (leading A run 4), `ACGT` (0/0), all-G (run = read_len) → assert
  `polyx3[G][3]`, `polyx5[G][3]`, `polyx5[A][4]` counts. Cap test: a run longer than read_len
  folds into the last bin.
- **Step 2:** run, expect FAIL (no module).
- **Step 3:** implement `measure_polyx(seq) -> (leading_base, leading_len, trailing_base,
  trailing_len)`; wire into `Accum::observe`.
- **Step 4:** run, expect PASS. Commit.

### Task 2: bounded repeated-sequence + k-mer counters
**Files:** new `tools/fastcheck/src/repeated.rs`, `accum.rs`.
- **Step 1 (test):** feed 1000 reads, 5 of one sequence, 2 of another, rest unique → assert
  the two are reported, singletons are not; assert the cap eviction keeps memory bounded (a
  synthetic run past the cap still reports the frequent sequence). Same for 12-mers (packed
  u64 encoding: `ACGT` → 0b00011011).
- **Step 2:** run, expect FAIL.
- **Step 3:** implement `RepeatedCounter` (insert + evict-singletons + report) and the
  12-mer packer.
- **Step 4:** run, expect PASS. Commit.

### Task 3: heatmap grid
**Files:** new `tools/fastcheck/src/heatmap.rs`, `accum.rs`.
- **Step 1 (test):** `TTCCGGAA` → x = #(T+C) = 4, y = #(T+G) = 4, index = 4*64+4; assert the
  cell increments. Assert the encoding matches fastmaps (`encode.rs`).
- **Step 2:** run, expect FAIL.
- **Step 3:** implement `encode(seq) -> (x, y)` + grid increment in `observe`.
- **Step 4:** run, expect PASS. Commit.

### Task 4: JSON report — metrics + blocks + positions columns
**Files:** `tools/fastcheck/src/report.rs`, `accum.rs` (expose new fields).
- **Step 1 (test):** on a synthetic Accum, assert the new `metrics` keys exist and are
  numeric; assert `polyx`, `repeated_sequences`, `kmer_spectrum`, `heatmap` blocks exist;
  assert the positions TSV gains `g_frac`, `polyx3_frac`, `polyx5_frac`.
- **Step 2:** run, expect FAIL.
- **Step 3:** implement the report additions + positions columns.
- **Step 4:** run, expect PASS. Commit.

### Task 5: HTML sections (polyX tabs, repeats, heatmap tabs)
**Files:** `tools/fastcheck/src/html.rs`.
- **Step 1 (test):** assert the HTML contains the new section headings, the 4 polyX tabs,
  the heatmap's two tabs, and the repeated-sequence table (string assertions).
- **Step 2:** run, expect FAIL.
- **Step 3:** implement the sections (polyX 4-tab charts, repeats table + k-mer chart,
  heatmap side-by-side/full tabs).
- **Step 4:** run, expect PASS. Commit.

### Task 6: fastcheck side-by-side (`--processed-1/--processed-2`)
**Files:** `tools/fastcheck/src/cli.rs`, `main.rs`, `report.rs`, `html.rs`.
- **Step 1 (test):** a synthetic raw + processed pair → assert the JSON has a `processed`
  block and the HTML has before/after views.
- **Step 2:** run, expect FAIL.
- **Step 3:** implement the second scan + `processed` block + before/after HTML.
- **Step 4:** run, expect PASS. Commit.

### Task 7: `fasthelper compare` verb
**Files:** `tools/fasthelper/src/cli.rs`, new `tools/fasthelper/src/compare.rs`, `report.rs`.
- **Step 1 (test):** two synthetic fastcheck JSONs → assert the compare TSV has
  `metric, before, after, delta` rows and the HTML has the side-by-side table.
- **Step 2:** run, expect FAIL.
- **Step 3:** implement the verb (reuse `collect`/`adapt_file`).
- **Step 4:** run, expect PASS. Commit.

### Task 8: concatenate compatibility end-to-end
**Files:** none (verification) + `docs/REPORT-SCHEMA.md`, `docs/CLI-SPEC.md`.
- **Step 1:** run `fastcheck --json r.json --positions-tsv p.tsv --html r.html` on a real
  sample (ERR10033173, 100 k pairs), and again with `--processed-*` on a fastprep output.
- **Step 2:** `fasthelper concatenate report r.json` → new metrics appear in the matrix.
- **Step 3:** `fasthelper concatenate --positions p.tsv` → new columns appear as metrics.
- **Step 4:** `fasthelper compare --before raw/ --after processed/ -o cmp/` → side-by-side
  report renders.
- **Step 5:** vision-check the HTML (tabs, heatmap, no overflow).
- **Step 6:** update `docs/REPORT-SCHEMA.md` and `docs/CLI-SPEC.md`. Commit.

---

## Files likely to change
- `tools/fastcheck/src/accum.rs` — new fields + per-read measurement.
- `tools/fastcheck/src/polyx.rs` (new) — run-length measurement.
- `tools/fastcheck/src/repeated.rs` (new) — bounded full-sequence + k-mer counters.
- `tools/fastcheck/src/heatmap.rs` (new) — grid + SVG renderer.
- `tools/fastcheck/src/report.rs` — new metrics, blocks, positions columns, `processed` block.
- `tools/fastcheck/src/html.rs` — new sections + before/after views.
- `tools/fastcheck/src/cli.rs`, `main.rs` — `--processed-*`, `mod` declarations.
- `tools/fasthelper/src/cli.rs`, `report.rs`, new `compare.rs` — the `compare` verb.
- `docs/REPORT-SCHEMA.md`, `docs/CLI-SPEC.md` — document the additions.

## Tests / validation
- Unit: polyX measurement, repeated + k-mer counters (incl. eviction bound), heatmap
  encoding, report keys/blocks, positions columns, HTML section/tab presence, `compare` TSV.
- Integration: real-sample run + `fasthelper concatenate` both verbs + `fasthelper compare`
  + vision check of the HTML.
- Existing fastcheck (29 + 10) and fasthelper tests must stay green.

## Risks / tradeoffs / open questions
- **Repeated-sequence memory:** the singleton-eviction counter is approximate (a sequence
  whose first occurrence is evicted loses 1 count). For over-representation thresholds
  (≥0.1%) this is negligible; document it. Exactness would need a two-pass or count-min
  sketch.
- **k-mer spectrum:** 12-mers packed into u64 (2 bits/base) keep memory small; the eviction
  counter bounds it. A different k (8 or 16) is a one-line change.
- **polyX cap at read_len:** runs longer than read_len fold into the last bin; fine for
  uniform-length reads.
- **`--processed-*` doubles the scan cost** (two full passes). Acceptable for a QC report;
  document that the processed pass is optional.
- **`fasthelper compare` scope:** it compares fastcheck reports (and any v1 reports). It
  does not re-run fastcheck; the user runs fastcheck on raw and processed separately (or
  fastcheck `--processed-*` in one go) and `compare` merges the two JSONs.
- **Open question:** should `fasthelper compare` also accept the fastcheck `--processed-*`
  single-JSON (which already contains both) and split it, or only two separate JSONs? Plan
  assumes two separate runs; the single-JSON case is a trivial extension.
