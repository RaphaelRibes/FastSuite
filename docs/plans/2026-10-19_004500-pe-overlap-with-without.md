# PE overlap detection — with/without test plan (rev. 2, after Phase 0)

> **For Hermes:** execute with subagent-driven-development, one task per
> subagent, parent verifies every number before it is committed.

**Goal:** decide on measurement whether restoring fastprep's paired-end overlap
detection earns its place — and if it does not, record that decision where the
paper and the code agree.

**Architecture: two cells.** **Cell A** (the existing corpus) is the
*spurious-fire guard*: adapter readthrough is configured in only 53 of the 102
corpus profiles and materialises in just 11,395 of 10.2M pairs (0.11%), so any
trimming the mode does outside those pairs is wrong by construction. **Cell B** (a synthetic cohort
with readthrough, inserts inside `2L-12`) is the *value test*: that is where a
benefit could exist. A mode that passes Cell A and wins Cell B is worth
restoring; a mode that fires on clean reads is a bug regardless of Cell B.

**Tech Stack:** Rust (fastprep CLI, `--json`, PXMT), Python bench harnesses
(`bench/scripts/closed_loop.py`, `trim_report.py` with `FASTPREP_VARIANTS`),
`fastsim` for ground truth, the desktop corpus at `/mnt/Fast/fastsuite-scratch`.

---

## What Phase 0 measured (do not re-derive)

102 datasets, 10.2M pairs, `read_len=150`, anchor 12, ceiling `2L-12 = 288`:

| domain | `I < L` | `L ≤ I ≤ 288` | censored | **mate carried adapter** |
|---|---|---|---|---|
| bacterial_isolate | 11.75% | 72.91% | 15.34% | **0.09%** |
| generalised | 24.37% | 75.63% | 0% | **0.03%** |
| human_gut | 0.83% | 35.06% | 64.11% | **0.01%** |
| metagenomics | 24.46% | 72.80% | 2.74% | **0.26%** |
| rna_seq | 25.56% | 74.44% | 0% | **0.24%** |
| virome | 8.49% | 23.53% | 67.97% | **0.05%** |

Two conclusions, in order:

1. **The overlap window is large** (73-76% in four domains). The mode *can* fire
   on most pairs — window size is not the constraint.
2. **Adapter readthrough is 0.01%-0.26% of pairs**, and the cause is structural
   but narrower than first written: **53 of the 102 profiles are configured with
   `adapter.mode: readthrough`** and 49 with `none`, so the corpus is *not*
   adapter-free — but readthrough only materialises in 11,395 of 10.2M pairs
   (0.112%), because it needs an insert shorter than the read and these libraries
   rarely have one. (The earlier revision of this plan said "adapter.mode: none in
   every profile" after sampling one file; that was wrong and is corrected here.)

**Therefore the corpus arm cannot be the value test.** Running with/without on
this corpus would measure a mode firing on ~0.1% of reads and produce noise. It
is kept — as the spurious-fire guard.

Script: `bench/scripts/insert_window_survey.py` (Task 0.3). It exits 0, prints
the table above, and takes `--root/--read-len/--anchor`.

---

## Context and assumptions

1. **The mode exists as a specification, not as code.** `grep -rn
   "detect_pe_overlap\|fn complement\|pe_overlap" tools/fastprep/src/` returns
   nothing. The paper §`sec:adp-pe` (marked "not implemented" by `21d9055`) is
   the spec: compare the 3′ suffix of R1 with the reverse-complement prefix of
   R2; when the insert is shorter than the read the overlap reveals the insert
   size and the adapter boundaries.
2. **The harness already has a variant slot.** `bench/scripts/trim_report.py`
   runs `FASTPREP_VARIANTS='[[tag,[args...]],...]'` alongside the two published
   configs, appending `--head-dropout 0` to each. A PE mode is exactly this
   shape: a new argv.
3. **The harness already has the mate-aware metric.** `closed_loop.py` gained
   `--per-mate` / length-matched columns in `7631292`, so a PE variant can be
   scored per mate rather than only pooled.
4. **Ground truth is the constraint.** fastsim's insert ceiling is
   `2L - 12 = 288 bp` at 150 bp reads; above it a pair is censored and no
   overlap rule can apply.
5. **All measurements run on the desktop** (`raphael@100.104.77.4`,
   `~/FastSuite` is a read-only mirror; `/mnt/Fast/fastsuite-scratch` is the
   corpus). Laptop edits only; scp scripts to `/tmp`.
6. **The suite's gate is the definition of done** for any code change:
   `pixi run ci` → 8/8 count assertions, 0 failures, GRAND as asserted.

### What this experiment must NOT become

It must not become "implement the whole PE feature and hope". The mode is
opt-in and default-off in every step until a measurement says otherwise. If the
measurement says no, the code is not written at all and the decision is
recorded. That asymmetry is the point.

---

## Phase 0 — DONE (see the table above)

**Task 0.3 — commit the survey.** Move `insert_window_survey.py` into
`bench/scripts/`, add its TSV under `bench/results/phase-d/`, and keep the table
in this file. Written and run this session; only the commit remains.

---

## Phase 1 — instrument (no fastprep code, ~1.5 h)

**Task 1.1 — per-mate overtrim/undertrim columns.** Extend the trim-report
variant table with `overtrim_R1`, `overtrim_R2`, `undertrim_R1`, `undertrim_R2`
(counts and rates), keeping the pooled columns.
- **Files:** `bench/scripts/trim_report.py` (+ the per-domain aggregation).
- **Verify:** `overtrim_R1 + overtrim_R2 == pooled overtrim` for every row, and
  the two published configs' numbers must be **identical** to the committed TSVs.
  A strict no-regression check on the instrument.

**Task 1.2 — a readthrough/window classifier in the ground truth.** From the
`gt.tsv` per-pair records compute, per pair: `readthrough = (r1_adapter>0 or
r2_adapter>0)`, `in_window = (L <= insert <= 2L-12)`, `below_L`, `censored`. The
new columns must be sliceable by all four.
- **Verify:** counts reproduced by an independent one-liner; totals match.

**Task 1.3 — a "fires at all" counter.** The variant table needs `pe_trims` =
pairs where the PE path actually changed the trim decision. On Cell A this must
be **~0**; a non-zero number is the guard failing and is a finding, not a
metric.
- **Verify:** on the corpus `pe_trims` is zero or explained; on Cell B it is a
  sane fraction of readthrough pairs.

---

## Phase 2 — implement, opt-in, default off (~1 day)

**Task 2.1 — failing test first.** `tools/fastprep/tests/pe_overlap.rs`: given a
pair whose insert is shorter than the read and whose 3' end carries a known
adapter, the detector returns the readthrough length on the correct mate; and it
**refuses** when the mates cannot overlap, when the alignment is ambiguous, or
when no anchor is found. Include the chance-match case (fastcheck has a
same-class regression test — reuse its discipline).
- **Verify:** fails against today's code; paste the failure.

**Task 2.2 — detection only.** `detect_pe_overlap(r1, r2) -> Option<(mate, len)>`
in a new `tools/fastprep/src/pe_overlap.rs`. Anchored suffix/prefix alignment,
reject ambiguity, respect the 12-base anchor.
- **Verify:** `cargo test` green; `cargo check --all-targets` warning-free.

**Task 2.3 — wire behind a flag, default off.** `--pe-overlap`, falling back to
the current behaviour when detection returns nothing. Documented in `--help` in
the neighbouring style and in `docs/CLI-SPEC.md`.
- **Verify:** a pair with no overlap runs **byte-identically** with and without
  the flag (the critical no-regression test); a pair with overlap trims
  differently; the gate stays 8/8 with the new count asserted in the same commit.

**Task 2.4 — expose to the harness.** `FASTPREP_VARIANTS='[["fp-pe",["--pe-overlap"]]]'`
produces an `fp-pe` row; confirm.

---

## Phase 3 — the two cells (desktop, ~2 h)

**Cell A — the guard (corpus, 102 datasets).**
- Run the two published configs plus `fp-pe` over the corpus.
- **Pass condition: `pe_trims ≈ 0`.** Any meaningful number means the mode trims
  clean reads — the same class of bug as the 4-base-overlap false trim that
  `--min-adapter-overlap 8` exists to prevent (~10% of clean reads, per its own
  help text). A failure here ends the experiment regardless of Cell B.
- Also record exact rate and retention with vs without; they should be
  indistinguishable. If the mode moves them on clean data, that is the bug.

**Cell B — the value test (synthetic, with readthrough).** — DONE, null result.
- Build a cohort whose mates carry adapter readthrough, using fastsim with an
  adapter-bearing profile (the corpus's `adapter.mode: none` must be overridden —
  this is the one place the profile is changed, and it must be recorded in the
  results).
- Run with and without `fp-pe`; score **per mate** and sliced by `readthrough`.
- **Primary metric:** exact rate at tolerance 0 on readthrough pairs.
- **Secondary:** retention, undertrim (a win bought with biology must show here),
  and the length-matched `g_frac` gate for the metagenomics domain.

**Correction to the Cell B spec (2026-10-19, after the run).** The spec above asked
for inserts "concentrated in `L ≤ I ≤ 2L-12`" *and* readthrough. **That set is
empty**: readthrough needs `I < L − t` (the read runs past the fragment) while the
overlap window needs `I ≥ L`. Verified on the built cohort — every `rt=1` row has
`win=0` and every `win=1` row has `rt=0`. The mode's real operating region is the
*below-window* readthrough mates, and the cohort was built for that instead.

**Task 3.2 — the decision rule, fixed before the numbers are read.**

| Cell A | Cell B | decision |
|---|---|---|
| `pe_trims ≈ 0` | exact-rate gain on readthrough pairs ≥ T | **restore**, default-on, paper updated to shipped |
| `pe_trims ≈ 0` | gain < T, or undertrim rises | **do not restore**; record both cells in `OVERTRIM-PAIRED.md` and put the numbers in the paper's banner |
| `pe_trims` non-zero | (any) | **do not restore**; the mode trims clean reads. Fix the detector, then re-run |

T = 0.5 pp is the proposal; the user sets it. Cell A's guard is not thresholded:
any firing on adapter-free data is disqualifying.

---

## Phase 4 — record

- `docs/OVERTRIM-PAIRED.md` gains a "with/without" section: both cells, the
  corpus's `adapter.mode: none` finding, and the numbers.
- The paper: banner gains the numbers, or comes down if restored.
- One commit per phase; manifest updated in the same commit as any new test;
  `pixi run ci` on both machines; `sync-fastsuite.sh`.

---

## Files likely to change

| file | change |
|---|---|
| `tools/fastprep/src/pe_overlap.rs` | new — detection (Phase 2) |
| `tools/fastprep/src/adapter.rs` / `pipeline.rs` / `cli.rs` | wire the flag |
| `tools/fastprep/tests/pe_overlap.rs` | new — the failing-first tests |
| `bench/scripts/trim_report.py` | per-mate columns + readthrough slicing + `pe_trims` |
| `bench/scripts/insert_window_survey.py` | new — the Phase 0 survey (Task 0.3) |
| `docs/OVERTRIM-PAIRED.md` | the with/without result |
| `papers/fastprep/paper.typ` | banner number or removal |
| `docs/CLI-SPEC.md` | the new flag |
| `scripts/expected-test-counts.txt` | only if a test is added (same commit) |

## Risks, open questions

- **Cell B is synthetic, so its result is conditional on the profile used.** The
  mode's measured value is "on a library with readthrough in the overlap window",
  not "on real Illumina data". The corpus cannot answer the latter because it has
  no adapter — state this limit in the write-up.
- **Readthrough is not adapter identity.** The mode finds where the fragment
  ends, not what the adapter was. On Cell B, measure how often readthrough
  implies adapter, and bound the claim with it.
- **Ambiguity rejection is the failure mode.** A permissive detector trims
  biology and looks like a win. Cell A is the cheap early warning; undertrim in
  Cell B is the expensive one.
- **Time.** Phase 2 is the only multi-hour block and is gated on Cell A's guard
  being defined first. Cheaper first cut: Phase 0 + 1 + a probe-only build
  (detection exposed, not wired) answers "does it fire correctly" for a fraction
  of the cost; the wiring decision then follows the numbers.

## Order of operations

Task 0.3 → Phase 1 → **write the decision rule** → Phase 2 → Phase 3 (Cell A
first; Cell B only if A passes) → Phase 4.
