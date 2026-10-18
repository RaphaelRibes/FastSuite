# PE overlap detection — with/without test plan

> **For Hermes:** execute with subagent-driven-development, one task per
> subagent, parent verifies every number before it is committed.

**Goal:** decide, on measurement rather than on the paper's claim, whether
restoring fastprep's paired-end overlap detection earns its place — and if it
does not, record that decision where the paper and the code agree.

**Architecture:** PE overlap is the *classical* adapter-trimming mode (mate
readthrough reveals the insert, hence the readthrough ends) and the *designed
headline* paired-end mode in fastprep's paper. It was removed in the CLI
unification (`90cb16c`) and the design is retained in the paper +
`docs/OVERTRIM-PAIRED.md`. A previous measurement (`OVERTRIM-PAIRED.md`) tested
the *disambiguation* question — can the mate resolve the artefact/biology G-run
ambiguity — and found the rule worth **~1% of overtrim events**. That is **not**
the same question as this one: does PE overlap detection improve trim accuracy
overall? This plan answers the second, which the 1% figure never addressed.

**Tech Stack:** Rust (fastprep CLI, `--json`, PXMT), Python bench harnesses
(`bench/scripts/closed_loop.py`, `trim_report.py` with `FASTPREP_VARIANTS`),
`fastsim` for ground truth, the desktop corpus at `/mnt/Fast/fastsuite-scratch`.

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
   scored per mate rather than only pooled. This matters because the mode's
   benefit should show at the mate that readthrough happened in.
4. **Ground truth is the constraint.** fastsim's insert ceiling is
   `2L - 12 = 288 bp` at 150 bp reads; above it a pair is censored and no
   overlap rule can apply. A with/without comparison is only meaningful on
   inserts **below** that ceiling, so the insert distribution is a first-class
   design variable of this experiment, not a detail.
5. **All measurements run on the desktop** (`raphelm@100.104.77.4`,
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

## Phase 0 — decide the insert distribution (no code, ~30 min)

The single biggest risk is measuring on a library where the mode cannot apply.

**Task 0.1 — measure the corpus's insert distribution.**
Use the existing closed loop's per-pair insert arithmetic (`insert <= 2L-12`
for measurability). Produce: the fraction of pairs with `I < L` (readthrough
possible at all), `L <= I <= 2L-12` (overlap measurable, the mode's window),
and `I > 2L-12` (censored, unmeasurable).
- **Verify:** a table over the 102-dataset corpus with the three buckets, plus
  the same buckets for the two profiles used for synthesis.
- **Exit criterion:** if the `L <= I <= 2L-12` bucket is tiny, the with/without
  test must be run on a **synthetic** library built for the window, not on the
  corpus — say so in the plan file and proceed to Phase 1 with a synthetic
  cohort.

**Task 0.2 — pick the two comparison libraries.** One "real" (corpus-derived
profiles) and one "in-window" (synthetic, inserts concentrated in the
measurable overlap window). Both must be reproducible from a command written
down in the plan file.

---

## Phase 1 — instrument the mate-aware comparison (no fastprep code, ~1 h)

**Task 1.1 — add a per-mate exact/over/undertrim breakdown to the trim report.**
The existing report gives aggregate exact/overtrim counts. The PE mode's effect
should be visible *per mate*, and specifically on the mate whose readthrough
was corrected. Extend the variant table with `overtrim_R1`, `overtrim_R2`,
`undertrim_R1`, `undertrim_R2` (counts and rates), and keep the pooled columns
for continuity.
- **Files:** `bench/scripts/trim_report.py` (and its per-domain aggregation if
  the columns must appear there too).
- **Verify:** rerun the existing two published configs and confirm the new
  columns are internally consistent — `overtrim_R1 + overtrim_R2` must equal the
  pooled overtrim count for every row, and the baseline numbers must be
  **identical** to the committed TSVs before the change. This is a strict
  no-regression check on the instrument, not on fastprep.

**Task 1.2 — add a "readthrough mate" classifier to the ground truth.**
From the PXMT/`gt.tsv` per-pair records, compute whether each mate ended in
adapter readthrough (`r1_adapter > 0` / `r2_adapter > 0`) and whether the pair
was in the overlap window. The report's new columns must be sliceable by this
classifier, because "PE overlap helps exactly the readthrough pairs" is the
hypothesis under test.
- **Verify:** the classifier's counts on one dataset are reproducible by an
  independent one-liner (a second script, or a pandas one-liner), and the
  totals match.

---

## Phase 2 — implement the PE mode (Rust, ~1 day, gated on a failing test)

**Task 2.1 — write the failing test first.**
A unit test on the overlap detector: given an R1/R2 pair whose insert is
shorter than the read length and whose 3′ end carries a known adapter, the
detector must return the readthrough length on the correct mate and must
**refuse** (return nothing) when the mates are too long to overlap, when the
overlap is ambiguous, or when the anchor is not found.
- **Files:** new `tools/fastprep/tests/pe_overlap.rs`; detection logic in
  `tools/fastprep/src/adapter.rs` (or a new `pe_overlap.rs` if it keeps
  `adapter.rs` focused — the plan's judgement, not a guess).
- **Verify:** the test fails to compile/pass against today's code. Paste the
  failure.

**Task 2.2 — implement detection only.**
`detect_pe_overlap(r1, r2) -> Option<(mate, readthrough_len)>`. Constraints
from the paper's spec: anchored at the suffix/prefix overlap, must reject an
ambiguous alignment (multiple admissible overlaps), must respect the 12-base
anchor used elsewhere in the suite, must not be fooled by a chance match
(there is an existing chance-hit regression test in fastcheck for the same
class of bug — reuse its discipline).
- **Verify:** `cargo test` green; `cargo check --all-targets` warning-free;
  the new test passes.

**Task 2.3 — wire it behind a flag, default off.**
`--pe-overlap` (or the CLI-SPEC-consistent spelling) that makes the pipeline
use the detected readthrough instead of internal matching **only** when it
returns a result; otherwise fall back to the current behaviour. The flag must
be documented in `--help` in the neighbouring style, and must appear in
`docs/CLI-SPEC.md`'s flag table.
- **Verify:** a pair with no overlap runs byte-identically with and without the
  flag (this is the critical no-regression test); a pair with overlap trims
  differently; the gate stays at 8/8 with the new test count asserted in the
  same commit.

**Task 2.4 — expose it to the harness.** No harness change needed: the variant
slot in `trim_report.py` takes argv. Add
`FASTPREP_VARIANTS='[["fp-pe",["--pe-overlap"]]]'` and confirm the harness
produces a `fp-pe` row.

---

## Phase 3 — measure with and without (desktop, ~2 h)

**Task 3.1 — the corpus arm.** Run the trim report over the 102-dataset corpus
with the two published configs plus the PE variant, on both the real and
in-window libraries from Phase 0.
- **Primary metric:** exact rate at tolerance 0, **per mate**, sliced by the
  readthrough classifier.
- **Secondary:** retention (kept %), overtrim, and the length-matched
  `g_frac` gate for the metagenomics domain.
- **Verify:** every row reproduced twice; the two published configs' numbers
  match the committed TSVs (proving the harness change did not move the
  baseline).

**Task 3.2 — the decision rule, written down before looking at the numbers.**
Pre-register it so the measurement cannot be rationalised afterwards:

| result | decision |
|---|---|
| PE variant raises exact rate by ≥ 0.5 pp on the corpus arm, with retention not worse | restore the mode, keep it default-on, update the paper to shipped |
| PE variant raises exact rate by < 0.5 pp, or costs retention | **do not restore**; record the measurement in `docs/OVERTRIM-PAIRED.md` and mark the paper's banner with the number |
| PE variant changes exact rate only on readthrough pairs and by < 0.5 pp overall | do not restore; record that the mode is real but immaterial, with the slice |

The 0.5 pp threshold is a proposal; the user sets it. What matters is that it
is fixed **before** the numbers are read.

**Task 3.3 — the closed-loop arm (optional but decisive).** If the corpus arm
says "restore", run the closed loop with the PE variant on the metagenomics
set to see whether the mate-aware gain reaches `g_frac`. A trim-accuracy win
that does not reach the closed loop is still a win for fastprep, but the two
must be reported separately.

---

## Phase 4 — record the outcome

**Task 4.1 — `docs/OVERTRIM-PAIRED.md` gains the with/without section.** The
existing doc answers the disambiguation question; this experiment answers the
accuracy question. Both belong there, with their numbers and their
distributions.
**Task 4.2 — the paper.** If the decision is "do not restore", the existing
banner gains the measured number in one sentence. If "restore", the banner
comes down and the section's status note is updated to shipped, with the
flag named.
**Task 4.3 — commit.** One commit per phase, manifest updated in the same
commit as any new test, `pixi run ci` on both machines, sync via
`sync-fastsuite.sh`.

---

## Files likely to change

| file | change |
|---|---|
| `tools/fastprep/src/pe_overlap.rs` | new — detection (Phase 2) |
| `tools/fastprep/src/adapter.rs` / `pipeline.rs` / `cli.rs` | wire the flag |
| `tools/fastprep/tests/pe_overlap.rs` | new — the failing-first tests |
| `bench/scripts/trim_report.py` | per-mate columns + readthrough slicing |
| `docs/OVERTRIM-PAIRED.md` | the with/without result |
| `tools/fastprep/paper/paper.typ` | banner number or removal |
| `docs/CLI-SPEC.md` | the new flag |
| `scripts/expected-test-counts.txt` | only if a test is added (same commit) |

## Risks and open questions

- **The 1% figure may generalise and mislead.** It was measured for the
  disambiguation question. Treating it as the answer to the accuracy question
  would be the same mistake as trusting the paper's prose. This plan exists to
  avoid that.
- **Readthrough detection is not the same as adapter identification.** PE
  overlap finds where the fragment ends; it does not name the adapter. If the
  corpus's readthrough is adapter-dominated, the two coincide and the test is
  clean. If it is not, the mode's benefit is bounded by how often readthrough
  implies adapter — measure that fraction explicitly.
- **Ambiguity rejection is where this mode can go wrong.** A permissive
  detector will trim biology and look like a win on reads where it is wrong.
  The test in Task 2.1 must include the ambiguous cases, and the report must
  show *undertrim* alongside overtrim so a win bought with biology is visible.
- **Time.** Phase 2 is the only multi-hour block, and it is gated on Phase 3's
  decision rule being written first. If the user wants a cheaper first cut:
  Phase 0 + Phase 1 + a **probe-only** implementation (detection exposed but
  not wired to trimming) answers "does it even fire correctly on this corpus"
  for a fraction of the cost, and the wiring decision then follows the numbers.

## Order of operations

Phase 0 → Phase 1 → (write the decision rule) → Phase 2 → Phase 3 → Phase 4.
Phase 2 is the only one that touches fastprep's behaviour, and it starts only
after the rule in Task 3.2 is fixed.
