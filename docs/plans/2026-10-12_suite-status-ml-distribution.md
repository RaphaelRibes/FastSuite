# FastSuite suite-level plan: status, ML campaign, distribution

> **For Hermes:** Use subagent-driven-development skill to implement this plan task-by-task.

**Goal:** Close the suite-level gaps the per-tool plans do not cover: a
living `docs/STATUS.md`, the PLAN-ML execution chain (QC→profile fitter,
corpus, pilot, sweep, ranker, `fastprep --auto-params`), and the distribution
layer (per-tool versions, bioconda recipes, CITATION.cff, RELEASES.md) that
decisions D-2/D-3/D-9 require before anything is publishable.

**Architecture:** No shared crate (D-1). The suite stays a monorepo (D-10) for
now. Distribution is per-tool: one semver, one `meta.yaml`, one container
recipe each, one release train listed in `docs/RELEASES.md`.

**Tech Stack:** Rust (tools), Python (bench harness), pixi + bioconda +
pixitainer (distribution), DuckDB/polars (aggregation).

---

## Context (verified 2026-10-12)

* All six tools build and pass tests: fastcheck 54, fastdedup 87, fastprep 110,
  fastsim 54, fastmaps 27, fasthelper 13 (272 tests, 0 failures).
* `examples/demo.sh` runs the chain end-to-end (fastsim → fastdedup →
  fastprep → fastmaps → aggregate).
* The QC→sim closed loop works in Rust: `fastcheck --insert-hist --profile
  p.json <fq>` writes a profile with a measured insert CDF that `fastsim
  --profile p.json` consumes — no Python fit step exists or is needed.
* Bugs fixed in this sitting (committed): `fasthelper compare` did not create
  its output dir; `pixi run test` omitted fastcheck; `doctor`/manifest omitted
  fastcheck; the uncommitted fastmaps `-j/--threads` change is now committed.
* Stale docs: `PLAN-ML-STEPS.md` lists A1 (`qc_to_profile.py`) as next, but
  `fastcheck --profile` supersedes it; `PLAN-QC.md` §0 does not mention the
  rich-report modules that shipped 2026-09-18; `REPLACEMENT-MAP.md` still
  lists k-mer/over-representation as gaps that now exist in bounded form.

---

## Phase 0 — Make the status honest (docs only, ~1 h)

### Task 1: Write docs/STATUS.md

**Objective:** One file that answers "what is done, what is next, per tool",
generated from verified facts, not aspiration.

**Files:** Create `docs/STATUS.md`.

Content skeleton (fill each row from the per-tool plans' "Context" sections):

```markdown
# FastSuite status

Verified 2026-10-12 against commit <sha>. Every claim here was executed:
builds, tests, demo chain, and the closed-loop fastcheck→fastsim run.

| tool | version | tests | state | next plan |
|---|---|---|---|---|
| fastdedup | 1.3.0 | 87 | shipped | plans/2026-10-12_fastdedup-streaming.md |
| fastprep | 0.1.0 | 110 | working | plans/2026-10-12_fastprep-cli-convergence.md |
| fastsim | 0.1.0 | 54 | working | plans/2026-10-12_fastsim-cli-convergence.md |
| fastmaps | 0.1.0 | 27 | working | plans/2026-10-12_fastmaps-cli-convergence.md |
| fastcheck | 0.1.0 | 54 | v0.1 | plans/2026-10-12_fastcheck-exactness.md |
| fasthelper | 0.1.0 | 13 | new | plans/2026-10-12_fasthelper-hardening.md |

## What works today (executed, not claimed)
- ...

## Known gaps (with the plan that closes each)
- ...

## What the suite refuses to do (D-12)
- ...
```

**Step 1:** create the file with the table and the verified facts above.
**Step 2:** `git add docs/STATUS.md && git commit -m "docs: status file, generated from verified runs"`.
**Verify:** `wc -l docs/STATUS.md` (non-trivial), every "next plan" path exists
(`ls docs/plans/2026-10-12_*.md | wc -l` → 6+1).

### Task 2: Reconcile the plan documents

**Objective:** Stop the docs contradicting each other.

**Files:** `docs/PLAN-ML-STEPS.md`, `docs/PLAN-QC.md`, `docs/REPLACEMENT-MAP.md`,
`docs/PLAN.md`.

1. **PLAN-ML-STEPS.md**: mark A1 `[done]` — `fastcheck --profile` is the
   fitter, in Rust, in the tool. Replace the "next three actions" list with
   A2/A3 (recovery + round-trip tests), which are genuinely next. Keep the
   rest (B–F) but note A2/A3 are free and local.
2. **PLAN-QC.md** §0: add a paragraph recording what shipped after the plan
   was written (rich report, polyX runs, bounded repeats/k-mers, heatmap,
   `--profile`, before/after, 54 tests, the Python differential harness), and
   point at `docs/plans/2026-10-12_fastcheck-exactness.md` for what remains
   (§5.9 overrep modes, `--max-mem`, adapter-position, tiles, `--only`,
   ground-truth recovery).
3. **REPLACEMENT-MAP.md**: k-mer content and over-representation move from
   "⚠ gap" to "partial — bounded approximate counters shipped; exact modes
   planned (fastcheck plan)". Per-tile stays a gap. The honest summary line
   updates accordingly.
4. **PLAN.md**: phase 1 items that are done (fastmaps `-j`, fastprep `-q`
   still missing — check the fastprep plan's status before marking), and
   phase 2 items now covered by the fasthelper plan, get a `→ see
   docs/plans/...` pointer instead of being re-listed.

**Verify:** `grep -c "docs/plans/" docs/PLAN*.md docs/REPLACEMENT-MAP.md` ≥ 4;
no doc claims a gap that a shipped feature contradicts.

---

## Phase A — The ML campaign chain (the research plan, ~local first)

> PLAN-ML-STEPS.md is the authoritative breakdown. This phase exists in the
> suite plan only to sequence it against everything else.

### Task A1: the recovery test (A2 in PLAN-ML-STEPS) — LOCAL, FREE, FIRST

**Objective:** prove the fitter (`fastcheck --profile`) recovers the
parameters of a profile it did not generate from. This is the go/no-go for
every downstream step.

**Files:** Create `bench/tests/recovery.rs`-style verification — but note:
the fitter lives in fastcheck, so the natural home is
`tools/fastcheck/tests/` (Rust) using the synthetic profile that ships:
`tools/fastsim/profiles/synthetic_novaseq.json`.

**Step 1 (test):** generate a known dataset:
```
./tools/fastsim/target/release/fastsim --pairs 20000 \
  --profile tools/fastsim/profiles/synthetic_novaseq.json -o /tmp/known
```
**Step 2:** fit it back:
```
./tools/fastcheck/target/release/fastcheck --insert-hist \
  --profile /tmp/fitted.json /tmp/known_R1.fastq.gz /tmp/known_R2.fastq.gz
```
**Step 3:** assert, in a test, against PLAN-ML-STEPS tolerances: head CDF
within 3 % per level; quality curve within 0.5 Q; GC within 0.005; tail run
lengths within one level.
**Step 4:** if it fails, the fix is in the fitter — stop and fix, do not
proceed to B.
**Commit:** `feat(fastcheck): ground-truth profile recovery test`.

**Note:** `synthetic_novaseq.json` has NO `insert` key (verified), so the
insert CDF part of the round-trip is measured-but-not-recoverable from that
profile; either add an insert CDF to a NEW synthetic profile (allowed — it is
hand-written, not data, per D-11) or scope the test to the recoverable
fields. Decide in the test, document the choice.

### Task A2: the round-trip test (A3) — LOCAL, FREE

**Objective:** simulate from the fitted profile, re-QC, compare to the
original QC — tolerances: per-position G fraction within 0.02, mean Q within
0.5, GC within 0.01, duplicate rate within 2×.

**Files:** same test file, second test.

**Step 1:** `fastsim --profile /tmp/fitted.json --pairs 20000 -o /tmp/resim`
**Step 2:** `fastcheck --json /tmp/resim.json /tmp/resim_R*.fastq.gz`
**Step 3:** assert the tolerances.
**Commit:** `feat(fastcheck): fitted-profile round-trip test`.

### Task A3: corpus + staging (PLAN-ML-STEPS Phase B) — GATE 1

**Objective:** `bench/scripts/ena_query.py` builds `manifests/real.tsv`
(100 datasets per domain, every row with `fastq_ftp` + `fastq_md5`),
`subsample.sh` produces deterministic 100 k-pair units with md5 verification
and **no raw file left on scratch** (D-11: never ship data — the bench
directory holds manifests and code, never reads).

**Gate:** confirm the pilot subject (30 datasets: 10 generalised, 10 human
gut, 10 virome proposed) before running anything on Io.

**Note:** `deploy/io/corpus-manifest.tsv` already exists with 19 rows — use
it as the shape reference, extend per domain. `deploy/io/` already has
`stage-to-io.sh`, `pilot.sbatch`, `sweep-cell.sbatch`, `dispatch.sbatch` —
reuse, do not reinvent.

### Task A4: pilot (Phase C) — GATE 2, ~2.5 core-h on Io

Full chain on 30 datasets × 60 configs × 3 variants, closed-loop check on
real data (the go/no-go), fastp comparison, decision memo `bench/PILOT.md`.

**Gate 2:** if the closed loop fails, Phase D does not run — the fix is in
the fitter.

### Task A5: full run (Phase D) — GATE 3, ~153 core-h on Io

300 000 exact evaluations. Per the standing rule: partition `cpu-dedicated`,
accounts `dedicated-cpu@cirad-normal` (cells) / `@cirad-long` (dispatcher),
`mem=256M`, throttled at 100 cores, dispatcher 7 days — all four values
confirmed before submission.

### Task A6: labels + models (Phases E–F) — LOCAL

Label table, leak test, study-level splits, per-domain rankers, model cards,
and finally `fastprep --auto-params <qc.json> --print-args` (Rust tree-walking,
no new dependency) + `fasthelper recommend runs/`.

**Note:** the 90-dataset corpus campaign already shipped the swept winner as
fastprep's default (winner score median 0.969, 90/90 datasets), so
`--auto-params` has nothing to switch on **for the current corpus**. Phase E
is what would justify it on a *larger, more diverse* corpus — record that
honestly in `docs/AUTO-PARAMS.md` rather than implying the shipped default is
the model's output.

---

## Phase D — Distribution (D-2, D-3, D-9) — ~1 day, parallelisable

### Task D1: per-tool versioning + CITATION.cff

**Files:** each `tools/*/Cargo.toml` `version`, each `tools/*/CITATION.cff`.

* fastdedup keeps 1.3.0 (it is the released one; bioRxiv citation exists —
  mirror it into CITATION.cff: Ribes, Mandier, Baniel 2026, doi:10.64898/
  2026.04.29.721745, with the ISDM/ISEM affiliations).
* fastprep gets 0.1.0 → its own CITATION.cff crediting the Mandrillus Project
  (D-9: fastprep and fastdedup carry project/ISDM/ISEM + Mandier/Baniel
  credits; fastsim/fastmaps/fastcheck/fasthelper credit Raphaël Ribes alone).
* The other four: 0.1.0, single-author CITATION.cff.

**Verify:** `ls tools/*/CITATION.cff | wc -l` → 6; each parses as YAML.

### Task D2: license text per tool

D-9 leaves the license text open (MIT proposed; fastdedup already MIT).
The suite root already has an MIT `LICENSE` (Copyright 2026 Raphaël Ribes).

* Copy the MIT text into `tools/fastprep/LICENSE` with the **Mandrillus
  Project** attribution line added to the copyright holder.
* `tools/fastdedup/Licence` already exists (MIT) — rename decision is
  cosmetic; keep both spellings working or normalise to `LICENSE` with a
  changelog note.
* The other four tools: `LICENSE` (MIT, Ribes).

**Verify:** `for f in tools/*/LICENSE tools/*/Licence; do head -1 $f; done`
shows a copyright line in all six.

### Task D3: bioconda recipes + pixitainer

**Files:** `recipes/fastdedup/meta.yaml`, `recipes/fastprep/meta.yaml`, …
(one per tool), plus pixitainer recipes.

Each `meta.yaml`: license file, `run:`/`test:` sections invoking
`--version` and a tiny **synthetic** fixture run (D-11: no data, ever; use
`profiles/test_small.json` for fastsim). No downloads at build time.

**Verify:** `ls recipes/*/meta.yaml | wc -l` → 6; each references the tool's
LICENSE and a fixture that exists in the repo.

### Task D4: docs/RELEASES.md

One file listing which tool versions made up each "FastSuite N" release with
each tool's maturity (stable/beta) — per D-2.

---

## Phase B — Benchmarks and the public claim (PLAN phase 6) — 2 days on Io

1. Refresh/regenerate the fastsim profile from the mandrill data (real
   profiles live OUTSIDE the repo, D-11:
   `~/hermes_scratch/fastsuite-local-artifacts/profiles/`).
2. Head-to-head on Io: fastdedup/fastprep vs fastp, fastmaps/fastcheck vs
   FastQC/Falco, fasthelper vs MultiQC — same data, same threads, wall + peak
   RSS, three reps.
3. Accuracy on ground truth: retention, retention of must-drop reads,
   adapter false-positive rate (use `fastsim`'s `PREFIX.gt.tsv`).
4. `docs/BENCHMARKS.md` with the raw matrix and caveats; store under `bench/`.

**Do not claim** "replaces FastQC" before the fastcheck exactness plan lands,
and "replaces MultiQC" before the fasthelper hardening plan lands —
`REPLACEMENT-MAP.md` is written to survive scrutiny.

---

## Suggested execution order

1. Phase 0 (docs, ~1 h) — makes every later claim honest.
2. Per-tool plans in value order: fasthelper hardening (fixes + wiring) →
   fastprep CLI convergence (v1 report is the most-used gap) → fastcheck
   exactness → fastmaps → fastsim → fastdedup streaming.
3. Phase A2/A3 (free, local, gates everything in ML).
4. Phase D (distribution) — unblocks publishing, needs no cluster.
5. Phase A3–A6 + Phase B — cluster work, gated.

## Risks

| risk | mitigation |
|---|---|
| docs drift again after this plan | STATUS.md is dated and generated; the per-tool plans are the source of truth |
| ML campaign proceeds on a broken fitter | A2/A3 run first and are free; GATE 2 blocks Phase D |
| distribution claims outrun reality | recipes test on synthetic fixtures only; RELEASES.md records maturity |
| the 90-dataset winner is mistaken for a model | AUTO-PARAMS.md states the default is the swept winner, not model output |
