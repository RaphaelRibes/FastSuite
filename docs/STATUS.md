# FastSuite status

Verified 2026-10-13 against the current `main` (commit pushed to
https://github.com/RaphaelRibes/FastSuite). Every claim here was executed:
builds, tests, the demo chain, the differential harness and the closed loop
(`fastcheck --profile` → `fastsim --profile`). **The suite is WIP: nothing
here is production-ready** (see the README banner).

## The tools

| tool | version | tests | state | next plan |
|---|---|---|---|---|
| fastdedup | 1.3.0 | 87 | shipped (published: bioRxiv 2026.04.29.721745) | `plans/2026-10-12_fastdedup-streaming.md` |
| fastprep | 0.1.0 | 117 | working — swept-winner default, `--force`/`--dry-run`/`-q` | v1 report wrapper + streaming (round 2) |
| fastsim | 0.1.0 | 57 | working — measured insert CDF, SE mode | clap migration (`plans/2026-10-12_fastsim-cli-convergence.md`) |
| fastmaps | 0.1.0 | 43 | **CLI-converged** (`--json` + v1 wrapper, `-j`, `--force`, `-q`, `FASTSUITE_THREADS`) | — |
| fastcheck | 0.1.0 | 77 | working — `--max-mem` + admission, **insert measurement fixed** (Task 10), harness green | `plans/2026-10-12_fastcheck-exactness.md` (tasks 5–9, 11–14) |
| fasthelper | 0.1.0 | 16 | working — fastcheck step, `check -c`, `--fail-on` | `plans/2026-10-12_fasthelper-hardening.md` (tasks 5–8) |

**397 tests, 0 failures.** `pixi run doctor` lists all six; `bash
examples/demo.sh 20000` runs the whole chain (fastsim → fastdedup → fastprep
→ fastmaps → aggregate).

## What works today (executed, not claimed)

* The canonical CLI option set: `-1/-2 -o/-p -j/--threads --gz-level --json
  --force --dry-run -q -v -V` — with the per-tool gaps closed this week
  (fastprep `--force/--dry-run/-q`, fastmaps `--json`/`--force`/`-q`,
  fastsim `--force`/`-q`, `FASTSUITE_THREADS` in fastmaps + fastsim).
* `fastsuite.report.v1` from every tool except fastprep (still legacy raw
  stats — the one wrapper left).
* The closed loop in Rust: `fastcheck --insert-hist --profile p.json <fq>`
  writes the profile `fastsim --profile p.json` consumes — no Python fit
  step.
* `fasthelper`: `run` (chain + manifest + intermediate cleanup),
  `concatenate` (both the report and the `--positions` matrix modes),
  `compare`, `check -c`, `doctor`, `schema`, `--fail-on` (exit 2).
* The Python differential harness in fastcheck is **green** (it was red: it
  read the positions TSV by position; now by column name).

## Measured this week (the benches behind the papers)

| claim | measurement |
|---|---|
| fastcheck memory is flat in dataset size | 133–140 MB from 100k to **20M pairs** (2.7 GB of reads) vs FastQC's 341–363 MB |
| fastmaps vs FastQC (5M pairs) | 12.25 s / 143 MB vs 72.22 s / 375 MB |
| fastprep vs fastp (5M pairs, pinned 1 core) | 32.9 s vs 84.7 s (2.58×); defaults: 33.9 s / 102 MB vs 24.1 s / 1143 MB |
| fastprep's shipped default is correct | 7.9% drops vs the 4.8% injected must-drop + polyG (fastp drops 0.0% — its length filter is off) |
| fasthelper vs MultiQC | zero orchestration overhead + manifest provenance; MultiQC v1.35 emits **no output** on `fastsuite.report.v1` (no module for it) |
| fastsim ground-truth recovery | 7/10 parameters exact; the failures are fastcheck defects (below) |

## Known defects (measured, not guessed)

1. **fastcheck's insert measurement** — FIXED 2026-10-13 (Task 10, commit
   `98296e3`): the anchor sat in R2's polyG tail (only ~25% of measurable
   inserts found) and the identity over-counted by the artefact lengths. Now:
   anchor at offset t2, `insert = b1 + b2 - o`, <= 3 mismatches tolerated.
   Verified on ground truth: 94.4% measured (seed 99), p50 201 vs injected
   200; on the strict window 97.3%, mean 206.62 vs GT 206.83. Residual known
   quirk: ~3% of sub-151 inserts leak out of `censored_below` (the
   `o >= min(l1,l2)` gate), flagged for Task 13.
2. **fastcheck's polyG detector** — reports 93.98% tail-free against a true
   66.01%, zero mass at t=3/t=15 (true 7.99%/6.99%). Reported in the fastsim
   paper §5.3 as a measured limitation.
3. **fastprep's `--json` is the legacy raw stats object** — every other tool
   emits the v1 wrapper; `fasthelper` adapts it with a note.

## Where each paper stands (local-only repos, never pushed)

`~/work/FastCheck-Paper` (11 p), `FastPrep-Paper` (19 p), `FastMaps-Paper`
(14 p), `FastHelper-Paper` (16 p), `FastSim-Paper` (18 p) — each with its
`-Bench` repo (committed `results.csv`, seeded generator, VmHWM peak-RSS
harness; reads gitignored). The fastsim paper's page 14 is the one page not
yet vision-QA'd.

## What the suite refuses to do (D-12)

PE overlap/merge, UMI handling, alignment-based QC, workflow-engine
semantics (arrays, resume), per-tile modules (until built), long reads (D-7).

## Where the ML campaign stands

Io is gone (2026-10-12); the plan is `docs/PLAN-ML-LOCAL.md` — all benchmark
compute on the desktop (16 cores, 48 GB), laptop interactive-only. The
pilot (~3 h) is the next real step, **gated on the fitter**: the recovery
run above is currently "fix the fitter", not "scale". The free A2/A3 tests
exist as bench scripts in `~/work/FastSim-Bench`; porting them into
`tools/fastcheck/tests/` is fastcheck plan Task 13.
