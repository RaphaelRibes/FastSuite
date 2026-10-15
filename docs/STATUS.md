# FastSuite status

Verified 2026-10-15 against the current `main` (commit `bff40c3`, 94 commits
ahead of `origin/main`; **nothing pushed**, per standing instruction). Every
number here was executed this session: the six cargo test suites, the
pathology-corpus integration test, the Phase D closed loop re-read from its
TSV, and `pixi run doctor`. **The suite is WIP: nothing here is
production-ready** (see the README banner).

## The tools

| tool | version | tests | state | next plan |
|---|---|---|---|---|
| fastdedup | 1.3.0 | 104 | shipped (published: bioRxiv 2026.04.29.721745) — **D-4 stdin/stdout streaming landed**: `-` for `-1/-2/-o/-p` with `--out-format`/`--out-compression`, and the documented refusal cases (two inputs as `-`, `--json -` with `-o -`, …) | — |
| fastprep | 0.1.0 | 128 | working — swept-winner default, `--force`/`--dry-run`/`-q`, emits `fastsuite.report.v1` with the raw stats under `legacy` | — |
| fastsim | 0.1.0 | 59 | working — measured insert CDF, SE mode, **clap-migrated** | — |
| fastmaps | 0.1.0 | 43 | **CLI-converged** — `--json` + v1 wrapper, `-j`, `--force`, `-q`, `FASTSUITE_THREADS`; `--corners` and `--preview-reads` | — |
| fastcheck | 0.1.0 | 107 | working — `--max-mem` + admission, insert measurement fixed (Task 10), **profile fitter fixed** (`430117a`: `BASE_POSITIONS` 60 → 512), harness green | `plans/2026-10-12_fastcheck-exactness.md` (tasks 5–9, 11–14) |
| fasthelper | 0.1.0 | 21 | working — fastcheck step, `check -c`, `--fail-on` | `plans/2026-10-12_fasthelper-hardening.md` — already done on disk: `run`, `compare`, `check -c`, `doctor`, `schema`, `--fail-on` (exit 2); the plan's open items are the *hardening* ones (5–8: `-P/--jobs`, resume/`--from`, `install --prefix`) |

**462 cargo tests, 0 failures**, plus the **pathology corpus: 9 assertions
passed, 0 failed** (`bash bench/tests/pathology_test.sh`, wired into
`scripts/test.sh`) — 471 assertions in total. Note that `scripts/test.sh`
counts the pathology suite as one aggregate line, so its own `GRAND` reads
463; the 9 is the number of individual assertions inside it. `pixi run
doctor` lists all six binaries; `bash examples/demo.sh 20000` runs the whole
chain end to end (re-run this week, report written to
`examples/runs/demo/report.html`).

## What works today (executed, not claimed)

* The canonical CLI option set: `-1/-2 -o/-p -j/--threads --gz-level --json
  --force --dry-run -q -v -V` — with the per-tool gaps closed this week
  (fastprep `--force/--dry-run/-q`, fastmaps `--json`/`--force`/`-q`,
  fastsim `--force`/`-q`, `FASTSUITE_THREADS` in fastmaps + fastsim).
* `fastsuite.report.v1` from every tool, **fastprep included** — verified by
  running the built binary in both its single and paired modes and reading the
  schema field. fastprep's `--dry-run` still correctly writes nothing, and the
  `--auto-params` path is a refusal stub today (no model trained yet), so there
  is no third report path to cover.
  * **Stale string found this session:** that `--auto-params` refusal text
    still claims "the shipped default IS the swept winner on the current
    corpus (90/90 datasets, score median 0.969)" — the source at
    `tools/fastprep/src/cli.rs:531`. `docs/AUTO-PARAMS.md` has the corrected
    figure (**1 of 101**, median winner score 0.9643) and says so explicitly.
    The string is user-facing and wrong; fixing it is a one-line edit.
* The closed loop in Rust: `fastcheck --insert-hist --profile p.json <fq>`
  writes the profile `fastsim --profile p.json` consumes — no Python fit
  step.
* `fasthelper`: `run` (chain + manifest + intermediate cleanup),
  `concatenate` (both the report and the `--positions` matrix modes),
  `compare`, `check -c`, `doctor`, `schema`, `--fail-on` (exit 2).
* The Python differential harness in fastcheck is **green** (it was red: it
  read the positions TSV by position; now by column name).
* The Phase D pathology corpus: five committed adversarial fixtures
  (`bench/fixtures/`, commit `aa955ef`) with a deterministic generator
  (`make_fixtures.py`) and a cross-binary integration test. See the caveat in
  "Known defects" — one of the five is a regression guard, not a
  reproduction.

## Measured this week (the benches behind the papers)

| claim | measurement |
|---|---|
| fastcheck memory is flat in dataset size | 133–140 MB from 100k to **20M pairs** (2.7 GB of reads) vs FastQC's 310–363 MB |
| fastmaps vs FastQC (5M pairs) | 9.17 s / 145 MB vs 16.1 s / 16 MB — FastQC is leaner and this is not a memory win; the honest claim is the composition plane plus thread scaling, not raw speed or RSS |
| fastprep vs fastp (5M pairs, 1 core, gzip level 6 both) | 34.37 s vs 28.50 s — fastp LEADS this cell; fastprep leads at 100k/500k/1M and on 2 cores (25.84 vs 39.61 s) |
| fastprep's shipped default is correct | 7.7–7.9% drops vs the 4.7–4.8% injected must-drop + polyG |
| fasthelper vs MultiQC | zero orchestration overhead (0.667 s vs 0.647 s) + manifest provenance; MultiQC v1.35 emits **no output** on `fastsuite.report.v1` (no module for it) — 0 report files at every N from 25 to 800 samples |
| fastsim ground-truth recovery | 8/10 parameters exact; the two failures are fastcheck defects (below) |
| fastsim insert resolution | 94–95% measured exactly at 200–288 bp (18 982 of 20 000 at insert 200, p50 200); the 288 bp ceiling is geometry, not a defect |
| trim accuracy vs simulator ground truth | **fastp-default wins** — 92.2% exact at tolerance 0 vs fastprep's 66.5%; 99.6% vs 98.4% at tolerance 2; see `docs/TRIM-ACCURACY.md` |

## Phase D — the closed loop is complete

The campaign ran the full corpus: **102 real ENA datasets × 55 configs ×
3 variants = 306 comparisons** (re-read this week from
`bench/results/phase-d/closed_loop.tsv`; the aggregation and campaign
scripts are under `bench/scripts/`, the fixtures under `bench/fixtures/`).

**PASS 167/306 = 54.6%**, 52 datasets passing every variant.

| domain | pass rate |
|---|---|
| human_gut | 92.2% (47/51) |
| bacterial_isolate | 76.5% (39/51) |
| generalised | 64.7% (33/51) |
| virome | 47.1% (24/51) |
| rna_seq | 41.2% (21/51) |
| metagenomics | **5.9% (3/51)** |

The pass rate moved through four bug fixes, and **not monotonically** — the
single most important result in the campaign:

| stage | pass rate |
|---|---|
| as first measured | 95/306 = 31.0% |
| after the 3′ boundary-rule fix alone | 86/306 = 28.1% (the fix *lowered* it — it exposed a GC under-model that had been partially cancelling) |
| after the merged dup_rate fix | 113/306 = 36.9% |
| after the duplicate-replay fix | 131/306 = 42.8% |
| after head-G deconvolution | **167/306 = 54.6%** |

Median residual over the final run: `g_frac_max_diff` 0.0140, `mean_q_diff`
0.1503, `gc_diff` 0.0002, `dup_ratio` 1.25, `polyx5_diff` 0.0215.

**The fitter is fixed**, and the pilot is long past — the "gated on the fitter"
framing in earlier revisions of this document no longer applies. The profile
fitter's `BASE_POSITIONS` was raised from 60 to 512 (commit `430117a`) so a
high-GC library's position-wise G ramp is no longer halved, alongside a
measured `dup_rate`; those were the two lossy terms behind the GATE-2 fitter
failures. fastsim's insert-size resolution went 0.5078 → 0.9491 at insert 200
as a result.

## Known defects (measured, not guessed)

1. **metagenomics closed-loop pass rate is 5.9% and has not moved across any
   of the four fixes.** Re-read at each commit in the
   chain — `efaf834` (first full-corpus run), `c3ad7a0` (+dup_rate),
   `0a19d06` (+duplicate replay), `bd31d9a` (+head-G deconvolution) — the
   count is **3/51 in all four**, while the other domains gained 25–60
   points over the same span. All 48 failures are `g_frac`, at a median
   over-tolerance error of 0.0529. Neither the dup model, nor the 3′
   boundary rule, nor the head-G deconvolution touched it: this is a
   **distinct, still-undiagnosed mechanism**, and it is the main remaining
   blocker to a publishable simulator accuracy number.
2. **fastcheck's polyG detector** — reports 93.98% tail-free against a true
   66.01%, zero mass at t=3/t=15 (true 7.99%/6.99%). Reported in the fastsim
   paper §5.3 as a measured limitation. **Not fixed.**
3. **fastcheck's worst 10-position quality block** — misses by 0.45 Phred
   (mean passes at 0.12). **Not fixed.** Together with item 2 these two are
   the whole of fastsim's 8/10 recovery score.
4. ~~**fastprep's `--json` was legacy raw stats** — every other tool emitted
   the v1 wrapper; `fasthelper` adapted it with a note.~~ **Fixed
   (PLAN.md §1.1):** fastprep now writes `fastsuite.report.v1` with the raw
   stats preserved under `legacy`; the legacy branch is still kept for
   pre-wrapper files. Verified by running the built binary in both its
   single and paired modes.
5. **The pathology corpus's k-mer fixture is a POST-FIX INVARIANT GUARD, not
   a reproduction.** The k-mer defect's 40-minute hang is **not reproducible
   by any small synthetic FASTQ** — verified directly: a faithful stream of
   1M distinct 12-mers (2M inserts, heavy count>1 retention) finishes in
   **< 0.2 s** even on the pre-fix every-insert eviction. The hang depends on
   real DRR509062's key distribution, which synthetic random DNA does not
   recreate. `kmer_saturation.fastq` saturates the 12-mer space (the necessary
   condition) and the test asserts the post-fix invariant — the scan
   completes within a generous bound and a k-mer spectrum is emitted. **A
   green line there does not mean the pathology was reproduced.**
6. ~~**pixitainer image build fails with `cargo: command not found` (exit
   127)**~~ **Retired 2026-10-15 — it builds.** The recipe, the committed
   `fastsuite.def` (with `keep-def = true`) and the pixi.toml
   `[tool.pixitainer]` block are in place, and two clean full builds were run
   this session (`--no-fakeroot` and the default fakeroot path), both with
   **exit 0** and an `✅ Success! Image built at: …` line. Inside each image
   all six binaries answer `--version` and `fasthelper schema` runs, confirming
   the `cargo: command not found` pitfall documented in `pixi.toml` (a bare
   `bash scripts/build.sh` outside an activated pixi env) is the one the
   recipe already avoids — the `%post` uses `pixi run --as-is build`, which
   activates the env itself. `pixi containerize -d` still emits a def whose
   only difference from the committed one is the version label. The
   RELEASES checklist item for the image is already ticked, and this session
   independently re-verified it end to end. Two minor drifts worth noting: the
   committed `fastsuite.def` carries `Pixitainer_Version 0.9.0` while the
   installed tool now reports 1.0.0, and `pixi.toml`'s `test` task comment
   still says "397 tests".

## Where each paper stands (local-only repos, never pushed)

All five papers are **written** and were sent to Discord on 2026-10-15. Page
counts measured from the built PDFs: FastCheck 17 p, FastPrep 24 p, FastMaps
14 p, FastHelper 16 p, FastSim 22 p. Each has its `-Bench` repo (committed
`results.csv`, seeded generator, VmHWM peak-RSS harness; reads gitignored).
The FastDedup paper is separate and older (28 p, bioRxiv
2026.04.29.721745).

**The FastPrep paper does not yet contain the trim-accuracy result.** That
section exists as an *uncommitted* addition in `sect-real-ena.typ`
(`@tbl-trim-acc`, 66.5% / 92.2% / 84.1% exact) — the edits are on disk but
not built into `paper.pdf` (mtime 2026-10-15 01:53) and not committed. It is
a separate in-flight task. Likewise, the Phase D closed-loop numbers
(167/306) are in `FastSim-Paper/benchmark/closed_loop.tsv` but are not yet
reflected in the papers' prose.

## What the suite refuses to do (D-12)

PE overlap/merge, UMI handling, alignment-based QC, workflow-engine
semantics (arrays, resume), per-tile modules (until built), long reads (D-7).

## Where the ML campaign stands

Io is gone (2026-10-12); the plan is `docs/PLAN-ML-LOCAL.md` — all benchmark
compute on the desktop (16 cores, 48 GB), laptop interactive-only. **The
pilot is done and the full campaign ran** — 306 comparisons over 102 real
ENA datasets, at 54.6% closed-loop pass. The sweep results stand independently:
`winners.tsv` (101 winning configs, 16 distinct of the 55-grid) and
`fastp-comparison.tsv` (303 rows at gz6, 4 threads) — fastprep-winner is
faster on 100/101 datasets (median 2.00×) and lighter on 101/101 (34.3 MB
median vs fastp's 1134.7 MB). The trim-accuracy benchmark
(`docs/TRIM-ACCURACY.md`, `bench/scripts/trim_accuracy.py`) is a separate
axis, and on that one fastp-default leads — stated plainly above. The free
A2/A3 tests exist as bench scripts in `~/work/FastSim-Bench`; porting them
into `tools/fastcheck/tests/` is fastcheck plan Task 13, not yet done.
