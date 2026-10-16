# FastSuite status

Verified 2026-10-16 against the current `main` (commit `50e8373`, 107 commits
ahead of `origin/main`; **nothing pushed**, per standing instruction). Every
number here was verified this session, either by running it or by re-reading the
cited artifact: the six cargo test suites, `scripts/test.sh` (pathology corpus
included), the Phase D closed loop re-read from its TSV, and the trim-accuracy
and `head_dropout` tables re-aggregated from `bench/results/phase-d/`. **The
suite is WIP: nothing here is production-ready** (see the README banner).

## The tools

| tool | version | tests | state | next plan |
|---|---|---|---|---|
| fastdedup | 1.3.0 | 104 | shipped (published: bioRxiv 2026.04.29.721745) — **D-4 stdin/stdout streaming landed**: `-` for `-1/-2/-o/-p` with `--out-format`/`--out-compression`, and the documented refusal cases (two inputs as `-`, `--json -` with `-o -`, …) | — |
| fastprep | 0.2.0 | 131 | working — `--force`/`--dry-run`/`-q`, emits `fastsuite.report.v1` with the raw stats under `legacy`. Shipped `--head-dropout` default is now **0** (was 4 in 0.1.0; `3d76a31`, version bumped to **0.2.0** — a behaviour change, so a minor bump) | — |
| fastsim | 0.1.0 | 59 | working — measured insert CDF, SE mode, **clap-migrated** | — |
| fastmaps | 0.1.0 | 43 | **CLI-converged** — `--json` + v1 wrapper, `-j`, `--force`, `-q`, `FASTSUITE_THREADS`; `--corners` and `--preview-reads` | — |
| fastcheck | 0.1.0 | 107 | working — `--max-mem` + admission, insert measurement fixed (Task 10), **profile fitter fixed** (`430117a`: `BASE_POSITIONS` 60 → 512), harness green | `plans/2026-10-12_fastcheck-exactness.md` (tasks 5–9, 11–14) |
| fasthelper | 0.1.0 | 21 | working — fastcheck step, `check -c`, `--fail-on` | `plans/2026-10-12_fasthelper-hardening.md` — already done on disk: `run`, `compare`, `check -c`, `doctor`, `schema`, `--fail-on` (exit 2); the plan's open items are the *hardening* ones (5–8: `-P/--jobs`, resume/`--from`, `install --prefix`) |

**465 cargo tests, 0 failures**, plus the **pathology corpus: 9 assertions
passed, 0 failed** — **474 assertions in total**. Note that `scripts/test.sh`
counts the pathology suite as one aggregate line, so its own `GRAND` reads
**466** (465 + the one aggregate line); the 9 is the number of individual
assertions inside it. `pixi run doctor` lists all six binaries; `bash
examples/demo.sh 20000` runs the whole chain end to end (report on disk at
`examples/runs/demo/report.html`, mtime 2026-10-15 21:40).

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
  * **Stale string found:** the `--auto-params` refusal doc-comment in
    `tools/fastprep/src/cli.rs` (`mod` ~line 139) still asserts the shipped
    default "is already the swept winner on the current corpus". It is wrong:
    `docs/AUTO-PARAMS.md` records that on Phase D's 101 real libraries the
    shipped default is the per-dataset winner on **1 of 101** (the 90/90,
    median-0.969 figure is from the older 90-dataset corpus). The string is a
    doc-comment, not a runtime message, but it is the source of the claim;
    fixing it is a one-line edit.
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
| fastprep's drop rate tracks the injected must-drop + polyG | 7.7–7.9% drops vs the 4.7–4.8% injected (a retention check, not a trim-accuracy claim) — carried from the paper campaign, not re-measured this session |
| fasthelper vs MultiQC | zero orchestration overhead (0.667 s vs 0.647 s) + manifest provenance; MultiQC v1.35 emits **no output** on `fastsuite.report.v1` (no module for it) — 0 report files at every N from 25 to 800 samples |
| fastsim ground-truth recovery | 8/10 parameters exact; the two failures are fastcheck defects (below) |
| fastsim insert resolution | 94–95% measured exactly at 200–288 bp (18 982 of 20 000 at insert 200, p50 200); the 288 bp ceiling is geometry, not a defect |
| trim accuracy vs simulator ground truth, tolerance 0 (101 real libraries, one harness) | **fastprep with `--head-dropout 0` leads**: 92.6% exact / 0.1% overtrim / 7.3% undertrim. fastp-default 92.2 / 0.2 / 7.6; fastp-tuned 84.1 / 8.6 / 7.3; **fastprep's shipped default (`--head-dropout 4`) is worst at 75.4 / 23.6 / 1.0**. Tolerance 2: 99.7 / 99.6 / 91.2 / 98.7. Aggregate of `bench/results/phase-d/overtrim-vs-fastp.tsv`; see `docs/TRIM-ACCURACY.md` |
| `head_dropout` axis (1515 runs, 3030 rows = 101 datasets × 3 variants × hd {0,1,2,3,4}) | `0` is the per-dataset winner on **100 of 101** datasets, `1` on the remaining one, `4` on **none**. Corpus total error at tolerance 0: **7.44% at `0` vs 24.64% at `4` (3.3×)**. `docs/OVERTRIM-SWEEP.md`, `bench/results/phase-d/head-dropout-axis.tsv` |
| drop asymmetry (the caveat under the trim result) | fastprep drops more reads than fastp at tolerance 0 (1.35 M vs 0.39 M), unevenly: virome **15.2% vs 0.0%**, RNA-seq 10.8% vs 2.1%, while on bacterial_isolate fastprep drops *fewer* (3.8% vs 5.1%). So the largest fastprep accuracy win (virome, 93.8% vs 92.1%) sits on a kept-set roughly one sixth smaller and is not apples-to-apples; the bacterial_isolate tie is. `docs/TRIM-ACCURACY.md` §Drop asymmetry |

## Phase D — the closed loop is complete

The campaign ran the full corpus: **102 real ENA datasets × 55 configs ×
3 variants = 306 comparisons** (re-read this week from
`bench/results/phase-d/closed_loop.tsv`; the aggregation and campaign
scripts are under `bench/scripts/`, the fixtures under `bench/fixtures/`).

**PASS 167/306 = 54.6%**, 52 datasets passing every variant.

The closed loop scores fastsim against each library's **fastcheck**-measured QC
profile. It is **fastcheck-only** and independent of fastprep's trimming and of
the sweep winners — the `head_dropout` work above does not move it.

| domain | pass rate |
|---|---|
| human_gut | 92.2% (47/51) |
| bacterial_isolate | 76.5% (39/51) |
| generalised | 64.7% (33/51) |
| virome | 47.1% (24/51) |
| rna_seq | 41.2% (21/51) |
| metagenomics | **5.9% (3/51)** |

The pass rate moved through four bug fixes, and **not monotonically** — the
single most important result in the campaign (stage counts as recorded during
the campaign; the final 167/306 and the domain split above are re-read from
`closed_loop.tsv`):

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

New and open (not defects of a measured component, but unresolved):

* **The trim-accuracy win is not apples-to-apples on every domain** — see the
  drop-asymmetry row above: fastprep drops more reads than fastp on RNA-seq and
  virome, and those are the domains where its accuracy lead is measured on the
  smaller kept-set. This bounds how the result may be stated, not whether it is
  real.
* **The `head_dropout` default is now 0** (committed `3d76a31`), which is the
  measured right default on this corpus (100/101 datasets; `1` on the remaining
  one; `4` on none). No version bump, so this is committed-but-unreleased: a
  binary built from `main` trims at 0, a released 0.1.0 tarball still trims
  at 4.

## Where each paper stands (local-only repos, never pushed)

All five papers are **written** and were sent to Discord on 2026-10-15. Page
counts measured from the built PDFs (as of 2026-10-15, not re-measured): FastCheck
17 p, FastPrep 24 p, FastMaps 14 p, FastHelper 16 p, FastSim 22 p. Each has its
`-Bench` repo (committed `results.csv`, seeded generator, VmHWM peak-RSS
harness; reads gitignored). The FastDedup paper is separate and older (28 p,
bioRxiv 2026.04.29.721745).

**The FastPrep paper now contains the reversed Section 5.6, and it is
committed.** `sect-real-ena.typ`'s trim-accuracy section (`@tbl-trim-acc`) was
rewritten to conclude that the *shipped default* was the problem, not the
algorithm: fastprep at `--head-dropout 4` is the worst of the four
configurations (75.4% exact at tolerance 0) and fastprep with the rule off
(92.6%) matches or beats fastp-default (92.2%). The reversal is commits
`b408bae` (reverse Section 5.6), `c089709` (drop asymmetry) and `f4a5425`
(corrected numbers, the `--head-dropout 0` default, a stale closed-loop figure
and a stale 90-of-90 claim) in the clean FastPrep-Paper repo. **`paper.pdf` is
rebuilt** — 28 pages, mtime 2026-10-16 14:35 — so the updated prose is in the
built PDF.
Separately, the Phase D closed-loop numbers (167/306) are in
`FastSim-Paper/benchmark/closed_loop.tsv` but are not yet reflected in the
papers' prose.

**Resolved (was flagged as a paper-side inconsistency):** the FastPrep paper's
§5.6 lead-in (`sect-real-ena.typ`) cited a stale closed-loop figure,
"94-of-305", together with the wrong dominant failure gate (it named
duplicate-ratio; the TSV shows the per-position G-fraction at 120 of 139
failing rows, against duplicate-ratio 18). Both were corrected in `f4a5425`.
The authoritative figure is `bench/results/phase-d/closed_loop.tsv` —
**167/306 = 54.6%**.

## What the suite refuses to do (D-12)

PE overlap/merge, UMI handling, alignment-based QC, workflow-engine
semantics (arrays, resume), per-tile modules (until built), long reads (D-7).

## Where the ML campaign stands

Io is gone (2026-10-12); the plan is `docs/PLAN-ML-LOCAL.md` — all benchmark
compute on the desktop (16 cores, 48 GB), laptop interactive-only. **The
pilot is done and the full campaign ran** — 306 comparisons over 102 real
ENA datasets, at 54.6% closed-loop pass (fastcheck-only; see above). The sweep
results stand independently: `winners.tsv` (101 winning configs, 16 distinct of
the 55-grid) and `fastp-comparison.tsv` (303 rows at gz6, 4 threads) —
fastprep-winner is faster on 100/101 datasets (median 2.00×) and lighter on
101/101 (34.3 MB median vs fastp's 1134.7 MB). The trim-accuracy benchmark
(`docs/TRIM-ACCURACY.md`, `bench/scripts/overtrim_vs_fastp.py`) is a separate
axis: there fastprep's **shipped** default is the worst of the four
configurations and fastprep with `--head-dropout 0` matches or edges
fastp-default. The earlier "fastp-default leads, 92.2% vs 66.5%" reading is
**superseded** — the 66.5% figure scored pre-regeneration sweep PXMT against a
regenerated simulator whose RNG stream had shifted, so it is not a fastprep
measurement at all (root cause: `docs/TRIM-ACCURACY-DISCREPANCY.md`; the
`fastprep` rows of `bench/results/phase-d/trim-accuracy.tsv` must not be cited —
its `fastp` rows remain valid). The free A2/A3 tests exist as bench scripts in
`~/work/FastSim-Bench`; porting them into `tools/fastcheck/tests/` is fastcheck
plan Task 13, not yet done.
