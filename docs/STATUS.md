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
| fastsim | 0.1.0 | 63 | working — measured insert CDF, SE mode, **clap-migrated**, **profile-driven ambiguous-base (N) emission** (`n_run`: read rate, run-length histogram, edge shares; recorded in the GT as `r1_n`/`r2_n`). This was the simulator's last-named composition gap: it emitted no `N` at all. See `docs/N-EMISSION.md` | — |
| fastmaps | 0.1.0 | 43 | **CLI-converged** — `--json` + v1 wrapper, `-j`, `--force`, `-q`, `FASTSUITE_THREADS`; `--corners` and `--preview-reads` | — |
| fastcheck | 0.1.0 | 108 | working — `--max-mem` + admission, insert measurement fixed (Task 10), **profile fitter fixed** (`430117a`: `BASE_POSITIONS` 60 → 512), **polyG head/tail and the per-position quality curve now MEASURED, not estimated** (see Known defects 2–3), harness green | `plans/2026-10-12_fastcheck-exactness.md` (tasks 5–9, 11–14) |
| fasthelper | 0.1.0 | 21 | working — fastcheck step, `check -c`, `--fail-on` | `plans/2026-10-12_fasthelper-hardening.md` — already done on disk: `run`, `compare`, `check -c`, `doctor`, `schema`, `--fail-on` (exit 2); the plan's open items are the *hardening* ones (5–8: `-P/--jobs`, resume/`--from`, `install --prefix`) |

**470 cargo tests, 0 failures**, plus the **pathology corpus: 9 assertions
passed, 0 failed**. Note that `scripts/test.sh` counts the pathology suite as one
aggregate line, so its own `GRAND` reads **471** (470 + the one aggregate line);
the 9 is the number of individual assertions inside it. (Measured this run;
the count moved up from the 465/466 recorded earlier because the tree gained
tests in the same session, not from the detector fix, which removed three
obsolete estimator tests and added three regression tests.) `pixi run doctor`
lists all six binaries; `bash
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
| fastsim ground-truth recovery | **10/10 parameters within tolerance** (was 8/10; the two failures were fastcheck fitter defects, now fixed — see "Known defects" 2–3) |
| fastsim insert resolution | 94–95% measured exactly at 200–288 bp (18 982 of 20 000 at insert 200, p50 200); the 288 bp ceiling is geometry, not a defect |
| trim accuracy vs simulator ground truth, tolerance 0 (101 real libraries, one harness) | **fastprep with `--head-dropout 0` leads**: 92.6% exact / 0.1% overtrim / 7.3% undertrim. fastp-default 92.2 / 0.2 / 7.6; fastp-tuned 84.1 / 8.6 / 7.3; **fastprep's shipped default (`--head-dropout 4`) is worst at 75.4 / 23.6 / 1.0**. Tolerance 2: 99.7 / 99.6 / 91.2 / 98.7. Aggregate of `bench/results/phase-d/overtrim-vs-fastp.tsv`; see `docs/TRIM-ACCURACY.md` |
| `head_dropout` axis (1515 runs, 3030 rows = 101 datasets × 3 variants × hd {0,1,2,3,4}) | `0` is the per-dataset winner on **100 of 101** datasets, `1` on the remaining one, `4` on **none**. Corpus total error at tolerance 0: **7.44% at `0` vs 24.64% at `4` (3.3×)**. `docs/OVERTRIM-SWEEP.md`, `bench/results/phase-d/head-dropout-axis.tsv` |
| drop asymmetry (the caveat under the trim result) | fastprep drops more reads than fastp at tolerance 0 (1.35 M vs 0.39 M), unevenly: virome **15.2% vs 0.0%**, RNA-seq 10.8% vs 2.1%, while on bacterial_isolate fastprep drops *fewer* (3.8% vs 5.1%). So the largest fastprep accuracy win (virome, 93.8% vs 92.1%) sits on a kept-set roughly one sixth smaller and is not apples-to-apples; the bacterial_isolate tie is. `docs/TRIM-ACCURACY.md` §Drop asymmetry |

## Phase D — the closed loop, re-measured

The campaign ran the full corpus: **102 real ENA datasets × 57 configs ×
3 variants = 306 comparisons**. Earlier revisions said 55 configs; the authored
grid is **57** (`sweep.py::GRID`, which now carries the `head_dropout` axis).

**PASS 142/306 = 46.4%**, 47 datasets passing every variant.

> **Re-measured 2026-10-16 — this replaces 167/306 = 54.6%.** That figure was
> inflated by a corpus artefact: `sim.v1.gt.tsv` and `sim.v2.gt.tsv` were
> **byte-identical in 100 of 102 datasets**, written in the same second, because
> `chain.sh` guards every artefact with `if [[ ! -s file ]]` and
> `run_campaign.sh` skips any `(dataset, variant)` with a status file — so when
> v0 was regenerated after an earlier fix, v1 and v2 were silently reused. All
> 18 duplicate-rate failures were that duplicated pair, and v0 had none.
> Rebuilt with `bench/scripts/rerun-campaign.sh`, which deletes each work dir
> before rebuilding: **`v1==v2` now 0/102, duplicate-rate failures 0.** Full
> analysis: `docs/CLOSED-LOOP-RERUN.md`. Evidence of the old state:
> `results/STALE-BEFORE-RERUN.txt`.

The closed loop scores fastsim against each library's **fastcheck**-measured QC
profile. It is **fastcheck-only** and independent of fastprep's trimming and of
the sweep winners.

| domain | pass rate | |
|---|---|---|
| human_gut | 94.1% (48/51) | |
| bacterial_isolate | 64.7% (33/51) | |
| generalised | 58.8% (30/51) | |
| virome | 35.3% (18/51) | 3′-biased |
| rna_seq | 19.6% (10/51) | 3′-biased |
| metagenomics | **5.9% (3/51)** | |

**`g_frac` is now the sole blocker**: it appears in **all 164 failing rows**.
Every other gate is noise by comparison (mean_q 8, gc 5), and the duplicate
gate is clear.

### The rate went *down* when the contamination was removed

On v0 — the variant that was never duplicated — the rebuild gives **14
pass→fail and 0 fail→pass**, with the `g_frac` mismatch 3–6× worse
(e.g. `DRR463958` 0.0033 → 0.0213, `SRR38668091` 0.0064 → 0.0400).

**The polyG fitter fix is the cause, not the N emission.** Isolated by
stripping `n_run` from the profile and re-running at the same seed, **12 of 14
still fail with the emission removed entirely**, at 0.0205–0.0395, with the
values barely moving. The emission tips 2 marginal cases.

The mechanism is the one `docs/OVERTRIM-CEILING.md` warned about: the fitter's
boundary deconvolution assumes terminal bases are i.i.d. at the measured
background G rate — *exact* for the simulator's random DNA (hence recovery
8/10 → 10/10), an *approximation* on real libraries with strong 3′ bias. The
closed loop fits from real data, and the failure lands on the read's last base.
`DRR001194` moves 0.0069 → 0.0198 from the fitter alone.

**The fitter became better at the metric it was tuned against and worse at the
metric the closed loop scores.** The fix direction is a per-position boundary
model, not a constant and not a wider tolerance.

### The path, for the record

Stages as recorded during the campaign. These were measured on the corpus that
turned out to be contaminated, so they are history, not a current baseline:

| stage | pass rate |
|---|---|
| as first measured | 95/306 = 31.0% |
| after the 3′ boundary-rule fix alone | 86/306 = 28.1% (the fix *lowered* it — it exposed a GC under-model that had been partially cancelling) |
| after the merged dup_rate fix | 113/306 = 36.9% |
| after the duplicate-replay fix | 131/306 = 42.8% |
| after head-G deconvolution | 167/306 = 54.6% (contaminated) |
| **after the corpus rebuild** | **142/306 = 46.4%** |

Median residual over the re-measured run: `g_frac_max_diff` **0.0218** (now
just *outside* the 0.02 bound), `mean_q_diff` 0.1735, `gc_diff` 0.0005,
`dup_ratio` 1.0154, `polyx5_diff` 0.0100. (The previous run's medians —
0.0140 / 0.1503 / 0.0002 / 1.25 / 0.0215 — are history, not a baseline.)

**The fitter is improved, but not unambiguously.** The "gated on the fitter"
framing in earlier revisions no longer applies: `BASE_POSITIONS` was raised
from 60 to 512 (`430117a`) so a high-GC library's position-wise G ramp is no
longer halved, alongside a measured `dup_rate`, and fastsim's insert-size
resolution went 0.5078 → 0.9491 at insert 200 as a result.

**But the later polyG run-length fix is a trade, not a win.** It made the
fitter's own recovery score 8/10 → 10/10 while *lowering* the closed loop,
because its boundary deconvolution assumes terminal bases are i.i.d. at the
measured background G rate — exact for the simulator's random DNA, approximate
on real 3′-biased libraries. That is the 54.6% → 46.4% drop above. Do not read
"the fitter is fixed" as a settled claim: it is better at the metric it was
tuned against and worse at the one the loop scores. See
`docs/CLOSED-LOOP-RERUN.md` §3.

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
   blocker to a publishable simulator accuracy number. (Checked 2026-10-16: the
   new ambiguous-base emission *does* move the metagenomics `g_frac` — e.g.
   DRR083326 0.0532 → 0.0396, DRR083353 0.0455 → 0.0293, DRR083390 0.0490 →
   0.0350 — but it does not close it (all still above the 0.02 bound) and the
   movement is confounded with the RNG-stream shift the emission causes, so it
   is not attributed to the N overlay. The metagenomics pass count is unchanged.
   `docs/N-EMISSION.md` §6.2.)
2. ~~**fastcheck's polyG detector** — reports 93.98% tail-free against a true
   66.01%, zero mass at t=3/t=15 (true 7.99%/6.99%).~~ **FIXED.** The profile
   fitter estimated the head/tail distributions from the per-position G
   fraction (`P(run >= k) ~= g_frac(k-1) - g_mean`), monotone-clamped and
   **capped at 12 bases**, while the tool already held the exact per-base
   run-length histograms (`polyx3`/`polyx5`) it publishes in its own report —
   so a 15-base run could not be expressed at all and whole modes of a
   multi-modal tail were lost to the clamp. It now reads the measured
   histogram and subtracts the biological boundary run (geometric,
   `(1-q) q^k`, q from the interior per-position G fraction). Measured against
   the simulator's ground truth the fitted R1 tail is now
   {0: 66.6%, 3: 7.98%, 5: 13.6%, 10: 4.75%, 15: 6.38%} against a true
   {66.0%, 7.99%, 14.0%, 5.0%, 7.0%}: **TV 0.280 → 0.013** (tolerance 0.20).
   The head row improves for the same reason, 0.025 → 0.005.
3. ~~**fastcheck's worst 10-position quality block** — misses by 0.45 Phred
   (mean passes at 0.12).~~ **FIXED (0.05 Phred).** The fitter stored the
   *measured* per-position mean as `pos_q`, while fastsim consumes `pos_q` as
   the **out-of-run** curve and applies the flat in-run mean to the run bases
   itself, so every tail base's collapse was counted twice and the curve came
   out 1.95 Phred low exactly in the block where the tails live. It now emits
   the deconvolved out-of-run curve (`out[i] = (measured[i] − p_i·tail_q −
   q_i·head_q) / (1 − p_i − q_i)`) plus the measured `run_q_mean`,
   `deep_head_q` and `pos_q_uncond`. Worst block **1.95 → 0.05** Phred
   (tolerance 1.5), curve mean 34.0 exactly.

   **Items 2+3 were the whole of fastsim's recovery score: the re-run
   recovers 10/10 parameters, up from 8/10.** Both fixes carry a regression
   test that fails on the pre-fix code (the tail fixture makes the old fitter
   emit `{"0": 7000, "12": 3000}` — no bin at 3, the 15-run folded into the
   12-base cap). Note the recovery fixture itself was stale: regenerating it
   with the current binaries moved the pre-fix score to 9/10 before the fix
   landed (the tail row passed at 0.126 by accident on the new sim), so the
   pre/post pair to quote is 0.280/1.95 against 0.013/0.05 on the same binary
   pair, and 8/10 → 10/10 on the papers' fixture lineage.
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
counts from the rebuilt PDFs: FastPrep **28 p**, FastCheck **18 p**, FastSim
**23 p**, FastMaps 14 p, FastHelper 16 p. Each has its
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
Separately, the Phase D closed-loop numbers have been re-measured to
**142/306 = 46.4%** (was 167/306 = 54.6%, inflated by the duplicated v1/v2
corpus rows — see `docs/CLOSED-LOOP-RERUN.md`). `bench/results/phase-d/closed_loop.tsv`
is the 46.4% file; `FastSim-Paper/benchmark/closed_loop.tsv` has been refreshed
from it and the FastSim and FastCheck prose updated.

**Resolved (was flagged as a paper-side inconsistency):** the FastPrep paper's
§5.6 lead-in (`sect-real-ena.typ`) cited a stale closed-loop figure,
"94-of-305", together with the wrong dominant failure gate (it named
duplicate-ratio; the TSV shows the per-position G-fraction at 120 of 139
failing rows, against duplicate-ratio 18). Both were corrected in `f4a5425`.
The authoritative figure is `bench/results/phase-d/closed_loop.tsv` —
**142/306 = 46.4%** (re-measured; see above).

## What the suite refuses to do (D-12)

PE overlap/merge, UMI handling, alignment-based QC, workflow-engine
semantics (arrays, resume), per-tile modules (until built), long reads (D-7).

## Where the ML campaign stands

Io is gone (2026-10-12); the plan is `docs/PLAN-ML-LOCAL.md` — all benchmark
compute on the desktop (16 cores, 48 GB), laptop interactive-only. **The
pilot is done and the full campaign ran** — 306 comparisons over 102 real
ENA datasets, at **46.4%** closed-loop pass (fastcheck-only; see above). The sweep
results stand independently: `winners.tsv` (102 rows, **11 distinct configs of
the 57-grid**, the winner being config 37 — the shipped default with
`head_dropout=0` — on 53 of them) and `fastp-comparison.tsv` (303 rows at gz6,
4 threads) —
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
