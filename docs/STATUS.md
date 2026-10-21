# FastSuite status

Verified 2026-10-19 against the current `main` (`3c0aa92`). Every number here was
verified by running it or by re-reading the cited artifact: the six cargo test
suites under `pixi run ci`, `scripts/test.sh` (pathology corpus included), the
Phase D closed loop re-read from its TSV, the trim-accuracy, `head_dropout` and
retention tables re-aggregated from `bench/results/phase-d/`, the `g_frac`
mechanism decomposition, and the paired-end overlap with/without experiment
(Cells A and B, both on the desktop corpus). **The suite is WIP: nothing here is
production-ready** (see the README banner).

## The tools

| tool | version | tests | state | next plan |
|---|---|---|---|---|
| fastdedup | 1.3.0 | 134 | shipped (published: bioRxiv 2026.04.29.721745) — **D-4 stdin/stdout streaming landed**: `-` for `-1/-2/-o/-p` with `--out-format`/`--out-compression`, and the documented refusal cases (two inputs as `-`, `--json -` with `-o -`, …) | — |
| fastprep | 0.3.0 | 203 | working — `--force`/`--dry-run`/`-q`, emits `fastsuite.report.v1` with the raw stats under `legacy`. Shipped `--head-dropout` default is now **0** (was 4 in 0.1.0; `3d76a31`, version bumped to **0.2.0** — a behaviour change, so a minor bump) | — |
| fastsim | 0.1.0 | 93 | working — measured insert CDF, SE mode, **clap-migrated**, **profile-driven ambiguous-base (N) emission** (`n_run`: read rate, run-length histogram, edge shares; recorded in the GT as `r1_n`/`r2_n`), and **all four (plane × artefact) pairs deconvolved on the consumer side** (`7631292`, the g_frac fix order's step 1(b)) — DRR083326's R1 terminal moved 0.5292 → 0.4140 against a real 0.4245. The insert ceiling (`2L − 12 = 288 bp`) is now stated in `--help` too | — |
| fastmaps | 0.1.0 | 96 | **CLI-converged** — `--json` + v1 wrapper, `-j`, `--force`, `-q`, `FASTSUITE_THREADS`; `--corners` and `--preview-reads`; **previews streamed instead of held resident** (`90573ea`: defaults 202.6 → 141.2 MB, `--buffer-mb 8` 100.5 → 36.1 MB, outputs verified byte-identical) | — |
| fastcheck | 0.1.0 | 143 | working — `--max-mem` + admission, insert measurement fixed (Task 10), **profile fitter fixed** (`430117a`: `BASE_POSITIONS` 60 → 512), **polyG head/tail and the per-position quality curve now MEASURED, not estimated** (see Known defects 2–3), **per-tile quality modules** (`--tiles`, opt-in, `93ed54a`), harness green | `plans/2026-10-12_fastcheck-exactness.md` — **tasks 5–14 all landed**; the ground-truth recovery claim is a gate now, not a session: `pixi run ground-truth` |
| fasthelper | 0.1.0 | 70 | working — fastcheck step, `check -c`, `--fail-on`, **post-trim residual-adapter guard**: a `run` step measures fastcheck's adapter content on the reads entering and leaving a trimming step and warns — or exits 2 under `--fail-on-adapter-residual` — when adapter survived (`tools/fasthelper/src/adapter_guard.rs`) | `plans/2026-10-12_fasthelper-hardening.md` — **every item is landed** (`run`, `compare`, `check -c`, `doctor`, `schema`, `--fail-on`, `-P/--jobs`, `install --prefix`, and `--from`/resume, which `runner.rs` implements via the per-step `ResumeEntry` re-run rule) |

**679 cargo tests, 0 failures** (143/122/187/89/77/60), plus the **pathology
corpus: 9 assertions passed, 0 failed**. Note that `scripts/test.sh` counts the
pathology suite as one aggregate line, so its own `GRAND` reads **680**
(679 + the one aggregate line);
the 9 is the number of individual assertions inside it. `pixi run ci` asserts
every per-tool count against `scripts/expected-test-counts.txt` and refuses a
rustc outside the pinned `>=1.99,<2`; **run it after any change** (the laptop's
ambient cargo is 1.95 while the desktop's is 1.99, and an unasserted total is
what let a stale test binary report 463 against 471 unnoticed). `pixi run doctor`
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
| fastmaps vs FastQC (5M pairs) | 9.17 s / 145 MB vs 16.1 s / 16 MB. The RSS is now decomposed term by term and the preview term is gone (`90573ea`: defaults **202.6 → 141.2 MB** on this batch's 5M-pair library, **100.5 → 36.1 MB** at `--buffer-mb 8`, outputs byte-identical). FastQC is still leaner, so the honest claim remains the composition plane plus thread scaling, not raw speed or RSS; that remaining lever was pulled on 2026-10-18 (`5e57224`: `--buffer-mb` now defaults to 8 MB, 48.6 → 22.7 MB on 2 M real reads at unchanged wall time), and the structural residue is the 8 bytes/read stream behind `heatmap.data` (PLAN-QC §1) |
| fastprep vs fastp (5M pairs, 1 and 2 cores pinned with `taskset`, gzip level 6 both) | **fastprep 73.99 s vs fastp 80.03 s on one core; 37.55 vs 40.55 on two.** fastprep leads at every size on both core counts (1.08–1.21× one core, 1.08–1.28× two). The earlier "34.37 vs 28.50, fastp leads" row came from a harness that set `-j N` without `taskset`: its `pinned` rows burned 2.4 cores. Source of truth: `FastPrep-Paper/benchmark/results.csv` (`cores` column, `cpu_s/wall_s` ≈ 1.00 on every row) |
| fastprep's drop rate tracks the injected must-drop + polyG | 7.7–7.9% drops vs the 4.7–4.8% injected (a retention check, not a trim-accuracy claim) — carried from the paper campaign, not re-measured this session |
| fasthelper vs MultiQC | zero orchestration overhead (0.667 s vs 0.647 s) + manifest provenance; MultiQC v1.35 emits **no output** on `fastsuite.report.v1` (no module for it) — 0 report files at every N from 25 to 800 samples |
| fastsim ground-truth recovery | **10/10 parameters within tolerance** (was 8/10; the two failures were fastcheck fitter defects, now fixed — see "Known defects" 2–3) |
| fastsim insert resolution | 94–95% measured exactly at 200–288 bp (18 982 of 20 000 at insert 200, p50 200); the 288 bp ceiling is geometry, not a defect |
| trim accuracy vs simulator ground truth, tolerance 0 (102 real libraries, one harness) | **fastprep with `--head-dropout 0` leads**: 89.1% exact / 0.6% overtrim / 10.3% undertrim. fastp-default 88.6 / 0.8 / 10.6; fastp-tuned 80.4 / 9.6 / 10.0; **fastprep's historical default (`--head-dropout 4`) is worst at 75.5 / 22.5 / 2.0**. Tolerance 2: 99.1 / 98.5 / 90.2 / 98.5. Re-measured on the regenerated corpus (see above); the ordering is unchanged but the level dropped ~3pp and undertrim roughly doubled, because the simulator now emits N-runs that no trimmer removes. Aggregate of `bench/results/phase-d/overtrim-vs-fastp.tsv`; see `docs/TRIM-ACCURACY.md` |
| `head_dropout` axis (1515 runs, 3030 rows = 101 datasets × 3 variants × hd {0,1,2,3,4}) | `0` is the per-dataset winner on **100 of 101** datasets, `1` on the remaining one, `4` on **none**. Corpus total error at tolerance 0: **7.44% at `0` vs 24.64% at `4` (3.3×)**. `docs/OVERTRIM-SWEEP.md`, `bench/results/phase-d/head-dropout-axis.tsv` |
| drop asymmetry (the caveat under the trim result) | **Re-attributed 2026-10-18: this is the dedup default, not the length floor.** The PXMT status byte shows 100% of the drops on the worst datasets are `Duplicate` — fastprep deduplicates by default, fastp does not: SRR20754397 167,162/167,162 drops, DRR015139 32,138/32,138, and `--min-len 15` leaves those histograms untouched. With `--no-dedup` the aggregate retention is **98.47% against fastp's 97.23%** (virome 99.74% vs 99.7%, rna_seq 98.90% vs 97.8%), i.e. the suite does not retain fewer reads, it deduplicates, and the trim harness charges it for that because a dropped read leaves the scored set. The earlier "virome 18.2% vs 0.3%" framing described the symptom and mis-named the cause; the corrected per-domain retention table is in `docs/TRIM-ACCURACY.md` §Read retention. **Read retention and exact rate together, never the rate alone**: `--min-len 0.9L` shows the highest fastprep exact *rate* (89.470% vs 89.088%) with 441k *fewer* exact mates and 2.8 pp lower retention, while `--no-dedup` raises both — so "dropping flatters the score" depends on the class dropped. |

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

## Fixed in the 2026-10-21 review round (full-suite review: six crates + six papers)

Every item below was reproduced on the release binary before the fix and re-run
after it. The row-level findings are in
`~/.hermes/skills/bioinformatics/fastsuite/references/review-2026-10-21.md`.

* **fastsim's ground truth contradicted itself.** `full_flag` was set per pair
  (either mate full) while the doc said "the pair is entirely artefact", so an
  R2-only full read marked a valid R1 `drop` (6.7% of rows on an R2-head-150
  profile). The GT now carries `r1_full`/`r2_full`, and `full_flag` is the AND
  of the two. The `rN` columns count the instrument's no-call overlay only; the
  no-adapter `N` filler is published as `adapter_pad_bases`. `is_dup=1` always
  has an earlier byte-identical record (the empty pool used to flag pair 0).
  `profile.gc` is now the authority for random DNA (`--gc` overrides it; a
  gc=0.30 profile used to emit 0.499). An empty `adapter.seq` is a no-adapter
  run instead of a modulo panic, and `--single` names no R2 and counts one read
  per pair. FASTQ golden pins are unchanged; the GT pin was re-based.
* **fastdedup's guards were bypassable and its PE writer mixed formats.**
  `-o ./x -p x` passed the string comparison and interleaved two writers into
  one mangled file (exit 0); formats are now compared on resolved paths, PE
  parallel refuses mismatched mate formats like the serial path, and a 0-byte
  output resumes (the documented crash-recovery case used to abort). Paired
  metrics gained `pairs_in`/`pairs_out`; the abundance temp file is removed on
  failure paths; `-j 0` reports the resolved thread count.
* **fastprep's adapter prefilter missed real matches.** The sound-k predicate is
  not monotone (false again at k=21,22,23,26,27,31 at eps=0.2), so a 21 bp
  TruSeq prefix with 5 mismatches was never examined (0 trimmed, 121 bp
  emitted, exit 0). The bound is now the monotone-safe `ms_safe` and the
  mismatch tolerance floors instead of ceilings. The dedup anchor no longer
  degenerates to the full read at shipped defaults, the report carries
  `assumptions` and `phred_offset` (D-6), the `--min-len` floor keeps a
  readthrough read instead of dropping it, the metadata viewer knows cause 8,
  and singleton salvage no longer trips the total-drop refusal. The inert
  `-C target-cpu=native` config and its AVX2 comments are gone.
* **fastcheck measured inserts short.** The overlap search took the leftmost
  anchor a 3-mismatch tolerance admits; on an exact 120 bp library 13.7% of
  pairs measured 1-10 bp short (mean 119.53, sd 1.22). It now searches from the
  largest overlap, so the fixture measures 120.0 with sd 0.0. `--only kmer`
  shares the fused rolling register (it was 2.6x slower than running everything),
  a truncated mate file is refused rather than silently rescaled, `.1/.2` SRA
  mate suffixes are accepted, `--overrep-mode second-pass` refuses stdin,
  `kmer_top_frac`/`base_composition` have documented denominators, and the
  fingerprint test pins the constant it claims to.
* **fastmaps' `--json` wrote the grid, not the report.** The v1 report now goes
  to `--json`'s PATH and the grid stays in `-o/heatmap.json`; the wrapper names
  the real grid file. A 1 bp read no longer underflows the corner window
  (`polyA_frac 1` for an all-T read, exit 0), qual!=seq is a refusal instead of
  an exit-101 panic, the output plane refuses a size that would OOM instead of
  silently allocating 3.8 GB, `--json` inside a fresh `-o` works, and the
  clobber guard covers reports outside `-o`.
* **fasthelper lost sample identity and could not run its own advertised chain.**
  Reports carry the sample name instead of re-deriving it from the filename
  (`fastcheck.json` used to resolve to sample "fastcheck"), `fastsim` works as a
  chain step (step-private prefix; its FASTQs become the next step's input), a
  step that exits 0 without writing its `--json` fails the run, `--from` is
  validated, the simulated control's name joins the duplicate check,
  `assumptions` and fastmaps' length scalars reach the aggregate, HTML/MD cells
  are escaped, and `-j` reaches every child (fastprep resolves the 0 sentinel).

Known open after this round: fastcheck's insert estimator is still a heuristic
about paired geometry (the fix removes the systematic bias, not the tail);
fastprep's clippy debt (~70 warnings) is reduced but not zero; the campaign TSVs
were produced with the pre-2026-10-21 GT (`full_flag` semantics changed), so
per-dataset reproduction of the closed loop needs a re-run before the numbers are
quoted against this build.

## Fixed in the 2026-10-20 review round (five commits, `a70e510..b0dc9a1`)

Each of these was reproduced on the built binary before the fix and re-run after
it; the commands and the before/after output are in
`~/.hermes/skills/bioinformatics/fastsuite/references/review-2026-10-20.md`.

* **fastdedup decoded only the first gzip member on the parallel path** — a
  200-member (BGZF-shaped) file of 20 000 reads gave 100 reads at `-j 4`, exit 0.
  Now `-j 1/4/8` all give 20 000. `--abundance` refuses an existing output, and
  `-o` == `-p` is refused. (`a70e510`)
* **fastcheck's `--profile` took `read_len` from the longest read** — 1000x151 bp
  + 1x300 bp produced a 300 bp profile, and the guard was keyed on a quantity
  that could not catch it. Now the mode, with a fit warning naming the outliers.
  The HTML insert table divided by the `--insert-pairs` cap (2000 pairs rendered
  as "20000", buckets summing to 20%); now the sampled count, buckets sum to
  100%. `--reads N --html` no longer claims "no read was skipped", and
  `scan_pair` checks mate headers instead of trusting ordinal pairing. (`3470a0b`)
* **fasthelper's aggregate ingested every `*.json` under the outdir** — a
  `[dedup,prep]` run into a directory that had held `fastmaps` re-published the
  fastmaps row. It now ingests only the paths this run's manifest records.
  Resume ignores the recorded exit code and could chain a crashed step's
  truncated output forward; it now requires exit 0 plus output signatures.
  `run` gained `--force` (the documented chain could not be re-run), unknown
  steps are refused like `check` refuses them, the tool precondition is derived
  from the enabled steps, and `--dry-run` no longer overwrites a real manifest.
  (`023cf9a`)
* **fastmaps' corner fractions are length-biased** — the window is anchored to
  `max_seq_len`, so adding one 200 bp read to ten 150 bp all-T reads took
  `polyT_frac` from 1.0 to 0.0 with no diagnostic. The grid and the wrapper now
  publish `min_seq_len`/`max_seq_len`/`corner_window` and warn when the spread
  exceeds the window. The corner definition is unchanged, deliberately: a
  length-relative plane is a design change, not a fix. clippy also builds the
  crate now (two `erasing_op` errors in a test had left it with no lint gate).
  (`7e05c5c`)
* **fastprep asserted Phred+64 from a preponderance of high bytes** — a Phred+33
  library at Q43 was read as +64, shifted down 31 and emptied by `--mean-q 25`
  (5000/5000 reads, exit 0), and `--phred-offset 33` was silently overridden
  because `--auto-phred` ships on. Now the auto path never asserts (+33, with
  `--phred-offset 64` for legacy data), an explicit offset disables detection,
  and a run that keeps nothing exits non-zero (`--allow-total-drop` opts out).
  `--min-len` now also filters: an absolute floor above the read length drops the
  read instead of being clamped to a no-op. (`24ac984`)

**Consequences for the campaign numbers.** The `--min-len` change alters integer
floors only (`30`, `50` on short-read datasets); the shipped `0.9L` rows behave as
before, and the shipped-default trim-accuracy rows are unaffected. The grid should
still be re-derived before the sweep numbers are quoted again.

## Known defects (measured, not guessed)

1. **metagenomics closed-loop pass rate is 5.9% and has not moved across any
   of the four fixes.** Re-read at each commit in the
   chain — `efaf834` (first full-corpus run), `c3ad7a0` (+dup_rate),
   `0a19d06` (+duplicate replay), `bd31d9a` (+head-G deconvolution) — the
   count is **3/51 in all four**, while the other domains gained 25–60
   points over the same span. All 48 failures are `g_frac`, at a median
   over-tolerance error of 0.0529. Neither the dup model, nor the 3′
   boundary rule, nor the head-G deconvolution touched it: this is a
   **distinct mechanism — DIAGNOSED 2026-10-18, and it is not in the fitter.**
   All 16 failing datasets are accounted for by three consumer-side mechanisms:
   a single scalar `read_len` applied to *both* mates while real libraries are
   asymmetric (DRR083073: 0.1029 pooled → 0.0358 length-matched), a
   double-counted injected run (the plane already holds the artefact Gs, and
   fastsim deconvolves one of four plane × artefact pairs: DRR083326 R1 terminal
   real 0.4245 vs sim 0.5270), and a frame shift (`build_read` prepends `[G;h]`,
   so read position p shows fragment p−h while the plane is a read-position
   measurement). The run-length round trip, the N overlay, the body zone and GC
   are ruled out with numbers. No single patch closes them: the 3′-plane fix
   moves a row's length-matched gate 0.1047 → 0.0288 while its *judged* pooled
   gate goes 0.0383 → 0.0969, because that row sits on two cancelling errors
   today. Full evidence and method: `docs/GFRAC-MECHANISM.md`. This remains the
   main blocker to a publishable simulator accuracy number; it is no longer
   undiagnosed. (Earlier check, 2026-10-16: the ambiguous-base emission *does*
   move the metagenomics `g_frac` — DRR083326 0.0532 → 0.0396 and others — but it
   does not close it, and the movement is confounded with the RNG-stream shift
   the emission causes. `docs/N-EMISSION.md` §6.2.)
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

New and open after the review round:

* **gzip level 6 costs ~7x the compressed output in peak memory.** The
  multi-threaded gzippy writer buffers the stream: on 5 M pairs (888 MB compressed)
  fastdedup peaks at 4.37 GB, against 153 MB at `--gz-level 1` on the single-threaded
  flate2 path — which is also *faster* there (18.3 s against 29.5 s) but produces
  1 681 MB instead of 888 MB. The footprint scales with input (~890 MB at 1 M pairs),
  so a 100 M-pair gzipped run needs ~80 GB. It is a trade, not a leak, but nothing
  documented it and the default is level 6. A bounded hand-off queue now caps the
  *other* contributor (see the thread table in `bench/profile/README.md`).
* **fastmaps' corner fractions cannot see a 3' tail.** `polyG_frac` counts reads
  that are *nearly all G* (the corner window is `max(size/50, 2)` at a grid sized
  for the longest read), so a tail is a displacement, never a corner. Measured on
  a 400-read R2 fixture: 200 reads with the last 20 bases G on a random body
  contribute 0, while 200 entirely-G reads give `r2_polyG_frac` 0.5. A NovaSeq
  polyG signature needs a positional metric (3'-end base composition), which does
  not exist yet; the paper states the limit, the tool does not warn about it.
* **fastsim ignores `profile.gc`** (`main.rs:221` uses `--gc`): a profile with
  `gc=0.30` emits GC 0.4881. Measured, open.
* **fastprep's `--adapter-error-rate` rounds the mismatch tolerance UP**
  (`ceil(k*eps)`), so 0.2 admits 3/14 = 21.4% at the shipped overlap. Switching
  to `floor` perturbs the pigeonhole seed-width bound and the differential test
  catches the prefilter missing a real short-overlap adapter — it needs its own
  test-first round, not a one-line edit. Open, deliberately.
* **`winners.tsv` and `fastp-comparison.tsv` disagree on the winning config for
  93 of the 101 datasets they share** (11 vs 16 distinct configs). Reconcile
  before quoting any sweep winner.
* The clone-vs-laptop checkouts differ (`6945219` before this round, desktop
  `57f3140`): run `fastsuite-local-artifacts/sync-fastsuite.sh` from the laptop.

New and open (not defects of a measured component, but unresolved):

* **The trim-accuracy win is not apples-to-apples on every domain** — see the
  drop-asymmetry row above. **Re-attributed 2026-10-18: the cause is the dedup
  default, not the length floor.** fastprep deduplicates and fastp does not, and
  the harness charges a dropped read to the kept-set; with `--no-dedup` the
  aggregate retention is 98.47% against fastp's 97.23%. So the win is not
  apples-to-apples *because the two tools are configured differently*, and the
  honest comparison quotes both configurations (decision D-13). The
  `--min-len 0.9L` default was measured and kept.
* **The `head_dropout` default is now 0** (committed `3d76a31`), which is the
  measured right default on this corpus (100/101 datasets; `1` on the remaining
  one; `4` on none). No version bump, so this is committed-but-unreleased: a
  binary built from `main` trims at 0, a released 0.1.0 tarball still trims
  at 4.
* **`--pe-overlap` is a measured null and ships opt-in** (2026-10-19). Exact
  matches identical (347,171) with and without it on a synthetic cohort built to
  be scorable, overtrim +106, +0.000 pp against a pre-registered 0.5 pp threshold.
  It is kept as an instrument because its guard is real — the boundary applies
  only when adapter starts exactly at it, so the mode can never trim what the
  ordinary adapter path could not itself confirm — and because it is the only way
  to re-run the question on a different corpus. Throughput cost measured: +3.7%
  single-threaded, +7.8% at 4 threads (best-of, 400k pairs).
* **The metagenomics `g_frac` blocker stands at 3/51.** Mechanisms A (harness
  pooling) and B (plane deconvolution) are fixed and measured; the residual is
  mechanism C (head-run frame geometry, `build_read` prepends its head run so
  read position *p* shows fragment *p−h*) and D (per-mate read length, one
  scalar applied to both mates). Both need corpus re-runs.

## Where each paper stands (local-only repos, never pushed)

All five papers are **written** and were sent to Discord on 2026-10-15. Page
counts from the rebuilt PDFs: FastPrep **28 p**, FastCheck **19 p**, FastSim
**24 p**, FastMaps 14 p, FastHelper 16 p. Each has its
`-Bench` repo (committed `results.csv`, seeded generator, VmHWM peak-RSS
harness; reads gitignored). The FastDedup paper is separate and older (28 p,
bioRxiv 2026.04.29.721745).

**The FastPrep paper now contains the reversed Section 5.6, and it is
committed.** `sect-real-ena.typ`'s trim-accuracy section (`@tbl-trim-acc`) was
rewritten to conclude that the *shipped default* was the problem, not the
algorithm: fastprep at `--head-dropout 4` is the worst of the four
configurations (75.5% exact at tolerance 0) and fastprep with the rule off
(89.1%) matches or beats fastp-default (88.6%). Those figures are the
re-measured ones — the pre-regeneration values were 75.4 and 92.6/92.2. The
reversal is commits
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
