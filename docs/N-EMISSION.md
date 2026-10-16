# Ambiguous bases (N) in fastsim — measurement, emission, and what it actually closes

> **Outcome in one paragraph.** `fastsim` now emits ambiguous bases (`N`),
> profile-driven, calibrated from the library's own measured no-call /
> padding structure rather than from a constant. On a synthetic profile whose
> N model matches a padding-heavy library the emission moves the simulated
> duplication level **0.0199 → 0.0397** and `dup_multiplicity_max` **5 → 4014**
> (real `bacterial_isolate/DRR000466` measures 2,035) — the padding-driven
> duplicate class now exists. **But it does not close the 18 duplicate-rate
> failures that motivated it, because those failures do not reproduce:** in the
> authoritative `closed_loop.tsv` the variant-1 and variant-2 rows are
> byte-identical in 100 of 102 datasets, which two independent seeds cannot
> produce, and every one of the 18 dup failures is in that duplicated
> v1/v2 pair. On freshly regenerated simulations the duplicate gate passes on
> every dataset rebuilt here. That is the honest result and it overturns the
> premise (see §6).

---

## 1. The gap

The closed loop compares a real library's QC against a simulation of that
library fitted from its own QC. Its duplicate-rate gate (ratio within 2×)
failed on 18 of 306 comparisons, and on all 18 the simulated rate was the
lower one. Two causes had been named: (a) `fastsim` emitted **no** ambiguous
bases, and (b) a duplicate was not sequence-identical to its original. (b) was
fixed in `781c396` (a duplicate now replays the pooled pair's exact bytes).
(a) is what this work addresses: `random_dna_bias` samples only ACGT, and the
adapter-padding path that could emit `N` is unreachable because a measured
profile's insert CDF sits above the read length.

## 2. What the real libraries actually contain

`bench/scripts/n_stats.py` measures, per FASTQ mate: the read-level N rate, the
run-length histogram, the run position (leading / interior / trailing), and the
concentration of N mass in duplicate reads. It is deliberately independent of
fastcheck, which reports only *that* a read carries a non-ACGT base
(`metrics.reads_with_non_acgt`) and nothing about runs.

Measured on rebuilt 100 k-pair units (deterministic subsample of the ENA runs;
`--pairs 100000`, the same stride the campaign used — verified by the real
report's `duplication_level` matching the corpus row to four decimals):

| library | R1 reads with N | runs/read | run p90 | run max | leading | trailing | interior | dup reads with N |
|---|---|---|---|---|---|---|---|---|
| human_gut/DRR425490 | 0.87% | 1.54 | 2 | 6 | 0.20 | 0.00 | 0.80 | **0.00** |
| human_gut/DRR425499 | 0.42% | 1.02 | 1 | 2 | 0.36 | 0.00 | 0.64 | **0.00** |
| metagenomics/DRR083326 | 0.15% | 1.06 | 1 | 2 | 0.22 | 0.00 | 0.78 | **0.00** |
| metagenomics/DRR083353 | 0.16% | 1.09 | 1 | 2 | 0.21 | 0.03 | 0.76 | **0.00** |
| metagenomics/DRR083390 | 0.22% | 1.01 | 1 | 1 | 0.44 | 0.00 | 0.56 | **0.00** |
| bacterial_isolate/DRR000466 R1 | 2.48% | 2.51 | 109 | 111 | 0.11 | 0.26 | 0.63 | **0.87** |
| bacterial_isolate/DRR000466 R2 | 2.30% | 2.36 | 111 | 111 | 0.37 | 0.09 | 0.54 | **0.95** |

Two regimes, and the split is real:

* **Sparse no-calls** — the human_gut and metagenomics libraries carry N on
  under 1% of reads, almost all as isolated single bases (max run 6), mostly
  interior or leading, and **not one duplicated read carries an N**. For these
  libraries the N term cannot move the duplicate rate at all.
* **Padding** — `bacterial_isolate/DRR000466` is the other extreme, and the
  rebuilt unit reproduces the corpus measurement exactly. 2.48% (R1) / 2.30%
  (R2) of reads carry an `N`; the run-length histogram is **bimodal**, ~2,800
  runs of length 1 (isolated no-calls) plus a hard spike at 108–111 bases
  (2035 runs at 111 on R1, 1759 at 111 on R2 — the full read), i.e. reads whose
  signal failed and were padded to the end. 87.3% of R1's duplicate reads and
  95.3% of R2's carry an `N` (pooled 92.6%, against the corpus record's 92.8%),
  and the top repeated sequences are `NNNN...`/`GNNNN...`. This is the class
  the simulation did not have.

The two regimes coexist inside a library, not only across libraries, so the
model has to allow both in one profile: the run-length histogram carries the
isolated no-call mass (length 1) and the padding spike in the same distribution,
and the placement terms decide which edge a run touches.

(The generalised libraries — `DRR000432`, `DRR013378` — sit in between; see the
prior corpus measurements in `docs/DUP-MECHANISM.md`: DRR000432 has 492 of its
duplicate reads carrying an N, 8.7%, against 5,164 genuine biological
duplicates.)

## 3. The mechanism, measured in isolation

A synthetic 150 bp profile whose only difference from a plain library is an
`n_run` block — `read_rate 0.02`, all runs length 150, all trailing — run
through the same fit → simulate → measure path, old vs new simulator:

| | duplication_level | dup_multiplicity_max | reads_with_non_acgt |
|---|---|---|---|
| before (no N emission) | 0.0199 | 5 | 0 |
| after (N emission) | **0.0397** | **4014** | 4014 |

The mechanism is exact, not approximate: with `h = t = 0` the biology window is
the whole read, a length-150 run fills it, and every padded read becomes the
identical string `NNNN...NN`, so the fingerprint table registers one sequence
with multiplicity ~4,000 and the duplication level doubles. This is the
`NNNN...` class the real libraries have and the simulation did not, and it also
fixes the other direction of the same defect — the simulated
`dup_multiplicity_max` of 3–4 against the real 2,035 was the same missing class.

## 4. The implementation

* **Profile schema** — a new optional object `n_run`:
  ```json
  "n_run": { "read_rate": 0.024, "len": {"1": 9000, "150": 400},
             "leading": 0.05, "trailing": 0.60 }
  ```
  `read_rate` = share of reads with ≥1 N; `len` = run-length histogram;
  `leading`/`trailing` = share of runs touching the read's first / last base
  (the rest interior). Loaded into `Profile::n_run: Option<NModel>`.
* **The fitter** (`fastcheck --profile`, `tools/fastcheck/src/profile.rs`) now
  measures it, over the **merged** R1+R2 accumulators exactly as it measures
  `dup_rate` and `reads_with_non_acgt`. This needed a real counter: `Accum`
  gained an ambiguous-run-length histogram (`n_run_hist`, `n_run_leading`,
  `n_run_trailing`), built in the existing one-pass `observe` loop — the
  per-position `base_pos[..][4]` plane cannot recover a run length, because a
  run is not a per-position event. A library with no ambiguous base emits **no
  `n_run` key**, and a zero `read_rate` parses back to `None`.
* **The emission** (`fastsim`, `quality.rs::inject_ambiguous`) — a per-read
  Bernoulli at `read_rate`; on fire, one run of length drawn from `len`,
  clamped to the biology window `[h, h + bio_len)`, placed at the window's 5′
  edge, its 3′ edge, or interior per `leading`/`trailing`. It is a sequencing
  no-call, so it lives in the sequencing stage next to the substitution model.
* **It never touches an injected G-run**, so `h`, `t` and the trim ground truth
  `[h, min(I, L−t))` are unchanged. This is deliberate: an N that overwrote a
  tail G-run would silently move the trim window and corrupt the sweep's GT.
* **Ground truth stays honest** — the GT gained `r1_n`/`r2_n` columns carrying
  the number of N bases the overlay put in each mate, so an injected N is a
  named artefact in the reference the whole trim benchmark is scored against,
  not a silent change to the read. The columns are appended, so every existing
  reader (sweep.py's fixed indices 1–9, the `r["r1_h"]`-by-name readers) still
  parses the file.
* **Report** — fastsim's `fastsuite.report.v1` now carries
  `metrics.ambiguous_bases`, `reads_with_ambiguous` and
  `reads_with_ambiguous_rate`.

## 5. Blast radius: the RNG stream shifts, and it is stated

`fastsim` draws artefacts, biology and quality from one seeded `ChaCha8Rng`,
and the N overlay draws from that same stream **only when the profile has an
`n_run` key**. Consequences, in order of what must be true:

* A profile **without** `n_run` (the hand-authored test profiles, any profile
  fitted before this change) consumes **zero** extra draws, so its output is
  byte-identical to before. Verified by
  `ambiguous_overlay_is_absent_without_the_profile_field`.
* A profile **with** `n_run` consumes extra draws, so every downstream draw
  shifts and the simulation becomes a **different random realisation of the
  same profile**. This is the same class of change as `324506d`, which shifted
  the stream and silently invalidated stored trim numbers
  (`docs/TRIM-ACCURACY-DISCREPANCY.md`).
* Therefore **any stored `sim.*` FASTQ/GT, sweep `c*.pxmt`, or `trim-accuracy`
  row produced from a profile that now carries `n_run` is stale.** The
  staleness guard and the input-identity resume key added in `73c2385` will
  refuse to score an old PXMT against a regenerated sim, so this surfaces as a
  refusal rather than a silent wrong number — but the corpus must be
  regenerated (fit → simulate → sweep → score) before its numbers are quoted
  again.

The g_frac movement reported in §6 is **confounded by this shift**: the
"after" runs change N emission *and* the stream at once, so the g_frac change
is **not** attributed to N emission.

## 6. The closed loop: what reproduced, and what did not

### 6.1 The 18 dup failures do not reproduce, and the corpus rows are duplicated

Rebuilt the deterministic 100 k-pair units for eight corpus datasets from ENA
(the same `read_count // 100000` stride, md5-verified against ENA) and re-ran
each through the exact `chain.sh` path. The **real** side reproduces the corpus
row: DRR425490 `dup_rate_real` 0.0103 (corpus 0.0103), DRR425499 0.0095
(0.0095), DRR083326 0.0005 (0.0005).

Before/after on those datasets, current (HEAD) simulator, all three variants:

| dataset | corpus dup (v0 / v1=v2) | rebuilt dup v0/v1/v2 | gate |
|---|---|---|---|
| human_gut/DRR425490 | 0.0103 / 0.0051 | 0.0103 / 0.0103 / 0.0103 | **PASS all 3** |
| human_gut/DRR425499 | 0.0095 / 0.0047 | 0.0096 / 0.0094 / 0.0102 | **PASS all 3** |
| metagenomics/DRR083326 | 0.0005 / 0.0001 | 0.0004 / 0.0006 / 0.0006 | **PASS on dup** (fails g_frac) |
| metagenomics/DRR083353 | 0.0001 / 0.0000 | 0.0002 / 0.0001 / 0.0001 | **PASS on dup** (fails g_frac) |
| metagenomics/DRR083390 | 0.0003 / 0.0001 | 0.0002 / 0.0003 / 0.0003 | **PASS on dup** (fails g_frac) |

The corpus rows' `dup_rate_sim` are, on the failing rows, about half the real
rate. Freshly generated simulations land *on* the real rate (ratio ≈ 1.0). Two
independent pieces of evidence say the failing rows are an artefact of stored
sims, not of the simulator:

1. **`v1` and `v2` are the same file.** In 100 of 102 datasets the variant-1
   and variant-2 rows of `closed_loop.tsv` are byte-identical across
   `g_frac_max_diff`, `mean_q_diff`, `gc_diff`, `dup_rate_real` and
   `dup_rate_sim` (e.g. `generalised/DRR000426` v1 = v2 =
   `0.0025 76 25.7103 25.6465 0.0638 … 0.0245 1.48571`). Two seeds (1042,
   2042) cannot produce identical statistics to six decimals; the stored v1
   and v2 simulations are one artefact measured twice.
2. **Every one of the 18 dup failures is a v1/v2 row.** Variant 0 has **zero**
   duplicate-rate failures across all 102 datasets. The regenerated
   simulations all pass the duplicate gate.

This is consistent with the campaign's own record that only `sim.v0` was
regenerated after the last simulator fix (`docs/TRIM-ACCURACY-DISCREPANCY.md`
§2, "15:54–15:58 all `sim.v0` FASTQ/GT regenerated") together with `chain.sh`
being idempotent (it re-simulates only when `sim.vN.gt.tsv` is missing), so
v1/v2 kept their pre-fix content while their reports were re-read.

### 6.2 The N emission does not move the pass rate

On the five rebuilt datasets, before vs after (N emission on):

| dataset | g_frac before | g_frac after | verdict before → after |
|---|---|---|---|
| human_gut/DRR425490 | 0.0026–0.0027 | 0.0030–0.0051 | PASS → PASS |
| human_gut/DRR425499 | 0.0022–0.0025 | 0.0025–0.0054 | PASS → PASS |
| metagenomics/DRR083326 | 0.0523–0.0532 | 0.0387–0.0396 | FAIL → FAIL (g_frac) |
| metagenomics/DRR083353 | 0.0426–0.0455 | 0.0293–0.0326 | FAIL → FAIL (g_frac) |
| metagenomics/DRR083390 | 0.0479–0.0490 | 0.0338–0.0353 | FAIL → FAIL (g_frac) |
| bacterial_isolate/DRR000466 | 0.0145–0.0159 | 0.0252–0.0255 | FAIL → FAIL (g_frac, mean_q) |

The duplicate-rate gate was already clear on fresh simulations and stays clear;
the g_frac gap narrows on the metagenomes but does not close, and that narrowing
is **not** attributable to N alone (§5). The pass rate does not move.

### 6.2a The one library that does move — and it overshoots

`bacterial_isolate/DRR000466` is the padding-heavy library, so it is the test
that matters. The emission behaves as the mechanism predicts — the simulated
duplicate rate **rises** — but it does not land on the real rate:

| DRR000466 | dup_real | dup_sim | ratio |
|---|---|---|---|
| before (HEAD simulator) | 0.0135 | 0.0135 / 0.0140 / 0.0130 | 1.00 / 1.04 / 1.04 |
| after (N emission) | 0.0135 | 0.0179 / 0.0188 / 0.0182 | 1.33 / 1.39 / 1.35 |

and its `g_frac_max_diff` rises from 0.0145–0.0159 (inside tolerance) to
0.0252–0.0255 (outside), turning a near-pass into a failure.

Two distinct effects, and they are worth separating:

* **The overshoot is structural, not a bug in the overlay.** The profile's
  `dup_rate` is the library's **total** measured duplication, and on this
  library that total is almost entirely padding (92.6% of duplicate reads carry
  an `N`). The simulator draws pool duplicates at `dup_rate` **and** now
  collapses the padded reads onto shared `NNNN...` strings, so the padding class
  is counted twice and the total lands ~1.35× high. Closing that properly needs
  the profile to separate the two duplicate sources (a padding-driven share and
  a pool draw that leaves room for it), which is a model change, not a
  calibration constant — recorded here as the honest next step rather than
  guessed at.
* **The g_frac movement is confounded by the RNG-stream shift** (§5), so it is
  not attributed to the overlay. It is, however, a real regression in outcome on
  the one dataset that has the padding class, and it is a reason not to claim
  the emission is a pass-rate win anywhere.

Note the *before* row is itself a correction: with the current (HEAD) simulator
the gate is 1.00–1.04 on this library, not the 1.57 the corpus records — more
evidence for §6.1.

### 6.3 What could not be rebuilt

Nine of the 102 corpus datasets are the duplicate-failure set; six of them were
rebuilt and re-run here (the two human_gut, three metagenomics, and
bacterial_isolate/DRR000466). Three were not affordable in this session:
`DRR000617` (10 GB per mate), `DRR000432` (2.6 GB per mate) and `DRR013378`
(3.2 GB per mate); and `DRR083337` failed ENA's md5 check mid-download
(`fail-md5-upstream` — ENA served different bytes than the manifest). The
finding in §6.1 is structural (the v1/v2 duplication holds across all 102
datasets), so it does not depend on the rebuilt subset; the per-library N
statistics for the three un-rebuilt generalised/isolate libraries are carried
from the corpus record (prior session) and marked as such in §2.

## 7. What this leaves open

* **The `dup_ratio` failures in `closed_loop.tsv` are a corpus artefact.**
  Regenerating all three variants and re-running the loop is the honest fix;
  the pass rate must be re-measured, not patched. Until then the 167/306 figure
  is a lower bound contaminated by a duplicated v1/v2 generation.
* **The per-position G fraction remains the real blocker** (38 of 102 datasets
  at variant 0, 41 at v1/v2), concentrated in metagenomics; the N emission is a
  genuine gap now closed but it is not that blocker.
* The ratio gate's tolerance caveat stands (`docs/DUP-MECHANISM.md`): a
  factor-of-2 ratio punishes sub-percent absolute differences at near-zero
  rates. With the duplicate gate clear on fresh sims, a companion absolute
  tolerance is still worth adding.
