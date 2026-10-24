# FastSuite Phase D — full failure analysis and fixes

> **Superseded (2026-10-16).** This is a point-in-time record from 2026-10-15.
> The closed loop later reached **167/306 = 54.6%** (52 of 102 datasets passing
> every variant), not the 113/306 recorded here, and the suite is now **852** on
> the `scripts/test.sh` GRAND line (not 461). Authoritative closed-loop figure:
> `bench/results/phase-d/closed_loop.tsv`. Current state: `docs/STATUS.md`.

**Corpus:** 102 real ENA datasets × 55 configs × 3 variants = 306 (dataset, variant) comparisons.
**Method (the closed loop):** fastcheck fits a profile from the real library → fastsim simulates from that profile → fastcheck re-measures the simulation → the two measurements are compared against tolerance. A comparison "closes" only when all four tolerances clear at once.

---

## 0. Headline

**PASS 113/306 = 36.9%**, 34/102 datasets passing every variant.

| domain | pass rate |
|---|---|
| human_gut | 90.2% (46/51) |
| virome | 49.0% (25/51) |
| bacterial_isolate | 39.2% (20/51) |
| generalised | 25.5% (13/51) |
| rna_seq | 11.8% (6/51) |
| metagenomics | 5.9% (3/51) |

Medians: `g_frac_max_diff` 0.0157 (tol 0.02) · `mean_q_diff` 0.1494 (tol 0.50) · `gc_diff` 0.0004 (tol 0.01) · `dup_ratio` 1.8205 (tol 2.0).

The pass rate did **not** move monotonically, and that is the most important result in this report:

| stage | pass rate | datasets passing all variants |
|---|---|---|
| as first measured | 95/306 = 31.0% | 27/102 |
| after fix 1 alone | 86/306 = 28.1% | 23/102 |
| after fixes 1 + 2 | **113/306 = 36.9%** | **34/102** |

---

## 1. The five defects, how each was found, and the fix

### 1.1 fastsim's 3′ boundary rule — `586a74a`

**How it failed.** The rule bounds an *injected* tail run: the base at index `L−t−1` must be non-G so the measured run is exactly `t`. The guard was `if t < read_len` — which is **true when t = 0**, and t = 0 is ~99.6% of reads. Every read *without* a tail artefact therefore had its last base forced non-G whenever biology put a G there.

**Why it hid.** The trailing-G histogram still matched the profile exactly (run = 0 on 99.58% of reads). Every artefact metric looked right. Only the per-position base composition diverged, at exactly one position.

**Measured on DRR000432** (read length 36, flat G curve ≈ 0.235):

| position | real G | sim G (before) | sim G (after) |
|---|---|---|---|
| 35 (last base) | 0.2353 | **0.0042** | 0.2352 |
| 34 | 0.2353 | 0.2328 | — |
| 0 | 0.2154 | 0.2173 | — |

A 56× shortfall on a value the profile supplies *directly* (`base.R1[35] = [0.2658, 0.2323, 0.2353, 0.2666]`). Every other position matched to within 0.006.

**The fix lowered the measured pass rate.** Four virome datasets had been *passing by accident*: the forced-non-G happened to correct a separate GC under-model in those libraries (sim GC 0.52 vs real 0.57). Fixing the boundary rule exposed that underlying gap — 9 virome rows flipped PASS→FAIL. Two defects had been partially cancelling.

### 1.2 The profile's duplicate rate — `6284829`, `ab7a081`

**How it failed.** `tools/fastcheck/src/profile.rs` computed the profile's `dup_rate` from `accs.first()` — **the first mate's accumulator only** — while `rich.rs` computes the report's `duplication_level` from the **merged** R1+R2 table. The two mates carry genuinely different duplicate rates, so the profile and the report disagreed *by construction*.

**Measured on DRR000466:**

| | reads | distinct | level |
|---|---|---|---|
| R1 only | 100,000 | 99,139 | **0.0086** ← what the profile stored |
| R2 only | 100,000 | 98,183 | 0.0182 |
| merged | 200,000 | 197,307 | **0.0135** ← what the report stored |

Only **15 sequences** are shared between the mates, so the merge was not collapsing across them — the mates simply differ. The closed loop then compared the simulator (fed the first-mate rate) against the report's merged rate, inflating `dup_ratio` by a median factor of **1.39** across 96 datasets and failing 144 of 306 comparisons.

**The fix** computes `dup_rate` over the merged fingerprint table, exactly as the report does. Verified: `profile.dup_rate == report.duplication_level` to full precision on DRR000466, DRR000978, DRR000432 and DRR425489 (ratio 1.000 each), against 1.570 / 1.234 / 1.491 / 1.980 before.

**A test had to be corrected too.** The regression test written for this fix expected a merged level of 0.125, but its input used `AAAA/GGGG/TTTT` in *both* mates, so the merge collapsed them into 4 distinct keys over 8 reads = 0.5. The production code was right; the test's arithmetic was wrong. The test now uses mate-disjoint sequences so it still proves the property it exists for.

### 1.3 `closed_loop.py` appended instead of overwriting — `ebcd370`

Re-running the loop silently duplicated rows (612 rows from a 306-row truth), because the script opened `--row` in append mode. Fixed to overwrite and refuse without `--force`. `fastp_compare.py` had the same pattern; `sweep.py` **keeps** its append deliberately, because concurrent chains legitimately share one rows file — documented at the write site.

### 1.4 fastcheck's k-mer counter — `a5005d5`

**How it failed.** `RepeatedCounter::insert` ran a full O(n) `retain()` eviction pass on **every insert** once the map exceeded its 1M cap. On real data the 12-mer space saturates (100k reads × ~140 windows against 4¹² = 16.7M keys), so the map sits at the cap permanently and every insert pays a ~1M-entry rehash: order 10¹³ operations.

**Symptom:** **40 minutes of CPU on one 100k-pair dataset with no output**, frozen RSS, ~200% CPU. Five datasets were spinning with 79 sweep pairs queued behind them — the campaign appeared frozen at 173/306 for 40 minutes.

**The fix** evicts lazily: one retain pass per `cap/4` keys past the last eviction. Amortized O(1) per insert, identical memory bound, identical counting semantics.

**Measured:** 40 min of spin → **2.48 s** on the exact dataset that stalled (DRR509062).

**Why this matters beyond the bug:** a simulated benchmark **structurally cannot** expose it, because the pathological case needs a real library's 12-mer saturation. It was found only by running the campaign on real data. It is written up as a section in the FastCheck paper.

### 1.5 The peak-RSS measurement — `5e60f49`

**How it failed.** `getrusage(RUSAGE_CHILDREN).ru_maxrss` is the maximum over **every reaped child the process ever spawned** and never resets. After one large child runs, every later child reports that same value regardless of its own footprint.

**Symptom:** all three tools in the fastp comparison reported an **identical 1209.8 MB** across 210 of 387 rows. The memory column was one process's historical maximum, not three per-run measurements.

**The fix** polls each child's `/proc/<pid>/status` VmHWM every 50 ms for its whole life and returns the last non-zero reading — a true per-run peak. Timings were never affected.

**The corrected result is far stronger than the bogus one:** fastprep-winner **34.3 MB** median vs fastp-default **1134.7 MB**, winning on 101/101 datasets (33×).

---

## 2. The remaining duplicate-rate gap — diagnosed, not yet fixed

> **Resolved (2026-10-16).** Both named causes are now addressed: the
> duplicate-identity defect was fixed in `781c396`, and the missing ambiguous
> base (`N`) emission is now implemented profile-driven (`docs/N-EMISSION.md`).
> **The 18 `dup_ratio` failures this section reports do not reproduce:** in the
> authoritative `closed_loop.tsv` the variant-1 and variant-2 rows are
> byte-identical in 100 of 102 datasets, and all 18 failures are in that
> duplicated v1/v2 pair; regenerated simulations pass the duplicate gate. So
> the residual failures below are a corpus artefact, not the emission gap this
> section names. The N gap was real and is closed; the failures were not the
> evidence for it. See `docs/N-EMISSION.md` §6.

122 rows still fail `dup_ratio`, and **the simulator is lower on all 122**, so this is a genuine model gap rather than residue from the fixed bugs. Two distinct causes, both verified by measurement.

### 2.1 The simulator emits no ambiguous bases at all

Measured: **0.00%** of simulated reads contain an `N`, against 4,772 of 200,000 in real DRR000466. `random_dna_bias` draws only ACGT, and the adapter-padding path that would emit `N` is unreachable because the profile's insert CDF sits above the read length.

### 2.2 Fresh sequencing errors break duplicate identity

A duplicate re-sequences a stored insert, but `qualify()` draws **fresh errors** for the copy. At mean Q ≈ 22.6 the per-base error rate is 0.0055, so a 150 bp read is error-free only ~43.5% of the time — most duplicates differ from their original somewhere, and fastcheck's exact sequence fingerprint therefore does not count them.

**Control measurement that isolates it:** a profile with `dup_rate = 0.05` measures 0.017, i.e. measured/p ≈ 0.35 — and that ratio is **independent of the ring-buffer size** (pool 1024 → 0.322, pool 1e6 → 0.346). The sim's `dup_multiplicity_max` is 3–4 against 2035 in the real data.

### 2.3 The N contribution is domain-specific — do not generalise it

An earlier analysis in this project claimed *"the real signal is N-padding, not PCR duplicates."* **That is true for one dataset and false for two:**

| dataset | duplicate reads containing N | genuine biological duplicates |
|---|---|---|
| bacterial_isolate/DRR000466 | 799 of 861 (**92.8%**) | 62 |
| generalised/DRR000432 | 492 (**8.7%**) | 5,164 |
| human_gut/DRR425489 | **0 (0.0%)** | 508 |

On DRR000466 the top "duplicates" are `NNNN…` ×276 and `GNNNN…` ×107 — padding, not biology. On DRR425489 **every** duplicate is a real sequence. Both terms exist and their mix varies by library, so a fix must address both.

### 2.4 Why no algebraic conversion was applied

An earlier hypothesis held that the sim's per-*pair* draw rate is measured over *reads* as p/2. **Tested and disproved** by emulating the exact emission algorithm: measured/p = 0.95–1.00 across p = 0.01/0.05/0.10, so there is no factor of 2 in the emission path.

With freshly drawn artefacts and errors there is no clean relation between a per-pair draw probability and the measured read-level rate. The fix belongs on the **emission side**: carry the artefact *and* the error realisation across a duplicate copy, so a copy is sequence-identical to its original — which is what a real optical duplicate actually is (same tile, same cycle). That is both more faithful and it removes the detection gap.

### 2.5 A caveat on the tolerance itself

**55 of the 122** failing rows have an *absolute* delta under 0.005 (sub-half-percent). A ratio gate of 2.0 punishes tiny absolute differences at low rates: at a real rate of 0.0011, a sim rate of 0.0005 fails the gate while differing by one duplicate read in a thousand. An absolute-delta alternative to the ratio gate is worth adding.

---

## 3. The remaining `g_frac` gap

150 rows fail `g_frac` at a median over-tolerance error of **0.0862** (p75 0.1143), concentrated in **metagenomics (48)** and **rna_seq (45)**. Only 13 are marginal. Worst cases: DRR014283 (0.1299) and DRR083295 (0.1261).

This is a genuine modelling gap in the per-position base composition, and unlike the duplicate gap it has **no single identified cause yet**. It is now the largest remaining blocker.

---

## 4. The harness defect that hid the true number — `f57f7c5`

`closed_loop.py` looked for the real QC report at the dataset root — the *old* shared-dir layout — but `run_campaign.sh` writes every artefact under per-variant dirs (`v0/`, `v1/`, `v2/`). On the Phase D tree that matched only the **23 pilot-era datasets**, so the closed loop reported 69–94 comparisons instead of 306 and a **52% pass rate instead of the true 31%**.

The wrong number was not a rounding difference: two whole domains (`rna_seq`, `metagenomics`) had **zero coverage**, and both gate almost entirely on `dup_ratio`/`g_frac`.

This is why the pilot-era verdict "SCALE WITH CAVEATS" (36/69) does not survive the full corpus — on the full corpus the loop closes on **36.9%**.

---

## 5. What these fixes do *not* affect

The sweep results — `winners.tsv` (101 winning configs, 16 distinct of the 55-grid; config 6 wins 33 datasets) and `fastp-comparison.tsv` (303 rows at gz6, 4 threads) — are **unaffected**. The sweep scores fastprep against the simulator's ground-truth *trim windows*, not against base composition.

Those numbers stand and were verified independently against the raw TSVs:

| | fastprep-winner | fastp-default | fastp-best |
|---|---|---|---|
| Time (median) | **0.20 s** | 0.50 s | 0.50 s |
| Memory (median) | **34.3 MB** | 1134.7 MB | 1134.9 MB |
| Keeps reads | 90.8% | 97.2% | 97.7% |

Faster on **100/101** datasets (median 2.00×), less memory on **101/101** (33×).

**The honest counterweight:** fastprep keeps fewer reads on 98/101 (median −5.1 pp) because **81 of 101** winners carry a *relative* length floor above 0.5L while fastp's default is an *absolute* 30 bp. That is a difference of contracts, not of trimming quality, and it is stated as such in the FastPrep paper rather than buried.

Also stated plainly in the FastPrep paper: the **shipped default wins on only 1 of 101** datasets (DRR455998), against the "90/90" the earlier corpus campaign suggested. That comparison was measured against the 0.4.x default, before D-18 set `--min-len` to 0.67L.

---

## 6. Corrections made during the analysis

Recorded because a report that only shows the final numbers hides how they were reached:

- A **p/2 duplicate mechanism** was hypothesised, then disproved by emulating the exact emission algorithm (measured/p = 0.95–1.00). A subagent was steered off that path before it shipped a conversion built on it.
- The **g_frac failure median** was briefed as 0.0971; it does not reproduce. The true value over the 150 failing rows is **0.0862**. The earlier figure had been computed over the pre-dup-fix file.
- A **ring-buffer explanation** of the duplicate gap was written into both papers and then disproved by the pool-size control (1024 → 0.322, 1e6 → 0.346). Both papers were corrected.
- A **regression test** for the merged-dup fix had wrong arithmetic in its expected value; the production code was correct.
- A subagent reported a pass rate of **131/306** and claimed the verified 113 was polluted by the append bug. Verified false: the authoritative file was unchanged (md5 `304abe02…`, mtime 09:40:06) and the run had deleted the output before executing. **113/306 stands.**

---

## 7. Suite status and provenance

**461 tests, 0 failures** across all six tools (fastcheck 107, fastdedup 104, fastprep 128, fastsim 58, fastmaps 43, fasthelper 21).

All work committed locally; **nothing pushed**, per standing instruction. 81 commits ahead of origin.

Key commits: `586a74a` (3′ rule) · `6284829`, `ab7a081` (profile dup) · `ebcd370` (append) · `a5005d5` (k-mer) · `5e60f49` (RSS) · `f57f7c5` (harness) · `c3ad7a0` (closed-loop re-run) · `6a77def`, `20b9d1d` (analysis docs).
