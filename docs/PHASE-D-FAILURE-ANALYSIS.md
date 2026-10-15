# FastSuite Phase D — failure analysis and fixes

**Corpus:** 102 real ENA datasets × 55 configs × 3 variants = 306 (dataset, variant) comparisons.
**Method:** closed loop — fastcheck fits a profile from the real library, fastsim simulates from that profile, fastcheck re-measures the simulation, and the two measurements are compared against tolerance.

---

## 1. What was broken, and what fixing it did

Three real bugs, each found by measurement rather than inspection. The pass rate did not move monotonically, and that is the most important result in this report.

| stage | pass rate | datasets passing all variants |
|---|---|---|
| as first measured | 95/306 = 31.0% | 27/102 |
| after fix 1 alone | 86/306 = 28.1% | 23/102 |
| after fixes 1 + 2 | **113/306 = 36.9%** | **34/102** |

### Fix 1 — fastsim's 3′ boundary rule (`586a74a`)

The rule bounds an *injected* tail run: the base at index `L−t−1` must be non-G so the measured run is exactly `t`. The guard was `if t < read_len`, which is **true when t = 0** — and t = 0 is ~99.6% of reads. Every read without a tail artefact had its last base forced non-G whenever biology put a G there.

Measured on DRR000432 (read length 36, flat G curve ≈ 0.235):

| position | real G | sim G (before) | sim G (after) |
|---|---|---|---|
| 35 (last base) | 0.2353 | **0.0042** | 0.2352 |
| 34 | 0.2353 | 0.2328 | — |
| 0 | 0.2154 | 0.2173 | — |

A 56× shortfall on a value the profile supplies directly. Every other position matched to within 0.006.

**This fix lowered the measured pass rate.** Four virome datasets had been *passing by accident*: the forced-non-G happened to correct a separate GC under-model in those libraries (sim GC 0.52 vs real 0.57). Fixing the boundary rule exposed that underlying gap, flipping 9 rows PASS→FAIL. Two bugs had been partially cancelling.

### Fix 2 — the profile's duplicate rate (`6284829`, `ab7a081`)

`tools/fastcheck/src/profile.rs` computed the profile's `dup_rate` from `accs.first()` — **the first mate's accumulator only** — while the report's `duplication_level` uses the **merged** R1+R2 table. The two mates carry genuinely different duplicate rates, so the profile and the report disagreed by construction.

Measured on DRR000466:

| | reads | distinct | level |
|---|---|---|---|
| R1 only | 100,000 | 99,139 | **0.0086** ← what the profile stored |
| R2 only | 100,000 | 98,183 | 0.0182 |
| merged | 200,000 | 197,307 | **0.0135** ← what the report stored |

Only **15 sequences** are shared between the mates, so the merge was not collapsing across them — the mates simply differ. The closed loop compared the simulator (fed the first-mate rate) against the report's merged rate, inflating `dup_ratio` by a median factor of **1.39** across 96 datasets and failing 144 of 306 comparisons.

After the fix, `profile.dup_rate == report.duplication_level` to full precision (ratio 1.000 on DRR000466, DRR000978, DRR000432, DRR425489; previously 1.570 / 1.234 / 1.491 / 1.980).

### Fix 3 — `closed_loop.py` appended instead of overwriting (`ebcd370`)

Re-running the loop silently duplicated rows (612 rows from a 306-row truth). Fixed to overwrite; `fastp_compare.py` had the same pattern.

---

## 2. Where the campaign stands now

**PASS 113/306 = 36.9%**, 34/102 datasets passing every variant.

| domain | pass rate |
|---|---|
| human_gut | 90.2% (46/51) |
| virome | 49.0% (25/51) |
| bacterial_isolate | 39.2% (20/51) |
| generalised | 25.5% (13/51) |
| rna_seq | 11.8% (6/51) |
| metagenomics | 5.9% (3/51) |

Medians: `g_frac_max_diff` 0.0157 (tol 0.02), `mean_q_diff` 0.1494 (tol 0.50), `gc_diff` 0.0004 (tol 0.01), `dup_ratio` 1.8205 (tol 2.0).

---

## 3. The remaining duplicate-rate gap — diagnosed, not yet fixed

122 rows still fail `dup_ratio`. **The simulator is lower on all 122**, so this is a genuine model gap, not residue from the fixed bugs. Isolating it on DRR000466 with a freshly fitted profile:

| step | value | ratio |
|---|---|---|
| profile rate (now exactly equals the real rate) | 0.0135 | — |
| duplicates the sim actually drew (ground truth) | 0.01073 | **0.795** |
| duplicates fastcheck detected in the output | 0.00905 | **0.851** |

Two distinct shortfalls, multiplying to the observed 0.670 (0.00905 / 0.0135):

**(a) The draw is short by ~20%.** The sim draws 1073 duplicate pairs when told 0.0135 → expected 1350.

**(b) ~15% of the duplicates the sim draws are not detected.** The ground truth flags 1073 duplicate pairs, but R1 contains only 913 exact-duplicate reads. A duplicate re-sequences a stored insert, and the artefact is deliberately drawn *fresh* for every pair (`main.rs`: "Drawn fresh even for duplicates: dropout is a per-cluster failure, so a copy gets its own polyG runs"). When a copy's polyG run differs from its original's, the two R1 sequences differ and fastcheck's exact fingerprint does not count them as duplicates.

That design is defensible for real optical duplicates — which *do* share their artefact, coming from the same tile and cycle — but it makes the simulator's duplicates systematically harder to detect than real ones. The fix is to carry the artefact across a duplicate copy rather than redrawing it, which is both more faithful and removes the detection gap.

**One caveat on the tolerance itself:** 55 of the 122 failing rows have an *absolute* delta under 0.005 (sub-half-percent). A ratio gate of 2.0 punishes tiny absolute differences at low rates — at a real rate of 0.0011, a sim rate of 0.0005 fails the gate while differing by one duplicate read in a thousand. An absolute-delta alternative to the ratio gate is worth adding.

---

## 4. The remaining `g_frac` gap

150 rows fail `g_frac` (median 0.0971, p75 0.1143), concentrated in **metagenomics (48)** and **rna_seq (45)**. Only 13 are marginal. Worst cases: DRR014283 (0.1299) and DRR083295 (0.1261).

This is a genuine modelling gap in the per-position base composition. It is the largest remaining blocker after the duplicate rate, and unlike the dup gap it has no single identified cause yet.

---

## 5. Other defects found and fixed during the campaign

| defect | symptom | fix |
|---|---|---|
| fastcheck's k-mer counter (`a5005d5`) | O(n) eviction on *every* insert once the map hit its 1M cap; **40 minutes of CPU on one 100k-pair dataset with no output**, stalling 79 sweep pairs | lazy amortized eviction → **2.48 s** |
| peak-RSS measurement (`5e60f49`) | `RUSAGE_CHILDREN.ru_maxrss` is a never-resetting high-water mark; three tools reported an identical 1209.8 MB across 210 rows | per-child `/proc/<pid>/status` VmHWM polling → true 34 MB vs 1135 MB |
| closed-loop harness coverage (`f57f7c5`) | looked for the QC report at the dataset root, matching only the 23 pilot-era datasets; reported a 52% pass rate | per-variant fallback → the true 31% |
| fastsim's 3′ rule (above) | last base forced non-G on ~100% of reads | `t > 0` guard |
| profile dup rate (above) | first mate instead of the merged accumulator | merged fingerprint table |

The k-mer defect is worth emphasising: a simulated benchmark **structurally cannot** expose it, because the pathological case needs a real library's 12-mer saturation. It was found only by running the campaign on real data.

---

## 6. Unaffected by these fixes

The sweep results — `winners.tsv` (101 winning configs) and `fastp-comparison.tsv` (303 rows at gz6, 4 threads) — are **not** affected. The sweep scores fastprep against the simulator's ground-truth *trim windows*, not against base composition. Those numbers stand: fastprep-winner faster on 100/101 datasets (median 2.00×), less memory on 101/101 (34.3 MB vs 1134.7 MB, 33×).

---

## 7. Suite status

**461 tests, 0 failures** across all six tools. All work committed locally; nothing pushed.
