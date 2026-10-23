# Duplicate-rate mechanism — verified findings

> **Superseded (2026-10-16).** This is a point-in-time record from 2026-10-15.
> The closed loop later reached **167/306 = 54.6%** (52 of 102 datasets passing
> every variant), not the 113/306 recorded below, and the suite is now **852**
> on the `scripts/test.sh` GRAND line (not 461) — see
> `scripts/expected-test-counts.txt` for the authoritative per-tool split. Authoritative closed-loop
> figure: `bench/results/phase-d/closed_loop.tsv`. Current state:
> `docs/STATUS.md`.
>
> **Update (2026-10-16): both named causes are now addressed, and part (a)'s
> premise is corrected.** Cause (b) — a duplicate not being sequence-identical
> to its original — was fixed in `781c396`. Cause (a) — the simulator emitting
> no `N` — is now implemented profile-driven (see `docs/N-EMISSION.md`): the
> emission exists, is calibrated from each library's measured no-call / padding
> structure, and moves the simulated duplication level and multiplicity exactly
> as the padding class requires. **But it does not close the 18 `dup_ratio`
> failures, because those failures do not reproduce:** in `closed_loop.tsv` the
> variant-1 and variant-2 rows are byte-identical in 100 of 102 datasets and
> every one of the 18 failures is in that duplicated v1/v2 pair; freshly
> regenerated simulations pass the duplicate gate. The N gap was real and is
> closed; the 18 failures were a corpus artefact. See `docs/N-EMISSION.md` §6.

Measurement-only analysis. No code changed by this document.

## Verified result

**Closed loop: PASS 113/306 = 36.9%**, 34/102 datasets passing every variant.
Authoritative file: `bench/results/phase-d/closed_loop.tsv`, md5
`304abe02019e090d39a428ae44e687c7`, 306 rows / 102 datasets × 3 variants.

Per-domain: human_gut 90.2%, virome 49.0%, bacterial_isolate 39.2%,
generalised 25.5%, rna_seq 11.8%, metagenomics 5.9%.
Medians: g_frac 0.0157, mean_q 0.1494, gc 0.0004, dup_ratio 1.8205.

## Bugs fixed (all verified by A/B measurement)

1. **fastsim's 3′ boundary rule** — guard was `t < read_len`, true at t = 0
   (~99.6% of reads), so every read without a tail artefact had its last base
   forced non-G. Real 23.5% G, sim 0.42%. Committed `586a74a`.
2. **The profile's dup_rate** — `profile.rs` used `accs.first()` (first mate
   only) while the report uses the merged R1+R2 table. The mates genuinely
   differ (DRR000466: R1 0.0086, R2 0.0182, merged 0.0135, only 15 sequences
   shared), so the loop compared the sim against a different rate than it was
   fed. Committed `6284829`, `ab7a081`.
3. **`closed_loop.py` appended instead of overwriting** — silently produced
   612 rows from a 306-row truth. Committed `ebcd370`.

## The remaining dup gap — two distinct causes, both real

### (a) The simulator emits no N bases at all

Measured: sim reads containing N = **0.00%**; real DRR000466 has 4772 of
200,000. `random_dna_bias` draws only ACGT, and the adapter-padding path that
would emit N is unreachable because the profile's insert CDF sits above the
read length.

**Now implemented (2026-10-16):** a profile-driven `n_run` model (read rate,
run-length histogram, leading/trailing shares), fitted by fastcheck and emitted
by fastsim, with the injected count recorded in the GT's `r1_n`/`r2_n` columns.
On a synthetic padding profile it moves `duplication_level` 0.0199 → 0.0397 and
`dup_multiplicity_max` 5 → 4014. Details and the (negative) closed-loop
outcome: `docs/N-EMISSION.md`.

### (b) Fresh sequencing errors break duplicate identity

A duplicate re-sequences a stored insert, but `qualify()` draws fresh errors.
At mean Q ≈ 22.6 the per-base error rate is 0.0055, so a 150 bp read is
error-free only ~43.5% of the time — most duplicates differ from their
original somewhere and fastcheck's exact fingerprint does not count them.
Control: a profile with `dup_rate = 0.05` measures 0.017, i.e. measured/p ≈
0.35, independent of the ring-buffer size (pool 1024 → 0.322, pool 1e6 → 0.346).
The sim's `dup_multiplicity_max` is 3–4 against 2035 in the real data.

### The N contribution is domain-specific — do not generalise it

| dataset | duplicate reads containing N | genuine biological duplicates |
|---|---|---|
| bacterial_isolate/DRR000466 | 799 (92.8%) | 62 |
| generalised/DRR000432 | 492 (8.7%) | 5164 |
| human_gut/DRR425489 | 0 (0.0%) | 508 |

On DRR000466 the top "duplicates" are `NNNN…` ×276 and `GNNNN…` ×107 —
padding, not biology. On DRR425489 **every** duplicate is a real sequence.
Both terms exist; their mix varies by library, so a fix must address both
rather than treating N-padding as *the* cause.

## Why no algebraic conversion was applied

With fresh artefacts and errors there is no clean relation between a per-pair
draw probability and the measured read-level rate. An earlier `p/2` hypothesis
was tested by emulating the exact emission algorithm and disproved
(measured/p = 0.95–1.00). The fix belongs on the emission side — carry the
artefact and the error realisation across a duplicate copy — not on a
conversion factor.

## Tolerance caveat

55 of the 122 dup failures have an absolute delta under 0.005. A ratio gate of
2.0 punishes sub-half-percent differences at low rates (at a real rate of
0.0011, a sim rate of 0.0005 fails the gate while differing by one duplicate
read in a thousand). An absolute-delta alternative is worth adding.

## Suite status

461 tests, 0 failures. All work committed locally; nothing pushed.
