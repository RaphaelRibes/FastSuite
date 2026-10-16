# The trim-accuracy discrepancy: 66.50% vs 75.59% — root cause

**Question.** `bench/results/phase-d/trim-accuracy.tsv` (published) scores the
Phase D sweep's winning fastprep config at **66.50% exact / 29.5% overtrim /
4.0% undertrim** (tol 0). `bench/results/phase-d/overtrim-vs-fastp.tsv` (fresh
harness) measures the same winner configs at **75.59% / 23.4% / 1.0%**. fastp
reproduces to the digit in both (92.19% default, 84.06% tuned). Why?

**Answer.** The published fastprep rows are not a measurement of fastprep at
all: they score trim decisions recorded on **sim generation 1** (the Phase D
sweep, 2026-10-14/15) against ground truth from **sim generation 2** (the
regenerated sims of 2026-10-15 15:54–15:58). The two generations are different
random realisations of the same fitted profile. fastp reproduces because both
tables run fastp fresh on the *current* sim and score against the *current* GT,
so both are generation-consistent; the published fastprep row is the only row
in either table built from stale artifacts.

**Authoritative table: `overtrim-vs-fastp.tsv`.** The published 66.50% must not
be cited as fastprep accuracy. `trim-accuracy.tsv` is left in place untouched;
see *Fix* below for the repair.

---

## 1. Both tables reproduce exactly from their sources

`bench/scripts/trim_discrepancy_probe.py` (run on the desktop, corpus under
`/mnt/Fast/fastsuite-scratch`) reconstructs each row from its actual input. On
31 datasets (spanning the whole gap distribution), every checked row
reproduces to the bucket — 62 unique rows per table (31 datasets × 2
tolerances), 124 comparisons, 0 mismatches:

| reconstruction | source | matches |
|---|---|---|
| A | stored sweep PXMT `work/<dom>/<acc>/v0/sweep.v0/c<cid>.pxmt` scored vs current `sim.v0.gt.tsv` | published `trim-accuracy.tsv` rows, 62/62 |
| B | fresh run of the current `fastprep` binary with `overtrim_vs_fastp.py`'s argv, scored vs current GT | fresh `overtrim-vs-fastp.tsv` rows, 62/62 |

So the discrepancy is fully localized: the two tables use the same GT, the same
scoring semantics and (except §4) the same flags — but the *published* fastprep
row is computed from a PXMT produced before the corpus on disk was regenerated.

Worked example, `bacterial_isolate/DRR014283` (winner config 22, `min_len=0.9`
so the flags are identical in both harnesses):

```
A (stored PXMT vs current GT)  tol0 exact= 92085 over=87834 under=20039 dropped=  42   <- published row
B (fresh run vs current GT)    tol0 exact=145806 over=51401 under=  113 dropped=2680   <- fresh row
A tol2 exact=198525  B tol2 exact=195831
```

## 2. The corpus changed between the two measurements

Timeline (desktop, mtimes + git):

| when | what |
|---|---|
| 10-13 01:36 | last change to `tools/fastprep` source (`029b41a`) |
| 10-14 23:17 – 10-15 01:06 | Phase D sweep writes `sweep.v0/c*.pxmt` (the campaign) |
| 10-15 01:30 | `winners.tsv` committed (`c76b9b4`) |
| 10-15 09:29 | fastsim fix `324506d` committed (3′ boundary rule needs the `t>0` guard) |
| 10-15 15:53 | `fastsim` rebuilt; **15:54–15:58 all `sim.v0` FASTQ/GT regenerated** |
| 10-15 16:55 | `trimwork/` fastp outputs written (post-regen) |
| 10-15 16:59 / 19:38 | `trim-accuracy.tsv` committed |
| 10-16 13:07 | `fastprep` rebuilt from unchanged source; fresh harness runs |

In **every** dataset/variant directory checked (102/102), `sweep.v0/c0.pxmt`
predates `sim.v0.gt.tsv`: the stored PXMT is from the pre-regen generation.

The regeneration is not a no-op. fastsim draws everything (artefacts, biology,
quality) from one shared seeded `ChaCha8Rng` (`tools/fastsim/src/main.rs`), and
`324506d` removes one RNG draw per `t == 0` read whose last base is G (the old
code forced that base non-G). Removing draws shifts the stream, so the
regenerated sim is a **different random realisation of the same profile**:
read lengths unchanged, artefact lengths and read content redrawn.

Per-read evidence (probe, 31 datasets):

* Stale vs fresh PXMT: mate cuts identical on only 29–74% of reads; where they
  differ the deltas are concentrated at ±1–2 bp (small-integer artefact
  lengths of two independent draws from the same profile).
* **Content check:** the bases a stale `t5>0` cut claims to remove are all G in
  the *current* read for only 14–65% of stale cuts, while fresh cuts are
  content-consistent to ~100% (e.g. DRR014283: 48065/109702 vs 95244/95247).
  Worked example, DRR014283 pair 39689 mate 2: stale cut 1 bp, but the current
  read has 0 leading G — the read under the stale PXMT is a different draw.
* Read lengths match everywhere (fixed `read_len` per dataset), which is why
  the mismatch is invisible at the file level.

Why tol 2 hides it: 76–100% of GT mates are `(h,t) = (0,0)` and artefact
lengths are small integers, so two independent realisations agree within ±2 by
chance on ~99% of mates (published tol-2: 98.42% exact) while exact agreement
at tol 0 collapses (66.50%). The published 29.5% overtrim + 4.0% undertrim is
the ±1–2 bp spread of two independent draws, not trimming behaviour.

## 3. Why fastp reproduces exactly

Both harnesses run fastp fresh on the current sim and score against the
current GT (the `trimwork/` cache was written 16:55, after the 15:54 regen).
All 404 fastp rows (2 modes × 2 tolerances × 101 datasets) are bucket-identical
between the two tables. The asymmetry is structural: fastprep's published row
is the only row in either table derived from pre-regen artifacts.

## 4. Secondary defect found: `--min-len` semantics differ between harnesses

`sweep.py` emits `--min-len` via `min_len_arg(scale_min_len(cfg, read_len))`:
the grid's integers are base counts authored against a 150 bp reference and
are converted to a *fraction of the actual read length*, clamped to ≤ 1.0
(the value is recorded as `min_len_effective` in `winners.tsv`). The fresh
harness's `_minlen` passes the literal integer as a base count.

24/101 winners carry integer `min_len`. For 22 of them the literal count
happens to resolve to the same threshold, but on the two clamped cases the two
harnesses ran genuinely different configurations:

| dataset | as-swept (`--min-len 1`) | fresh harness (`--min-len 50`) |
|---|---|---|
| `virome/SRR19649201` | 75.9% exact, 95,052/200,000 mates dropped | 97.7% exact, 129,626 dropped |
| `generalised/DRR000432` | 78.0% exact, 16,920 dropped | 99.7% exact, 56,820 dropped |

Re-aggregated with as-swept semantics on those datasets, the fresh table's
headline becomes **75.35%** (vs 75.59%) — the tool comparison is unchanged in
direction and magnitude, but `overtrim_vs_fastp.py`'s `_minlen` should
reproduce `scale_min_len` if "the winner config" is meant literally.

## 5. Ruled out

* **(a) different fastprep binary** — no. The desktop's `tools/fastprep` source
  last changed 10-13 01:36, before both runs; both harnesses ran the same
  binary path. The behavioural difference tracks the sim, not the binary
  (fresh cuts are content-consistent with the current reads to ~100%).
* **(b) flags** — only `--min-len` (§4). `--head-dropout`: the campaign's
  `sweep.py` never emitted the flag (0 occurrences in the desktop's copy), so
  fastprep's CLI default applied; the fresh harness passes `4` explicitly,
  which is the same value for all winner configs (`min_run ≥ 8`, cap is a
  no-op). `--gz-level 6` and all other flags match.
* **(c) GT/sim change** — confirmed, this is the cause (§2).
* **(d) scoring logic** — no. `trim_accuracy.py::score_pxmt` and
  `overtrim_vs_fastp.py::score_fastprep` have identical bucket semantics
  (exact: `|t5-h|≤tol ∧ |t3-t|≤tol`; overtrim: `t5>h ∨ t3>t`; else undertrim;
  `status != 0` → dropped; PXMT mate byte +1 normalisation identical). Each
  scorer reproduces its own table to the bucket from the same PXMT/GT.

## 6. Why the stale PXMT survived (harness mechanics)

* `chain.sh` is idempotent: it only simulates when `sim.v0.gt.tsv` is missing,
  and `sweep.py` skips any config whose `c<cid>.row` exists ("resumable"). The
  post-regen chain runs therefore **replayed cached rows and never re-ran
  fastprep on the new sim** — the stale `sweep.v0/*.pxmt` were never refreshed.
* The recurring `chain: sweep failed on <acc> v0` (in every `regen.v0.log`) is
  spurious: `sweep.py` exits `0 only if n_done == args.configs`, and the chain
  passes `CONFIGS=60` while the grid has 55 configs — a full 55/55 sweep still
  exits 1.
* `trim_accuracy.py` then scored `sweep.v0/c<cid>.pxmt` against the new GT with
  no staleness check (no mtime/hash comparison against the sim).

Note also that `winners.tsv` (winner selection) is itself from the old
generation's sweep rows; both tables inherit this, so both measure "the
old-generation winner config" and the fastp-vs-fastprep comparison stands.

## 7. Fix (applied in `73c2385`)

`trim-accuracy.tsv` is **proven corrupt as a fastprep measurement**: its
fastprep-winner rows pair generation-1 trim decisions with generation-2 ground
truth. Its fastp rows are valid.

All four repairs below landed in commit `73c2385` ("bench: fix the three
trim-accuracy defects; re-run the four-tool comparison"). `trim-accuracy.tsv`
is left in place untouched as superseded; the authoritative table is
`overtrim-vs-fastp.tsv`, regenerated over all 101 datasets at both tolerances.
Its corrected tolerance-0 headline: **shipped default 75.35%** (was 75.59%
before the `_minlen` repair), **fastprep-hd0 92.60%**, **fastp-default 92.19%**,
**fastp-tuned 84.06%**.

1. **Staleness guard — applied.** `trim_accuracy.py` now refuses to score a
   PXMT older than the sim it is scored against (mtime guard), so a regenerated
   corpus can no longer be scored against stored trim decisions from a
   different draw.
2. **Resume-key identity — applied.** `sweep.py`'s resume cache is keyed on the
   inputs' identity (path + size + mtime of `r1`/`r2`/`gt`) instead of
   `c<cid>.row` existence alone, and invalidates cached rows when the corpus
   changes — a regenerated corpus no longer silently replays stale metrics.
3. **Exit code — applied.** `sweep.py` now exits 0 on a full sweep by comparing
   against `min(--configs, len(GRID))`, so the chain's `CONFIGS=60` over a
   57-config grid no longer reports `sweep failed`.
4. **`_minlen` semantics — applied.** `overtrim_vs_fastp.py::_minlen` reads
   `min_len_effective` from `winners.tsv` instead of passing the grid's raw
   integer as a base count, so it runs the same config the sweep ran (this is
   what changed the tol-0 headline from 75.59% to 75.35%: it affected the two
   clamped winners, `virome/SRR19649201` and `generalised/DRR000432`).

Reproduce: `python3 bench/scripts/trim_discrepancy_probe.py <domain>/<acc> ...`
on the desktop (paths and binary locations are module constants).
