# The head_dropout axis: the sweep could not express 0, and 0 wins everywhere

**Question.** The shipped `fastprep` default at the time (4; it was briefly 0
after `3d76a31` and is 4 again since 0.4.1, D-16) trims a fixed 5' dropout of
`head_dropout=4` bases before the polyX scan. The config sweep never tested
whether that default was right — so this phase measured it. Two bugs guaranteed
the answer was invisible:

1. `effective_head_dropout()` in `bench/scripts/sweep.py` floored the value at
   1 (`max(1, min(hd, min_run - 1))`), so the grid could **never** express
   `head_dropout 0` even when a config authored it. The floor is now **0**
   (the `min(hd, min_run - 1)` upper clamp stays — fastprep refuses
   `head_dropout >= min_run`, and it accepts 0, which disables the 5' hard
   trim).
2. `build_argv()` never emitted `--head-dropout`, so every swept config ran at
   fastprep's CLI default 4 regardless of what the grid authored. The flag is
   now emitted explicitly on every run and routed through
   `effective_head_dropout()` — the same function the results row records, so
   the emitted argv and the recorded `head_dropout` column cannot disagree.

The grid now authors the axis as 0, 1, 4, 6 (0/1 as single-axis moves off the
shipped default at `min_run=12`; 4/6 at lower `min_run`), 57 configs total.
`bench/scripts/test_sweep_scorer.py` covers: 0 is expressible, 1 passes, the
cap still applies (6 at `min_run=2` → 1), 0 survives the cap, every grid
config emits the flag, and the emitted value equals the recorded value.

## What was measured

Measured 2026-10-16, before D-18: every run used the winning config's own
`min_len` (the then-shipped `0.9L` floor on most datasets; the short-read
configs carry explicit 30/50 bp or `0.5882L` floors) and the pre-D-18 mercy
clause, so the tables below are measurements of that build, not of the `0.67L`
default.

Not a re-run of the 57-config grid. For each of the **101 datasets with a
winning config** (`bench/results/phase-d/winners.tsv`; 102 corpus dirs exist —
`bacterial_isolate/DRR000466` has sim data but no winner row), the winning
config was re-run with **only `head_dropout` overridden to 0, 1, 2, 3, 4** and
scored the same way the sweep scores: PXMT trim decisions against the sim
ground truth, per mate, exact / overtrim / undertrim at tolerance 0 and 2.
All three sim variants (v0/v1/v2) exist for every dataset and were scored —
1515 fastprep runs, 3030 rows. Harness: `bench/scripts/head_dropout_axis.py`;
data: `bench/results/phase-d/head-dropout-axis.tsv`.

Since the pre-fix sweep never emitted the flag, **the current winner's value
is 4 on every dataset** — that is the value the winners were selected at.

Faithfulness anchors:

- hd=4 reproduces the established `overtrim-head-dropout.tsv` baseline row for
  row on 98/101 datasets (e.g. DRR000978 v0: 157138 / 28895 / 97 exactly).
- The 3 divergences are an argv-fidelity artifact of the *older* harness
  (`overtrim_real.py` passed the authored `min_len`, this run passes the
  recorded `min_len_effective` — what the winning run actually used):
  `DRR000432` (36 bp reads) and `SRR19649201` (40 bp) carry `min_len=50`,
  which the sweep clamps to `"1"`; the old harness passed `"50"`, an
  impossible floor that dropped 20 pp more pairs. `DRR000616` (51 bp) differs
  by one mate on the `30` vs `0.5882L` floor boundary. This run is the
  faithful one; the established headline (75.6 / 23.4 / 1.0 at hd=4) comes
  from that older file.

## Corpus-wide result (all 3 sim variants pooled, 56.96 M scored mates)

tolerance 0 (% of kept mates):

| head_dropout | exact | overtrim | undertrim | total error | dropped |
|---|---|---|---|---|---|
| 0 | **92.56 %** | **0.07 %** | 7.37 % | **7.44 %** | 3,634,937 |
| 1 | 80.27 % | 18.44 % | 1.29 % | 19.73 % | 3,638,828 |
| 2 | 76.39 % | 22.54 % | 1.07 % | 23.61 % | 3,641,460 |
| 3 | 75.53 % | 23.44 % | 1.03 % | 24.47 % | 3,642,980 |
| 4 (then-shipped) | 75.36 % | 23.64 % | 1.00 % | 24.64 % | 3,643,718 |

tolerance 2:

| head_dropout | exact | overtrim | undertrim | total error |
|---|---|---|---|---|
| 0 | **99.69 %** | 0.02 % | 0.28 % | **0.31 %** |
| 1 | 99.69 % | 0.05 % | 0.25 % | 0.31 % |
| 2 | 99.69 % | 0.06 % | 0.25 % | 0.31 % |
| 3 | 98.82 % | 0.97 % | 0.21 % | 1.18 % |
| 4 (then-shipped) | 98.65 % | 1.17 % | 0.18 % | 1.35 % |

Per sim variant the numbers are within 0.05 pp of each other (v0 tol 0:
hd=0 92.60 / 0.07 / 7.33, hd=4 75.35 / 23.64 / 1.00 — v1/v2 identical to two
decimals), so nothing here is a v0 seed artifact.

## The trade: which error is cheaper?

hd=0 trades overtrim for undertrim. At tolerance 0 the trade is **3.3x
cheaper**: total error 7.44 % vs 24.64 % (17.2 pp). It is not close even
before weighting: hd=4 destroys 23.64 % of kept mate windows (overtrim kills
real biology — the sweep's own score penalises it 2x) while hd=0 leaves
7.37 % undertouched. Weighting overtrim at 2x as the sweep score does:
hd=0 7.51 % vs hd=4 48.3 %. At tolerance 2 total error is 0.31 % vs 1.35 %.
hd=0 also drops fewer pairs (3.635 M vs 3.644 M mates), so it does not buy
its accuracy by discarding data.

## Per-dataset winners (max exact @ tol 0, ties → min total error → lower hd)

| best head_dropout | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| datasets (tol 0) | **100** | 1 | 0 | 0 | **0** |
| datasets (tol 2) | 88 | 7 | 4 | 2 | 0 |

- **The best head_dropout differs from the current winner's value (4) on
  101/101 datasets.** hd=4 is never the strict max-exact and never even ties
  for it; hd=0 has at least as many exact trims as hd=4 on every dataset
  (strictly more on all but the one below). The result is the same if the
  criterion is minimal total error (hd=0 on 100, hd=1 on 1) or if each sim
  variant is counted separately (100/1 in v0, v1 and v2 alike).
- **hd=0 wins 100, hd=1 wins 1, hd=4 wins 0** — so yes, low dropout wins far
  more often than the shipped 4.
- The single hd=1 dataset is `bacterial_isolate/DRR014283` (250 bp reads,
  config 22): hd=0 leaves its long 5' artefacts undertouched (459,840 exact /
  132,620 under) where hd=1 removes them exactly (464,753 exact / 127,381
  over / 358 under) — hd=1 wins on exact *and* total error (127,739 vs
  132,652). This is the one regime (long reads, artefact exactly 1 base deep)
  where the shipped rule helps at all.

## What propagates

- **The default question is settled by the corpus numbers, and now also by
  per-dataset selection:** `head_dropout 0` wins on this corpus (100/101
  datasets, +17.2 pp exact at tol 0, 3.3x lower total error) *at a 1:1 error
  price*. A paper reporting that price should record `head_dropout=0` for
  fastprep; the default it shipped at the time (4) is badly wrong for these
  artefact regimes — it eats 23.4 % of kept windows. **The shipped default did
  not end up here:** `0` held only from `3d76a31` to 0.4.1; D-16 re-priced
  undertrim at 7:1 and put `4` back, which is the default today.
- The closed loop must treat `head_dropout` as a real axis from 0:
  `effective_head_dropout` now floors at 0, and any future sweep grid that
  authors 0 will actually run it. hd=1 is worth keeping as a candidate for
  long-read regimes (DRR014283), hd=2/3 only ever win tolerance-2 ties and
  are dominated at tolerance 0.
- One caveat worth carrying into the paper: hd=0's residual error is
  *undertrim* (7.37 % at tol 0, ~0.25-0.28 % at tol 2). If a downstream task
  needs hard 5' cleanliness more than raw window accuracy, hd=1 is the
  compromise point — hd=2/3/4 buy almost nothing over hd=1 (80.3 % → 75.4 %
  exact) for the same overtrim regime.

## Files

- `bench/scripts/sweep.py` — floor 0 in `effective_head_dropout`, grid
  authors 0/1/4/6, `--head-dropout` emitted via the same function the row
  records.
- `bench/scripts/test_sweep_scorer.py` — coverage for the above.
- `bench/scripts/head_dropout_axis.py` — the measurement harness (winning
  config + hd override, per-mate PXMT vs GT, tol 0/2, v0/v1/v2, sharded).
- `bench/results/phase-d/head-dropout-axis.tsv` — 3030 measured rows
  (101 datasets x 3 variants x 5 values x 2 tolerances).
