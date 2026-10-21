# fastprep — undertrimting mechanisms and ranked fixes (raw dossier, 2026-10-21)

Evidence artifacts (all on the built release binary, tree `30e8fb3` + `af268d4`):
`bench/results/undertrim-ab.json` (two-axis corpus sweep),
`undertrim-bands.json` (residual by overhang band),
`undertrim-probe.json` (arms A-D: ladder / errors / contamination / clean),
`undertrim-probe-pe.json` (paired geometries x pe-overlap modes),
`undertrim-probe-pe-guard.json` (decisive random-tail-mate A/B),
`undertrim-probe-pe-dissim.json` (confounded control, see M1).
Reproduce: `bench/scripts/undertrim_{probe,probe_pe,probe_detect,ab}.py`.
The curated write-up is `docs/UNDERTRIM-REDUCTION.md`; this file is the
evidence-level record. Calibration: the corpus is synthetic at Q37; "residual"
counts a partial cut (0 < cut < t) as undertrim; the carrier rate is 29.9 %.

## Measured baseline (shipped: min_overlap 14, eps 0.2, illumina preset)

200 000 reads, 59 843 carriers: exact 174 875 / **residual 25 105 (41.95 % of
carriers)** / past-truth 20, of which 10 false trims on 140 157 clean reads
(0.007 %). Bands: 1-5 bp 9 026 under (36.0 % of residual), 6-9 bp 7 309
(29.1 %), 10-13 bp 7 374 (29.4 %), 14-21 bp 752, 22-33 bp 644.
**94.44 % of residual is strictly below the 14-base floor.**

## Mechanism map

| mechanism | code path | trigger | prevalence evidence |
|---|---|---|---|
| M1 sub-floor overhangs + PE-guard circularity | `adapter.rs` `match_adapters` (`k < min_overlap -> continue`); `pe_overlap.rs:275-288` `readthrough_is_adapter` calls `match_adapters(tail)` re-imposing the floor | any 3' overhang t <= 13; paired: insert >= L-13 | 94.44 % of residual (`undertrim-bands.json`); PE arms t=3/6/10 are 0/500 exact in auto/always/never (`undertrim-probe-pe.json`); decisive: random-tail mate trimmed 200/200 at t=20/mo14 and 200/200 at t=10/**mo10**, 0/200 at t=10/mo14 (`undertrim-probe-pe-guard.json`) — the detector pins the boundary, the guard's floor blocks it |
| M2 terminal G-run overprint (NovaSeq no-signal) | `match_adapters` is 3'-anchored (`k = min(a, l-p)` must reach the terminus); G's sit inside the window and spend the `floor(k*eps)` budget | 20 bp tail + terminal G-run | Arm C: g=0/2/4 -> 2000 exact, **g=6/8 -> 0 exact / 2000 none**, g=10/12 -> partial (`undertrim-probe.json`); 3 673 of 25 105 residual are G-class (`undertrim-ab.json`) |
| M3 `--pe-overlap auto` off for auto-detected sets | `cli.rs` `PeOverlapMode::Auto` runs only for `explicit` sources (FP-2) | default runs without `--adapter-preset`/`--adapters` | not in the 41.95 % (corpus used an explicit preset); affects real default runs |
| M4 mismatch-budget quantization / N as full mismatch | `floor(k*eps)`; N consumes budget | e >= 3 within the floor; N-bearing tails | Arm B: t=14 e=3 -> 27.2 %, e=4 -> 16.2 %, e=5 -> 6.2 %; t=20 e=5 -> 40.2 %; 1 185 n-class + 851 sub-class residual |

Note on stage order: `trim()` runs adapter detection BEFORE polyX on the
trimmed slice (`algorithm.rs:1304-1336`) — "poly-G runs first and hides the
adapter" is NOT the mechanism; M2 is anchoring + mismatch budget.

## Why the naive fix fails (both axes, `undertrim-ab.json`)

| config | under/carrier | false/clean |
|---|---|---|
| shipped (mo=14, eps=0.2) | 41.95 % | 0.007 % |
| mo=12 | 36.32 % | 0.022 % |
| mo=10 | 30.74 % | 0.140 % |
| mo=8 | 25.63 % | 0.243 % |
| mo=6 | 20.56 % | 1.466 % |
| mo=10, eps=0.3 | 27.49 % | 1.448 % |
| mo=14, eps=0.3 | 39.55 % | 0.373 % |

Arm D (20 000 clean reads) confirms the clean axis: mo=14 -> 0, mo=12 -> 5,
mo=10 -> 22, mo=8 -> 45, mo=6 -> 258 (1.290 %), mo=5 -> 861 (4.305 %).
Every global loosening costs 20-50x false trims for <= 16 pp residual, and
even mo=6 leaves the 1-5 band (36 % of residual) untouched.

## Ranked fixes (both axes named; full text in docs/UNDERTRIM-REDUCTION.md)

| ID | fix | expected undertrim reduction | overtrim risk | validation A/B |
|---|---|---|---|---|
| UT-1 | PE boundary as primary evidence: on a UNIQUE insert pin, verify the tail against `adapter[..t]` directly (t >= 2 exact, t >= 3 with `floor(t*eps)`), dual-mate corroboration, bypass `min_overlap_se` | M1 population: bands 6-13 (58.5 %) + most of 1-5 under dual-mate -> **41.95 % to 12-18 % on PE libraries** [E] | false pins refused on ambiguity + anchor ~1e-5; t<=5 gated on both mates; keep Cell A `pe_trims ~= 0` | `undertrim_ab.py` paired mode; `undertrim_probe_pe.py` t=3/6/10 must go 0/500 -> ~500/500 |
| UT-2 | G-overprint retry: on match failure with terminal G-run >= 4, re-match on the pre-G slice, cut at match start (adapter + G together) | M2: Arm C g=6..8 -> ~100 % exact; 3 673 reads (14.6 % of residual) | G-run is cut by polyX anyway; full k >= 14 bar on pre-G bases | Arm C ladder + `g under` column + Cell A |
| UT-3 | Two-tier corroborated short tails (SE, t=6..13): exact 3' prefix + 3' quality drop or dataset-level 3'-prefix enrichment | SE half of bands 6-13 (58.5 %) | uncorroborated form = Arm D mo=6 (1.29 %); corroboration must hold < 0.02 % [E] | new arm in `undertrim_ab.py` |
| UT-4 | Quality/N-aware scoring (N = 0.5 mismatch if >= 12 unambiguous bases agree; Q-weighted) | 1 185 n-class + part of 851 sub-class + Arm B e=3 gap | small [E]; no eps change | AB n/sub columns + clean corpus |
| UT-5 | Adapter-only reads (insert 0) come out PARTIAL (200/200, Arm C0) — trim-to-empty vs drop semantics | correctness wart | none | Arm C0 |

Deliberately rejected: global mo/eps lowering (table above); `--head-dropout`
revival (it overtrimmed real sequence; default correctly 0).

## Standing honesty rules

Quote CPU-seconds + busy cores, never the broken `pinned` wall column
(optim-dossier item 1). `winners.tsv` vs `fastp-comparison.tsv` disagree on
93/101 datasets — name both or neither. The 41.95 % is conditional on this
synthetic corpus (29.9 % carriers, Q37); real-data curves will move.
