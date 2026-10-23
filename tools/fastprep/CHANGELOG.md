# fastprep changelog

Newest first. Default-behaviour changes are marked **DEFAULT** and always carry
the measurement that justified them.

## Unreleased

* `--pe-correct [GATE]` (opt-in; paired-end FASTQ only): posterior consensus
  correction of base errors in the doubly-retained mate overlap. At overlap
  disagreements a two-hypothesis posterior over the mates' quality strings
  picks the call (`|log10 odds| >= GATE`, default 2.0; ties never resolved by
  fiat) and rewrites only the losing mate's base and its Q, rescored to the
  capped posterior Q (max 40). The insert is inferred from the pair and
  accepted only when exactly one candidate verifies — a repeated fragment
  fitting several inserts is refused, not guessed; `--pe-correct-min-overlap`
  (default 30) matches `fastp --correction`'s own overlap floor. Trim
  decisions run before correction and are never revisited, so the flag cannot
  move a trim boundary (flag-off output is byte-identical to 0.4.1 on a
  50k-pair corpus, both mates, sha256-checked). Measured on the 20k-pair
  truth corpus at matched retention (40k reads out, fragment-only
  denominators; `bench/results/consensus-headroom*.tsv`): `--pe-correct=2.0`
  leaves **5.2% fewer residual base errors than `fastp 1.4.0 --correction`**
  (22,190 vs 23,405 fragment errors; 3.93 vs 4.15 per kb) at comparable
  corruption (13 vs 10 correct bases broken); `=0.5` leaves 13.0% fewer
  (20,355) at 104 corrupted. The margin holds when quality strings are
  miscalibrated 2× in either direction (5.3% and 5.0% fewer, corruption 5 vs 3
  and 15 vs 15). Both axes reported; this is not a default change.

## 0.4.1

* **DEFAULT — `--head-dropout` is back to 4.** The rule trims a 1..=4 bp
  leading G-run exactly, with no quality gate. It is a prior, not a detector:
  trimming it is wrong on roughly three reads in four that it touches. Whether
  that is a mistake depends on the price of the two errors, and the price was
  never written down before this round.
  Measured at the 7:1 anchor ("prefer less data to wrong data": an artefact
  left in the output costs seven deleted bases), 102 libraries, desktop only,
  18.96 M mates (`bench/results/four-tool/`): `4` scores **+43.6** of mates,
  `0` scores **+17.1** — a 26.5-point margin for the rule being on. The
  break-even is a cost ratio of **2.4**, measured four ways: the quality
  posterior for a 1..4 bp leading run is 0.27–0.29 artefact, so
  P(biology)/P(artefact) ≈ 2.5–2.7; the 102-library four-tool table gives
  19.97/8.31 = 2.40; the twelve-arm campaign gives 2.39; the 101-library
  head-dropout axis gives 3.70. At 1:1 the ordering inverts (total error 22.2 %
  with the rule, 10.5 % without), which is the reading that made `0` the default
  on 2026-10-16 (`3d76a31`).
  So the two defaults are both defensible and encode different prices; this
  release encodes the 7:1 one, because that is the project's stated preference.
  `--head-dropout 0` remains the setting for a pipeline that prices a lost base
  and a kept artefact equally. Full table: `docs/plans/2026-10-21_fastprep-exactness.md`
  §7 and `docs/TRIM-ACCURACY.md`; decision record: `docs/DECISIONS.md` D-16.
  Verified behaviourally, not by inspection: on `bacterial_isolate/DRR014197`
  (config 45, a dataset where the knob bites) the bare invocation and an
  explicit `--head-dropout 4` produce the same PXMT md5
  (`3afed057…`) and `--head-dropout 0` produces a different one (`cfc0186f…`) —
  `bench/scripts/prove_hd_default.py`, which refuses a dataset where the two
  agree because that would prove nothing.

## 0.4.0

* **DEFAULT — `--boundary-anchored` is now on.** The polyX end cut trims exactly
  the matched run extent instead of the mercy walk's end, so the ~2 bp flank
  overshoot the walk absorbs is no longer cut. `--no-boundary-anchored` (or
  `--boundary-anchored=false`) restores the pre-0.4.0 cut; the bare
  `--boundary-anchored` spelling still parses as a no-op so existing bench
  invocations keep working.
  Measured, 102 datasets at the 7:1 trim cost anchor (desktop only, 1 224 rows,
  `bench/results/exactness/`): +0.0121 per mate, 0.32 pp of mates moved from
  overtrim to exact, undertrim unchanged (1.87 → 1.88 %), clean-read cuts
  unchanged. Full table: `docs/plans/2026-10-21_fastprep-exactness.md` §7.
  The same measurement prices the leading-G rule: the campaign's per-dataset
  winner configuration runs `--head-dropout 4`, and the arm that reaches the
  highest strict exactness (`0`, 82.80 %) pays for it with 9.58 % undertrim and
  scores −0.4947 at the anchor. (The shipped default at the time of 0.4.0 was
  already `0`, from `3d76a31`; 0.4.1 changes that, above.)

## 0.3.0 (unreleased baseline)

* `--overtrim-margin N`, `--short-tail-longest`, `--pe-ambiguous-smallest`, all
  default-off. Measured 0.0000 at the 7:1 anchor on the real-library GT corpus,
  because their classes (3′ polyX boundary, adapter readthrough, PE insert
  ambiguity) are 0.003 % / 0.01 % / absent there. Kept as instruments for
  corpora that contain them; see plan §7.
