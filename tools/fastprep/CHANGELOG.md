# fastprep changelog

Newest first. Default-behaviour changes are marked **DEFAULT** and always carry
the measurement that justified them.

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
  The same measurement says the shipped `--head-dropout 4` is *not* to be
  touched: hd0 reaches 82.80 % strict exactness but buys it with 9.58 % undertrim
  and scores −0.4947 at the anchor.

## 0.3.0 (unreleased baseline)

* `--overtrim-margin N`, `--short-tail-longest`, `--pe-ambiguous-smallest`, all
  default-off. Measured 0.0000 at the 7:1 anchor on the real-library GT corpus,
  because their classes (3′ polyX boundary, adapter readthrough, PE insert
  ambiguity) are 0.003 % / 0.01 % / absent there. Kept as instruments for
  corpora that contain them; see plan §7.
