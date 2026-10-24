# Reducing the strict-tolerance overtrim — what works and what does not

Bench reference: `fastprep-winner` over 101 real libraries, strict tolerance 0
= 66.5% exact / 29.5% overtrim / 4.0% undertrim; tolerance ±2 = 98.4 / 1.3 /
0.3. `fastp` default: 92.2 / 0.2 / 7.6 at tolerance 0. (The 66.5% is the
published cross-generation number, later shown invalid; the 0.9L-era
shipped-default figure is **75.4%** — see `docs/TRIM-ACCURACY-DISCREPANCY.md`.
The D-18 scan-change re-measure (15 datasets, `--min-len` pinned) puts the
`head-dropout 4` arm at **79.19%** exact / **19.39%** overtrim; the `0.67L`
floor change is not in those numbers.)

Four independent investigations were run (ceiling, quality, paired-end,
implementation). This is the consolidated answer to "fix the overtrim".

All measurements below were taken 2026-10-16, before D-18: `0.9L` floor, the
pre-D-18 mercy clause, and `--boundary-anchored` not yet the default. They
describe that build; the shipped default today is `0.67L` with the anchored cut
on.

## The mechanism, corrected

`docs/OVERTRIM-DIAGNOSIS.md` attributed the overtrim to artefact and biology
G's merging into one run. That mechanism **cannot occur in fastsim output**:
`tools/fastsim/src/artefact.rs`'s boundary rule forces a non-G base adjacent to
every injected run, so the observed run is exactly the artefact length and there
is nothing to merge. Verified independently by the ceiling and implementation
workers.

What the data actually shows is two separate errors:

1. **A mercy-overshoot** (~2 bp). The scan walks a clean homopolymer run, then
   keeps walking into flanking non-run biology until the gap mass trips `mercy`,
   so a 12/16/20 bp run is cut at 14/18/22. This is a boundary-placement bug
   and it is fixable.
2. **An unconditional leading-G hard trim.** `head_dropout` (the shipped
   default 4 at the time; briefly 0 after `3d76a31`; 4 again since 0.4.1, D-16)
   trims any leading G-run of length 1..=4 exactly, with no quality gate and no
   question of whether the G is artefact. On a library whose first base is often
   a real G, this removes real sequence.

## Finding 1 — the then-default `head-dropout 4` is dominated by `1`

Measured per mate against fastsim ground truth, tolerance 0 (30k pairs, the
`bench/fixtures/overtrim/` profiles):

| profile | mate | hd=4 (default) | hd=1 | hd=2 | hd=0 |
|---|---|---|---|---|---|
| mid (R1: no artefact) | R1 | 74.16 / 24.33 / 1.52 | **79.86** / 18.50 / 1.64 | 75.48 / 22.97 / 1.54 | **97.92** / 0.08 / 2.00 |
| mid (R2: 18.5% h=1) | R2 | 77.80 / 20.49 / 1.71 | **82.51** / 15.66 / 1.83 | 78.87 / 19.40 / 1.73 | 78.96 / 0.53 / **20.51** |
| low (no artefact) | R1 | 75.05 / 24.60 / 0.35 | **81.05** / 18.57 / 0.38 | 76.43 / 23.20 / 0.36 | **99.43** / 0.07 / 0.49 |
| low (no artefact) | R2 | 75.07 / 24.56 / 0.37 | **81.18** / 18.42 / 0.40 | 76.40 / 23.23 / 0.37 | **99.43** / 0.08 / 0.50 |

Columns are `exact / overtrim / undertrim`, percentages.

**`head-dropout 1` strictly dominates `4`: +5 to +6 pp exact, with undertrim
unchanged to within 0.1 pp, on every mate of every profile.** `2` is worse than
`1` (it trims two G's instead of one, and the second is usually biology). The
default is the worst of the four settings tested.

`0` (rule disabled) is nearly perfect where there is no head artefact (99.4%
exact) and badly undertrims where there is (20.5% undertrim on mid's R2). So it
is the right setting per library, not a universal one.

## Finding 2 — the rule could not be switched off at all

`--head-dropout 0` was rejected: *"must be ≥ 1"*. The validation is now relaxed
to allow `0`, because the rule is unconditional and there is no other way to
avoid trimming a leading G that is real biology. The old comment claiming a fix
made the rule "monotone in `min_run`" does not survive inspection: with the
CLI-enforced `head_dropout < min_run`, `min_run.min(n+1) == n+1`, so runs of
length `1..=n` are hard-trimmed exactly as before. The fix only moved the
`>= min_run` boundary, which the scan already owned.

## Finding 3 — `--boundary-anchored` removes the overshoot, at zero cost

Shipped opt-in at `fd52adb`; on by default since 0.4.0 (`--no-boundary-anchored`
restores the old cut). Cuts the homopolymer run extent (the last
matching base) instead of the mercy walk's `k_star`, dropping the ~2 bp flanking
overshoot. Keeps error bases inside the run (the last match sits past them).
Measured: heavy-tail overtrim 0.2444 → 0.2233, exact 0.2181 → 0.2392, undertrim
and drops unchanged. **Speed and peak RSS identical to HEAD** (0.345 s per 100k
pairs at `-j4`, 34.3 MB VmHWM).

It does not touch the leading-G overtrim — that is finding 1's job.

## What does NOT work

- **Quality at the 5' boundary.** fastsim gives the first head base normal
  quality by construction (`deep_head_q` only applies from position 2), so
  there is no step to detect at the dominant h=1 case. Measured P(step) = 23.2%
  = rounding noise. A quality rule cannot help there.
- **Quality at the 3' boundary.** Here it *does* work — fastsim's own model
  writes artefact runs 3.0 Phred below biology (`tail_q_collapse`), 942/942
  boundaries separable, and a `q_cut=36` rule lifts tolerance-0 exact 81.2% →
  90.5% (mid) at −1.4% throughput. **But this is a property of the simulator's
  quality model.** Real libraries have more noise than fastsim's sd 0.5 and
  measured collapse varies 1.7–11.7 Phred, so it is unverified on real reads and
  should not be shipped on this evidence alone.
- **Paired-end / insert-length.** Mate corroboration resolves the ambiguity
  perfectly where the mate covers the boundary (60.2% → 99.9% exact, zero
  breakage), but the error sits at the fragment's *outer* ends, which the other
  mate covers only when `I <= L` — 4.7% of the fitted insert distribution. Net
  effect: +0.1 pp pooled (304 fixed, 86 broken). Not a cure.
- **Posterior / bounded-h cuts on fastsim data.** Moot: the run-continuation
  mechanism they target does not exist in fastsim output.

## The ceiling

Any rule reading only the observed run length `r` can reach **90.15%** strict
exact on the mid profile (current `f(r)=r` scores 60.59% on the same synthetic
arm). 75% of the error is reachable from `r` alone; the residual 25% is the 5'
`r=1` cell, a genuine 36/64 artefact/biology coin flip that base identity cannot
settle. Total error is **not conserved**: the Bayes rule converts 29.4 pp
overtrim into exact while creating only 9.9 pp undertrim (39.40% → 9.88% total
error). See `docs/OVERTRIM-CEILING.md`.

## Recommendation

1. Ship `--boundary-anchored` (done; on by default since 0.4.0) — fixes the
   overshoot, free.
2. Change the default `--head-dropout 4` → `1`, or make it per-mate. At a 1:1
   error price this is a strictly-better default: +5–6 pp exact, no cost.
   **Superseded** — the default was set to `0` (`3d76a31`), then reverted to
   `4` in 0.4.1 on D-16's 7:1 price, which is the shipped default today;
   `1` remains the compromise for long-read regimes and `0` is still reachable
   as `--head-dropout 0`.
3. Allow `--head-dropout 0` (done) and document that the right value is
   library-dependent.
4. Do not ship the 3' quality rule on fastsim evidence alone.
5. Report the residual honestly: the 5' `r=1` ambiguity is 36/64 and no
   base-identity rule can settle it.

## Reproduce

```bash
python3 bench/scripts/overtrim_head_dropout.py --profile bench/fixtures/overtrim/mid.json --out /tmp/hd.tsv
python3 bench/scripts/overtrim_permate.py        # the per-mate table above
python3 bench/scripts/overtrim_speed.py          # timing + RSS
python3 bench/scripts/overtrim_ceiling.py --pairs 200000 --regimes mid
python3 bench/scripts/overtrim_quality.py --part all
python3 bench/scripts/overtrim_paired.py sim
```
