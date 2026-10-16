# The Bayes-optimal ceiling for run-length trimming

**Question.** fastprep trims a poly-G artefact by cutting the entire observed run,
`f(r) = r`. On the Phase D bench that rule overtrimmed 29.5% of mates at strict
tolerance 0 (1.3% at tolerance ±2 — the errors are 1-2 bp). Is that overtrim
fixable by *any* decision rule that reads only the observed run length `r`, or is
`f(r) = r` already optimal given base identity alone?

**Short answer.** Fixable in large part, not in full — and the current rule is
*not* near the optimum. On the synthetic profiles that reproduce the mechanism
(fastsim arm, `mid` regime, 200k pairs / 400k mates):

| mate-joint, tolerance 0 | exact | overtrim | undertrim |
|---|---|---|---|
| current `f(r) = r` | 60.59% | 39.40% | 0.00% |
| Bayes-optimal rule of `r` | 90.12% | 0.02% | 9.86% |
| **ceiling (any function of `r`)** | **90.15%** | — | — |

75% of the strict error (29.5 of the 39.4 pp) is reachable from `r` alone; the
residual 9.9 pp is not. Overtrim is *not* conserved against undertrim: the Bayes
rule converts 29.5 pp of overtrim into exact calls and only 9.9 pp into
undertrim. At tolerance ±2 the current rule already scores 97.38% and the
ceiling is 99.93%.

Everything below is reproduced by `bench/scripts/overtrim_ceiling.py`.

## What is computed

Reads are simulated with `tools/fastsim` from the synthetic profiles in
`bench/fixtures/overtrim/` (`profiles_overtrim.py`), scored against fastsim's
exact per-mate ground truth (`sim.v0.gt.tsv`: `pair_id insert r1_h r1_t r2_h
r2_t r1_adapter r2_adapter full_flag is_dup`). For every mate end we build the
empirical conditional **P(true artefact | observed run r)** and then score a
family of rules `f(r) -> cut` (constant, `r`, `r-c`, `min(r,k)`, `round(a*r)`,
posterior-thresholded mode) at tolerance 0 and 2. The ceiling is
`sum_r max_cut P(|h-cut| <= tol | r) P(r)` — the best any function of `r` can do.
Classification follows `overtrim_sensitivity.py` (mate-joint: exact iff both
ends within tolerance; overtrim iff either end cuts past the artefact).

## The conditional table (fastsim arm, `mid`, 400k mates)

5' end — `h` is the true head artefact:

| r | P(r) | P(h=0\|r) | P(h=1\|r) | argmax |
|---|---|---|---|---|
| 0 | 0.686 | 1.000 | 0.000 | 0 |
| 1 | 0.259 | 0.644 | 0.356 | 0 |
| 2 | 0.042 | 1.000 | 0.000 | 0 |
| 3 | 0.010 | 1.000 | 0.000 | 0 |
| 4 | 0.003 | 0.999 | 0.001 | 0 |
| ≥5 | 0.001 | 1.000 | 0.000 | 0 |

3' end — `t` is the true tail artefact (columns are `P(t | r)`):

| r | P(r) | t=0 | t=r | argmax |
|---|---|---|---|---|
| 0 | 0.755 | 1.000 | — | 0 |
| 1 | 0.172 | 0.983 | 0.017 | 0 |
| 2 | 0.040 | 0.983 | 0.017 | 0 |
| 3 | 0.011 | 0.901 | 0.099 | 0 |
| 4 | 0.004 | 0.504 | 0.495 | 0 |
| 5 | 0.001 | 0.543 | 0.453 (t=6, error-truncated) | 0 |
| 6 | 0.005 | 0.030 | 0.970 | 6 |
| 7 | 0.002 | 0.007 | 0.991 | 7 |
| 8 | 0.001 | 0.011 | 0.989 | 8 |
| ≥9 | 0.008 | ≈0 | ≈1 | r |

Two structural facts drive everything:

1. **The table has support {0, r} only.** fastsim's boundary rule
   (`artefact.rs`) forces the base next to a non-zero artefact run to be non-G,
   so a simulated run is *exactly* the artefact unless the artefact is zero — in
   which case the run is pure biology. Substitution errors (0.55%/base) add the
   thin `h > r` truncation tail. The probe's merging mechanism (biology
   continuing the run) cannot appear in this ground truth.
2. **Given `r`, the only question is artefact-or-biology for the whole run.**
   `r = 1` at 5' is a near coin flip (36% artefact pooled; 55% on R2 alone, whose
   prior carries the 1 bp dropout; 0% on R1). `r >= 2` at 5' is pure biology with
   certainty. At 3', `r <= 3` is ~pure biology and `r >= 6` is ~pure artefact.

## The four questions (fastsim arm, `mid`)

**Maximum strict exact rate for any function of `r` alone.** Per end:
5' **90.79%**, 3' **99.24%**. Mate-joint (both ends decided together):
**90.15%** (per-end composition gives 90.12% — the other end's run carries no
extra signal). If the rule may additionally know which mate it reads (R1 and R2
have different priors), the ceiling rises to 5' 92.45% (R1 100%, R2 84.89%),
3' 99.28%.

**How close is `f(r) = r`?** Not close. Per end: 5' 77.81% (gap 13.0 pp),
3' 77.86% (gap 21.4 pp); mate-joint 60.59% (gap 29.6 pp). `f(r) = r` sits at the
zero-undertrim corner of the achievable set and pays 39.40% overtrim for it. The
real tool (fastprep default config, scored on its kept 357,645 mates) gets
76.02% exact / 22.34% over / 1.64% under; its empirical policy is `cut = r` for
5' runs of 1-4 bp (0 for ≥5) and almost no 3' trimming (only runs ≥12 bp).

**Can overtrim drop without an equal undertrim increase?** Yes — the exchange
rate is far from 1:1, and total error is *not* conserved.

- Class flow `f(r)=r ->` Bayes rule, tolerance 0: at 5', 22.19 pp
  overtrim→exact at the cost of 9.21 pp exact→undertrim (2.4:1); at 3',
  22.12 pp overtrim→exact for 0.74 pp exact→undertrim (30:1). Mate-joint:
  39.38 pp overtrim removed, 29.52 pp of it landing exact, 9.86 pp landing
  undertrim. Total error falls 39.40% → 9.88%.
- Part of the gain is strictly free: any rule that abstains on runs that are
  certain biology cannot undertrim there. That covers 5.5 pp of mates at 5'
  (the `r >= 2` mass) and ~21 pp at 3' (`r = 1,2,3,5` are 90-98% pure biology).
- The one place where the trade is strict is the 5' `r = 1` cell: trimming it
  overtrims 16.7 pp of mates, abstaining undertrims 9.2 pp. From `r` alone that
  cell cannot be told apart; it is the entire irreducible residue.

**Does the answer change at tolerance ±2?** Yes — it mostly disappears.
Current `f(r)=r` mate-joint: 97.38% exact / 2.62% over. Ceiling: 99.93%.
Reachable: 2.55 pp of the 2.62 pp; unreachable from `r` alone: 0.07 pp. The
`cut = max(0, r-2)` rule alone reaches 99.83%.

## The merging mechanism (analytic merge arm — NOT fastsim GT)

The probe (`/tmp/probe_overtrim.py`) mechanism — biology continuing the artefact
run, `r = h + b` — is structurally different and worse for `r`-only rules. Using
the same priors with `b ~ Geom(0.235)` (iid DNA, GC 0.47), `mid`:

| mate-joint, tol 0 | exact | over | under | ceiling |
|---|---|---|---|---|
| `f(r) = r` | 58.52% | 41.48% | 0.00% | — |
| Bayes rule of `r` | 80.43% | 0.33% | 19.24% | **80.43%** |

Here the task's hypothesis is confirmed: with `h ∈ {0, 1}`, `f(r)=r` is never
exactly right once the biology adds ≥1 G, and the optimal cut is `min(r, 1)`-like.
But the 5' posterior is 49/51 at *every* run length (the geometric boundary run
makes `P(h=1|r)` constant), so the cut/no-cut decision is a near coin flip and
the strict ceiling from `r` alone is only 81.54% at 5'. At tolerance ±2 the
problem again mostly vanishes (`f(r)=r` scores 97.42%, ceiling 99.98%), because
the overshoot is exactly +1..2.

## What this does NOT claim

Only the observed run length is in scope. Base identity is not the only
information: the quality collapse inside the run (3.0 Phred, the two-colour
signature, `tail_q_collapse`), mate agreement, insert length and cycle position
all carry signal that no rule of `r` uses. The residual 9.9 pp (tol 0) is
irreducible *from `r` alone* — it is not claimed irreducible in general. The
most promising target is the 5' `r = 1` cell (16.7% of mates): the R1/R2 prior
split already moves the ceiling 2.3 pp, and in-run quality is the obvious next
feature. The low/high regimes behave like `mid` on the fastsim arm (ceiling
99.98% / 99.65% mate-joint; current 58.45% / 71.28%), so the conclusion is not
an artefact of one profile.

## Reproduce

```
python3 bench/scripts/overtrim_ceiling.py --pairs 200000 --regimes mid low high
```

Needs `tools/fastsim/target/release/fastsim` (and, optionally,
`tools/fastprep/target/release/fastprep` for the reference row). Regenerates the
fixtures via `profiles_overtrim.py` if missing, writes the full report to stdout
and a machine-readable copy to `~/.hermes/cache/scratch/overtrim_ceiling/ceiling.json`.
