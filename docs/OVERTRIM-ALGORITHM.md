# Overtrim algorithm work — what was proposed, measured, and recommended

Task: reduce fastprep's strict-tolerance (tol 0) polyX overtrim with a genuinely
new trim algorithm, without regressing speed or memory. This document records
what was tested, what the synthetic ground truth actually shows, and the honest
outcome. Measurement harness: `bench/scripts/overtrim_impl.py`. Candidate-rule
model: `bench/scripts/overtrim_rulesim.py`. Diagnostics: `bench/scripts/overtrim_poly.py`
(poly-only scoring), `bench/scripts/overtrim_diagnose.py` (per-read geometry).

## TL;DR

* The task framed the overtrim as **run continuation** (biology continuing with
  G's past an artefact run, so cutting the whole observed run eats biology).
  **fastsim never produces that mechanism.** `tools/fastsim/src/artefact.rs`
  ("Boundary rule") forces a non-G base immediately adjacent to every injected
  run, so on fastsim data the observed run equals the true artefact length
  `r == h` and run-continuation cannot occur. Verified empirically (below).
* What fastsim's ground truth *does* show is two different effects:
  1. **Mercy-overshoot (fixable, this change).** A polyX scan walks a clean
     homopolymer run and then keeps walking into the flanking non-run biology
     until the gap mass trips `mercy` — cutting **~2 bp past the run boundary**
     (measured on controlled clean runs). This is a boundary-placement error,
     and it is the entire strict overtrim in the heavy-tail regime.
  2. **1 bp leading-G ambiguity (not fixable by boundary placement).** The
     `head_dropout` rule trims a 1 bp leading G on reads whose true head
     artefact `h = 0`, i.e. a *biological* G. This dominates the overtrim in the
     low/mid regimes. On a 1 bp G-lead, artefact and biology are
     indistinguishable (same base, same quality, same boundary context) — this
     is `docs/OVERTRIM-DIAGNOSIS.md`'s irreducible ambiguity, and no boundary
     rule can recover it. It is prior-driven (which is why the sweep's
     `head_dropout` knob moves it).

## The change: `--boundary-anchored` (opt-in, default off)

A new `Params::boundary_anchored` (CLI `--boundary-anchored`). When set, a polyX
scan that terminates by mercy cuts the **homopolymer run extent** — from the scan
anchor to the last matching base — instead of `k_star`, the full mercy walk.

Why this is correct: `k_star` counts every step the mercy walk took, which
includes a couple of flanking non-run bases absorbed before `Δ_m` trips. The run
extent `run_span = last_match_step + 2` stops exactly at the run/flanking
boundary. It still keeps substitution-error bases *inside* the run (the last
match sits past an internal mismatch), so error-containing runs are not
short-changed — verified by a unit test.

Default (`boundary_anchored = false`) is byte-for-byte the shipped whole-run/mercy
cut. The `run_span` bookkeeping is gated on the flag, so the default path does no
extra work.

## Measurements (fastsim profiles, exact GT, poly-only scoring)

Scoring compares the **polyX** trim (`t5p`, `t3p`) against the GT head/tail
artefact lengths, per kept mate. Adapter trims are excluded (they are correct and
are not part of the polyX window; folding them in inflates 3′ "overtrim").

`bench/scripts/overtrim_impl.py accuracy`:

| regime | cut               | exact0 | over0  | und0   | exact2 | over2  | und2   | drops |
|--------|-------------------|--------|--------|--------|--------|--------|--------|-------|
| low    | whole-run         | 0.8116 | 0.1847 | 0.0037 | 0.9954 | 0.0010 | 0.0037 | 2329  |
| low    | boundary-anchored | 0.8116 | 0.1847 | 0.0037 | 0.9954 | 0.0010 | 0.0037 | 2329  |
| mid    | whole-run         | 0.8138 | 0.1693 | 0.0169 | 0.9818 | 0.0039 | 0.0143 | 2437  |
| mid    | boundary-anchored | 0.8145 | 0.1686 | 0.0169 | 0.9822 | 0.0035 | 0.0143 | 2437  |
| high   | whole-run         | 0.2181 | 0.2444 | 0.5375 | 0.2917 | 0.1708 | 0.5375 | 7162  |
| high   | boundary-anchored | 0.2392 | 0.2233 | 0.5375 | 0.3036 | 0.1589 | 0.5375 | 7162  |

* **high regime:** overtrim 0.2444 → 0.2233 (−2.1 pp), converted to exact
  (0.2181 → 0.2392, +2.1 pp) with **undertrim and drops unchanged**. This is the
  mercy-overshoot being removed. The residual undertrim (53.75%) is `min_run 12`
  missing the profile's 4–11 bp tails — a *config* gap, not this algorithm.
* **low/mid:** unchanged. The overshoot only appears when a run reaches `min_run`
  (≥12 bp); low/mid tails are mostly 0 or short, so there is nothing to fix — and
  their overtrim is the 1 bp `head_dropout` lead (mechanism 2), which boundary
  placement cannot touch.
* No regime regresses: overtrim only ever moves toward exact, never toward
  undertrim, and drop counts are identical (the length budget is unaffected).

Controlled clean-run probe (`tools/fastprep/examples/probe_clean_run.rs`), run
length 12/16/20, `min_run 12`, `mercy 3`:

| run length | whole-run cut (err) | boundary-anchored cut (err) |
|------------|---------------------|-----------------------------|
| 12         | 14 (+2)             | 12 (0)                      |
| 16         | 18 (+2)             | 16 (0)                      |
| 20         | 22 (+2)             | 20 (0)                      |
| 14 (1 internal error) | 14 (0)   | 14 (0) — error base kept    |

## Speed and memory (100k-pair file, `-j 4`)

`bench/scripts/overtrim_impl.py speed --pairs 100000` (peak RSS = the child's
own `/proc` VmHWM polled every 20 ms — `ru_maxrss` cannot isolate one child, see
`bench/scripts/fastp_compare.py`):

| cut                        | wall (median) | throughput       | peak RSS | drops |
|----------------------------|---------------|------------------|----------|-------|
| original (git HEAD)        | 0.363 s       | ~275k pairs/s    | 34.3 MB  | 21133 |
| new, default (whole-run)   | 0.363 s       | ~275k pairs/s    | 34.3 MB  | 21133 |
| new, `--boundary-anchored` | 0.363 s       | ~275k pairs/s    | 34.3 MB  | 21133 |

All three are identical: **no throughput regression, no memory regression.** The
default path is unchanged from `git HEAD` (the `run_span` bookkeeping is gated on
the flag). The flag itself is free. Drops are identical across all three.

## Candidates evaluated (the four proposed ideas)

1. **Posterior-expected cut of the artefact length** (`bench/scripts/overtrim_rulesim.py`):
   on a *run-continuation* generative model (artefact `h` + geometric biological
   continuation `k`, observed run `r = h+k`), posterior-mean/MAP cut raises exact
   from 0.352 to 0.383 and halves overtrim on a heavy-tail (high) prior, at no
   undertrim cost. **But run-continuation does not occur on fastsim data** (the
   boundary rule removes it), so on the actual ground-truth measurement this rule
   is a no-op. It also needs the artefact-length prior as input. Not shipped — no
   measurable gain on real data.
2. **Quality-gated boundary stop**: the quality signal exists (the tail run has a
   collapsed `tail_q` ~3 Phred below flanking biology) but the scan's gate
   `g(q) = min(1, q/phred_sat)` **saturates at `phred_sat = 20`**, so both the
   collapsed tail (~34) and normal biology (~37) read as `g = 1`. The gate is
   blind to the 3-Phred split. The *structural* boundary (run → non-run) is the
   reliable signal, and that is exactly what `--boundary-anchored` cuts on — so
   this idea is realized structurally rather than by quality.
3. **Interior-scan extension**: the boundary is already locatable from the run
   itself (a run followed by non-run); the interior scan is off and stays off.
   Subsumed by the boundary-anchored cut.
4. **Bounded-h rule (`min(r, h_max)`)**: fragile in the model — +0.007 exact at
   `h_max = 20` but −0.08 at `h_max = 12` (too-aggressive caps undertrim genuine
   long runs). Not robust enough to ship.

## Honest limitations

* The change **does not** reduce the low/mid overtrim (17–18%). That overtrim is
  the 1 bp `head_dropout` lead — an artefact/biology ambiguity, not a boundary
  error. It is prior-driven (which is why the `head_dropout` knob moves it) and,
  per `docs/OVERTRIM-DIAGNOSIS.md`, no base/quality feature separates the two.
* The pre-existing `head_dropout` "FIX" in the working tree (the comment claims it
  fixes the 29.5% overtrim) is a **no-op**: under the CLI's enforced
  `head_dropout < min_run`, `min_run.min(n+1) == n+1`, so the rule fires on
  `1..=n` exactly as before. Verified by reasoning over the code and by the fact
  that the 1 bp lead overtrim is unchanged before/after. Its comment overstates
  the effect.
* The huge undertrim in the high regime (53.75%) is `min_run 12` sitting above
  the profile's 4–11 bp tails; lowering `min_run` recovers it (that is knob
  tuning, out of scope here).

## Recommendation

Ship `--boundary-anchored` as an opt-in boundary-placement fix: it removes the
mercy-overshoot (measured −2.1 pp overtrim / +2.1 pp exact on the heavy-tail
regime, no undertrim/drop cost) at zero speed and memory cost, and is default-off.
Report the 1 bp `head_dropout` lead as the known, irreducible overtrim with the
±2 tolerance number (overtrim falls to ~0.4–1.6% at tol 2), which is what the
paper already does. Do **not** expect a boundary rule to recover the 1 bp
ambiguity — the information is not in the read.
