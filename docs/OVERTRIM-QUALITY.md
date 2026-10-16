# Can base quality discriminate an artefact G-run from biology at the boundary?

Short-tolerance overtrim root cause: fastprep trims the whole *observed* G-run
(artefact + any biological G's continuing at the boundary). Base identity
cannot split it (see `OVERTRIM-DIAGNOSIS.md`). This note settles whether
**quality** can — with probes that survive fastprep's length budget (the
earlier `/tmp/probe_q_split.py` could not: 28/30 reads were budget-dropped
before any boundary decision), and with the rule evaluated in **two**
settings: quality drawn from **fastsim's own model output, unmodified**, and
quality drawn from **chosen constants** (my probe worlds).

Reproduce: `python3 bench/scripts/overtrim_quality.py --part all`
(parts: `model`, `probe`, `rule`, `cost`).

## Verdict

**Split by boundary class — and the split is a property of fastsim's quality
model, not of the probe.**

| boundary class | signal in fastsim's model (unmodified) | can a quality rule help? |
|---|---|---|
| **3' tail run + biology ending in G's** | **yes: 3.0 Phred step** (tail run Q34.0 vs biology Q37.0, sd 0.55; P(artefact base < boundary base) = 942/942 = 100%) | **yes** — tol-0 exact 81.2→90.5% (mid) and 21.6→98.9% (high) |
| **5' head run, h = 1** (the dominant 5' case) | **none** — the model gives run position 0 NORMAL quality (37.00 = biology; P(step) = 23.2% = rounding noise) | **no** — and the rule's entire false-positive cost lives here |
| **5' head run, h >= 2** | only if the profile carries `deep_head_q`/`run_q_mean.head` (bench fixtures do not: 37.00 vs 37.00 — no signal at all) | yes where the profile carries it (positions 2..h at 25.0 vs biology 37.0 measured) |

So: the improvement on the tail is **real in fastsim's own output** (the
as-is world applies the rule to qualities fastsim wrote; no constants were
chosen). The head h=1 case — which is most of the 5' strict overtrim — has
**no signal in the model by construction**, and on fastsim data the rule
cannot fix it.

**Transferability to real libraries is UNVERIFIED.** The 3-Phred tail step
in fastsim comes from the fitted profile's `tail_q_collapse = 3.0` (measured
real values vary: DRR014283 R1 11.7, R2 1.7), and fastsim's per-base noise
is small (`qual_noise = 0.5`), so the distributions barely overlap. Real
reads have much larger noise and 3'-quality decay; the `overlap` probe world
(noisy, decaying) shows the gain degrading into a threshold tradeoff (mid:
88.6% exact at best, or 82.9% with 10% overtrim). The Phase-D real-profile
data is on the unreachable desktop; **the rule's real-data gain is not
measured anywhere in this note.**

## 1. What fastsim actually writes for artefact quality (source + measured)

Source (`tools/fastsim/src/quality.rs` `qualify()`, `profile.rs`
`build_profile`):

* **Head run `[0, h)`**: position 0 *always* gets the normal per-position
  curve (`if i == 0 { pos_q[i] + offset }` — a hardcoded exception);
  positions 2..h (1-based) get `head_q_mean`, which is `deep_head_q` or
  `run_q_mean.head` when the profile carries them, else **`mean(pos_q)` —
  i.e. no collapse at all**. The bench fixtures (`bench/fixtures/overtrim`)
  carry neither field.
* **Tail run `[L-t, L)`**: *all* positions get `tail_q_mean =
  mean(pos_q) − tail_q_collapse` (fixtures: 37 − 3 = 34). No first-base
  exception, **no decay through the run**, no run-length dependence: one
  flat run mean plus per-base Gaussian noise (sd 0.5).
* **First biology base after/before the run** (`pos_q`): normal quality.

Measured (20k pairs from `mid.json`, adapter-free, GT-labeled):

```
head_pos0 (h>=1)        mean Q=37.00   (== biology: no signal at h=1)
tail_run (all t>=1)     mean Q=33.99   sd 0.57   (biology 37.00)
tail by run offset      off0 34.00  off1 33.97  off2 34.01  off3 33.99  (flat)
bio base before tail    mean Q=36.98   -> step = -3.01 Phred
P(last artefact < first biology base):
    tail boundary       942/942  = 100.0%
    head boundary h=1   863/3714 =  23.2%   (ties dominate: pure noise)
head-probe h=2..6, profile without deep_head_q: pos0 37.00, positions 2..h 37.00
head-probe h=2..6, deep_head_q=25:              pos0 37.00, positions 2..h 25.00
```

The diagnosis's claim is therefore verified but needs precision: the model
gives **the first head-run base** (position 0) normal quality — so for h=1
(the dominant 5' artefact) the artefact base itself is quality-identical to
biology, and the boundary is invisible. Deeper head bases collapse only
when the profile carries measured deep-head quality; the bench fixtures
don't, so in fixture-simulated data the **whole head artefact** is
quality-flat with biology.

## 2. Does fastprep's cut vary with quality at all? (probe grid)

Design (`--part probe`): 4800 reads, 128 bp, runs 3–6 bp, k ∈ {0,1,2}
biological G's at the boundary, q_art/q_bio ∈ {2,10,20,30,40}, interior
q ∈ {15,35}, both ends. `--min-run 3 --min-len 30`: trim budget 98 bp,
never fires (PXMT status bytes: **kept = 4800/4800**, zero drops — the
earlier probe's failure mode is gone).

Results (binary built from this tree):

* **The cut never lands inside the observed run: 0/4800** at `--mercy 1`
  and `--mercy 3`. The run is always trimmed whole (plus overshoot), so
  quality cannot currently split artefact from biology at the boundary —
  the walk's in-run steps are matches with gap `mu(1-W)` <= 0.1, far below
  the mercy mass 3.0, so it cannot stop inside a 3–6 bp run at any quality.
* `--mercy 1`: cut == observed run exactly, 4800/4800; cut − truth = the
  k boundary G's exactly (+0:1600, +1:1600, +2:1600).
* `--mercy 3`: cut = observed run + 1..20 bp overshoot into *non-G*
  biology (mean +2.56, sd 1.32). Quality moves only this overshoot, and
  perversely: low-q artefact runs overshoot *less* (+1.6..+2.3 at
  q_art <= 10) than high-q runs (+2.8..+3.3 at q_art >= 20), because
  low-quality matches accumulate mercy mass slightly faster. The boundary
  base's quality has no systematic effect. (The overshoot is a separate
  overtrim source; the base-identity `--boundary-anchored` cut in this tree
  removes it — orthogonal to the quality question.)

## 3. The quality-boundary rule

Trim at each end the **maximal G-run segment whose bases are all below a
quality threshold `q_cut`** (walk inward while base == G and q < q_cut).
Two semantics: `replace` (the rule owns the end cut) and `refine` (fall
back to the scan's cut when the rule finds nothing). Implemented in
fastprep as `--q-boundary <Q>` (default off, `replace` semantics,
budget-guarded; `head_dropout` subsumed) and offline in the analysis
script. Needs **no look-ahead** beyond the run: the walk stops at the first
non-G or high-quality base.

## 4. Evaluation — both settings

Scoring is the trim-accuracy classifier per mate (exact / overtrim /
undertrim at tol 0 and 2) against fastsim GT, adapter-free mates, full-flag
pairs excluded, only fastprep-kept mates scored. Baseline = fastprep with
the sensitivity BASE config (`min_run 12, mercy 3, phred_sat 20, min_len
0.9L, mean_q 25, head-dropout 1`).

### Setting A — quality from FASTSIM'S MODEL, unmodified (the settled answer)

Rule best `replace` q_cut = 36 (tol 0 / tol 2):

| dataset | fastprep now tol0 (exact/over/under) | rule tol0 | fastprep tol2 | rule tol2 |
|---|---|---|---|---|
| low  | 81.1 / 18.5 / 0.4 | **99.9 / 0.1 / 0.0** | 99.4 | **100.0** |
| mid  | 81.2 / 17.0 / 1.8 | **90.5 / 0.1 / 9.4** | 98.1 | **100.0** |
| high | 21.6 / 24.6 / 53.8 | **98.9 / 0.0 / 1.1** | 29.1 | **99.2** |

Sweeps show a threshold plateau at q_cut 35–36 (between the model's Q34
artefact and Q37 biology), not a knife-edge fit. The mid undertrim 9.4% is
**entirely the h = 1 head artefacts** the model endows with normal quality
(11.3% of baseline-exact mates become 1-bp undertrim — see §5). `refine`
semantics is weaker (mid tol0 83.1 / 16.8 / 0.0): it only removes undertrim
where the rule finds nothing.

Note the high-regime gain is inflated by that regime's own pathology: ~90%
of mates carry 4–20 bp tails which `min_run 12` mostly leaves (undertrim
53.8%, 4.9 bp/read) — the rule trims short low-quality runs too, so it
fixes both boundary overtrim and min-run undertrim there.

### Setting B — quality from CHOSEN constants (probe worlds)

| world | best q_cut | tol0 (exact/over/under) | tol2 exact | what it is |
|---|---|---|---|---|
| strong-Q2 (artefact run at Q2) | 5 | mid/low 100.0 / 0 / 0; high 99.7 / 0 / 0.3 | 100.0 / 99.8 | the 2-colour "no signal" world |
| decay (biology 37→22, artefact local−12) | 15 | mid/low 100.0 / 0 / 0; high 99.7 / 0 / 0.3 | 100.0 / 99.8 | decay stress |
| overlap (both noisy/decaying, distributions overlap) | 20–25 | mid 88.6 / 0 / 11.4 (q_cut 20); 82.9 / 10.0 / 7.1 (q_cut 25); high 96.9 / 2.5 / 0.6 (q_cut 25) | 98.1–99.8 | realistic-noise stress |

These numbers **do depend on my constants** — they show what the rule would
give *if* real artefact quality behaves like the construction. The as-is
table in Setting A is the one that requires no belief in my construction.

## 5. False-positive cost (the rule's downside)

At the best q_cut, on mates fastprep had EXACT at tol 0 (Setting A):

| dataset | -> undertrim (artefact left) | -> overtrim (biology cut) |
|---|---|---|
| low  | 0 (0.0%) | 11 (0.0%), 1 bp each |
| mid  | **3353 (11.3%), 1.00 bp each** | 11 (0.0%), 1 bp each |
| high | 2 (0.0%), 6 bp each | 3 (0.0%), 1 bp each |

The mid cost is structural: those are exactly the h = 1 head artefacts whose
single base the model gives normal quality — a quality rule must leave them
(artefact left behind) or give up the whole 5' win. In the overlap world the
other direction appears (high: 675 mates = 9.8% of baseline-exact get 1 bp
of biology cut at q_cut 25), i.e. when qualities overlap the threshold pays
for overtrim removal with undertrim and vice versa.

## 6. Why the local baseline is not 66.5 / 29.5 / 4.0

The Phase-D reference (desktop, 101 fitted real-library profiles,
fastprep-winner: tol 0 = 66.5/29.5/4.0, tol 2 = 98.4/1.3/0.3) is quoted as
external context only; those datasets are on the unreachable desktop. The
local fixtures are regime-shifted: `high` puts 4–20 bp tails on ~90% of
mates against `min_run 12` (hence 53.8% undertrim), `low`/`mid` land at
81% exact / 17–18.5% overtrim. All before/after numbers above are measured
against the **same-regime baseline**, so the deltas are internally valid;
the levels are not comparable to Phase-D.

## 7. Throughput cost

`--part cost`: synthetic 100k pairs x 2 mates, 150 bp, with low-quality end
runs on 35% of mates (the rule's maximal walk), gzip level 6, 2 threads,
best of 3:

```
fastprep (rule off)        0.56 s   355,224 mates/s
fastprep (--q-boundary 25) 0.57 s   350,377 mates/s   ->  -1.4%
```

Per read the rule adds at most two comparisons per run base walked and
stops at the first non-G or high-quality base (typically 0–2 steps per
end; worst case the run length). No look-ahead. The single streaming pass
structure is unchanged.

## 8. Bottom line

* fastsim's model **does** carry a discriminating quality step at the 3'
  tail boundary (3 Phred, from the fitted `tail_q_collapse`), and a quality
  rule evaluated on **fastsim's own output** converts most tail-driven
  overtrim into exact cuts (tol-0 exact +9 pp on mid, +77 pp on high).
* fastsim's model **does not** carry a signal at the dominant 5' boundary
  (h = 1: first head base is normal quality by construction); the rule
  cannot help there and its whole false-positive cost (11.3% of
  baseline-exact mates leaving 1 bp of artefact on mid) is exactly that
  case.
* Whether real libraries carry a tail-boundary step strong enough for a
  threshold rule is **unverified** — real noise is larger than the
  simulator's and the measured collapse varies 1.7–11.7 Phred by run. A
  real-data test (Phase-D reads + GT) is the one missing experiment.
* Recommendation: the rule is cheap (-1.4%) and correct in direction, but
  do not claim the fastsim-measured gains for real data. If adopted, gate
  it behind a flag (as `--q-boundary` is) and settle `q_cut` on real data.
