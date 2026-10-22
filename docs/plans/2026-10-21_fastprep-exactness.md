# FastPrep exactness: exact where provable, overtrim where not

Status: PROPOSAL — nothing implemented, no default changed. 2026-10-21.

## 0. The ask

Make FastPrep's trimming as exact as possible; where a boundary decision is
genuinely ambiguous, resolve it toward **overtrimming** (the user's standing
cost model: "we prefer less data than wrong data").

## 1. Diagnosis — where the inexactness actually lives

**Baseline corrected 2026-10-22.** The round this section was first drafted
from (`bench/results/phase-d/trim-accuracy.tsv`: fastprep-winner 66.5 / 29.5 /
4.0) is **invalid**: `docs/TRIM-ACCURACY-DISCREPANCY.md` shows it scored
sim-generation-1 trim decisions against sim-generation-2 ground truth, so its
29.5% / 4.0% split is the ±1-2 bp spread of two independent draws, not trimming
behaviour. The authoritative, generation-consistent re-measurement
(`bench/results/phase-d/overtrim-vs-fastp.tsv`, 102 datasets / 20.4 M mates),
tol 0:

| tool | tol 0 exact |
|---|---|
| fastp-default | 92.19% |
| fastp-tuned | 84.06% |
| fastprep-hd0 (head-dropout 0) | 92.60% |
| fastprep winner (shipped default) | 75.59% (75.35% after the `_minlen` repair) |

fastprep-winner's tol-0 over/under split: **23.4% / 1.0%**. fastp-default's
0.2% / 7.6% split reproduces identically in both rounds, so it was never
affected. Every bucket at both tolerances is in the TSV.

What survives the correction: at tol 0 fastprep puts ~23% of mates off-boundary
while losing only ~1.0% to artefact left behind, so overtrim still carries
essentially the whole ~17 pp exactness gap against fastp, and the tie-break
toward cutting is still the lever. What does NOT survive is the specific jitter
claim ("29.5% collapses to 1.3% at ±2"): that decay was measured on the invalid
generation. The jitter-vs-deliberate split must be re-measured on the desktop
(Step 4) from tolerance decay **and** the margin-on/off delta before any number
in this section is quoted again.

**The thesis: the algorithm's tie-breaks currently resolve ambiguity toward
KEEPING bases, while the user's cost model resolves toward CUTTING them.**
Three concrete decision points encode the wrong bias:

1. **The mercy walk overshoots into the flank and the default cut keeps it.**
   `algorithm.rs::scan` accumulates gap mass until `gamma > mercy` (Δ_m = 3.0);
   `resolve_trim` on `StopMode::Mercy` cuts `k_star`, the full walk, which
   absorbs ~2 bp of flanking non-run biology past the true run boundary
   (measured). `Params::boundary_anchored` already exists and cuts `run_span`
   (the run extent, keeping error bases *inside* the run) instead — but it is
   **default `false`** and has never been measured on both axes as a default.

2. **Nested corroborated adapter prefixes resolve shortest-wins.**
   `adapter.rs::short_tail_trim` returns the *shortest* corroborated exact
   prefix ("the more conservative claim about how much of the read is
   adapter"). Under the user's cost model the longest consistent candidate is
   the defensible cut: a nested longer exact prefix is either more adapter or
   an exponentially unlikely coincidence.

3. **PE-overlap ambiguity refuses = no trim.**
   `pe_overlap.rs::detect_with_rc2` returns `None` when more than one insert
   verifies ("refuse rather than guess"). When readthrough is real, refusing is
   the **undertrim** choice — exactly the more expensive error.

The mismatch budget edge (`match_adapters`, `error_rate` 0.2 / `min_overlap_se`
14) is NOT in scope: the floor/eps sweep already showed the residual-optimal
point buys 365× the false-trim rate. The UT-2/UT-3/UT-4 retry chain is the
instrument for that class and stays.

**Premise-census guardrails that bound the whole design** (`premise_census.py`,
20.4 M mates): 3′ A/C/T terminal runs are artefact once in ~2,000 and N runs
never — any rule below is **per-base (G) or adapter-specific, never generic
polyX**. The 3′ terminal-run band's closure verdict ("the undertrim direction
is closed", +0.094 pp recoverable vs 236 k overtrim) was computed on the
*symmetric* exact-rate. Under an asymmetric cost it must be **re-derived**, not
reused: the measured trade (−35,196 undertrim bought +236,779 overtrim) turns
net-positive at w_under/w_over ≳ 6.7. That re-scoring is Step 5 below.

## 2. The design: three tiers

### Tier 1 — exact by construction (no guessing involved)

* PE-overlap pinning already wins over sequence matching (`trim()` Step 1:
  `params.pe_adapter` short-circuits the matcher). Keep it decisive. Verify the
  detector/verifier threshold-cancellation class is closed (the verifier must
  test prefix identity at the pinned position, never re-route through
  `min_overlap_se` — see `references/measurement-pitfalls.md`).
* UT-3 corroborated exact prefix (quality-cliff or dataset-level
  `prefix_enrichment` corroboration): the cut lands exactly where the exact
  prefix ends. Keep.

### Tier 2 — exact by measurement: ship `boundary_anchored` as default

`resolve_trim` / `resolve_subscan` (algorithm.rs):

```rust
} else if boundary_anchored {
    res.run_span            // run extent: start ..= last matching base
} else if res.mode == StopMode::Mercy {
    res.k_star              // today's default: includes the flank overshoot
}
```

Flipping the default converts the mercy overshoot (~2 bp) from cut-signal into
nothing-to-cut. Expected effect: a large share of the valid 23.4% tol-0
overtrim → exact. That is a **hypothesis until Step 4 measures it** — the
invalid round's tolerance decay cannot stand in for the measurement. This is a
**user-visible default change**: measure it, document the trade at the flag,
hold the flip for the user's OK with the version bump (0.1.0 → 0.2.0).

### Tier 3 — the ambiguity band, resolved toward overtrim

Wherever the algorithm has a genuine either/or and no evidence settles it, take
the larger cut. Each rule is **bounded, per-base or adapter-specific, and
recorded in PXMT metadata** so the bias stays auditable per read.

**(a) Run-flank extension** (`resolve_trim`, both cut paths). After the
run-extent cut, extend through *immediately adjacent* bases that (i) match the
run base (case-insensitive) and (ii) carry `q < phred_sat` (an unconfident base
is more likely artefact continuation than a biological boundary base — the
simulator's own boundary rule forces the base at the true boundary to be
non-run, i.e. the only ambiguous class on real data). Bounded at +2 bp. Never
generic: the run base is whatever the scan matched (G at 3′ on NovaSeq two-
colour), so this is base-specific by construction. Applies to the forward and
reverse scan cuts at both ends. **Excluded by construction from the
`--head-dropout` cut** (implementation note, 2026-10-22): head-dropout trims the
*maximal* case-sensitive leading G-run, so the base immediately past that cut can
never equal the run base — the only continuation the margin could match is a
lowercase `g`, which real FASTQ does not carry. That boundary is exact by
construction, so there is no ambiguity for a margin to resolve; the margin stays
on the scan cuts, where the boundary genuinely is ambiguous.

**(b) Nested adapter prefixes: longest-wins.** In `short_tail_trim`, iterate
candidates longest-first and return the first corroborated one instead of the
shortest. One-line semantic change behind a flag (below).

**(c) PE-overlap multi-verifier: smallest verifiable insert.** In
`detect_with_rc2`, when >1 `o` verifies, do not blanket-refuse: pick the
smallest admissible `o` (= shortest insert = largest readthrough trim) **only
if** `readthrough_is_adapter` confirms adapter behind that boundary on the
pinned mate; otherwise keep today's refusal. The adapter check is what keeps a
repeated-fragment false pin from trimming biology.

**(d) Keep the existing +1 at `StopMode::Boundary`** (`k_star + 1`) — already
overtrim-leaning, unchanged.

**Guardrails (non-negotiable):**

* The margin consumes the trim budget (`extra_committed` / `MinLen::floor`):
  a margin that pushes a read under `--min-len` must drop it via the normal
  `DropReason::Length` path, never silently over-cut past the floor.
* False-trim rate on **clean reads** (reads whose GT window is the whole read)
  is a pre-registered ceiling on every rule, measured per rule (per the
  two-axis discipline: residual on carrier reads AND false trims on clean
  reads; a one-axis win is rejected).
* No rule touches A/C/T terminal runs or N runs (census verdict).
* Every extension is counted and reported, but **not through PXMT**: the v3
  record layout is unpacked at fixed offsets by 25+ bench scripts, and a field
  insertion would silently misalign every one of them. The margin is attributed
  from the config (the margin flag is known per run) and aggregated in the JSON
  report (`margin_bases_5p`/`margin_bases_3p` + `reads_with_margin`); PXMT
  stays byte-stable. The ≤2 bp class is split jitter-vs-margin from the
  margin-on/margin-off delta per dataset — tolerance decay alone cannot tell
  them apart, because both decay at ±2. A per-read PXMT v4 sidecar stays a
  separate, explicitly scheduled change if the report ever needs it.

**CLI surface** (all default-off until the sweep clears them; the shipped
default is the swept winner and defaults change only with the user's OK):

* `--boundary-anchored` already exists (flip candidate, Tier 2).
* `--overtrim-margin N` (Tier 3a; N ≤ 2).
* `--short-tail-longest` (Tier 3b).
* `--pe-ambiguous-smallest` (Tier 3c).
* **Decided 2026-10-21 (implementation time): three separate flags, no
  umbrella.** The sweep must arm each rule alone and in combination, so the
  flags double as the measurement knobs; a `--prefer-overtrim` umbrella can be
  added later as sugar once the swept N is known.

## 3. Scoring: the cost model must be explicit

`score = (exact + drop_correct − Σ wᵢ·errorᵢ)/pairs` with
`w_under > w_over` (the standing preference; `references/trim-scoring.md`).

**Anchor decided (user, 2026-10-21): w_under : w_over = 7 : 1.** The sweep is
reported on the ladder r ∈ {1, 3, 5, 7, 10} plus a per-rule break-even
r* = −Δover/Δunder ("this rule pays off above r*"), with 7 as the
ship/no-ship anchor. Weights live in `bench/weights/trim-weights-v1.json`
(written before the first run so the margin cannot be tuned to its result);
every result row records `weights_version = v1`. Scoring is per-read (primary,
comparable to the existing trim-accuracy numbers); a per-base variant is a
secondary column only. The per-axis raw counts never depend on the weights. **Scale warning:**
`score` divides the numerator by **pairs** (the pre-registered model), while the
numerator counts **mates**, so a perfect paired run reads 2.0, not 1.0. Ladder
rankings, every Δ and every `r*` are invariant to that factor, but never compare
`score` to a per-mate percentage (75.59%) — compare the per-axis columns.

Known gap to close before the first weighted row: `trim_accuracy.py`'s
`classify()` emits five buckets and leaves `dropped_or_missing` unweighted,
while this formula needs `drop_correct` vs `drop_wrong` (from the GT
`full_flag` and the min-len floor). `bench/scripts/score_weighted.py` reads
`results/trim-accuracy.tsv` + the GT and emits the weighted score at each r —
that is what makes r = 7 actionable. **The sweep runs over CONFIGS once and
scores at all five r values post-hoc; never multiply the run matrix by five.**

## 4. Verification plan (all on the desktop; ground truth = fastsim corpus)

1. **Boundary-sweep fixture** (regression gate, cheap, local): synthetic reads
   with the artefact boundary at *every* position (adapter slid across every
   offset; G-run boundary at every index; run merged with a biological
   boundary G vs a non-G boundary), assert the recovered cut equals the GT cut
   at every offset for the unambiguous cases and equals `true + margin` for the
   ambiguous class. Catches off-by-ones before any corpus run.
2. **Two-axis corpus sweep** on the GT corpus as it exists on the desktop
   (306 datasets = 102 accessions × 3 variants, `/mnt/Fast/fastsuite-scratch/
   work-gt2026`): each rule
   measured alone and combined; report tol 0 and tol ±2 exact/over/under,
   false-trim on clean reads, and the tolerance-decay split (jitter vs
   deliberate). Compare against fastp-default (the current exactness leader,
   92.19% tol-0 in the valid round) head-to-head.
3. **Re-score the closed 3′ band** (the −35,196/+236,779 trade) under the
   user's weights before deciding the band stays closed.
4. Acceptance (pre-registered):
   * tol-0 exact rate up substantially from the valid 75.59% baseline;
   * undertrim (tol 0 and ±2) ≤ the re-derived current 1.0% (no regression),
     with the residual adapter-readthrough class (scored by `undertrim_ab.py`)
     down;
   * clean-read false trims under the per-rule ceiling;
   * overtrim surviving at tol ±2 = only the deliberate margin, reported as
     such.

## 5. Cost / time estimate (order of magnitude, confirm from measured unit cost)

* Per-config evaluation: ~15–17 µcore-s/pair (fastprep + tier-1 score) →
  100 k pairs ≈ 1.6 core-s; scoring itself ~0.05 core-s/100 k.
* Boundary-sweep fixture: minutes (unit tests with the suite on the laptop;
  re-confirmed on the desktop before any of its numbers are quoted).
* **Venue (user rule, 2026-10-21): every bench row is produced on the desktop,
  exclusively.** The laptop is edit + `cargo test` only. Each bench row records
  the  machine identity, and the acceptance baseline is re-derived on the desktop
  before it is used as a yardstick: the only valid repository numbers are
  75.59% / 23.4% / 1.0% from `overtrim-vs-fastp.tsv` (see §1), and even those
  are re-checked against `work-gt2026` before use. The retired 66.5/29.5/4.0
  round is not an acceptance value at all.
* The corpus sweep re-scores existing sims (GT files exist) — a handful of
  configs × 306 datasets is **low core-hours; overnight on the desktop with
  3–4 parallel subagents is the honest continuous estimate**. The sims need
  regenerating only if we also change fastsim (we do not).
* Implementation + TDD: one session per tier, committed per task.

## 6. Explicit non-goals

* No lowering of `--min-adapter-overlap` / raising of `--adapter-error-rate`
  (the floor/eps sweep's residual-optimal point is a 365× false-trim
  regression).
* No generic polyX threshold changes (census: A/C/T/N closed).
* No hard trim for exactness (a mandatory cut makes `exact` structurally
  unreachable — the scorer's own warning).
* No default flips without the user's OK + version bump.

## 7. Measured result — 2026-10-22, desktop only

`bench/scripts/exactness_sweep.py`, 102 datasets (v0) × 12 arms = 1 224 rows, one
host (`cachyos-x8664`), one binary (md5 `defb04d4`), winner config per dataset,
scored by `score_weighted` at the 7:1 anchor. Raw:
`bench/results/exactness/sweep-raw.jsonl`; rows:
`sweep.tol0.tsv` / `sweep.tol2.tsv`.

| tol 0, per mate | exact | overtrim | undertrim | clean-read cuts | score @7 | Δ vs base |
|---|---|---|---|---|---|---|
| base (hd 4, as shipped) | 72.04% | 19.03% | 1.87% | 20.37% | 0.6628 | — |
| `--boundary-anchored` | 72.35% | 18.71% | 1.88% | 20.37% | 0.6749 | **+0.0121** |
| `--overtrim-margin 1` | 72.04% | 19.03% | 1.87% | 20.37% | 0.6628 | 0.0000 |
| `--overtrim-margin 2` | 72.04% | 19.03% | 1.87% | 20.37% | 0.6628 | 0.0000 |
| `--short-tail-longest` (vs `st6`) | 72.04% | 19.03% | 1.87% | 20.37% | 0.6628 | 0.0000 |
| `--pe-ambiguous-smallest` (vs `peov`) | 72.03% | 19.04% | 1.87% | 20.38% | 0.6626 | 0.0000 |
| `--head-dropout 2` | 72.50% | 18.20% | 2.24% | 19.45% | 0.6373 | −0.0255 |
| `--head-dropout 1` | 74.13% | 15.12% | 3.69% | 16.05% | 0.5290 | −0.1338 |
| `--head-dropout 0` | **82.80%** | 0.57% | 9.58% | 0.02% | 0.1680 | −0.4947 |

At tol ±2 the whole spread collapses: base 91.58% exact, hd0 92.06%, every
other arm within 0.01 pp.

**What this says, in order of how much it matters:**

1. **The exactness question on this corpus is the 5′ G-run policy, not boundary
   resolution.** A decomposition over 1.2 M mates put **100% of the overtrim at
   the 5′ end and 80% of it at exactly 1 bp**, with a 1-2-3-4 bp histogram that is
   `--head-dropout`'s own 1..=4 window. The 3′ polyX boundary — the target of
   Tier 2 and Tier 3a — produces **40** overtrimmed mates in 1.2 M (0.003%).
2. **The user's 7:1 anchor already picks the shipped policy.** Turning
   head-dropout off buys the highest strict exactness available (82.80%, +10.8 pp)
   by converting 18.5 pp of overtrim into 7.7 pp of undertrim — and at 7× that is
   −0.49 per mate, a decisive loss. hd 2/1 are also negative. `hd = 4` is the
   measured optimum of the shipped knob at this anchor; no change is needed.
3. **Tier 2 pays, barely: +0.0121 at r = 7** (0.32 pp of mates move from
   overtrim to exact, no undertrim cost, clean-read cuts unchanged). It is free
   and it is real, but it is 40× smaller than the head-dropout trade above.
4. **Tier 3a/3b/3c are 0.0000 on this corpus, and correctly so.** The margin
   never fires because its target class is 0.003% here; the short-tail and
   PE-ambiguity arms cannot fire because adapter readthrough is 0.01%. They are
   not wrong, they are simply unpriced by this corpus — their classes live in the
   synthetic readthrough cohorts (`OVERTRIM-PAIRED.md`, `undertrim_ab.py`).
5. **A quality gate cannot rescue the 1 bp class.** The leading G-run's mean Q is
   statistically identical for the reads the GT calls biology (median 36–37,
   p10 5) and the reads inside artefact windows (median 36–37, p10 5) — the
   simulator gives artefact bases normal quality, by design. There is no
   per-read evidence separating "dark cycle" from "biology" here, which makes the
   1 bp class a genuine either/or, and the 7:1 anchor is the only thing that
   should decide it. That is exactly what point 2 says it does.

**Recommendation, therefore:** flip `--boundary-anchored` to default (needs the
user's OK + the version bump, per §2), keep `--overtrim-margin`,
`--short-tail-longest` and `--pe-ambiguous-smallest` default-off as instruments
for corpora that contain their classes, and do not touch `--head-dropout`: at
the 7:1 anchor the shipped value is already the optimal one, and the honest
statement of the exactness ceiling on real-library data is "82.8% strict
exactness, 0.02% clean-read damage, paid for with 9.6% undertrim" — a trade this
cost model rejects.
