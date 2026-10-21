# FastPrep exactness: exact where provable, overtrim where not

Status: PROPOSAL — nothing implemented, no default changed. 2026-10-21.

## 0. The ask

Make FastPrep's trimming as exact as possible; where a boundary decision is
genuinely ambiguous, resolve it toward **overtrimming** (the user's standing
cost model: "we prefer less data than wrong data").

## 1. Diagnosis — where the inexactness actually lives

Measured (101-dataset corpus, exact simulator ground truth, gzip 6 / -j 4,
`bench/scripts/trim_accuracy.py`; see `references/measurement-pitfalls.md`
"Trim accuracy"):

| tool | tol 0 exact / over / under | tol ±2 exact |
|---|---|---|
| fastp-default | 92.2% / 0.2% / 7.6% | 99.6% |
| fastprep-winner | 66.5% / **29.5%** / 4.0% | 98.4% |

The 29.5% overtrim collapses to 1.3% at ±2 → it is **boundary jitter**, not
biology loss. So "make it exact" is mostly "kill the jitter", and the residual
undertrim (4.0%, mostly adapter readthrough the matcher never sees) is where the
overtrim bias must bite.

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
nothing-to-cut. Expected effect: most of the 29.5% tol-0 overtrim → exact. This
is a **user-visible default change**: measure it (Step 4), document the trade at
the flag, hold the flip for the user's OK with the version bump (0.1.0 → 0.2.0).

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
colour), so this is base-specific by construction. Applies symmetrically to the
5′ head (head_dropout and the forward scan).

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
secondary column only. The per-axis raw counts never depend on the weights.

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
2. **Two-axis corpus sweep** on the existing 101-dataset GT corpus: each rule
   measured alone and combined; report tol 0 and tol ±2 exact/over/under,
   false-trim on clean reads, and the tolerance-decay split (jitter vs
   deliberate). Compare against fastp-default (the current exactness leader,
   92.2% tol-0) head-to-head.
3. **Re-score the closed 3′ band** (the −35,196/+236,779 trade) under the
   user's weights before deciding the band stays closed.
4. Acceptance (pre-registered):
   * tol-0 exact rate up substantially (the jitter class must convert);
   * undertrim (tol 0 and ±2) ≤ current 4.0%, and the residual adapter-
     readthrough class (scored by `undertrim_ab.py`) down;
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
  the machine identity, and the 66.5/29.5/4.0 baseline is re-derived on the
  desktop before it is used as an acceptance yardstick.
* The corpus sweep re-scores existing sims (GT files exist) — a handful of
  configs × 101 datasets is **low core-hours; overnight on the desktop with
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
