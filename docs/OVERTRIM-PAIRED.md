# Can paired-end structure resolve the artefact/biology G-run ambiguity?

Question. fastprep trims the whole terminal G-run it observes; when the
biology next to the run is itself G's, those G's are trimmed too. Base
identity cannot tell `h` artefact G's + `k` biological G's from `h+k`
artefact G's (docs/OVERTRIM-DIAGNOSIS.md). This note tests whether the
**mate** can: (a) a run position the other mate corroborates as G is
biology, not artefact; (b) knowing the insert `I` pins the run beyond the
fragment.

Reproduce: `python3 bench/scripts/overtrim_paired.py all --out summary.json`
(probe / sim / bench subcommands). Findings below are its measured output.

## Verdict

**The rule works when it applies, but it barely applies.** On the controlled
probe the mate-overlap rule removes the ambiguity *completely*; on the real
(fitted-profile) insert distributions it reaches about **1% of the overtrim
events**, because the overtrim sits at the mates' 5' ends — the fragment's
*outer* boundaries — and a mate only covers the fragment's outer boundary
when the insert is no longer than the read. On the four fastsim datasets the
strict-tolerance exact rate moves **61.3% -> 61.4%** (tol 0, pooled).
Paired structure is a targeted fix for **short-insert libraries**, not a
general cure.

## 1. The controlled probe: the rule resolves the ambiguity

Synthetic pairs, known truth, forced ambiguous boundaries (the biology edge
next to each artefact run is a G- or C-run of 0..2 bases — the exact case
`/tmp/probe_overtrim.py` established as the error source). Scoring follows
trim_accuracy semantics (correct cut = exactly `h` at 5', `t` at 3').

| cohort (geometry) | fastprep tol 0 | + mate-overlap add-back | fixed | broke |
|---|---|---|---|---|
| head_full (I = 150−h ≤ L, full overlap) | exact 60.2% / over 39.7% / under 0.1% | **exact 99.9% / over 0.0% / under 0.1%** | 476 | 0 |
| head_partial (L < I = 220 < 2L) | 62.3 / 37.5 / 0.2 | 62.3 / 37.5 / 0.2 (abstains) | 0 | 0 |
| head_none (I = 380 ≥ 2L) | 61.5 / 38.4 / 0.1 | 61.5 / 38.4 / 0.1 (no overlap) | 0 | 0 |
| tail_abuts (I = 138 = L−t, fragment-end regime) | 0.0 / 100 / 0 | 0.0 / 100 / 0 at tol 0; tol 2: 27.9% -> 50.3% | 0 (see §4) | 0 |
| tail_under (I = 220, fragment under the artefact) | 0.0 / 100 / 0 | 1.5 / 98.2 / 0.3 (tol 2: 55.6%) | 18 | 4 |

Where the mate sees the boundary (`head_full`: both mates carry the whole
fragment), the add-back rule is exact: it recovers every biological G and
trims only artefact, at both tolerances, with **zero breakage**. Where the
mate cannot see the boundary it abstains and nothing changes. The one real
error mode is `tail_under`: an artefact slot that hides a genuine G in the
fragment is indistinguishable from a biological G even with the mate
(see §4).

## 2. Why the gain on real libraries is tiny

The measured overtrim is a 5'-end phenomenon. On the mid dataset the 5' end
carries 8474 of 8568 overtrim events (98.9%; repro: `repro_base.py` pattern,
production defaults), and 93% of the excess is 1-2 bases — exactly the
biology-begins-with-G class of the diagnosis.

A read's 5' end is the fragment's **outer** boundary (R1: fragment 5' end,
R2: fragment 3' end). The other mate covers that position only when it reads
the whole fragment, i.e. **I ≤ L** (minus its own artefact run). The insert
distributions the profiles carry (bench/fixtures/overtrim/*.json):

* P(I < L = 150) = **4.7%**, P(I < 2L = 300) = 82.2%.
* For L ≤ I < 2L the mate covers only the fragment *interior* — never the
  outer edge where the ambiguity sits.
* For I ≥ 2L there is no overlap at all and no per-pair insert measurement.

Measured overtrim events by insert class (sim, per regime, tol 0):

| regime | I < L | L ≤ I < 2L | I ≥ 2L | total |
|---|---|---|---|---|
| low | 81 | 7650 | 1749 | 9480 |
| mid | 83 | 6805 | 1664 | 8552 |
| high | 102 | 7927 | 1816 | 9845 |
| mandril | 6 | 4587 | 7713 | 12306 |

Only 0.05-1.1% of overtrim events fall in the class the mate can fully
verify. Consequently (fastsim datasets, 20k pairs each, production defaults,
--no-dedup; scored exactly like bench/scripts/overtrim_sensitivity.py):

| regime | before tol 0 (exact/over/under) | after | before tol 2 | after | fixed | broke |
|---|---|---|---|---|---|---|
| low | 75.0 / 24.7 / 0.3 | 75.3 / 24.3 / 0.3 | 98.1 / 1.6 / 0.3 | 98.1 / 1.5 / 0.3 | 123 | 0 |
| mid | 76.1 / 22.3 / 1.6 | 76.4 / 22.0 / 1.6 | 97.0 / 1.7 / 1.3 | 97.0 / 1.7 / 1.3 | 116 | 27 |
| high | 20.4 / 29.6 / 50.1 | 20.5 / 29.3 / 50.2 | 28.7 / 21.2 / 50.1 | 29.1 / 20.7 / 50.2 | 45 | 50 |
| mandril | 67.8 / 30.8 / 1.4 | 67.9 / 30.7 / 1.4 | 97.2 / 1.7 / 1.1 | 97.2 / 1.7 / 1.1 | 20 | 9 |
| **pooled** | **61.3 / 26.8 / 12.0** | **61.4 / 26.5 / 12.0** | **82.2 / 6.0 / 11.8** | **82.3 / 5.9 / 11.9** | **304** | **86** |

(The published figure in docs/OVERTRIM-DIAGNOSIS.md — 66.5 / 29.5 / 4.0 at
tol 0 — is the winner config's number; it is the *published* cross-generation
figure, since shown invalid — the corrected shipped-default figure is
75.4 / 23.4 / 1.0, see `docs/TRIM-ACCURACY-DISCREPANCY.md`; the closest reproduction here is
mandril with production defaults: 67.8 / 30.8 / 1.4. Configs differ in
head-dropout and dedup; the before/after comparison above is internally
consistent on one config.)

304 of ~40 000 overtrim events are fixed and 86 mates are pushed into
undertrim — a net strict-tolerance gain of **+0.1 pp exact**, i.e. within
noise of doing nothing. The rule fires on 18 219 mates (11.4%) but almost
never at the boundary that matters: it corroborates interior run positions
whose trims were already correct.

**Transfer to real libraries.** The gain scales with P(I ≤ L):

* short-insert libraries (small-RNA, cfDNA, heavy readthrough: most pairs
  I ≤ L): the applicable fraction approaches the whole class and the probe
  says the fix is exact (39.7% overtrim -> 0.0%);
* the fitted profiles here (P(I<L) = 4.7%): ~1% of overtrim events;
* long-insert libraries (mandril-like, P(I<L) ≈ 0): nothing — mandril fixes
  20 of 12 306 events.

## 3. Rule (b), insert arithmetic: necessary but not sufficient

Knowing I in the fragment-end geometry gives `h + t = L - I` per mate. That
bounds the split but does **not** determine it: every `h` in
[max(0, L−I−r3), min(r5, L−I)] is self-consistent with some fragment — the
same ambiguity, one dimension smaller. The insert resolves the boundary only
when one end carries **no** run at all (that end's artefact is 0, so the
other end is pinned). Measured (probe): the insert rule fixes 151 mates in
head_full and 400 in tail_abuts with 0 breakage; it fires nowhere else.

fastcheck already measures I per pair from the read overlap
(`tools/fastcheck/src/insert.rs`: `I = b1 + b2 - overlap`, `b = L - h - t`
after terminal homopolymer trimming, mismatch tolerance 3, censored below/
above reported). A faithful re-implementation of that identity
(`measure_insert` in the bench script) validates on the probe:

* fragment-end regime (tail_abuts): measured on 579/600 pairs, **exact on
  577/579**;
* merged-run regime (head_full): measurable on 402/600, **exact on only
  195/600** — the homopolymer trimming counts artefact+biology G's as one
  run, so `b` is undercounted by the merged G's and the identity is biased.

So the rule is implementable with existing machinery (fastcheck for I,
fastprep's PXMT for the cuts), no new measurement required — but its reach
is capped by the overlap, not by the tooling.

## 4. Costs and pitfalls of the rule (measured)

1. **Collapse-regime undertrim (inherent).** When the fragment continues
   under the 3' artefact slots (fastsim's `bio = min(I, L−h−t)` case), an
   artefact slot hiding a real fragment G is corroborated by the mate
   exactly like a biological G. P ≈ 25% per slot; the add-back then keeps
   it and undertrims (probe tail_under: 18 fixed vs 4 broke; the 3'-scan
   overshoot of §4.3 happens to absorb most of it). No rule can be right in
   both readings — the input genuinely does not determine the answer there.
   In the fragment-end regime (I = L−h−t) the slots have no fragment base
   and the rule is exact.
2. **Adapter pads falsely corroborate head artefact slots.** A head slot's
   mate-side position can land in the mate's adapter pad; the pad's first
   base (adapter `CTGT...`, C) complements to G and "corroborates". On the
   sim data breakage concentrates on pad-bearing mates (mid: 21 of 27 broke
   mates have pads; high: 41 of 50). A production rule must void evidence
   inside adapter windows — fastprep already detects adapters, so the
   information is available. The probe cohorts are pad-free by construction
   for this reason.
3. **The 3' scan overshoots runs ≥ min_run (12) by 2-5 bases**, even on a
   clean boundary, and the overshoot depends on flank G content (uniform
   Q40: a pure-C flank + 20 G run is not trimmed at all; a random flank +
   12 G run is trimmed 14-17). This is *not* the boundary ambiguity — it is
   the quality-weighted scan's own boundary placement — and it is what makes
   `tail_abuts` 0% exact at tol 0. It is largely invisible in the sim stats
   because profile tails are mostly < 12. The diagnosis doc's "exact when
   the boundary is not a G" was established on the 5' head-dropout path
   (runs 1..4); the 3' scan path behaves differently and deserves its own
   investigation.
4. **Scorer convention.** Adapter trims count in the observed t3
   (trim_accuracy: correct cut is `h` and `t`), so a handful of the 3'
   "overtrim" events on the sim data (magnitudes 2-15) are detected adapter
   trims, not G-run errors. Both before/after columns carry this equally.

## 5. Cost / throughput

100 000 synthetic pairs (bench subcommand):

* fastprep paired path: 0.20 s = **493 345 pairs/s** wall (4 threads);
* fastprep single-mate path: 0.18 s = 563 912 reads/s;
* add-back rule, this Python prototype (1 thread, exhaustive shift search
  for the mate correspondence): 146.2 s = **1461 µs/pair**.

The prototype number is dominated by the O(L²) brute-force correspondence
search — an implementation choice, not the algorithm. The algorithmic work
is one shift search (fastcheck's 12-base anchor + verify is O(L)) plus an
O(run) walk: a native implementation is well under 1 µs/pair against
fastprep's ~2 µs/pair wall (4 threads) / ~8 µs/pair single-thread scan
cost, i.e. **low single-digit percent on the paired path and exactly zero
on the single-mate path**, which never runs the rule (it needs both mates).

## 6. Recommendation

Do not sell the paired rule as a fix for the 29.5%. On realistic insert
distributions it removes <1% of the strict-tolerance overtrim and adds a
small undertrim tail (adapter pads, collapse slots). It *is* a complete,
zero-breakage fix for the short-insert regime (I ≤ L on most pairs), where
the whole overtrim class is reachable — worth having as a paired-mode option
for small-RNA/cfDNA libraries if such libraries enter the corpus. The
dominant class — a 1-2 bp biology-begins-with-G overtrim on I ≥ L pairs —
remains information-theoretically ambiguous with the data these reads
carry, which is the diagnosis doc's conclusion, now with the paired
dimension measured too.

---

# Addendum (2026-10-19): the mode was rebuilt and measured with/without

The sections above answer the **disambiguation** question: can the mate resolve
the artefact/biology G-run ambiguity? This addendum answers the different
**accuracy** question — does a paired-end overlap mode improve trimming? — by
building the mode behind a flag and measuring it, rather than reasoning about it.

## What was built

`fastprep --pe-overlap` (default off). Detection is anchored and verified
(`tools/fastprep/src/pe_overlap.rs`, 12-base anchor, full-window verification,
ambiguity refused) and the boundary is applied only when the adapter **starts
exactly at** it. That last condition is the whole safety property:

> **The PE mode can never trim something the ordinary adapter path could not
> itself confirm.**

An overlap proves where the fragment ends; it does not prove that anything was
ligated past it. The first implementation got this wrong and trimmed biology —
see "The defect" below.

## Cell A — the spurious-fire guard (real corpus, 102 datasets, 10.2M pairs)

The guard: on data with no adapter readthrough, the mode must not change any
trim decision.

| | pre-fix | post-fix |
|---|---|---|
| `pe_trims` (pairs whose trim decision changed) | 50 | 18,436 |
| new overtrims vs `fastprep-hd0` on a 200k slice | 46 | — |
| lost exact mates on that slice | 48 | — |

The pre-fix detector applied the overlap boundary as an adapter trim
unconditionally, so any mate whose fragment was shorter than its read was
trimmed to the insert **even with nothing ligated past it**. That is the same
class of failure as the 4-base-overlap false trim that `--min-adapter-overlap 8`
exists to prevent (~10% of clean reads, per its own help text), which is why the
guard is not thresholded: any firing on adapter-free data is disqualifying.

**The corpus cannot answer the value question**, for two independent reasons:

1. Adapter readthrough is configured in **53 of 102 profiles** but materialises
   in only **11,395 of 10.2M pairs (0.112%)** — it needs an insert shorter than
   the read, which these libraries rarely have. A 0.5 pp threshold on 0.11% of
   pairs measures noise.
2. The scorer compares the tool's trim positions against the ground truth's
   **G-run** lengths. An adapter trim is not a G-run, so on a cohort whose GT
   carries no runs, every adapter trim is labelled `overtrim`. Cell A's
   readthrough slice shows 8,533 of 11,299 pairs "overtrimmed" — a scorer
   artefact, not a detector verdict.

`pe_trims = 18,436` against an 11,395-pair readthrough population is therefore a
*ceiling* on legitimate firings, not a count of errors: the column measures a
difference from the `fastprep-hd0` baseline, which already trims readthrough
pairs its own way.

## The defect, in case it matters later

Two bugs, both found by measuring rather than reasoning:

1. **The boundary was applied unconditionally.** Fixed by requiring adapter to
   start exactly at the boundary. Verified on synthetic cohorts, flag off vs on,
   both mates: adapter present → trims to the insert; adapter absent → refuses;
   unrelated mates → does not fire; adapter tail shorter than the ordinary
   14-base minimum → refuses (the matcher will not confirm it).
2. **"Contains adapter" was too weak.** A 38-base random tail hits a 14-base
   adapter prefix by chance, so the first guard version fired on ordinary
   sequence. The check now requires the match to begin at the boundary.

## Status

The mode ships **default off** and is safe by construction. Whether it earns its
place is answered by Cell B (a synthetic cohort whose ground truth carries both
G-runs and readthrough, so the scorer is not confounded) — see
`docs/plans/2026-10-19_004500-pe-overlap-with-without.md`.

## Cell B — the value test (synthetic cohort, 2026-10-19)

Cell A could not answer the accuracy question, for two reasons: readthrough is
0.11% of corpus pairs, and the scorer compares trim positions against the
ground truth's **G-run** lengths — so an adapter trim is labelled `overtrim`
wherever the GT has no runs. Cell B is a `fastsim` cohort whose ground truth
carries **both** (489,546 of 800,000 mates with an injected G-run; 112,658 of
400,000 pairs with readthrough; 74,916 readthrough mates also carrying a run),
so the G-run boundary is a real reference.

| tool | exact | overtrim | und-rec | und-irr | dropped |
|---|---|---|---|---|---|
| `fastprep-hd0` | 347,171 | 18,021 | 228,567 | 74,968 | 131,273 |
| `fp-pe` | **347,171** | **18,127** | 228,567 | 74,809 | 131,326 |

**Exact is identical; overtrim rises by 106.** On the readthrough slice that is
+0.094 pp, and +0.000 pp overall — against the pre-registered 0.5 pp threshold,
**do not restore**.

The mechanism is worth recording because it explains the null without hand-waving:

* All 5,328 changed decisions are on readthrough mates, but **4,866 are
  `dropped → dropped`** — the scorer cannot see a dropped mate at all.
* Of the 462 that move a class, the effect is **net-negative** (+211 overtrim,
  −251 dropped): the mode trims 1–2 bp *less* on mates whose adapter the baseline
  had already found. It corrects an over-trim; it does not discover one.
* The addressable set is **1,328 of 109,013 readthrough mates (1.2%)**. Readthrough
  length shrinks as the 3′ run grows, so the same 3′-run mass that makes the
  G-run boundary scorable starves the adapter.

**A structural result the plan got wrong.** The spec asked for inserts inside the
overlap window (`L ≤ I ≤ 2L−12`) *and* readthrough on the same pairs. That set is
empty: readthrough needs `I < L − t`, the window needs `I ≥ L`. Verified on the
cohort — every `rt=1` row has `win=0`. The mode's real operating region is the
below-window readthrough mates.

## Verdict

**`--pe-overlap` does not improve trimming and is not restored as a default.** It
ships opt-in, and it is safe by construction — the boundary applies only when
adapter starts exactly at it, so it can never trim what the ordinary adapter path
could not itself confirm. That property is worth keeping even though the feature
is not: the guard is what turned a silent biology-trimming bug into a measured
null.

Cohort: `bench/work/cellB/` (profile, deterministic build, README with the full
spec and its limits). Plan:
`docs/plans/2026-10-19_004500-pe-overlap-with-without.md`.
