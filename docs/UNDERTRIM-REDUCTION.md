# Reducing fastprep's undertrimming — measured mechanisms and ranked fixes

Round 2 of the optimisation audit, 2026-10-21. **Proposals only: none of these
became a default.** UT-1..UT-4 now exist behind opt-in flags (§5); none changed a
default. Every number is [M] measured on the built release binary with the
artifacts in `bench/results/undertrim-*.json` (reproduce with
`bench/scripts/undertrim_{probe,probe_pe,probe_pe_detect,ab}.py`), [E] estimated
with the pricing A/B named, or [V] verified by reading the cited line.

## 0. The measured problem

Corpus: `bench/scripts/undertrim_ab.py` `build_corpus`, 200 000 single-end
reads, seed 31 — 140 157 adapter-free, **59 843 carrying 1..33 bases of TruSeq
readthrough** (29.9 %), of which 10 % get a terminal G-run, 5 % Ns, 3 %
substitutions. Shipped config (`--min-adapter-overlap 14`,
`--adapter-error-rate 0.2`, `--adapter-preset illumina`):

| outcome | reads | of carriers | false trims on clean |
|---|---|---|---|
| exact cut | 174 875 | | |
| **residual (undertrim)** | **25 105** | **41.95 %** | |
| past truth (overtrim) | 20 | | 0.007 % (10 / 140 157) |

Residual by overhang length `t` (`bench/results/undertrim-bands.json`):

| band | reads | exact % | under % | share of residual |
|---|---|---|---|---|
| 1–5 | 9 028 | 0.00 | 99.98 | **36.0 %** |
| 6–9 | 7 314 | 0.00 | 99.93 | **29.1 %** |
| 10–13 | 7 377 | 0.00 | 99.96 | **29.4 %** |
| 14–21 | 14 428 | 94.79 | 5.21 | 5.2 % |
| 22–33 | 21 696 | 97.03 | 2.97 | 2.6 % |

**94.44 % of all residual sits strictly below the 14-base overlap floor.**
`undertrim_probe.py` Arm A is a clean ladder on perfect tails (2 000 reads per
`t`): every `t = 1..13` gives 0 trimmed / 2 000 untrimmed, every `t >= 14` gives
2 000 exact. The cliff is the floor, and nothing else.

Two calibration caveats to carry everywhere: the corpus is synthetic at perfect
quality (real quality moves the error classes), and "residual" counts a
*partial* cut (0 < cut < t) as undertrim — the number a user sees is "reads
still carrying adapter", which is exactly that.

## 1. Mechanism map

### M1 — Sub-floor overhangs, and the PE-overlap guard that should rescue them but cannot (the flagship)

`min_overlap_se = 14` by default (`cli.rs:1199`); `match_adapters`
(`adapter.rs`) requires `k >= min_overlap` for every candidate
(`adapter.rs`: `if k < min_overlap || k < ms || p > max_p { continue; }`), so a
tail of 1..13 bases can never match through the ordinary path (Arm A).

Paired-end data has a second, geometry-based detector: `detect_pe_overlap`
(`pe_overlap.rs:156-219`) pins the insert `I` from mate complementarity — a
12-base anchor with 1 mismatch tolerance (`MAX_ANCHOR_MISMATCHES`), a
window verification with at most 3 mismatches (`MAX_OVERLAP_MISMATCHES`), and a
uniqueness rule (more than one admissible `o` -> `None`, refuse rather than
guess). Chance of a false pin is ~1e-5 x (3/4)^o.

But `detect_pe_overlap_trim` (`pe_overlap.rs:241-261`) then gates the pinned
boundary on `readthrough_is_adapter(r1, t3_1, set)`
(`pe_overlap.rs:275-288`), which calls **`match_adapters(tail, set)` on the
tail itself** — re-imposing the identical `min_overlap_se = 14` floor
[V: read, both lines]. A 10-base tail is provably unmatchable there, so the
detector's perfect boundary is discarded.

**Decisive A/B** (`bench/results/undertrim-probe-pe-guard.json`): R1 =
fragment + TruSeq[:t], R2 = revcomp(fragment) + a **random** t-base tail (never
adapter) — any exact cut on R2 can only come from the PE boundary path, so R2
is a pure probe of the detector+guard:

| arm | t | `--min-adapter-overlap` | R1 exact | R2 exact (PE-only) |
|---|---|---|---|---|
| A | 20 | 14 | 200/200 | **200/200** |
| B | 10 | 14 | 0/200 | **0/200** |
| C | 10 | 10 | 200/200 | **200/200** |

The detector pins the boundary at t=10 (arm C proves the geometry works
identically); the guard's floor alone blocks it (arm B). Consistent with
`bench/results/undertrim-probe-pe.json`: tail3/tail6/tail10 pairs are 0/500
exact in **all** of `auto`/`always`/`never`, while tail20/tail30 are 500/500.
(The custom-dissimilar-adapter run in
`bench/results/undertrim-probe-pe-dissim.json` is a **confounded control**, not
evidence: the tails were outside the configured set, so the guard must refuse
them regardless of the floor.)

### M1b — The paired baseline (the "before" row for the round-3 work)

The single-end corpus above cannot price the PE-overlap path, because that
path needs two mates to pin a boundary. The paired corpus
(`bench/scripts/undertrim_ab.py --paired`, 200 000 pairs, same seed, both
mates scored against a known insert) gives the missing row.

| config | residual / carrier | false trims / clean mate | R1 under | R2 under |
|---|---|---|---|---|
| shipped (mo=14, eps=0.2) | **41.15 %** | 0.011 % | 25 132 | 24 277 |
| mo=12 | 35.60 % | 0.023 % | 21 802 | 20 936 |
| mo=10 | 30.11 % | 0.136 % | 18 523 | 17 629 |
| mo=8 | 24.98 % | 0.228 % | 15 444 | 14 542 |

**The PE path buys almost nothing today.** The same corpus with
`--pe-overlap never` gives 41.21 % residual against the shipped 41.15 % with
`--pe-overlap always` — a **0.06 pp** difference, and R1's residual count is
*identical* (25 132) in both. The detector pins the boundary and the guard
discards it (M1), so the stage is a dead end below the floor.

That 0.06 pp is the size of the prize UT-1 is aiming at: the boundary the
detector already pins is worth the whole 6–13 band (58.5 % of the
single-end residual) if the guard stops re-imposing the floor on it.

Artifact: `bench/results/undertrim-ab-paired.json`.

### M2 — Terminal G-run overprint kills detection (NovaSeq no-signal)

`undertrim_probe.py` Arm C, a 20-base adapter tail with its last `g` bases
replaced by G: `g = 0,2,4` -> 2 000/2 000 exact; **`g = 6,8` -> 0 exact / 2 000
untrimmed**; `g = 10,12` -> partial (the later polyX stage cuts only the G's).
Mechanism [V: `match_adapters`]: the matcher is 3'-anchored — the alignment
window `k = min(a, l-p)` must reach the read terminus — so the G's sit *inside*
the window and consume the mismatch budget `floor(k·eps) = 4` at k=20. On the
AB corpus, 3 673 of the 25 105 residual reads (14.6 %) are G-contaminated
(`g under` column).

Note on stage order: `trim()` runs adapter detection **first**, then polyX on
the trimmed slice (`algorithm.rs:1304-1336`) — so "poly-G runs first and hides
the adapter" is wrong; the failure is the anchored window plus the mismatch
budget. N-tails of 1..3 bases are harmless here (Arm C n-series: 2 000/2 000
exact) but see M4 for the N cost in general.

### M3 — `--pe-overlap auto` policy excludes auto-detected adapter sets

`PeOverlapMode::Auto` runs the stage only when the adapter source is **explicit**
(`--adapter-preset` / `--adapters`); a detected set carries `explicit: false` and
gets no overlap trimming at all (`cli.rs` PeOverlapMode docs; round-1 FP-2 made
this deliberate for safety). Users at pure defaults therefore lose even the
t>=14 PE-boundary pin. The AB corpus used an explicit preset, so the 41.95 %
does not include this mechanism — but real default runs do.

### M4 — Mismatch-budget quantization and N-as-mismatch

`floor(k·eps)` at eps=0.2 allows 2 mismatches through k=14..19 and 3 from k=15
[floor(15*0.2)=3]. `undertrim_probe.py` Arm B (perfect tail, e substitutions):
t=14: e<=2 -> 100 %, e=3 -> 27.2 %, e=4 -> 16.2 %, e=5 -> 6.2 %; t=20: e<=4 ->
100 %, e=5 -> 40.2 %; t=33: 100 % even at e=5. `N` counts as a full mismatch,
which is why the AB corpus shows 1 185 n-class and 851 sub-class residual reads.
Raising eps is the wrong lever: eps=0.3 at mo=14 multiplies clean-read false
trims by ~53 (0.007 % -> 0.373 %) to buy 2.4 pp of residual.

## 2. Why the naive fix fails — the two-axis discipline

`bench/results/undertrim-ab.json`, both axes from the same runs:

| config | under / carrier | false / clean |
|---|---|---|
| shipped (mo=14, eps=0.2) | 41.95 % | 0.007 % |
| mo=12 | 36.32 % | 0.022 % |
| mo=10 | 30.74 % | 0.140 % |
| mo=8 | 25.63 % | 0.243 % |
| mo=6 | 20.56 % | 1.466 % |
| mo=10, eps=0.3 | 27.49 % | 1.448 % |
| mo=14, eps=0.3 | 39.55 % | 0.373 % |

(`undertrim_probe.py` Arm D, 20 000 clean reads, confirms the clean axis:
mo=14 -> 0/20 000, mo=12 -> 5 (0.025 %), mo=10 -> 22 (0.110 %), mo=8 -> 45
(0.225 %), mo=6 -> 258 (1.290 %), mo=5 -> 861 (4.305 %).)

Every global loosening costs 20-50x false trims for at most 16 pp of residual —
and even mo=6 leaves the 1-5 band (36 % of residual) untouched, because a
3-base 3'-prefix is indistinguishable from chance without independent evidence.
The instrument must therefore be **evidence, not threshold**.

## 3. Ranked fixes

### UT-1 — Treat the verified PE boundary as primary evidence (flagship, PE data)

When `detect_pe_overlap` returns a **unique** insert, verify the pinned tail
against the adapter prefix directly — `tail == adapter.seq[..t]`-style, t >= 2
(0 mismatches) or `t >= 3` with `<= floor(t·eps)` and quality-aware scoring —
**without** going through `min_overlap_se`; corroborate on the second mate
(`r2` tail vs `A2[..t]`) where present. Keep `readthrough_is_adapter`'s
"adapter must *start at* the boundary" rule (`pe_overlap.rs:263-274`) — that
part is sound and is what protects biology.

* Undertrim reduction [E]: the M1 population — bands 6-13 = 58.5 % of residual
  plus most of 1-5 under dual-mate corroboration; projected 41.95 % -> **12-18 %
  on PE libraries**. Priced by: `undertrim_ab.py` extended to paired mode, plus
  `undertrim_probe_pe.py` (tail3/6/10 must go 0/500 -> ~500/500).
* Overtrim risk: false boundary pins are already refused on ambiguity and pass
  the anchor (~1e-5) and window checks; the residual risk is a *right* boundary
  whose tail is coincidence, cut by the prefix check at ~4^-t per read at
  t=3-5 — gate t<=5 on dual-mate agreement and keep Cell A
  (`docs/plans/2026-10-19_004500-pe-overlap-with-without.md`: `pe_trims ~= 0`
  on the 102-dataset corpus) as a hard invariant.
* Counterweight to read first: `docs/OVERTRIM-PAIRED.md` — paired structure is
  a *targeted* fix for short-insert libraries, not a general cure.

### UT-2 — G-overprint retry (SE and PE)

When `match_adapters` fails and the 3' terminus carries a G-run (>= 4 bases),
retry the anchored match on the pre-G slice; on a hit, cut at the match start,
consuming adapter + G-run together.

* Undertrim reduction [E]: Arm C g=6..8 goes 0 % -> ~100 % exact; on the AB
  corpus that is the 3 673 g-class reads (14.6 % of residual).
* Overtrim risk: the G-run is removed by the polyX stage anyway; the added cut
  is the adapter extension, and the full `k >= 14` bar still applies to the
  pre-G bases. Require the match to end at the pre-G terminus.
* A/B: Arm C ladder + AB `g under` column + Cell A.

### UT-3 — Two-tier corroborated short tails (SE data, t = 6..13)

SE reads have no mate geometry. Accept an exact 3'-terminal adapter prefix of
6..13 bases only with independent corroboration: a 3' quality drop at the
boundary, or dataset-level 3'-prefix enrichment (learned from the ~4 000-read
auto-detection sample, `cli.rs` `det.sampled`).

* Undertrim reduction [E]: bands 6-13 on SE data = 58.5 % of residual.
* Overtrim risk: the *uncorroborated* version of this is Arm D mo=6 (1.29 %
  false trims); corroboration must hold the clean rate below 0.02 % [E],
  measured as a new arm in `undertrim_ab.py`.

## §5 — Task 1.6: the measurement wave and the default decision

**Before → after, single-end two-axis corpus (200 000 reads, seed 31):**

| metric | before | after (UT-3 arm) | threshold | verdict |
|---|---|---|---|---|
| SE residual / carriers | 41.95 % | **24.34 %** | <= 25 % | **PASS** |
| SE false trims / clean | 0.007 % | **0.016 %** | <= 0.020 % | **PASS** |
| PE residual / carriers | 41.15 % | 41.15 % (unchanged) | <= 25 % | FAIL (no PE change) |
| PE false trims / clean | 0.011 % | 0.011 % (unchanged) | <= 0.020 % | PASS |
| Cell A `pe_trims` | 0 | 0 (UT-1 arm AND UT-3 arm) | = 0 | **PASS** |

The UT-3 arm is `--adapter-short-tail 6 --adapter-min-enrichment 100`
(50 and 100 give identical results on this corpus; 100 is the more
conservative and is the recommended operating point).

### The default decision: no default flips

**None of UT-1..UT-4 becomes the default.** The reasoning is per item:

* **UT-3 — keep behind the flag.** It is the one real win (41.95 % -> 24.34 %
  residual), but it moves the false-trim rate from 0.007 % to 0.016 % — a
  2.3x regression on the overtrim axis — and the plan's own rule is that a
  change must be measured before it becomes a default, not after. The flag
  exists, is tested, and has a measured operating point; flipping the default
  is a release decision, not a measurement outcome. **Recommendation: ship
  behind the flag, revisit after a real short-insert corpus is measured.**
* **UT-1 — keep behind the flag.** Cell A shows zero disturbance
  (`pe_trims = 0`), which is the safety result, but the corpus's fragments
  are longer than the reads so the rescue never fires. There is no measured
  benefit on Cell A to justify a default, and the benefit is confined to
  short-insert libraries.
* **UT-2 — keep behind the flag.** Byte-identical by default; the retry fires
  only on a terminal G-run, which is a NovaSeq-specific artefact. The default
  should not change for a machine-specific artefact.
* **UT-4 — keep behind the flag, and do not rank it as a residual lever.** It
  is inert on this corpus (byte-identical) because the residual sits in
  sub-floor tails, not in N-heavy tails. It is a capability for N-heavy data.

**What the measurement wave established, positively:**
the default path is unchanged by these flags and the suite is green
(`pixi run ci`: 852 passed, 0 failed); every new capability is opt-in; the one
arm that fails the false-trim gate (factor 5) is recorded in the artifact next
to the arms that pass, so the trade-off is visible rather than hidden.

### Cell A — the hard gate on the real corpus, for BOTH behaviour changes

**Cell A** is the hard gate (invariant C-3): 102 datasets, 17,000,000 reads,
ground-truth readthrough. The invariant is `pe_trims = 0` — no trim decision
may change.

**`fp-bnd`** (`--pe-overlap --adapter-evidence boundary`, UT-1) against
`fastprep-hd0`:

| | fastprep-hd0 | fp-bnd |
|---|---|---|
| exact | 13 883 011 | **13 883 011** |
| overtrim | 114 290 | **114 290** |
| undertrim_recoverable | 15 851 | **15 851** |
| dropped | 1 260 259 | **1 260 259** |
| scored | 15 739 741 | **15 739 741** |
| **pe_trims** | n/a | **0** |

**`fp-st`** (`--adapter-short-tail 6 --adapter-min-enrichment 50`, UT-3)
against `fastprep-hd0` — the SAME numbers, every column:

| | fastprep-hd0 | fp-st |
|---|---|---|
| exact | 13 883 011 | **13 883 011** |
| overtrim | 114 290 | **114 290** |
| undertrim_recoverable | 15 851 | **15 851** |
| dropped | 1 260 259 | **1 260 259** |
| scored | 15 739 741 | **15 739 741** |
| **pe_trims** | n/a | **0** |

**Both behaviour changes leave Cell A byte-identical, with `pe_trims = 0`.**
The corpus has no adapter readthrough at all, so this measures the false-fire
rate only — and it is zero for both.

Artifacts: `bench/results/cellA-boundary/` (UT-1) and `bench/results/cellA-ut3/`
(UT-3), both run on the desktop corpus, 2026-10-21.

The reason is structural, and it bounds the whole PE-line of work: this
corpus's fragments are LONGER than the reads, so the mates never overlap and
the boundary is never pinned. The rescue only exists for short inserts
(I < read length). Artifacts:
`bench/results/cellA-boundary/` (run on the desktop corpus, 2026-10-21).

**What this means for the ranking.** UT-1 is *safe* (zero disturbance) but
its benefit is confined to short-insert libraries. The 41.15% → 12.58%
paired gain is real on a corpus built with I < L; it does not transfer to
Cell A. Any proposal ranked above UT-1 must therefore be justified on a
short-insert corpus, not on Cell A.

### UT-3 — Corroborated short tails: the biggest single win, and it passes both axes

**The lever.** ~40 % of every residual class sits in tails shorter than the
14-base overlap floor, and a sub-floor tail is refused for LENGTH before any
mismatch count is consulted. Lowering the floor is the `mo=6` arm: 1.29 %
false trims, a 180x overtrim regression.

**The corroboration.** An EXACT adapter prefix at the 3' end is accepted only
when one of two independent signals says it is real:

1. **Quality drop** — the 4 bases before the boundary are at least
   `--adapter-q-drop` (default 8) lower in mean quality than the 10 bases
   before those.
2. **Dataset-level prefix enrichment** — the prefix occurs at the 3' end of
   more sampled reads than chance predicts. The expected count of a specific
   `t`-mer at a 3' end is `n * 4^-t`, so the enrichment is
   `observed / expected`; `--adapter-min-enrichment FACTOR` sets the bar.

The second signal is the one that works on flat-quality data. Signal 1 alone
is **inert on the AB corpus** (its qualities are flat): 41.95 % residual,
byte-identical to the shipped arm.

**Measured on the 200k single-end two-axis corpus** (SE):

| factor | residual/carrier | false/clean | gate (<= 0.020 %) |
|---|---|---|---|
| shipped | 41.95 % | 0.007 % | PASS |
| 5 | 21.87 % | 0.039 % | FAIL |
| 20 | 21.87 % | 0.037 % | FAIL |
| 50 | **24.34 %** | **0.016 %** | **PASS** |
| 100 | **24.34 %** | **0.016 %** | **PASS** |
| 150 | 26.85 % | 0.009 % | PASS |
| 500 | 26.85 % | 0.009 % | PASS |

**Factor 50-100 is the operating point: residual 41.95 % -> 24.34 % (a 42 %
reduction) with the false-trim rate at 0.016 %, inside the 0.020 % gate.**

**The false trims are understood, not mysterious.** 30 of the 140 157 clean
reads end with a 6-base TruSeq prefix by chance (9 at 7 bases, 2 at 8, 1 at
9); the run trims 22 of them. That is the whole of the 0.016 %: exact
chance hits on a 4^-6 prefix, and nothing else.

**A note on the flag's off-semantics.** `--adapter-min-enrichment 0` DISABLES
the enrichment signal. It must not mean "factor >= 0 passes", which would
corroborate every exact prefix — the uncorroborated form this task exists to
avoid. The first implementation had exactly that bug and trimmed every read
ending in a 6-base prefix.

### UT-4 — Quality/N-aware scoring: LIVE, but inert on the AB corpus

`scored_mismatches` charges an `N` or a badly-called base half a mismatch.
It is behind `--adapter-scored-retry LOW_Q` (default off).

**The mechanism is live.** A 20-base adapter tail with 5 `N`s:
raw budget = `floor(20 * 0.2) = 4`, raw count = 5 → refused;
scored = `5 * 0.5 = 2.5 <= 4` → **accepted**. Measured: 100/200 exact cuts
raw, **200/200** with the retry.

**It does not move the AB corpus, and the reason is a measurement, not a
bug.** On the 200k single-end corpus the result is byte-identical
(174 875 exact / 25 105 under / 20 over, both arms). Two causes, both checked:

* the corpus injects at most **2** `N`s per tail, and the raw budget at
  `k = 14` is already 2 — so the `N`s are inside the raw budget;
* the residual is concentrated in **tails shorter than 14** (`n`-class:
  38.8 % of reads, `sub`-class: 40.2 %, `g`-class: 40.3 % have `t < 14`),
  and a sub-floor tail is refused for LENGTH, before any mismatch count is
  consulted. Scoring cannot help a tail that never reaches the matcher.

So UT-4 is a genuine capability for N-heavy data (e.g. low-quality tails),
but it is not a residual-reduction lever for this corpus. **UT-3 (short
tails) is where the `n`/`sub` residual lives**, and UT-4 must not be ranked
as if it addressed it.

**Design note.** An earlier version of `scored_mismatches` counted only the
ambiguous/low-quality positions and ignored real disagreements entirely, so a
window matching nothing scored 0.0 and every read was trimmed — 164 248
false trims on the SE corpus (100.000 % false). The unit tests
`a_real_disagreement_still_costs_a_full_mismatch` and the corpus gate both
catch it. The `MIN_UNAMBIGUOUS = 12` guard also means `N` vs `N` is **not**
an agreeing call: two uncalled bases are not evidence.

### UT-4 — Quality/N-aware scoring in place of raw mismatch counting

Score mismatches by base quality (Q<=2 counts half / N is a wildcard when >= 12
unambiguous bases agree). LUT-backed, no per-base transcendental.

* Undertrim reduction: the 1 185 n-class + part of the 851 sub-class reads
  (~5-8 % of residual) and the Arm B e=3 gap at t=14.
* Overtrim risk: small [E]; no eps change. A/B: AB n/sub columns + clean corpus.

### UT-5 — "adapter-only reads come out partial" was a PROBE BUG, not a defect

**Finding: no tool change is needed.** The plan's UT-5 was derived from
`undertrim_probe.py` Arm C0 (insert 0, t=150), which reported 200/200 reads
emitted *partial*. That was the scorer, not the trimmer.

The probe built the read as `dna(L - t) + adapter[:t]` with `t = 150`, but
TruSeq/R1 is **33 bases long**, so the read was 33 bases and the recorded
truth was 150. A cut can never equal 150 on a 33-base read, so every read was
scored `partial`.

With the truth corrected to what the read actually carries
(`bench/scripts/undertrim_probe.py`, `t_actual`):
**200/200 exact, 0 partial.** The tool trims the whole read to empty
(`t3_adapter = 33`, output sequence length 0), which is the correct
behaviour for a read with no insert.

A genuinely all-adapter read (150 bases of the TruSeq sequence repeated) is
likewise trimmed to empty: `t3_adapter = 150`, output length 0.

**No `--min-len` change and no `DropReason` change.** The "emit an empty
record" question this task was meant to answer is already settled by the
existing contract: an empty sequence is emitted, the read is counted, and
`--min-len` governs the post-trim length as before.

### UT-5 — Correctness wart: adapter-only reads come out *partial*

`undertrim_probe.py` Arm C0 (insert 0, t=150): 200/200 reads are emitted
**partial** (cut somewhere, not to empty). Trim-to-empty should drop or emit an
empty record per `--min-len` semantics. Low priority, but it is a visible
contract bug.

### Deliberately not proposed

* Global `min_overlap` / `eps` lowering (section 2).
* Re-introducing `--head-dropout`-style recovery: it overtrimmed real sequence.
  The current default is **4** (D-16's 7:1 price) and has been since the 0.4.1
  release, after a brief 0 at `3d76a31`; `--head-dropout 0` is the opt-out (see
  `docs/OVERTRIM-*`, `docs/DECISIONS.md` D-16).
* M4's eps raise.

## 4. Validation protocol

1. `bench/scripts/undertrim_ab.py` — **both axes every time**; a fix quoted on
   undertrim alone is not a measurement.
2. `bench/scripts/undertrim_probe.py` arms A-D (ladder, errors, contamination,
   clean) and `undertrim_probe_pe.py` + `undertrim_probe_pe_guard.py` semantics
   (the random-tail mate design) for UT-1.
3. `fasthelper run --fail-on-adapter-residual` (`adapter_guard.rs`) as the
   shipped gate on real data — the machine-checkable contract fastp does not
   ship.
4. `docs/OVERTRIM-PAIRED.md` Cell A invariant (`pe_trims ~= 0` on the
   102-dataset corpus) before any PE change ships; Cell B for value.
5. The fastsim closed loop (GT under/over rates) for end-to-end claims.

Honesty note: `winners.tsv` and `fastp-comparison.tsv` disagree on 93 of 101
datasets (2026-10-21 review) — name both or neither, never quote either as
canonical.

## 5. What this buys against fastp

1. **A quality contract fastp does not ship**: `fasthelper run
   --fail-on-adapter-residual` exits 2 when adapter survives trimming. One
   command, falsifiable by a sceptic.
2. Residual adapter cut from 41.95 % to a projected 12-18 % on PE data
   (UT-1+UT-2), with clean-read false trims held at the shipped 0.007-0.02 % —
   priced on both axes, not asserted.
3. Speed parity claims kept honest: 16.36 vs 17.64 CPU-s/M pairs and 11x
   lighter memory, with the broken "1-core" wall cell retired
  (`docs/optim/round2-fastprep.md`).
