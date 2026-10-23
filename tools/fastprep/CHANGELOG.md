# fastprep changelog

Newest first. Default-behaviour changes are marked **DEFAULT** and always carry
the measurement that justified them.

## 0.5.0 — 2026-10-23

### Added — `--auto-params`: the model recommendation is real

`fastprep --auto-params qc.json [--print-args] [--model PATH]` recommends the
trimming configuration for a library from a `fastcheck --json` report of its
input, and `--print-args` prints the argv (stdout) plus the rule that fired,
its training coverage and the measured acceptance (stderr).

* **What it is:** a decision list walked in-tree — ordered `(field, op, value)`
  rules over fastcheck's input-side metrics, each producing a patch of knobs.
  First match wins; no inference runtime, no new dependency (PLAN-ML-STEPS F1).
  An explicitly-passed flag always wins over the model (`ValueSource` is
  checked, so "typed 12" is distinguishable from "default 12").
* **Where it comes from:** `bench/scripts/train_auto_params.py`, trained on
  Phase D's 101 real ENA libraries (6 domains) and 39 expressible configs of
  the 55-config grid; 13 configs are excluded because they need the bench-only
  `r2_reverse` transform, which is not a fastprep flag.
* **Measured (group 5-fold by dataset, stratified by domain):** mean regret vs
  the per-dataset oracle **0.0331** against the shipped default's **0.0721**;
  delta vs the default **+0.0389**, beating it on **70.9 %** of held-out
  datasets; top-5 agreement **0.533** against the default's **0.020**. Leave-
  one-domain-out stays positive on all six domains (+0.0004 metagenomics …
  +0.130 generalised) — see `docs/MODEL-CARD-auto-params.md`, which also states
  the one fold that loses to the default and the regimes the corpus does not
  cover (150 bp reads, 2-colour chemistry, long 3′ tails).
* **The model is embedded in the binary** (`models/auto-params.v1.json`, 8 KB,
  `include_str!`) and still shipped as a file, so `--auto-params` works with no
  setup and the rules can be read without a build. `--model` or
  `$FASTPREP_MODEL` override the discovery order.
* **Refusals that remain:** a `--print-args` without a report, a report whose
  `metrics` are not numeric, a model with an unknown schema tag, an unknown
  feature name (a typo would otherwise be a rule that never fires), an
  unconditional rule before the last one, and any out-of-range knob value.
* Verified by walking all 102 real QC reports in Rust and in the trainer's
  Python reference: **102/102 byte-identical argv** (`/tmp/parity.py`).

**DEFAULT — `--min-len` is `0.67L`, not `0.9L`.** The floor is also the trim
budget (a read needing more trim than the floor allows is dropped, not trimmed),
so `0.9L` let a read lose only 10 % of its length. Measured on a 150 bp read
carrying a 30 bp 5′ polyG artefact: **dropped** at `0.9L`, trimmed to 117 bp at
`0.67L`; in paired mode the same library gave **0 pairs out of 2000** at
`0.9L` and 500/500 pairs at `0.67L`. `0.67L` keeps 100 bp of a 150 bp read — a
third of the read is trimmable — and is a pipeline-owner decision, NOT a sweep
result: the 90-dataset campaign and Phase D both ran at `0.9L`, so the trim-accuracy
tables in `docs/TRIM-ACCURACY.md` describe the old floor. Re-benchmark before
comparing across the change. `--min-len 0.9L` pins the old behaviour, and every
run now prints the resolved floor in bases ("100 bp of a 150 bp mean read;
50 bp of trim allowed") so the budget is never invisible.

**DEFAULT — the mercy clause no longer charges a run's own bases.** Δ_m bounds
how far the cut walks *past* a run's end (the overshoot window), not how long
the run may be. Charging matches turned the budget into a length cap because a
match costs `μ(1−W)`: at Q2 that is 0.09, so 30 bases exhausted Δ_m = 3 and a
**low-quality (4-colour) polyG tail — the normal case for that chemistry — was
trimmed only 37 bases deep however long it was** (measured: 50, 70 and 100 bp
tails all cut at 37; Q5 at 43, Q10 at 64, Q20 uncapped). A match now costs
nothing until the first mismatch; after it the walk is charged as before, so the
overshoot window is unchanged (still the ~2 flanking bases the mercy walk
absorbs). Measured after: every tail length cuts tail + 2 at Q0…Q37. The corpus
cannot see the tail case (its longest injected 3′ tail is 20 bp and
`undertrim-decomp.tsv` attributes 0.95 % of undertrim to 3′ runs ≥ 6), so it
needs a long-tail arm. It does see the rest: on a 15-dataset subset of the GT
corpus, each at its own swept winner with `--min-len` pinned, the shipped arm
moved **76.92 % → 79.19 % exact** and **21.70 % → 19.39 % overtrim**, with the
`scored` denominator within 0.1 % (D-18 has the table); `--head-dropout 0` moved
−0.05 pp. Desktop, 2026-10-23, binary md5 `7b27b2af…`.

**The report closes at base level.** `bases_trimmed` counted only reads that
were emitted, so a run that removed 1 M bases from reads it then dropped
reported "0 bases trimmed (0.00 %)". Trims on dropped reads count; `bases_out`
(a recommended shared key in REPORT-SCHEMA.md that no tool emitted until now) and
`bases_dropped` are reported, and `bases_in == bases_out + bases_trimmed +
bases_dropped` holds on every run. `reads_out` counts records written, not trim
decisions: a kept mate discarded for want of `--singletons1` no longer inflates
it. `bases_interior_mass` is in the v1 `metrics` block. `by_name` accumulates
when a FASTA reuses an adapter ID.

**Correctness fixes** (each pinned by a test in `tests/safety_2026_10_23.rs`,
`tests/pe_overlap.rs` or the module's own suite):

* `--adapter-scored-retry` set `t5` to the window offset as well as
  `t3 = l − p`, so `t5 + t3 == l` and every read the retry rescued was emitted
  at **length 0**, exit 0 (measured 50/50). `t5` is 0 (a 3′-anchored window
  trims nothing from the 5′ end) and the adapter-choice rule is the general
  matcher's: leftmost window, ties by adapter index.
* `--json PATH` naming `-o PATH` overwrote the FASTQ with the report and exited
  0; `-o x -p x`, `-o x --metadata x` and `-o x` with `-1 x` wrote two streams
  into one `.part` and failed at the second rename, after the first commit. Any
  pair of colliding paths is refused before anything is created.
* `.bz2` / `.xz` outputs fell through to the plain-text arm: an uncompressed
  file under a compressed name, exit 0. Refused now.
* the final flush of every FASTQ/singleton writer was `flush().ok()`: a write
  error discarded the last staged block, the truncated file was renamed into
  place and the run exited 0 (reproduced with `ulimit -f 1`: 1024 of 6440
  bytes). All five sites propagate.
* `infer_insert` computed `l1 + l2 − min_ovl` unchecked: a `--pe-correct-min-overlap`
  above the read sum wrapped to ~1.8e19 and the first pair spun forever.
* `--pe-correct` never corrected a Q0 base (p = 1 → degenerate posterior): the
  most decisive disagreement was the one skipped. Decoded Q is clamped at 1.
* `--pe-correct`'s overlap floor applies to the doubly-retained window, not the
  geometric one.
* `--mean-q` was measured on the post-adapter insert, so the same read was kept
  with `--adapters` and dropped with `--no-adapter`; it is the whole-read mean.
* `--target` (anchored) lost its modal trim: the tier-1 prefilter hard-codes
  `b[i] == b[i+1]` and rejected windows the anchored scan accepts (an 11 bp
  terminal G-run trimmed 0 bases). No prefilter in non-homopolymer rules.
* `--no-5 --interior` disabled the interior scan (interior mass 12 → 0): the
  walk now runs whenever any scan is enabled, and an interior-only pass does not
  become the 5′ cut.
* interior runs were measured from an anchor one base past the pair the
  sub-scan re-examines, so every interior run's start, length and mass was off
  by one; an interior-only budget drop was reported as Poly5/Poly3.
* `--head-dropout` matched uppercase `G` only (soft-masked reads lost the trim).
* a trim that consumed the whole read emitted a zero-length FASTQ record when
  `--min-len 0` made both sides of the floor comparison 0.
* adapter matching is case-insensitive end to end (a lowercase adapter tail was
  seeded, then failed every verification compare: 0 trims).
* `--adapter-error-rate 1.0` is refused; short-tail candidates are sorted
  shortest-first as documented; a FASTA with > 255 adapters is refused.

**Flags and report**: `--no-mean-q`/`--mean-q` and `--no-auto-phred`/`--auto-phred`
are last-wins like `--dedup`; `--no-adapter` conflicts with an explicit adapter
set; `--dedup-max-mem` refuses at admission as its help promises; `-c` is
accepted for `--gz-level`; `--dry-run --json` emits the report with
`"dry_run": true` instead of leaving a stale one; `--metadata` help says one
record per read; the dedup block reports the window actually used.

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
