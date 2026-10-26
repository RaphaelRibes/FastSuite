# Trim accuracy against simulator ground truth

## Reading guide — the conclusion has flipped, twice

The earliest revision of this document concluded that **fastp-default is the
strongest tool** (92.2% exact at tolerance 0 against fastprep's 66.5% — the
66.5% being an invalid cross-generation figure, see §4 of
`docs/TRIM-ACCURACY-DISCREPANCY.md`). That conclusion was **reversed** on a
single harness: **fastprep's then-shipped default was the worst configuration at
strict tolerance** (75.5% exact with 22.5% overtrim) while **fastprep with the
leading-G rule disabled (`--head-dropout 0`) is the strongest tool** (89.1%
exact at tolerance 0, 99.1% at tolerance 2), marginally ahead of fastp-default
(88.6% / 98.5%). The reversal stands.

**Everything here was then re-measured again** on the regenerated corpus
(`docs/CLOSED-LOOP-RERUN.md`): fastcheck's polyG fitter fix and fastsim's new
ambiguous-base emission both changed the simulated reads. The **ordering did
not move** — not one of the four tools changed rank, at either tolerance — but
the *level* dropped roughly 3 points and **undertrim roughly doubled**
(fastprep-hd-off 7.3% → 10.3%): no trimmer removes an N-run, so reads now end
up under-trimmed rather than over-trimmed. Pre-regeneration values are quoted
alongside the current ones in the tables below.

This is a **configuration finding, not an algorithmic one**. The loss is
entirely attributable to one shipped default: `--head-dropout 4`, a rule that
trims any leading G-run of length 1..=4 exactly, with no quality gate and no
test of artefact versus biology. Turn the rule off and fastprep matches or
beats fastp on every aggregate axis. That default was **0** between `3d76a31`
and the 0.4.1 release and is **4** again since D-16 (the 7:1 price), so for the
`--head-dropout` axis a build from current `main` runs the `fastprep-hd4`
column — the `fastprep-hd-off` column is `--head-dropout 0`; the floor is a
separate change (next paragraph).

**Every number below was measured at `--min-len 0.9L`.** As of D-18 the shipped
floor is **0.67L** (100 bp retained of a 150 bp read), which raises the trim
budget from 15 bp to 50 bp on that read — three and a third times as much — and
therefore drops fewer reads. The tables are NOT comparable across that change:
the `dropped` column is scaled by a different budget and the exact-rate
denominator moves with it. Re-run
`bench/scripts/overtrim_vs_fastp.py` on the desktop before quoting any of these
numbers against a 0.67L build; the ordering of the tools at a fixed floor is
what this document established, and that ordering is expected to hold.

**Provenance caveat (2026-10-26).** The four-tool table comes from
`bench/results/phase-d/overtrim-vs-fastp.tsv` (817 lines: header + 816 rows,
102 datasets, four tools, two tolerances), but the script this document used to
name as its producer — `bench/scripts/overtrim_vs_fastp.py` — **is not in this
repository and never has been** (`git log --all -- bench/scripts` has no such
path). The fastp columns therefore cannot be regenerated from the committed
tree; only the fastprep rows can, with `bench/scripts/trim_score.py` (committed
2026-10-26 with `bench/results/fix2026/`), which re-derives them from the
committed GT. Treat the fastp numbers as campaign records, not as something the
repo can reproduce today.

## What is measured

The simulator knows the exact composition of every read it emits. For a read

```
[h Gs][biology][adapter?][t Gs]      length L
```

the correct trim is exactly `h` from the 5' end and `t` from the 3' end. A read
with `h = t = 0` must be left untouched. So for each mate, with observed cuts
`t5` (front) and `t3` (back):

| bucket      | condition                          | meaning                  |
|-------------|------------------------------------|--------------------------|
| `exact`     | `|t5-h| <= tol` and `|t3-t| <= tol` | removed exactly the artefact |
| `overtrim`  | `t5 > h` or `t3 > t`               | removed real biology     |
| `undertrim` | `t5 < h` or `t3 < t`               | left artefact behind     |

Scored at **tolerance 0** (strict) and **tolerance 2** (forgiving). Reads a tool
dropped are counted separately in `dropped` and are excluded from every rate.

Tools, gzip level 6 and 4 threads each, all four under one harness, one ground
truth, one read-keying scheme:

- **fastprep-hd-off** — fastprep with `--head-dropout 0`, which disables the
  leading-G hard trim entirely. Cuts come from fastprep's PXMT metadata, which
  records the exact trim applied per read.
- **fastp-default** — fastp 1.4.0 with no trimming flags.
- **fastp-tuned** — fastp with `--cut_front`/`--cut_tail` (window 5, mean
  quality 20), `--qualified_quality_phred 20`, `--length_required 30`,
  `--low_complexity_filter`.
- **fastprep-shipped** — fastprep with `--head-dropout 4`, the shipped value on
  that axis (it was 0 for a while after `3d76a31` and D-16 re-priced it back to
  4; 0.4.1 shipped 4 and 0.5.0 still does): any leading G-run of length 1..=4 is
  trimmed exactly, unconditionally. The `0.9L` floor applies here too, so this
  column is the shipped `--head-dropout` value, not the shipped build.

fastp reports no trim positions, so its cuts are recovered by aligning the
trimmed read back onto the original. A pure trim is a substring of the
original, so the alignment is exact.

## Read identity

fastsim emits headers of the form

```
@LH00916:1:SIMFLOW1:1:1101:<x>:<y> <mate>:N:0:<index>
```

where `x = pair_idx % 10 000` and `y = pair_idx / 10 000`, so the pair id is
recovered **exactly** as `x + 10 000 * y` (header fields 5 and 6). No
occurrence counting, and no assumption that the GT shares emission order with
the FASTQ.

Verified: **100 000 / 100 000** headers match their GT `pair_id` on every
dataset tested, including DRR000978, DRR001780 and DRR002643.

Three earlier revisions keyed reads wrongly — by occurrence count, by the wrong
colon fields, and positionally. Each produced a phantom `unaligned` bucket by
comparing a read against a *different* read's ground truth. In particular the
occurrence-counted key is only valid when the GT and FASTQ agree on order,
which held for some datasets and not others; the resulting unaligned rate was
bimodal (exactly 0% on seven datasets, 27–50% on five), which is the signature
of mispairing rather than trimming behaviour. The positional scheme fails
because fastp drops reads (90 711 of 100 000 survive on DRR000978), so
everything after the first drop pairs with an unrelated read.

`unaligned` is now exactly **0** for every tool and dataset — earned, not
circular: each read is matched to its own original by header-derived pair id,
and the surviving counts match fastp's real drop counts (9 289 on DRR000978).

## Results

102 datasets, scored by `bench/scripts/overtrim_vs_fastp.py`; 18.9–19.8 M
mates per tool per tolerance. Re-measured on the regenerated corpus.

### Tolerance 0 (strict)

| tool             | mates      | exact     | overtrim | undertrim | dropped   |
|------------------|------------|-----------|----------|-----------|-----------|
| fastprep-hd-off  | 18,945,051 | **89.1%** | 0.6%     | 10.3%     | 1,454,949 |
| fastp-default    | 19,834,186 | 88.6%     | 0.8%     | 10.6%     | 565,814   |
| fastp-tuned      | 19,849,862 | 80.4%     | 9.6%     | 10.0%     | 550,138   |
| fastprep-shipped | 18,928,152 | 75.5%     | 22.5%    | 2.0%      | 1,471,848 |

Read the `undertrim` column here only next to the decomposition below (full
argument in §"What the undertrim figure actually contains"): **almost all of it
is irreducible, not a defect.** Of fastprep-hd-off's 10.3%, **0.094 pp is
recoverable** and **10.21 pp is not.**

### Tolerance 0, shipped default re-derived after the 2026-10-26 fix wave

`bench/scripts/trim_score.py` on `bench/results/fix2026/` (the same 102 datasets,
the fixed fastprep at its shipped defaults — `0.67L`, `--min-run 12`,
`--mercy 3`, `--head-dropout 4`, boundary-anchored, dedup on):

| | mates | of scored | of all |
|---|---|---|---|
| exact | 13,913,273 | **77.87%** | 68.20% |
| overtrim | 3,573,808 | 20.00% | 17.52% |
| undertrim | 380,327 | 2.13% | 1.86% |
| dropped | 2,532,592 | — | 12.41% |
| scored | 17,867,408 | | |

That is 75.5 / 22.5 / 2.0 in this section's terms, moved ~2.4 pp from overtrim to
exact: the fix closes the case where a 3′ polyX run whose trim landed at the
trim-budget boundary was emitted **whole** (a 48-base artefact on a 151 bp read
at `0.67L`). The row is measured on the fixed build, not spliced in.

### The same measurement after the second fix wave (2026-10-26, `bench/results/wave2/`)

Same 102 datasets, same scorer, same shipped defaults, on the build that fixes
the **two-sided** polyX case (a 3′ artefact was emitted whole when a 5′ artefact
had already spent the trim budget):

| | mates | of scored | of all | Δ vs the row above |
|---|---|---|---|---|
| exact | 13,914,963 | **77.88%** | 68.20% | +1,690 |
| overtrim | 3,574,241 | 20.00% | 17.52% | +433 |
| undertrim | 378,184 | **2.12%** | 1.85% | **−2,143** |
| dropped | 2,532,612 | — | 12.41% | +20 |
| scored | 17,867,388 | | | −20 |

The movement is entirely in the **undertrim** row, which is the defect: 2,143
mates had kept an artefact at the 3′ end, and 20 of them now carry both ends
past the floor and are dropped instead of emitted with a trimmed-head,
untrimmed-tail read. 1,690 of the 2,143 land on the exact cut, 433 on an
overcut. Nothing else in the table moves: the fix is confined to the
(both-ends-artefact × budget-exhausted) corner, which this corpus hits in
0.012% of scored mates.

The closed loop is **unchanged** by this wave — 183/306, 60/102, median
G-fraction 0.0092, byte-identical to `bench/results/fix2026/` over the scored
columns. That is expected and it is not evidence of a thorough fix: every
scored quantity is produced by the same `fastcheck` on both sides of the
comparison, so a report-side error cancels. The wave's real gains are in the
report, the GT and the output paths, which this loop does not measure.

### Tolerance 0 (strict), undertrim decomposed

`bench/scripts/trim_report.py`, output `bench/results/phase-d/trim-report-20261018.tsv`.
"kept mates" = `scored`; shares are of `scored`. **Recoverable** = a failing
mate whose ground-truth artefact is a clean G-run (no adapter, no N) whose
*observed* run length has `P(artefact | r) > 0.5` in the measured posterior;
**irreducible** = every other failing mate (P ≤ 0.5, or an adapter/N-run the
posterior does not model). The totals reproduce the table above to the mate.

| tool             | kept mates | exact      | overtrim | und-recoverable  | und-irreducible      | recoverable share |
|------------------|------------|------------|----------|------------------|----------------------|-------------------|
| fastprep-hd-off  | 18,945,051 | **89.09%** | 0.61%    | 17,793 (0.094%)  | 1,933,465 (10.21%)   | 0.91% of undertrim |
| fastp-default    | 19,834,186 | 88.61%     | 0.76%    | 16,124 (0.081%)  | 2,091,973 (10.55%)   | 0.76%             |
| fastp-tuned      | 19,849,862 | 80.37%     | 9.59%    | 13,538 (0.068%)  | 1,978,415 (9.97%)    | 0.68%             |
| fastprep-shipped | 18,928,152 | 75.48%     | 22.55%   | 14,570 (0.077%)  | 358,848 (1.90%)      | 3.90%             |

fastprep-shipped's recoverable *share* is the largest (3.90%) only because its
unconditional leading-G rule already eats most short G-runs, leaving a
disproportionately long-run residual — it is the worst tool overall, not the one
with the most to fix. For every tool the recoverable class is under 1% of the
undertrim bucket except fastprep-shipped's 3.9% residual.

### Tolerance 2 (forgiving)

| tool             | exact     | overtrim | undertrim |
|------------------|-----------|----------|-----------|
| fastprep-hd-off  | **99.1%** | 0.1%     | 0.9%      |
| fastp-default    | 98.5%     | 0.7%     | 0.9%      |
| fastp-tuned      | 90.2%     | 9.1%     | 0.7%      |
| fastprep-shipped | 98.5%     | 1.1%     | 0.4%      |

### Per-domain exact rate (tol 0 / tol 2)

| domain            | fastprep-hd-off | fastp-default  | fastp-tuned    | fastprep-shipped |
|-------------------|-----------------|----------------|----------------|------------------|
| bacterial_isolate | 93.4% / 99.2%   | 93.3% / 99.2%  | 77.6% / 82.7%  | 77.2% / 98.4%    |
| generalised       | 94.2% / 98.5%   | 94.2% / 98.6%  | 74.0% / 78.7%  | 77.8% / 97.9%    |
| human_gut         | 98.9% / 99.9%   | 98.8% / 99.8%  | 98.8% / 99.8%  | 76.7% / 98.8%    |
| metagenomics      | 79.4% / 99.2%   | 79.2% / 99.1%  | 79.2% / 99.1%  | 73.9% / 99.1%    |
| rna_seq           | 81.0% / 98.4%   | 80.7% / 98.4%  | 80.4% / 98.1%  | 75.4% / 98.6%    |
| virome            | 87.3% / 99.1%   | 86.0% / 95.8%  | 71.7% / 81.6%  | 71.6% / 97.9%    |

The ordering is identical to the pre-regeneration run in all six domains;
only the levels moved.

### Harness validation

Under the same harness, fastp's numbers reproduce to the digit from its own
published table: **fastp-default 92.19% -> 92.19%** and **fastp-tuned 84.06%
-> 84.06%** — both on the pre-regeneration corpus; on the regenerated corpus the
same two configurations read 88.6% and 80.4% above, the ~3-point level shift.
fastprep's published baseline does **not** reproduce: **66.50%
published vs 75.59% fresh**. That discrepancy was a separate defect in the
earlier harness — root-caused in `50e8373` and repaired in `73c2385` (see
`docs/TRIM-ACCURACY-DISCREPANCY.md`) — so the current measurement supersedes
the published baseline.

## Interpretation

**The then-shipped default loses; the algorithm does not.** fastprep at the
then-shipped default (`--head-dropout 4`) is the worst of the four at tolerance
0 (75.5% exact,
22.5% overtrim). With `--head-dropout 0` it is the best (89.1%, 0.6%
overtrim). The single rule accounts for the entire gap: `--head-dropout 4`
trims any leading G-run of length 1..=4 exactly, with no quality gate and no
test of artefact versus biology. On a library where roughly a quarter of reads
start with a real G, that removes real sequence — 22.5% of kept mates are
overtrim at tolerance 0, and even at tolerance 2 (which forgives off-by-one
and off-by-two boundaries) 1.2% remain, because a run of 3–4 genuine leading
G's is trimmed wholesale. The project's own sweep clamped `head_dropout` to
`>= 1` at the time, so the disabled alternative was inexpressible in the sweep
and the default was never calibrated against it — the knob accepts 0 now
(`c53dbd6`), but `--head-dropout 0` is **not** the shipped default: D-16
restored 4 in the 0.4.1 release and 0.5.0 still ships 4.

**fastprep-hd-off is the strongest tool by this measurement**, but the margin
over fastp-default is within noise at tolerance 2 (99.1% vs 98.5%) and slight
at tolerance 0 (89.1% vs 88.6%). Both are effectively at the boundary-jitter
ceiling; the honest claim is parity, not dominance. fastp-default's remaining
weakness is undertrim (10.6% at tolerance 0, falling to 0.9% at tolerance 2 —
largely boundary jitter), and fastprep-hd-off is barely better on that axis
(10.3% / 0.9%). Neither tool meaningfully removes real biology (0.6–0.8%
overtrim).

**fastp's tuned mode is worse than fastp's default and the only configuration
where tuning hurt.** Its aggressive `--cut_front`/`--cut_tail` raised overtrim
from 0.8% to 9.6% and lowered exactness from 88.6% to 80.4%. The overtrim does
*not* decay with tolerance (9.6% at tolerance 0 against 9.1% at 2), so unlike
boundary jitter
these are genuine over-trims of several bases. Note the domain split:
fastp-tuned matches fastp-default exactly on human_gut (98.8%) but is far
worse on bacterial_isolate (77.6% vs 93.3%) and virome (71.7% vs 86.0%) —
polyG-heavy libraries.

**fastprep-shipped's one axis of strength is undertrim, and it is bought.**
2.0% undertrim against fastp-default's 10.6% — but at the cost of 22.5%
overtrim, because trimming leading G-runs unconditionally removes artefact
that a quality-gated rule would leave behind along with real biology. It is
not a favourable trade: at tolerance 0 the configuration destroys 22.5% of
kept mates' real sequence to recover 6.6 points of undertrim.

## Drop asymmetry (caveat)

Dropped reads are never scored, so tools with different drop rates are scored
on different kept-sets. Share of each domain's reads dropped at tolerance 0:

| domain            | fastprep-hd-off | fastp-default |
|-------------------|-----------------|---------------|
| bacterial_isolate | 5.6%            | 8.7%          |
| generalised       | 7.3%            | 4.5%          |
| human_gut         | 2.5%            | 0.0%          |
| metagenomics      | 0.4%            | 0.0%          |
| rna_seq           | 12.2%           | 2.3%          |
| virome            | 18.2%           | 0.3%          |

The bacterial_isolate tie (93.4% vs 93.3%) is meaningful — comparable
kept-sets, and fastp drops slightly *more* there. The largest apparent
fastprep-hd-off win (virome, 87.3% vs 86.0%) is measured on a kept-set roughly
substantially smaller (18.2% of virome reads dropped against fastp's 0.3%), so it
is not an apples-to-apples comparison.

### Read retention (kept vs total)

`bench/scripts/trim_report.py --mode retention` reads the same TSV and reports
`kept = scored`, `total = scored + dropped`, tolerance 0
(`bench/results/phase-d/trim-retention-20261018.tsv`; the drop-share table above
uses a different denominator, these are shares of *all* mates):

| domain            | fastprep-hd-off kept/total | fastp-default kept/total | kept-set gap |
|-------------------|----------------------------|--------------------------|--------------|
| bacterial_isolate | 94.7% (3,218,614/3,400,000) | 92.0% (3,128,564)        | +2.6 pp      |
| generalised       | 92.1% (3,131,038)          | 93.8% (3,189,934)        | −1.7 pp      |
| human_gut         | 97.2% (3,305,267)          | 100.0% (3,399,744)       | −2.8 pp      |
| metagenomics      | 99.5% (3,384,568)          | 100.0% (3,399,580)       | −0.4 pp      |
| rna_seq           | 89.1% (3,029,509)          | 97.8% (3,325,108)        | −8.7 pp      |
| virome            | 84.6% (2,876,055)          | 99.7% (3,391,256)        | −15.2 pp     |
| **ALL**           | **92.87%** (18,945,051)    | **97.23%** (19,834,186)  | **−4.4 pp**  |

fastprep-hd-off's exact rate is earned on **92.9%** of mates against fastp's
**97.2%**. The asymmetry is concentrated, not uniform. The widest per-dataset
gaps are all virome — **SRR20754397 retains 16.4%** against fastp's 100.0%
(−83.6 pp), **SRR19649201 52.0%** vs 99.1%, **SRR17327626 57.3%** vs 100.0% — so
the virome per-domain win (87.3% vs 86.0%) is the least apples-to-apples number
in the table. The bacterial_isolate domain is the opposite case: fastprep keeps
*more* (+2.6 pp) and still ties, which is the comparable-kept-set comparison the
caveat is meant to isolate.

**Re-attributed 2026-10-18: this gap is the dedup default, not `--min-len`.**
The PXMT status byte shows that on those same datasets **100% of the drops are
`Duplicate`** (SRR20754397: 167,162 of 167,162; DRR015139: 32,138 of 32,138), and
`--min-len 15` leaves the histogram untouched to the drop — a ladder of absolute
floors and fractions moves retention by ±0.15 pp at best and costs 2.7–2.9 pp at
worst. With dedup off the aggregate is **98.47% against fastp's 97.23%**:

| config (tol 0) | kept | kept % | exact | exact/kept | overtrim |
|---|---|---|---|---|---|
| fastprep-hd-off (shipped) | 18,945,051 | 92.87% | 16,877,693 | 89.088% | 116,100 |
| fastprep-hd-off `--no-dedup` | 20,087,297 | **98.47%** | 17,814,876 | 88.687% | 138,037 |
| fastp-default (no dedup) | 19,834,186 | 97.23% | 17,575,859 | 88.614% | 150,230 |

So the suite does not keep fewer reads than fastp; it deduplicates, and the
harness charges it for that because a dropped read leaves the scored set. The
ship default stays dedup-on (decision D-13 in `docs/DECISIONS.md`, which also
records the rule that follows: quote both configurations whenever fastp is in the
table, and never read the exact rate without the retention column beside it).

## The two opt-in flags: measured, and neither pays

Two flags were added on the hypothesis that the recoverable class above could be
taken (`--min-run-3p`, a separate 3′ threshold; `--trim-terminal-n`). Both were
then measured over the whole corpus, 102 datasets / 20.4 M mates, against the
same ground truth (`bench/results/phase-d/trim-report-flags-20261018.tsv`,
`fastprep-flags-analysis-20261018.tsv`; the baseline rows reproduce the table
above to the mate). **Neither converts undertrim into exact — both convert it
into overtrim, and both lose exact mates.**

| config (tolerance 0) | exact | overtrim | und-recoverable | und-irreducible | kept | Δ exact |
|---|---|---|---|---|---|---|
| `--head-dropout 0` (baseline) | 16,877,693 | 116,100 | 17,793 | 1,933,465 | 92.868% | — |
| `+ --min-run-3p 6` | 16,674,168 | 352,879 | 3,352 | 1,912,710 | 92.858% | **−1.074 pp** |
| `+ --trim-terminal-n` | 16,776,348 | 232,594 | 17,793 | 1,917,211 | 92.863% | **−0.535 pp** |
| `+ both` | 16,575,521 | 466,384 | 3,352 | 1,896,784 | 92.853% | **−1.595 pp** |

Trade made explicit: `--min-run-3p 6` buys −35,196 undertrim with +236,779
overtrim — one exact mate lost per ~1.4 "gained". Both flags are off by default
and should stay off. Their `--help` text carries this verdict.

### Why — a census of the premise, not a trimmer comparison

`bench/scripts/premise_census.py` reads the raw reads and the ground truth only
(no trimmer involved) and asks, per terminal homopolymer, whether the ground
truth moves that boundary at all. That is the ceiling for any run-length rule:

| 3′ terminal run | mates | P(GT moves that end) |
|---|---|---|
| **G** r=6 | 12,111 | **0.5361** |
| **G** r=7 | 7,435 | 0.6975 |
| **G** r=8 | 4,409 | 0.7521 |
| A / C / T r=6 | 5,418 / 3,582 / 5,319 | **0.0002 / 0.0006 / 0.0002** |
| N, r=1..>15, either end | 66,583 mate-ends | **0.0000** |

Two conclusions, both fatal to the flags as built:

1. **The posterior is a polyG posterior, but `--min-run` is a generic polyX
   threshold.** Short 3′ A/C/T runs are artefacts once in ~2,000; lowering the
   threshold for the run length alone therefore trims biology. G-runs are the
   only end class where the 6–11 band is majority-artefact, and the whole of it
   is worth +0.094 pp.
2. **`--trim-terminal-n` cannot win on this corpus, structurally.** `fastsim`
   places N *inside* the biology window, so no terminal N-run is ever an
   artefact to remove: 66,583 mate-ends carry one, the ground truth moves the
   boundary for none. The earlier "2.9% of undertrim is N-runs" observation is
   correct, and `trim_report.py` already classed those mates **irreducible** —
   the inference that trimming them recovers mates was the error.

A second, independent effect compounds (1): a lower 3′ threshold shrinks the 3′
rejection window, so the reverse scan engages on reads it previously rejected
and overshoots boundaries that were already correct (on one dataset 1,771 of
3,046 moved mates had been exact). Any future attempt at the 3′ band must keep
the window sized to the global `min_run` and lower only the terminal sub-scan
threshold, and it must be G-specific.

**The undertrim direction is closed**: its recoverable ceiling is +0.094 pp, and
every instrument tried against it so far costs at least 11× that in overtrim.



## Caveats

- This measures agreement with the **simulator's artefact model**, not wet-lab
  truth. It is exact ground truth for trimming, but it is the same model the
  closed loop tests, so a result favourable to fastprep would partly reflect
  fastprep being tuned to fastsim's conventions. When fastp appeared to win,
  that bias could not explain the result, since fastp was not fitted to the
  simulator. The polarity has now flipped: with `--head-dropout 0` fastprep
  matches or edges fastp, so the circularity caveat applies again and cuts
  against the new headline. The mitigating fact is that the *reversal* comes
  from disabling one fastsim-unrelated rule (the leading-G trim was a
  heuristic, not a fitted parameter), and fastp's own numbers reproduce to the
  digit under the same harness.
- fastp's cut positions are recovered by alignment, not reported; the
  `unaligned = 0` result is the evidence that the recovery is sound.
- fastprep's numbers come from its own PXMT metadata (exact per-mate trim
  positions), while fastp's are reconstructed — the two sides of the
  comparison are not equally direct.
- The old published baseline (66.5% for fastprep-winner) does not reproduce
  under the current harness (75.59% pre-regeneration, 75.5% on the regenerated
  corpus, for the same shipped configuration). This was resolved in `73c2385`
  — the stored-PXMT staleness bug, see
  `docs/TRIM-ACCURACY-DISCREPANCY.md`; the numbers here come from the
  single-harness run and should be cited in preference to the older table.

## What the undertrim figure actually contains

**Read this before quoting the undertrim rate.** It is 10.3% for fastprep-hd-off
at tolerance 0 — but almost none of it is a trimmer failure, and the raw number
overstates the deficiency by roughly two orders of magnitude.

`bench/scripts/trim_attribution.py` buckets each failing mate by cause, using the
ground truth's per-mate columns (`r{1,2}_adapter`, `r{1,2}_n`). Over the corpus:

| cause of undertrim | share of 1,951,258 |
|---|---|
| **G-run** (no adapter, no N) | **97.0%** |
| N-run | 2.9% |
| adapter | 0.1% |

| end | share |
|---|---|
| 5′ only | 76.7% |
| 3′ only | 19.9% |
| both | 3.3% |

So the question becomes whether those G-runs are *detectable*. The GT records
`r1_h` / `r1_t` — the run length the simulator **injected** — which a real
trimmer cannot see. `bench/scripts/trim_posterior.py` computes
`P(artefact | observed run length)` from the reads plus that column, which is
the ceiling for any run-length rule. On the regenerated corpus:

```
 5'   r    artefact    biology   P(artefact|r)    3'   P(artefact|r)
      1   1,329,574  3,538,296      0.273              0.071
      2     322,766    765,543      0.297              0.107
      3      66,903    168,272      0.284              0.156
      4      15,636     38,409      0.289              0.236
      5       5,527      9,154      0.376              0.384
      6       1,673      2,345      0.416              0.536  <- crosses
      7       1,037        615      0.628              0.698
      9+       ~900        ~60      0.92-0.98          0.79-0.91

  Bayes-optimal recovery from run length alone:
    5'  0.3% of injected runs        3'  34.0% of injected runs
```

**The 5′ dropout is not detectable, only decidable.** `P(artefact)` sits at
~0.27–0.29 for r=1–4, so trimming a short leading G-run is wrong about **73%** of
the time, and it crosses 0.5 only at r≥7 (1.5% of cases). The artefact and the
biology are the same base with the same quality; the simulator knows which is
which only because it injected one. Keeping those bases is the **correct**
decision, and the loop counts it as error.

Decomposed, the 10.3% is (shares re-derived from
`bench/results/exactness/undertrim-decomp.tsv` on 2026-10-26; this table said
~54% / ~40%, which the committed file contradicts by 20–25 points):

| | share of undertrim | verdict |
|---|---|---|
| 5′ r=1–4 | **79.58%** | irreducible — trimming would be wrong ~73% of the time |
| 3′ r=1–5 | **19.01%** | irreducible — P = 0.07–0.38 |
| **3′ r≥6 that the scan missed** | **~0.9%** (~17,700 mates) | **genuinely recoverable** |

The default `min_run` is 12, so 3′ runs of 6–11 fall below it even though
P > 0.5 from r=6. Recovering all of them would move the tolerance-0 exact rate
by about **0.09 pp**.

**Consequence for the paper:** quote the undertrim rate with this decomposition
attached, or it misrepresents the trimmer by ~100×. The defensible claim is
that the trimmer keeps those bases because they are majority biology, and that
the ground truth can label them artefact only because a simulator injected them.

`bench/scripts/trim_report.py` now produces this split mechanically
(`bench/results/phase-d/trim-report-20261018.tsv`) and it reproduces the
hand-derived figures above. Measured across all four tools, the recoverable
class is **0.68–0.91%** of each tool's undertrim bucket (3.90% for
fastprep-shipped's long-run residual). For fastprep-hd-off it is **17,793 mates
— 0.91% of the 1,951,258 undertrimmed mates — worth +0.094 pp of tolerance-0
exact rate**; the other 1,933,465 are irreducible. The recoverable set is
essentially invariant to tolerance (17,793 at tolerance 0, 17,721 at 2), the
signature of a long missed run rather than boundary jitter.

## Reproducing

```bash
# The fastprep rows, re-derived from the committed GT (one configuration per
# invocation, the binary's defaults unless wrapped):
python3 bench/scripts/trim_score.py \
  --work <worktree with sim.v0.*> \
  --fastprep tools/fastprep/target/release/fastprep \
  --out bench/results/fix2026/trim-accuracy-shipped.tsv
```
`bench/scripts/overtrim_vs_fastp.py` (cited by older revisions of this file)
does not exist in the repository; the four-tool table's fastp columns can only
be taken from `bench/results/phase-d/overtrim-vs-fastp.tsv`.

The tables above are aggregates of that file (tool x domain x accession x
tolerance). Two further scripts back the undertrim decomposition:

```bash
python3 bench/scripts/trim_attribution.py    # bucket each failing mate by cause
python3 bench/scripts/trim_posterior.py --end 5   # P(artefact | observed run length)
python3 bench/scripts/trim_posterior.py --end 3
```

`bench/scripts/trim_report.py` merges both into the decomposed table and the
retention table. It re-runs the same four tools, so its `exact`/`overtrim`/
`dropped` totals reproduce the TSV above to the mate, and it reuses
`trim_posterior.py`'s `leading_g`/`trailing_g` (replicating its tabulation):

```bash
python3 bench/scripts/trim_report.py          # corpus + TSV -> both tables
#   -> <results>/trim-report-<YYYYMMDD>.tsv, <results>/trim-retention-<YYYYMMDD>.tsv
python3 bench/scripts/trim_report.py --mode retention   # TSV only; no corpus needed
```

It never overwrites an existing `.tsv` unless `--force` is passed.

Both re-run fastprep and re-parse its PXMT, so their totals reproduce the TSV
above. `overtrim_vs_fastp.py` **refuses to overwrite** an existing `--out`
unless `--force` is given (it appends, so a second run would mix two corpora in
one file); the other two only read. `bench/scripts/aggregate_trim.py` still reproduces the
**superseded** 3-tool table from `bench/results/phase-d/trim-accuracy.tsv`
(606 rows, fastprep-winner baseline) — it is kept for the historical record
and its output must not be cited as current. `bench/scripts/trim_accuracy.py`
is likewise the original single-side harness behind that older table.
