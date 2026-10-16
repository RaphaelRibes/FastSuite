# Trim accuracy against simulator ground truth

## Reading guide — the conclusion has flipped

The earlier revision of this document concluded that **fastp-default is the
strongest tool** (92.2% exact at tolerance 0 against fastprep's 66.5%). That
conclusion is **reversed**. Re-measured under a single harness, **fastprep's
shipped default is the worst configuration at strict tolerance** — 75.6% exact
with 23.4% overtrim — while **fastprep with the leading-G rule disabled
(`--head-dropout 0`) is the strongest tool** (92.6% exact at tolerance 0,
99.7% at tolerance 2), marginally ahead of fastp-default (92.2% / 99.6%).

This is a **configuration finding, not an algorithmic one**. The loss is
entirely attributable to one shipped default: `--head-dropout 4`, a rule that
trims any leading G-run of length 1..=4 exactly, with no quality gate and no
test of artefact versus biology. Turn the rule off and fastprep matches or
beats fastp on every aggregate axis. The shipped default is still 4 as of this
writing; a default change is pending.

`bench/scripts/overtrim_vs_fastp.py` produces these numbers from
`bench/results/phase-d/overtrim-vs-fastp.tsv` (809 lines: header + 808 rows,
101 datasets, four tools, two tolerances).

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
- **fastprep-shipped** — fastprep with its shipped default
  `--head-dropout 4`: any leading G-run of length 1..=4 is trimmed exactly,
  unconditionally.

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

101 datasets, scored by `bench/scripts/overtrim_vs_fastp.py`; 18.8–19.8 M
mates per tool per tolerance.

### Tolerance 0 (strict)

| tool             | mates      | exact     | overtrim | undertrim | dropped   |
|------------------|------------|-----------|----------|-----------|-----------|
| fastprep-hd-off  | 18 846 688 | **92.6%** | 0.1%     | 7.3%      | 1 353 312 |
| fastp-default    | 19 810 648 | 92.2%     | 0.2%     | 7.6%      | 389 352   |
| fastp-tuned      | 19 810 632 | 84.1%     | 8.6%     | 7.3%      | 389 368   |
| fastprep-shipped | 18 769 306 | 75.6%     | 23.4%    | 1.0%      | 1 430 694 |

### Tolerance 2 (forgiving)

| tool             | exact     | overtrim | undertrim |
|------------------|-----------|----------|-----------|
| fastprep-hd-off  | **99.7%** | 0.0%     | 0.3%      |
| fastp-default    | 99.6%     | 0.1%     | 0.3%      |
| fastp-tuned      | 91.2%     | 8.6%     | 0.2%      |
| fastprep-shipped | 98.7%     | 1.2%     | 0.2%      |

### Per-domain exact rate (tol 0 / tol 2)

| domain            | fastprep-hd-off | fastp-default  | fastp-tuned    | fastprep-shipped |
|-------------------|-----------------|----------------|----------------|------------------|
| bacterial_isolate | 96.3% / 99.8%   | 96.2% / 99.8%  | 80.8% / 83.5%  | 77.7% / 98.9%    |
| generalised       | 96.5% / 99.7%   | 96.5% / 99.8%  | 76.1% / 79.2%  | 78.9% / 98.7%    |
| human_gut         | 99.8% / 99.9%   | 99.8% / 99.9%  | 99.8% / 99.9%  | 76.8% / 98.7%    |
| metagenomics      | 83.7% / 99.8%   | 83.6% / 99.8%  | 83.6% / 99.8%  | 72.9% / 99.1%    |
| rna_seq           | 85.7% / 99.4%   | 85.3% / 99.3%  | 84.9% / 98.9%  | 73.0% / 98.5%    |
| virome            | 93.8% / 99.4%   | 92.1% / 99.1%  | 78.4% / 84.7%  | 74.1% / 98.0%    |

### Harness validation

Under the same harness, fastp's numbers reproduce to the digit from its own
published table: **fastp-default 92.19% -> 92.19%** and **fastp-tuned 84.06%
-> 84.06%**. fastprep's published baseline does **not** reproduce: **66.50%
published vs 75.59% fresh**. That discrepancy is a separate, open
investigation; the current measurement supersedes the published baseline
rather than explaining it.

## Interpretation

**The shipped default loses; the algorithm does not.** fastprep with
`--head-dropout 4` is the worst of the four at tolerance 0 (75.6% exact,
23.4% overtrim). With `--head-dropout 0` it is the best (92.6%, 0.1%
overtrim). The single rule accounts for the entire gap: `--head-dropout 4`
trims any leading G-run of length 1..=4 exactly, with no quality gate and no
test of artefact versus biology. On a library where roughly a quarter of reads
start with a real G, that removes real sequence — 23.4% of kept mates are
overtrim at tolerance 0, and even at tolerance 2 (which forgives off-by-one
and off-by-two boundaries) 1.2% remain, because a run of 3–4 genuine leading
G's is trimmed wholesale. The project's own sweep clamps `head_dropout` to
`>= 1`, so the disabled alternative was inexpressible in the sweep and the
default was never calibrated against it.

**fastprep-hd-off is the strongest tool by this measurement**, but the margin
over fastp-default is within noise at tolerance 2 (99.7% vs 99.6%) and slight
at tolerance 0 (92.6% vs 92.2%). Both are effectively at the boundary-jitter
ceiling; the honest claim is parity, not dominance. fastp-default's remaining
weakness is undertrim (7.6% at tolerance 0, falling to 0.3% at tolerance 2 —
largely boundary jitter), and fastprep-hd-off is barely better on that axis
(7.3% / 0.3%). Neither tool meaningfully removes real biology (0.1–0.2%
overtrim).

**fastp's tuned mode is worse than fastp's default and the only configuration
where tuning hurt.** Its aggressive `--cut_front`/`--cut_tail` raised overtrim
from 0.2% to 8.6% and lowered exactness from 92.2% to 84.1%. The overtrim does
*not* decay with tolerance (8.6% at both 0 and 2), so unlike boundary jitter
these are genuine over-trims of several bases. Note the domain split:
fastp-tuned matches fastp-default exactly on human_gut (99.8%) but is far
worse on bacterial_isolate (80.8% vs 96.2%) and virome (78.4% vs 92.1%) —
polyG-heavy libraries.

**fastprep-shipped's one axis of strength is undertrim, and it is bought.**
1.0% undertrim against fastp-default's 7.6% — but at the cost of 23.4%
overtrim, because trimming leading G-runs unconditionally removes artefact
that a quality-gated rule would leave behind along with real biology. It is
not a favourable trade: at tolerance 0 the configuration destroys 23.4% of
kept mates' real sequence to recover 6.6 points of undertrim.

## Drop asymmetry (caveat)

Dropped reads are never scored, so tools with different drop rates are scored
on different kept-sets. Share of each domain's reads dropped at tolerance 0:

| domain            | fastprep-hd-off | fastp-default |
|-------------------|-----------------|---------------|
| bacterial_isolate | 3.8%            | 5.1%          |
| generalised       | 7.3%            | 4.5%          |
| human_gut         | 2.5%            | 0.0%          |
| metagenomics      | 0.4%            | 0.0%          |
| rna_seq           | 10.8%           | 2.1%          |
| virome            | 15.2%           | 0.0%          |

The bacterial_isolate tie (96.3% vs 96.2%) is meaningful — comparable
kept-sets, and fastp drops slightly *more* there. The largest apparent
fastprep-hd-off win (virome, 93.8% vs 92.1%) is measured on a kept-set roughly
one sixth smaller (15.2% of virome reads dropped against fastp's 0.0%), so it
is not an apples-to-apples comparison.

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
  under the current harness (75.59% for the same shipped configuration). This
  is unresolved; the numbers here come from the single-harness run and should
  be cited in preference to the older table.

## Reproducing

```bash
python3 bench/scripts/overtrim_vs_fastp.py   # writes bench/results/phase-d/overtrim-vs-fastp.tsv
```

The tables above are aggregates of that file (tool x domain x accession x
tolerance). `bench/scripts/aggregate_trim.py` still reproduces the
**superseded** 3-tool table from `bench/results/phase-d/trim-accuracy.tsv`
(606 rows, fastprep-winner baseline) — it is kept for the historical record
and its output must not be cited as current. `bench/scripts/trim_accuracy.py`
is likewise the original single-side harness behind that older table.
