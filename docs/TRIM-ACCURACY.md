# Trim accuracy against simulator ground truth

`bench/scripts/trim_accuracy.py` scores overtrim / undertrim / exact for
three tools over 101 simulated datasets. `bench/scripts/aggregate_trim.py`
reproduces every table below from
`bench/results/phase-d/trim-accuracy.tsv`.

## What is measured

The simulator knows the exact composition of every read it emits. For a read

```
[h Gs][biology][adapter?][t Gs]      length L
```

the correct trim is exactly `h` from the 5' end and `t` from the 3' end. A
read with `h = t = 0` must be left untouched. So for each mate, with
observed cuts `t5` (front) and `t3` (back):

| bucket     | condition                        | meaning                     |
|------------|----------------------------------|-----------------------------|
| `exact`    | `|t5-h| <= tol` and `|t3-t| <= tol` | removed exactly the artefact |
| `overtrim` | `t5 > h` or `t3 > t`             | removed real biology        |
| `undertrim`| `t5 < h` or `t3 < t`             | left artefact behind        |

Scored at **tolerance 0** (strict) and **tolerance 2** (forgiving). Reads a
tool dropped are counted separately in `dropped_or_missing` and are not in
any rate.

Tools, gzip level 6 and 4 threads each:

- **fastprep-winner** — the per-dataset winning config from the Phase D
  sweep, scored from its PXMT metadata, which records the exact trim
  applied per read.
- **fastp-default** — fastp 1.4.0 with no trimming flags.
- **fastp-tuned** — fastp with `--cut_front`/`--cut_tail` (window 5, mean
  quality 20), `--qualified_quality_phred 20`, `--length_required 30`,
  `--low_complexity_filter`.

fastp reports no trim positions, so its cuts are recovered by aligning the
trimmed read back onto the original. A pure trim is a substring of the
original, so the alignment is exact.

## Read identity — the bug this benchmark had

fastsim emits illumina headers of the form

```
@LH00916:1:SIMFLOW1:1:1101:<tile>:<cycle> <mate>:N:0:<tag>
```

where `<tile>` takes only **10** values (0..9) and each repeats **10 000**
times per mate file. No header field is a stable read identifier:
`(tile, cycle, mate)` collides 10 000 times per file.

An earlier revision keyed reads as `tile + 10000 * occurrence` and so
compared nearly every read against a *different* read's ground truth. That
produced a phantom **5.1% "unaligned"** bucket for fastp — about 1M mates
whose trimmed output was not a substring of the original they had been
scored against. The stored fastp outputs from those runs were also found to
be internally inconsistent, so the current scorer re-runs fastp from
scratch and scores only freshly generated output.

Reads are now paired **positionally**: both fastp and fastprep preserve
input order and only drop reads, so the i-th surviving trimmed read belongs
to the i-th surviving original. The GT file is in the same emission order
and is replayed positionally. With correct pairing, `unaligned` is exactly
**0** for every tool and every dataset — every trimmed read really is a
substring of its own original, which is the sanity check the previous
numbers failed.

## Results

101 datasets, 19.5–19.8M mates per tool per tolerance. Percentages are of
`exact + overtrim + undertrim + unaligned` (i.e. excluding dropped reads).

### Tolerance 0 (strict)

| tool             | mates     | exact  | overtrim | undertrim | unaligned | dropped |
|------------------|-----------|--------|----------|-----------|-----------|---------|
| fastprep-winner  | 19 477 749| 66.5%  | 29.5%    | 4.0%      | 0.0%      | 722 251 |
| fastp-default    | 19 810 648| 92.2%  | 0.2%     | 7.6%      | 0.0%      | 389 352 |
| fastp-tuned      | 19 810 632| 84.1%  | 8.6%     | 7.3%      | 0.0%      | 389 368 |

### Tolerance 2 (forgiving)

| tool             | exact  | overtrim | undertrim | unaligned |
|------------------|--------|----------|-----------|-----------|
| fastprep-winner  | 98.4%  | 1.3%     | 0.3%      | 0.0%      |
| fastp-default    | 99.6%  | 0.1%     | 0.3%      | 0.0%      |
| fastp-tuned      | 91.2%  | 8.6%     | 0.2%      | 0.0%      |

### Per-domain exact rate

| domain            | fastprep tol0 / tol2 | fastp-default tol0 / tol2 | fastp-tuned tol0 / tol2 |
|-------------------|----------------------|---------------------------|-------------------------|
| bacterial_isolate | 72.8% / 98.7%        | 96.2% / 99.8%             | 80.8% / 83.5%           |
| generalised       | 73.3% / 98.6%        | 96.5% / 99.8%             | 76.1% / 79.2%           |
| human_gut         | 76.9% / 98.6%        | 99.8% / 99.9%             | 99.8% / 99.9%           |
| metagenomics      | 56.8% / 98.7%        | 83.6% / 99.8%             | 83.6% / 99.8%           |
| rna_seq           | 54.6% / 98.2%        | 85.3% / 99.3%             | 84.9% / 98.9%           |
| virome            | 64.9% / 97.7%        | 92.1% / 99.1%             | 78.4% / 84.7%           |

## Interpretation

**fastprep does not win at strict tolerance.** At tolerance 0 it is the
worst of the three on exactness: 66.5% against fastp-default's 92.2%. Its
29.5% strict overtrim is not wrong trimming — it collapses to 1.3% at
tolerance 2, so those are overwhelmingly off-by-one-or-two boundary
decisions about where a polyG run ends. At tolerance 2 fastprep is the most
accurate tool at 98.4%.

**fastp-default is the most conservative.** 0.2% overtrim at either
tolerance — it essentially never removes real biology — but it leaves 7.6%
of artefact behind at tolerance 0. Its undertrim does fall to 0.3% at
tolerance 2, meaning most of what it leaves behind is again a base or two
at a boundary.

**fastp's tuned mode is the worst of the three, and the only case where
tuning hurt.** Its aggressive `cut_front`/`cut_tail` raised overtrim from
0.2% to 8.6% and cut exactness from 92.2% to 84.1% at tolerance 0, and at
tolerance 2 it is still at 8.6% overtrim against fastp-default's 0.1%. That
is a genuine finding about fastp's tuned defaults on polyG-heavy data: the
quality-window trimming does not merely move the boundary, it eats real
sequence.

**Neither tool dominates.** fastprep trades a little overtrimming for
near-complete artefact removal; fastp-default trades undertrimming for never
overtrimming. Which is preferable depends on whether leaving a few adapter
bases or losing a few biology bases matters more downstream.

Note that fastprep drops far more mates than fastp (722k vs 389k), and that
its winner configs were selected per dataset by the Phase D sweep, whereas
both fastp runs use one global setting.

## Caveats

**This measures agreement with the simulator's artefact model, not wet-lab
truth.** The ground truth is exact *for trimming* — the simulator knows
precisely which bases it injected as polyG and adapter. But it is the same
artefact model the closed loop is currently testing, so the benchmark
cannot independently validate that model. A trimmer that happened to share
the simulator's assumptions would score perfectly here and could still be
wrong on real data.

**fastprep's numbers come from its own metadata; fastp's are inferred.**
fastprep's PXMT records the exact trim applied per read, so its scores are
direct. fastp reports nothing, so its cuts are recovered by substring
alignment, which is exact only because a pure trim is a substring. Any
non-trimming edit fastp made would show up as `unaligned` — and since
`unaligned` is now exactly 0, fastp made no such edits on this data.

**Tolerance is a modelling choice, not a measurement.** The gap between
tolerance 0 and tolerance 2 is where most of the disagreement lives. Any
single number quoted from this benchmark should name its tolerance.

## Reproducing

```sh
python3 bench/scripts/trim_accuracy.py     # writes trim-accuracy.tsv
python3 bench/scripts/aggregate_trim.py    # prints the tables above
```

The scorer honours `FASTP`, `TRIM_WORK` and `TRIM_RESULTS` environment
overrides. It re-runs fastp from scratch on every invocation and writes
trimmed output to `trimwork2/`, never reusing stored output.
