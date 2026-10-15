# Trim accuracy against simulator ground truth

`bench/scripts/trim_accuracy.py` scores overtrim / undertrim / exact for three
tools over 101 simulated datasets. `bench/scripts/aggregate_trim.py`
reproduces every table below from
`bench/results/phase-d/trim-accuracy.tsv` (606 rows).

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
dropped are counted separately in `dropped_or_missing` and are excluded from
every rate.

Tools, gzip level 6 and 4 threads each:

- **fastprep-winner** — the per-dataset winning config from the Phase D sweep,
  scored from its PXMT metadata, which records the exact trim applied per read.
- **fastp-default** — fastp 1.4.0 with no trimming flags.
- **fastp-tuned** — fastp with `--cut_front`/`--cut_tail` (window 5, mean
  quality 20), `--qualified_quality_phred 20`, `--length_required 30`,
  `--low_complexity_filter`.

fastp reports no trim positions, so its cuts are recovered by aligning the
trimmed read back onto the original. A pure trim is a substring of the
original, so the alignment is exact.

## Read identity

fastsim emits headers of the form

```
@LH00916:1:SIMFLOW1:1:1101:<x>:<y> <mate>:N:0:<index>
```

where `x = pair_idx % 10 000` and `y = pair_idx / 10 000`, so the pair id is
recovered **exactly** as `x + 10 000 * y`. No occurrence counting, and no
assumption that the GT shares emission order with the FASTQ.

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

101 datasets; 19.5–19.8 M mates per tool per tolerance.

### Tolerance 0 (strict)

| tool            | mates      | exact  | overtrim | undertrim | dropped |
|-----------------|------------|--------|----------|-----------|---------|
| fastprep-winner | 19 477 749 | 66.5%  | 29.5%    | 4.0%      | 722 251 |
| fastp-default   | 19 810 648 | **92.2%** | 0.2%  | 7.6%      | 389 352 |
| fastp-tuned     | 19 810 632 | 84.1%  | 8.6%     | 7.3%      | 389 368 |

### Tolerance 2 (forgiving)

| tool            | exact | overtrim | undertrim |
|-----------------|-------|----------|-----------|
| fastprep-winner | 98.4% | 1.3%     | 0.3%      |
| fastp-default   | **99.6%** | 0.1% | 0.3%     |
| fastp-tuned     | 91.2% | 8.6%     | 0.2%      |

### Per-domain exact rate (tol 0 / tol 2)

| domain            | fastprep     | fastp-default | fastp-tuned |
|-------------------|--------------|---------------|-------------|
| bacterial_isolate | 72.8% / 98.7%| 96.2% / 99.8% | 80.8% / 83.5% |
| generalised       | 73.3% / 98.6%| 96.5% / 99.8% | 76.1% / 79.2% |
| human_gut         | 76.9% / 98.6%| 99.8% / 99.9% | 99.8% / 99.9% |
| metagenomics      | 56.8% / 98.7%| 83.6% / 99.8% | 83.6% / 99.8% |
| rna_seq           | 54.6% / 98.2%| 85.3% / 99.3% | 84.9% / 98.9% |
| virome            | 64.9% / 97.7%| 92.1% / 99.1% | 78.4% / 84.7% |

## Interpretation

**fastprep does not win at strict tolerance.** At tolerance 0 it is the worst
of the three on exactness: 66.5% against fastp-default's 92.2%. Its 29.5%
strict overtrim is not wrong trimming — it collapses to 1.3% at tolerance 2,
so those are overwhelmingly off-by-one-or-two boundary decisions about where a
polyG run ends. At tolerance 2 fastprep is competitive at 98.4% but still
behind fastp-default's 99.6%.

**fastp-default is the strongest tool by this measurement.** 92.2% exact at
tolerance 0 and 99.6% at tolerance 2, with 0.2% overtrim at either tolerance —
it essentially never removes real biology. Its weakness is undertrim: 7.6% at
tolerance 0, though that too is largely boundary jitter since it falls to 0.3%
at tolerance 2.

**fastp's tuned mode is the worst of the three and the only case where tuning
hurt.** Its aggressive `--cut_front`/`--cut_tail` raised overtrim from 0.2% to
8.6% and lowered exactness from 92.2% to 84.1%. The overtrim does *not* decay
with tolerance (8.6% at both 0 and 2), so unlike fastprep's these are genuine
over-trims of several bases, not boundary jitter. Note the domain split:
fastp-tuned matches fastp-default exactly on human_gut and metagenomics
(99.8%) but is far worse on bacterial_isolate (80.8% vs 96.2%) and virome
(78.4% vs 92.1%) — polyG-heavy libraries.

**fastprep's undertrim is the one axis where it clearly leads.** 4.0% against
fastp-default's 7.6% at tolerance 0, and 0.3% against 0.3% at tolerance 2 —
so the edge at tolerance 2 is gone, but at strict tolerance fastprep removes
artefact more completely.

## Caveats

- This measures agreement with the **simulator's artefact model**, not wet-lab
  truth. It is exact ground truth for trimming, but it is the same model the
  closed loop tests, so a result favourable to fastprep would partly reflect
  fastprep being tuned to fastsim's conventions. fastp's win here is not
  subject to that bias, since fastp was not fitted to the simulator.
- fastp's cut positions are recovered by alignment, not reported; the
  `unaligned = 0` result is the evidence that the recovery is sound.
- fastprep's numbers come from its own PXMT metadata (exact per-mate trim
  positions), while fastp's are reconstructed — the two sides of the
  comparison are not equally direct.

## Reproducing

```bash
python3 bench/scripts/trim_accuracy.py     # writes results/trim-accuracy.tsv
python3 bench/scripts/aggregate_trim.py    # prints the tables above
```

Both default to `/mnt/Fast/fastsuite-scratch/{work,results}`; override with
the `TRIM_WORK` and `TRIM_RESULTS` environment variables.
