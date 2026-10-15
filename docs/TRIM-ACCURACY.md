# Trim accuracy — STATUS: NOT TRUSTWORTHY

**Do not cite the numbers in `bench/results/phase-d/trim-accuracy.tsv` or the
tables below. Both pairing schemes tried so far are provably wrong.**

## Why

The benchmark needs to pair each trimmed read with its own ground-truth row.
Two schemes were tried; both fail, in different ways.

### The read id carries no stable identifier

fastsim emits headers of the form

```
@LH00916:1:SIMFLOW1:1:1101:<tile>:<cycle> <mate>:N:0:<tag>
```

Measured on DRR000978 R1: `<tile>` takes **10 000 distinct values (0..9999)**,
each with multiplicity exactly 10. No header field is unique per read.

### Scheme A — decode `pair_id = tile + 10000 * occurrence` (my original)

Valid only when the GT file and the FASTQ are in identical emission order.
That holds for some datasets and not others. When it does not hold, every
read is paired with a *different* read's ground truth.

Measured with this scheme over 12 bacterial_isolate datasets:

| dataset | unaligned |
|---|---|
| DRR000978 | 32.97% |
| DRR001194 | 40.65% |
| DRR001780 | 50.24% |
| DRR001787 | **0.00%** |
| DRR001885 | **0.00%** |
| DRR002002 | **0.00%** |
| DRR002022 | **0.00%** |
| DRR002497 | 27.52% |
| DRR002643 | 50.32% |
| DRR003693 | **0.00%** |
| DRR013910 | **0.00%** |
| DRR014039 | **0.00%** |
| **overall** | **15.36%** |

The bimodality (exactly 0% on some, 27-50% on others) is the signature of
mispairing, not of trimming behaviour. fastp's output quality does not swing
between perfectly aligned and half garbage by dataset. The 5.1% figure in the
committed results is real but inflated by this artifact.

### Scheme B — positional pairing (subagent revision 0cb765e, reverted)

Pair the i-th surviving trimmed read with the i-th surviving original. This
is wrong because **fastp drops reads**: on DRR000978, 90 711 of 100 000
survive. Everything after the first drop is paired with an unrelated read.

Measured: 28.09% of the first 60k originals pair with a trimmed read that is
a substring of **0 out of 5000** other originals — two completely unrelated
reads. That revision reported `unaligned = 0.0%` everywhere because each
read was scored against its own position in a reordered stream, which is
circular.

That revision also claimed `<tile>` has 10 values; it has 10 000. It was
committed with its own validation crashing (a malformed path
`/mnt/Fastsuite-scratch/...`) and was never verified.

## What would fix it

**Match by read content, not by id or position.** fastp preserves read ids,
so the trimmed file's *headers* are a reliable key. The correct approach:

1. Build the map from **trimmed-file header → sequence**.
2. For each original read, look up by header (not by decoded pair_id).
3. Only then attach the GT row, whose ordering must be established
   independently — and if the GT order cannot be established, the GT must be
   re-emitted keyed by the read id fastsim actually wrote.

Step 3 is the open question: whether `sim.v<N>.gt.tsv` rows are keyed by
emission order at all, or whether they carry a field that ties them to the
read id. **That has not been checked and is the prerequisite.**

## What is still valid

- fastprep's PXMT metadata is unaffected: it carries exact per-mate trim
  positions keyed by read id, and all 200 000 entries resolve on DRR000432.
  fastprep's numbers are trustworthy *if* the GT rows are correctly paired.
- The measurement design (overtrim / undertrim / exact at tol 0 and tol 2
  against exact simulator ground truth) is sound.

## The fastp-tuned finding is probably robust

The tuned-mode regression (overtrim 0.2% -> 8.5%) is large enough that it is
unlikely to be a pairing artifact — mispairing would scatter errors rather
than produce a systematic one-directional bias. Worth re-verifying, not worth
discarding.
