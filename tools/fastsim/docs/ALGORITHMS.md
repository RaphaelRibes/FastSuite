# Algorithms

## CDF sampling (`randutil.rs`)

Histograms arrive as `{len: count}` (the profiler's native format). At load
time they become monotonic cumulative distributions `[(len, cumfrac)]`. A
sample draws `u ~ U(0,1)` and binary-searches the first bin with
`cumfrac >= u` — O(log n) per sample, no per-read allocation.

## Joint R1×R2 tail sampling

The profiler emits a joint histogram capped at 64 bp (`JOINT_CAP`). This is
flattened to a 65×65 cumulative grid (row-major, R1 outer). One uniform draw
+ binary search yields an `(r1, r2)` pair whose marginal R1/R2 distributions
exactly match the joint table — reproducing the measured correlation between
the mates' tail lengths (the old simulator drew each mate independently,
which our 82B-read profile showed to be wrong).

## Read construction (`artefact.rs`)

Window arithmetic, per mate with read length L:

```
window  = L − h − t            (bases for biology + adapter)
bio     = min(window, I)       (biology available from the insert)
adapter = window − bio         (> 0 only when I < window: readthrough)
```

R1 biology = `insert[0..bio]`; R2 biology = reverse complement of
`insert[I−bio..I]` (the fragment's far end). Adapter bases cycle the adapter
sequence; tail Gs and head Gs are appended/prepended. `bio < 1` → that mate
is near-full → entirely artefact (`r1_full` / `r2_full`); the pair-level
`full_flag` is 1 only when both mates are.

## Quality model (`quality.rs`)

Per-read baseline variance is sampled from the measured per-read mean-Q
histogram (`read_meanq`). The offset shifts the whole position curve, which
reproduces both the per-position means and the between-read variance (Q30
reads and Q20 reads coexist). Run collapse subtracts `tail_q_collapse` inside
G-runs; this is what makes R1 tails read at ~Q26 while the rest of the read
is ~Q38 — the two-colour chemistry signature.

Errors: `p(err) = 10^(−q/10)` per base, substitution uniform over the 3 other
bases (N untouched). The quality string is the *true* quality of the emitted
base, so downstream tools see a consistent FASTQ.

## Adapter readthrough

When the insert is shorter than the read window, real NovaSeq reads continue
into the adapter. We reproduce this exactly: `read = [h Gs][insert (I
bases)][adapter (window−I bases)][t Gs]`. The GT insert field lets a scorer
distinguish biology from adapter without any per-base annotation.
