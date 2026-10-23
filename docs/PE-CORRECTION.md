# PE-overlap consensus correction (`--pe-correct`)

Opt-in correction of residual base errors where both mates read the same
fragment. Paired-end FASTQ only (the flag is refused in single-end mode);
the kernel is `tools/fastprep/src/correct.rs`, wired in `pipeline.rs`.

## The rule

For every position in the **doubly-retained** overlap (both mates kept the
base after trimming) where the mates disagree:

1. Convert both mates' quality to error probabilities, `p = 10^(-Q/10)`.
2. Two-hypothesis posterior over which mate carries the true base:
   `odds = (1-p1)·p2 / ((1-p2)·p1)` (wrong-base-uniform factors cancel).
3. Act only when `|log10 odds| >= GATE` (default **2.0**). A tie
   (`odds == 1`) is never resolved by fiat.
4. Rewrite **only the losing mate**: its base becomes the consensus call and
   its Q becomes the posterior Q of that call, capped at 40 — an
   evidence-based rescore, not an inflation.

fastp 1.4.0's `--correction` uses a rigid rule instead: fire only when one
mate is ≥Q30 and the other ≤Q14. The posterior uses the full quality pair, so
it also fires on (and only on) the cases the rigid rule misses where the
evidence is still decisive.

## Geometry

The insert is inferred from the pair (`infer_insert`) and accepted only when
**exactly one** candidate verifies within a mismatch budget over the shared
window. Ambiguous pairs (a repeated fragment that fits several inserts) and
pairs without a usable overlap are left untouched — censoring beats
fabrication. `--pe-correct-min-overlap` (default 30) bounds the shared
window, matching fastp's own overlap floor so comparisons are fair. Both
overlap regimes are covered: `I ≤ L` (full overlap + readthrough) and
`L < I < 2L` (partial overlap).

Correction runs **after** all trim decisions and never revisits them: the
trim path is byte-identical with the flag on or off (reads_out and bases_out
are identical across flag states in every benchmark arm).

## Measured behaviour (20k-pair truth corpus, matched retention)

`bench/results/consensus-headroom*.tsv` — fragment-only denominators, both
axes reported:

| arm | frag err/kb | corrupted | edits |
|---|---|---|---|
| fastp 1.4.0 default | 6.35 | 0 | 0 |
| fastp 1.4.0 `--correction` | 4.15 | 10 | 12,426 |
| fastprep default | 6.35 | 0 | 0 |
| fastprep `--pe-correct=2.0` | **3.93** | 13 | 13,652 |
| fastprep `--pe-correct=0.5` | **3.61** | 104 | 15,692 |

Gate 2.0 (the default arm) leaves 5.2% fewer residual base errors than
`fastp --correction` at comparable corruption; gate 0.5 is the aggressive arm
(13.0% fewer errors, ~10× the corruption, still 99.3% precision). The margin
holds under 2× quality miscalibration in either direction. Consensus can only
reach error mass inside the overlap: with the 30 bp floor that is ~41% of
injected errors (`consensus-headroom-sweep.tsv`).

## Scope and caveats

* Opt-in; default behaviour is unchanged and byte-identical to 0.4.1.
* Paired-end FASTQ only (the flag is refused in single-end mode; on FASTA
  input the pipeline silently no-ops — nothing reads qualities there). The Q written on edited bases is a posterior
  estimate, not a measurement: at gate 2.0 the mean written Q is 33.9
  (claiming ≈0.04% error) while at least 13 of 13,652 rewritten bases are
  wrong (≈0.1%, ≈Q30) — the written Q is somewhat optimistic. Judge the
  feature by the empirical error columns above, not the Q columns; the same
  caution applies to fastp's corrected-Q output.
* Pairs without a unique insert inference are left alone, so on data with
  flat, uniformly-high quality the flag is a no-op (it correctly refuses to
  act on ties).
* The comparison axis is **in-place correction at matched single-read
  output**; fastp's `-m/--merge` produces a different output object and is
  not compared here.

## Throughput

Insert inference scans candidate inserts and is the flag's cost: on a 50k-pair
corpus at `-j 1`, a run goes from ≈0.17 s (flag off) to ≈2.1 s (flag on).
Candidates fail fast on mismatch, but the scan is O(candidates × window);
a seed-and-extend inference would remove most of it.
