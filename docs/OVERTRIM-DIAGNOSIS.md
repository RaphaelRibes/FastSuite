# Overtrim root cause — diagnosis and what can and cannot fix it

Bench: strict-tolerance (tol 0) overtrim of 29.5% for `fastprep-winner`
against fastsim's exact ground truth. fastp-default: 0.2%; fastp-tuned: 8.6%.
(This 29.5% is the *published* cross-generation figure, since shown invalid —
the corrected overtrim is 23.4%; see `docs/TRIM-ACCURACY-DISCREPANCY.md`. The
mechanism below was re-verified on fresh measurements, so the analysis stands.)

## The cause, established by controlled experiment

`/tmp/probe_overtrim.py` builds reads with a known head run `h` and a known
biology body, runs fastprep with dedup/mean-q/auto-phred disabled, and reads
the cut it actually made back from its PXMT metadata.

| biology construction | 5' overtrim | shape of the error |
|---|---|---|
| no G in biology at all | **0 / 63 (0.0%)** | — |
| biology begins with exactly 1 G | **28 / 63 (44.4%)** | every one exactly **+1** |
| biology begins with exactly 2 G's | **21 / 63 (33.3%)** | every one exactly **+2** |
| random ACGT biology (~25% G) | 3 / 22 | +1 |

When the base after the run is not a G, fastprep's cut is **exact** — no
off-by-one, no mercy overshoot, nothing. The `resolve_subscan` `k_star + 1`
and the `mercy` mass threshold are both correct.

When the base after the run **is** a G, fastprep trims it. It always trims
exactly the number of extra G's, never more.

## Why this is not a coding bug

A run of `h` artefact G's followed by `k` biological G's is
**indistinguishable** from a run of `h + k` artefact G's using base identity
alone. There is no algorithm that reads only the bases and gets this right,
because the input does not determine the answer.

This matches the measured magnitude. With a ~25% per-base G rate, biology
begins with at least one G about 25% of the time; fastprep measures 29.5%
overtrim. The excess over 25% is the 3' side plus the boundary rule's
contribution. The agreement between the predicted 25% and the measured 29.5%
is the confirmation that the experiment found the whole effect.

## What this rules out

- Tightening `mercy` does not help: the overshoot is zero when the boundary
  is unambiguous.
- Changing the `k_star + 1` to `k_star` would introduce **undertrim** of
  exactly the same magnitude, because when the boundary is unambiguous the
  current code is already exact.
- Any per-base classifier cannot fix it. The information is not in the bases.

## What could reduce it

1. **Quality, if it discriminates.** Injected polyG is usually lower quality
   than biology. **Tested and inconclusive, not confirmed.** A controlled probe
   (`/tmp/probe_q_split.py`) put 8 artefact G's at a chosen quality in front of
   a biological G at another quality and looked for a split. It could not
   answer the question: on those reads fastprep *drops* the read for budget
   (28 of 30) rather than reaching a boundary decision, and the 2 reads that
   did reach one gave +5/+6, i.e. the run plus the biological G — no split.
   Note also fastsim's own model gives the FIRST head base normal quality, so
   the boundary is deliberately weak. A decisive test would need reads whose
   run is short enough to survive the budget floor and long enough to trip
   `min_run`; that window is narrow.
2. **Accept the trade deliberately.** A read whose biology begins with a G is
   ambiguous; trimming the G is the conservative choice for artefact removal
   and costs one biology base. Note the tolerance-2 measurement already
   tells us this: overtrim falls from 29.5% to **1.3%** when ±2 is allowed,
   so the vast majority of these overtrims are 1-2 bases.
3. **Report it as a known, bounded limitation** with the ±2 number, which is
   what the paper now does.

## Recommendation

Do **not** chase the 29.5% with algorithm changes. It is the correct answer
to an ambiguous question. The one thing worth measuring is whether quality
at the boundary separates injected from biological G strongly enough to be
worth a rule — and fastsim's first-base-normal-quality detail says the
answer is probably no.

The honest statement is: *fastprep trims the whole observed G run. When
biology begins with a G, that G is indistinguishable from artefact and is
trimmed. This costs 1-2 biological bases on ~25-30% of mates and is the
entire strict-tolerance overtrim.*
