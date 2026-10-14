# `fastprep --auto-params`: what it would do, what it cannot do yet

**Status: the flag exists and refuses.** `fastprep --auto-params <qc.json>
--print-args` is implemented as a stub: it validates the path, prints the
refusal below, and exits non-zero. No model is trained, no model file exists,
and nothing in the repo produces the argv `--print-args` would emit. This
document is the honest record of what would have to be true before it does.

---

## What `--auto-params` would do

Given a **fastcheck QC report of the real input** (`fastcheck --json qc.json
R1 R2`), recommend the fastprep configuration for that dataset, then print the
argv a run would take:

```
fastprep --auto-params qc.json --print-args
# -> --min-run 12 --mercy 3 --phred-sat 20 --gate-shape 1.0 --min-len 0.9L \
#    --mean-q 25 --adapter-error-rate 0.2 --min-adapter-overlap 14 ...
```

The input is exactly the fields `fastcheck` measures on the reads **before**
trimming (its `metrics` block: `mean_len`, `gc_fraction`,
`mean_read_mean_q`, `reads_q20_frac`, `polyx5_frac_gt_1`, `head_g_frac`,
`duplication_level`, …). The output is the ten config knobs the sweep
explored. The model in between is the Phase E ranker: one decision list per
domain, distilled from the label table `bench/scripts/labels.py` builds.

**The design constraint (Phase F): Rust tree-walking, no new dependency.**
The trained model is exported to a compact per-domain decision list — a small
JSON file of ordered `(field, threshold) → params` rules — and fastprep walks
it in-tree. No ONNX runtime, no tract, no inference crate: the whole
recommendation is a sequence of comparisons against the QC report's measured
fields. The shipped default is the fallback for every dataset whose QC report
is missing a field or falls outside every rule.

## What it CANNOT do today — the honest part

**The shipped default is already the swept winner on the current corpus.**
The 90-dataset corpus campaign (the `FastSuite-bench` repo, the earlier
campaign) swept the config grid over 90 datasets and found **one
configuration best everywhere** (winner score median 0.969, 90/90 datasets).
That configuration **is fastprep's default** — it is locked by the
`shipped_defaults_are_the_swept_winner` test in `tools/fastprep/src/cli.rs`
(min-run 12, mercy 3, phred-sat 20, gate-shape 1.0, min-len 0.9L, mean-q 25,
adapter-error-rate 0.2, min-adapter-overlap 14, dedup on, interior scan off).

So `--auto-params` on the current corpus has **nothing to switch on**: the
model's recommendation and the default would agree on every dataset the
default was already validated on. Shipping the flag today, wired to the
default, would be **dishonest** — it would imply the default is a model's
output when it is the output of an exhaustive sweep over the same datasets.
The flag therefore refuses until a model trained on a *larger, more diverse*
corpus (Phase D's full run) exists and can be shown to beat the default on
**held-out studies** — the acceptance criterion in `docs/PLAN-ML-STEPS.md`
E5: regret vs the sweep's best, top-5 agreement, and a fastp comparison on
held-out real datasets.

**The pilot does not change that verdict.** The 30-dataset pilot
(`bench/PILOT.md`, GATE 2 verdict "scale with caveats") measured per-dataset
winners, but the pilot's corpus is a subset of the same ENA space the 90-dataset
campaign swept — a ranker trained on it would still be validated against the
same default the sweep already chose. Phase D's 500-dataset corpus is what
makes "a model beats the default on held-out studies" a question with an
answer.

## What Phase D's results would change

The label table (`bench/scripts/labels.py`) turns sweep rows into the training
matrix, and the leak test (`bench/scripts/leak_test.py`) guarantees the
held-out evaluation is honest. What Phase D adds over the pilot:

* **Scale**: ~500 datasets × 200 configs × 3 variants of exact GT-scored
  evaluations — enough rows that the per-domain ranker is trained on studies
  that are not the ones the 90-dataset campaign swept.
* **Diversity**: the full corpus spans more of the QC space (read lengths,
  GC extremes, duplication regimes), so a decision list has rules to
  generalise over — the pilot's 28 scored studies (28/30 complete) are too
  few and too narrow.
* **The honest acceptance test**: on Phase D's held-out studies, does the
  model's pick beat the shipped default? If yes, the flag ships a model file
  and this document records which domains it was trained on and where it must
  not be trusted (E5/F3). If no, `--auto-params` stays a refusal — the
  default is the winner, and that is the correct engineering outcome, not a
  failure.

## The label side of the machinery (already built, tested on the pilot)

Phase E's components exist and run green on the pilot's results — they are
waiting on Phase D's sweep rows, which have the identical schema:

| file | what it does | pilot status |
|---|---|---|
| `bench/scripts/labels.py` | sweep_rows.tsv + pilot.tsv + fastcheck reports → `labels.tsv` (one row per dataset × config; 22 QC-side features + the winner/quantile label) | 4614 rows / 28 datasets, green |
| `bench/scripts/leak_test.py` | (1) no study in two splits, (2) no feature is an outcome column or uncomputable at predict time | green; verified to FAIL on a broken split, a leaked `score` feature, and a ghost feature |
| `bench/scripts/splits.py` | study-level group k-fold by `run_accession`, stratified by domain, seeded (default seed 7), byte-deterministic | 28 studies, k=5, green; same seed → identical output |

The features are exactly the QC fields `fastcheck` measures on the input plus
the config knobs. The label is the sweep's own winner (from pilot.tsv) or the
score quantile — never a re-derivation. No GT-derived field is a feature;
`seconds` is not one either (it is a property of the run, not the input).

## The stub's contract

`fastprep --auto-params <qc.json> --print-args` today:

* validates that `<qc.json>` exists and parses as a fastcheck report;
* prints, on stderr, the refusal: no model trained yet, pointing here;
* exits non-zero (CLI-SPEC §1 rule 4: a correct refusal is non-zero with a
  message saying what to do).

It never emits an argv, never falls back to the default silently, and never
claims to recommend anything. When Phase D's model ships, the stub is replaced
by the tree-walk over the model file — the CLI surface and the refusal message
are the only things that change.
