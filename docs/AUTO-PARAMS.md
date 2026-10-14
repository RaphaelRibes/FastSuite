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

**The shipped default is the swept winner on the old corpus, and no longer on
the real one.** The 90-dataset corpus campaign (the `FastSuite-bench` repo,
the earlier campaign) swept the config grid over 90 datasets and found **one
configuration best everywhere** (winner score median 0.969, 90/90 datasets).
That configuration **is fastprep's default** — it is locked by the
`shipped_defaults_are_the_swept_winner` test in `tools/fastprep/src/cli.rs`
(min-run 12, mercy 3, phred-sat 20, gate-shape 1.0, min-len 0.9L, mean-q 25,
adapter-error-rate 0.2, min-adapter-overlap 14, dedup on, interior scan off).

Phase D's real-library sweep does **not** reproduce that verdict, and this
document records the correction rather than the original framing. Over 101
real ENA datasets in six domains the shipped default (config 0) is the
per-dataset winner on **1 of 101**; a per-dataset winner from the 55-config
grid beats it on **100 of 101** (median winner score 0.9643, median gap vs
default 0.0049 — 66 of 101 datasets sit inside a 0.01 margin, 12 above 0.1).
So the default is a good *central* configuration, not a universal one.

What Phase D settles, and what it does not:

* It **does** give `--auto-params` something to switch on. The winners are not
  all the same configuration: 16 distinct configs win somewhere, led by
  config 6 (33 datasets), 44 (16), 52 (10) and 45 (9). The knobs that vary
  most are `mercy` (76 of the 100 non-default winners differ), `min_run` (45)
  and `min_len` (46); the effective length floor is *relative* in 99 of 101
  winners (0.9L on 55, 0.85L on 22), against fastprep's shipped relative 0.9L
  and `fastp`'s absolute 30 bp floor. `r2_reverse` is on in 32 winners, so
  read orientation is genuinely dataset-dependent.
* It **confirms two defaults**: `adapter-error-rate 0.2` and
  `min-adapter-overlap 14` hold in **all 101** winners, and the interior scan
  stays off in 100 of 101 — the one point where the old campaign's verdict
  survives intact.
* `head-dropout` is **not a swept knob** at all — it is absent from the Phase D
  grid, so this campaign says nothing about it either way. It remains a
  default carried from the earlier campaign, not a validated one.

Shipping the flag wired to the default would still be **dishonest**: on this
corpus the default is not the winner, so a ranker that merely reproduces it
would be recommending a losing configuration on 100 of 101 libraries. The
flag keeps refusing until a model trained on a *larger, more diverse* corpus
exists and beats the default on **held-out studies** — the acceptance
criterion in `docs/PLAN-ML-STEPS.md` E5: regret vs the sweep's best, top-5
agreement, and a fastp comparison on held-out real datasets.

**Phase D's real-library head-to-head is the measured context for that test.**
Each of the 101 datasets was also run through `fastp` at its shipped defaults
and at a tuned best (`--cut_front --cut_tail` window 5 / mean-quality 20,
`--qualified_quality_phred 20 --length_required 30 --low_complexity_filter`),
all tools at 4 threads with gzip level 6 output. The per-dataset winners are
faster than default `fastp` on 100 of 101 libraries (median wall-time ratio
0.500, i.e. 2.00x) and use 33x less peak memory on 101 of 101; tuning `fastp`
moves neither axis (0.444 against the tuned configuration, same memory). The
one axis where they lose is retention — a median −5.1 pp — and that is an
artefact of the floors being compared, not of trimming quality: 83 of the 101
winners carry a relative floor above 0.5L while `fastp`'s default floor is an
absolute 30 bp. Any ranker must therefore be scored on the same
floor-normalised axis, or it will learn to recommend looser floors purely to
win a retention number.

**Phase D's corpus is 101 datasets, not 500 — the sweep as run is narrower
than planned.** The committed Phase D results are 101 real ENA datasets × 55
configs × 3 variants (`bench/results/phase-d/winners.tsv`, 101 rows;
`closed_loop.tsv`, 276 rows over 23 datasets), with the closed loop at a
144-of-276 pass rate and the duplicate-ratio gate the dominant failure (120
of 132). That is enough to establish the correction above — the default is
not the universal winner, and per-domain signal exists — but it is **not**
enough to train and honestly evaluate a per-domain ranker: 101 datasets
across six domains is roughly 17 per domain, and a group k-fold over that
leaves a handful of held-out studies per fold. The 500-dataset corpus
originally scoped for Phase D remains the prerequisite for the acceptance
test, not something Phase D has already delivered.

## What Phase D's results changed

The label table (`bench/scripts/labels.py`) turns sweep rows into the training
matrix, and the leak test (`bench/scripts/leak_test.py`) guarantees the
held-out evaluation is honest. What Phase D delivered, and what it did not:

* **Delivered — the correction**: 101 real datasets × 55 configs, each scored
  against simulator ground truth for the same run. This is what overturned the
  "one configuration best everywhere" framing (1 of 101 for the default) and
  what makes the ranker a question worth asking.
* **Delivered — the head-to-head**: the 101-dataset `fastp` comparison
  (`bench/results/phase-d/fastp-comparison.tsv`, 303 rows) is the E5
  "fastp comparison on held-out real datasets", already measured for the
  *winners*. It is not yet measured for a *model's picks* — that is the part
  that still needs the ranker.
* **Not delivered — training scale**: ~500 datasets × 200 configs was the
  original scope; 101 × 55 × 3 variants is what ran. Six domains at ~17
  datasets each is too thin for a group k-fold to leave a meaningful held-out
  set per fold, so the ranker cannot be evaluated honestly on this corpus
  alone.
* **The honest acceptance test, still open**: on a corpus large enough for it,
  does the model's pick beat the shipped default? If yes, the flag ships a
  model file and this document records which domains it was trained on and
  where it must not be trusted (E5/F3). If no, `--auto-params` stays a refusal
  — and on the current corpus that is the *likely* outcome, since the default
  already loses to the per-dataset winners by only a median 0.0049, a margin
  a ranker must clear on held-out studies to be worth shipping at all.

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
