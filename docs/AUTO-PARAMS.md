# `fastprep --auto-params`: the recommendation, and where it can be trusted

**Status: implemented, trained and evaluated (fastprep 0.5.0).**
`fastprep --auto-params qc.json --print-args` prints the recommended argv and the
rule that produced it; without `--print-args` the recommendation is applied to
the run. The model is a decision list walked in-tree, embedded in the binary and
shipped as a readable file at `tools/fastprep/models/auto-params.v1.json`. Its
measured card is `docs/MODEL-CARD-auto-params.md`.

```
$ fastcheck --json qc.json R1.fq.gz R2.fq.gz
$ fastprep --auto-params qc.json --print-args
--min-run 12 --mercy 1 --phred-sat 20 --gate-shape 1 --min-len 0.9L \
  --mean-q 25 --adapter-error-rate 0.2 --min-adapter-overlap 14
auto-params: model <built-in>
  trained on {"corpus":"phase-d (real ENA libraries, GT-free features)", ...}
  profile   12 of 12 feature(s) present
  rule 4/7: polyx5_frac_gt_1 > 0.25125 && kmer_top_frac <= 0.0433 &&
            read_len > 95 [n=12, regret 0.0001]
  measured  CV delta vs the default +0.0389 (beats it on 71 % of held-out
            datasets), top-5 0.533
```

The argv goes to **stdout** (pipe it, or paste it); the explanation goes to
stderr, so `--print-args` composes with a shell. In apply mode nothing goes to
stdout. An **explicit flag always wins** over the model: `--auto-params qc.json
--min-run 14` keeps 14 (`ValueSource` is checked, so "typed 12" is
distinguishable from "default 12").

## What the model is

A decision list: ordered rules, each a conjunction of `(field, op, value)` tests
over the fields `fastcheck` measures **on the input before trimming**, each
producing a patch of trimming knobs. First matching rule wins; if none matches,
the model's `default` patch applies; a knob no patch sets keeps its shipped
default. Walking it is a handful of float comparisons — no ONNX, no tract, no
inference crate (PLAN-ML-STEPS F1) — and the whole model is 8 KB of JSON a user
can read.

Why a decision list and not a boosted ranker: a recommendation you cannot audit
is not usable on a trimming tool whose failure mode is silently wrong reads.
`--print-args` names the rule, its training coverage and its mean regret, so the
user can see the evidence and disagree with it.

The features are `QC_FEATURES` from `bench/scripts/labels.py` — the same list the
training matrix uses, and nothing measured on the trimmed output. Two are
**derived**, and the derivation lives in `auto_params::Profile` so training and
serving cannot drift: `read_len ← mean_len`, `bases_non_acgt_frac ←
bases_non_acgt / bases`. Without that alias the model's `read_len` rules would
parse, never fire, and fall through to the default without saying so.

## Why the default is not good enough (the Phase D correction)

The 90-dataset corpus campaign (the earlier `FastSuite-bench` campaign) swept the
config grid over 90 datasets and found **one configuration best everywhere**
(winner score median 0.969, 90/90 datasets). That configuration **was**
fastprep's default. `shipped_defaults_are_the_swept_winner` in
`tools/fastprep/src/cli.rs` still pins the knobs it shares, but two moved after
the sweep — D-18 took `--min-len` to **0.67L** and D-16 took `--head-dropout` to
**4** — so config 0 is the default *the sweep measured*, not the one 0.5.0
ships.

Phase D's real-library sweep does **not** reproduce that verdict, and this
document records the correction rather than the original framing. Over 101 real
ENA datasets in six domains the sweep's default (config 0) is the per-dataset
winner on **1 of 101**; a per-dataset winner from the grid beats it on **100 of
101** (median winner score 0.9643, median gap vs default 0.0049 — 66 of 101
datasets sit inside a 0.01 margin, 12 above 0.1). The default is a good *central*
configuration, not a universal one, which is why the model has to beat it on
held-out libraries before shipping.

What Phase D settles: the winners are not all one configuration (16 distinct
configs win somewhere, led by config 6 on 33 datasets, 44 on 16, 52 on 10, 45 on
9); the knobs that vary most are `mercy` (76 of the 100 non-default winners
differ), `min_run` (45) and `min_len` (46); `r2_reverse` is on in 32 winners.
It also **confirms two defaults**: `adapter-error-rate 0.2` and
`min-adapter-overlap 14` hold in **all 101** winners, and the interior scan stays
off in 100 of 101.

Phase D's head-to-head is the measured context: each dataset also ran through
`fastp` at its shipped defaults and at a tuned best, all tools at 4 threads with
gzip level 6. The per-dataset winners are faster than default `fastp` on 100 of
101 libraries (median wall-time ratio 0.500, i.e. 2.00x) and use 33x less peak
memory on 101 of 101; tuning `fastp` moves neither axis. The one axis where they
lose is retention — a median −5.1 pp — and that is an artefact of the floors
being compared (83 of the 101 winners carry a relative floor above 0.5L while
`fastp`'s default floor is an absolute 30 bp), not of trimming quality.

## How it was trained

`bench/scripts/train_auto_params.py` reads the label table (`labels.py` output,
15 847 rows over 101 datasets) and grows a depth-3 tree whose leaves are configs,
then flattens it into the ordered list. Leaves hold ≥ 10 datasets and a split
must gain ≥ 0.0005 in mean score, which is what keeps a 101-dataset corpus from
producing 55 one-dataset rules.

* **Action space:** 39 of the sweep's 52 distinct knob sets. The 13 excluded
  configs set `r2_reverse`, a **bench-side** transform (the sweep revcomp'd R2 to
  simulate the other strand convention); it is not a fastprep flag, so
  recommending one would recommend an effect the tool cannot produce. **23
  datasets' oracle winner needed one of those** — a permanent ceiling on this
  model's regret, not a tuning problem.
* **Score:** the sweep's score (keep precision/recall, trim MAE, over/under-trim
  and drop accuracy combined). A config's score for a dataset is the mean over
  variants, the same rule `aggregate.py` uses to name winners.
* **No GT in the features:** the label side needs ground truth; the feature side
  must be available before running, which `labels.py` enforces and
  `bench/scripts/leak_test.py` checks.

## How it was evaluated

Group 5-fold **by dataset**, stratified by domain, so every score is a library
the model has never seen, and all rows of a dataset stay in one fold.

| metric | model | shipped default |
|---|---|---|
| mean regret vs the per-dataset oracle | **0.0331** | 0.0721 |
| delta vs the default (positive = model better) | **+0.0389** | — |
| datasets where it beats the default | **70.9 %** | — |
| top-5 agreement | **0.533** | 0.020 |

Per-domain top-5 (CV): metagenomics 0.765, bacterial_isolate 0.600, rna_seq
0.588, human_gut 0.529, generalised 0.375, virome 0.294.

Leave-one-domain-out — train on five study types, test on the sixth — stays
positive everywhere: generalised **+0.130**, bacterial_isolate **+0.088**,
virome **+0.053**, human_gut **+0.005**, rna_seq **+0.005**, metagenomics
**+0.0004** (no transfer at all there).

**The failure it does have:** one of the five CV folds loses to the default —
fold 1, regret 0.1180 against the default's mean 0.0721, delta −0.0462. The fold
is dominated by the `generalised` domain, where a leaf recommends `min_len 30` on
libraries the default would have left alone. A mean of +0.0389 with one fold at
−0.046 is a real improvement on average and not a guarantee per library, which is
why `--print-args` exists.

The corpus is 101 datasets, not the 500 originally scoped. Six domains at ~17
datasets each leaves 15–17 held-out studies per fold, and the fold-to-fold spread
above (0.0010 to 0.1180) is the width of the uncertainty. The closed loop that
produced the corpus has been re-measured on the current build: **183/306 =
59.8 %** (62 of the loop's 102 libraries passing every variant, median
per-position G-fraction difference 0.0093 against the 0.02 bound). The
167/306 = 54.6 % this document first carried is retracted — it was inflated by
duplicated v1/v2 artefacts (`docs/CLOSED-LOOP-RERUN.md`) — and on the re-run the
per-position G-fraction gate appears in **121 of 123** failing rows, against
**zero** duplicate-rate failures. That is a limit on the *labels*, and it does
not go away with more modelling.

## Where it must not be trusted

These are properties of the corpus, not of the model file, and none of them is
fixed by training harder on the same data:

* **Read length.** The sweep ran on **101 bp reads**. No corpus dataset is 150
  bp, so `read_len` splits were fitted on a range that excludes common modern
  libraries; two of the seven rules test `read_len` and extrapolate there.
* **Chemistry.** No 2-colour (NovaSeq/NextSeq) library is in the corpus, so the
  high-quality 5′ polyG head — the case `--head-dropout` exists for — is not
  represented at all. `--head-dropout` is not an axis of this model.
* **Long low-quality 3′ tails.** The longest injected 3′ tail in the corpus is
  20 bp; D-18's unbounded low-Q run trimming is exactly the regime the sweep
  cannot score.
* **The floor.** Every rule's `min_len` comes from the Phase D grid, which ran at
  0.9L-era floors. D-18 moved the shipped floor to **0.67L** after the sweep, so
  applying a recommendation **replaces** that floor. Pass `--min-len` explicitly
  (explicit flags win) to keep the 0.67L budget.
* **Domain coverage.** Six study types. The LODO table above is the measurement
  of how well it transfers to a seventh — from none (metagenomics) to a lot
  (generalised) — and it is the number to look at before using this on a study
  type the corpus has never seen.

## Retraining

```
bench/scripts/labels.py --rows bench/results/phase-d/sweep_rows.tsv \
    --pilot bench/results/phase-d/pilot.tsv --qc-dir qc/ \
    --out labels.tsv --features features.txt
bench/scripts/train_auto_params.py --labels labels.tsv \
    --out tools/fastprep/models/auto-params.v1.json \
    --card docs/MODEL-CARD-auto-params.md --leave-one-domain-out
```

The QC directory is `<domain>/<acc>.json`, each the `fastcheck --json` report of
the **real** reads the dataset's profile was fitted from (on the desktop corpus:
`work-gt2026/<domain>/<acc>/v0/real.report.json`). The trainer **exits non-zero
if the model does not beat the default** on held-out datasets; `--force` writes
it anyway and the card has to say so. Re-running on the same inputs is
deterministic (`--seed`), which is how the committed model was verified.

## How the walk is verified

* `the_embedded_model_is_the_file_on_disk` — the bytes compiled in are the bytes
  committed, and every rule is satisfiable and not shadowed by a later rule.
* Rust unit tests cover each operator, first-match order, "a missing field is not
  agreement", `mean_q: null` vs absent, range validation, and the fact that an
  explicitly-passed knob survives `apply_recommendation`.
* **Cross-language parity:** walking all 102 real QC reports in Rust and in a
  Python re-implementation of the same rules gives **102/102 byte-identical
  argv**, so the trainer's view of the model and the binary's cannot drift.

## What is still missing

* `fasthelper recommend runs/` — the batch case (one recommendation per QC report
  in a directory). The walk is a library function; only the CLI front is missing.
* The held-out **fastp** comparison on the model's *picks* (E5's third
  criterion). The Phase D head-to-head measured the per-dataset winners, not the
  model's choices, and the two differ on the 23 datasets above.
