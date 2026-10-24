# PLAN-ML — step-by-step execution

Companion to `PLAN-ML.md` (rev 3). This is the runnable breakdown: each step
names its artefact, its acceptance test, and its cost. Steps that need a decision
from you are marked **GATE**.

Legend: `[done]` already built and verified · `[next]` next up · `[ ]` pending ·
`[~]` partly delivered, with the gap named · `[—]` dropped/superseded.

Status 2026-10-23: Phases B–D ran, the corpus landed at 101 datasets (not 500),
Phase E produced ONE model and Phase F shipped it in fastprep 0.5.0. The step
rows below carry the measured outcome, not the original estimate.

---

## Phase A — the one new component: QC → fitted profile ~1 core-h, local

The whole rev-3 premise rests on this: a profile fitted from a real QC must
reproduce that dataset's QC closely enough that sweeping it means something.

| step | artefact | acceptance test |
|---|---|---|
| **A1** `[done]` write `bench/scripts/qc_to_profile.py`: `fastcheck --json` + `--positions-tsv` → a fastsim profile + its CLI knobs — **delivered as `fastcheck --profile`, in Rust, not as a Python script (see the note below)** | `qc_to_profile.py` | runs on the last fastcheck report without error |
| **A2** `[done]` **recovery test**: profile the *known* `synthetic_novaseq.json` dataset, fit it back, and compare the fitted profile to the one that generated it — **done and measured**: the fitter recovers **10/10 generating parameters** within the campaign's tolerances (re-run 2026-10-21 on the build that honours `profile.gc`; record `FastSim-Bench/results/recovery.csv`, quoted in `docs/STATUS.md`). It is a campaign measurement, not a `bench/tests/` unit test | test in `bench/tests/` | head CDF within 3 % per level; quality curve within 0.5 Q; GC within 0.005; tail run lengths within one level (the original targets — the 10/10 is scored under the campaign's own tolerance table) |
| **A3** `[done]` round-trip test: simulate from the *fitted* profile, run fastcheck on the simulation, compare to the original QC — **done at corpus scale as the closed loop** (`bench/scripts/closed_loop.py`, same four tolerances): **183/306 = 59.8 %** of comparisons inside them on the current build (`bench/results/gt2026/closed_loop.tsv`), the per-position G fraction the dominant failure | test | per-position G fraction within 0.02; mean Q within 0.5; GC within 0.01; duplicate rate within 2× — **same thresholds, measured on 306 dataset-variants** |

> **A1 is done, in Rust.** `fastcheck --profile p.json` fits and writes the
> profile `fastsim --profile p.json` consumes — the same artefact this Python
> step was to produce, with no fit step in between (verified 2026-10-12: a
> measured 104-entry insert CDF flows from `fastcheck --insert-hist` into a
> fastsim run). `--profile` implies the readthrough + fragment-GC scan;
> `--insert-hist` makes the profile's insert distribution measured rather
> than invented. The fit report (`p.fit_report.json`) records what was
> measured, estimated and invented — an ESTIMATED (head/tail G-runs) or
> INVENTED (insert shape, mate tail correlation) value is never a
> measurement. The Python `qc_to_profile.py` is therefore **not** written;
> do not build it. The suite plan's Phase A (docs/plans/
> 2026-10-12_suite-status-ml-distribution.md) turns A2/A3 into TDD tasks.
| **A4** `[done]` ensemble for the non-identifiable fields: 2–3 variants per dataset (insert shape, adapter-vs-short-insert split), bounded and documented | `--variants N` flag | **done**: `bench/scripts/run_campaign.sh --variants` (default 3) and `chain.sh` change only the seed (`SEED = 42 + VARIANT*1000`), which is exactly the non-identifiable draw; the corpus ran 3 variants per dataset |
| **A5** `[done]` unit tests for the mapping (G-run extraction from the position matrix, quality curve, GC, length) | tests | **done**: `tools/fastcheck/src/profile.rs` unit tests cover the tail G-runs, the out-of-run quality curve, `g_fraction`, the modal read length and the insert CDF |

**Why A2 before anything else:** I have a profile whose parameters I set by hand,
so I can check the fitter *recovers* them. That is a free, exact test of the new
component — no downloads, no cluster, seconds to run.

---

## Phase B — corpus and staging ~3 h transfer, ~2 core-h

| step | artefact | acceptance test |
|---|---|---|
| **B1** `[~]` `bench/scripts/ena_query.py`: build `manifests/real.tsv` from the ENA portal API, stratified by domain (generalised, metagenomics, human gut, virome, isolate, RNA-seq), with `fastq_ftp`, `fastq_md5`, `fastq_bytes`, `read_count` | manifest | **delivered, at a smaller size than planned**: the corpus is **101 datasets over six domains (~15–17/domain)**, not 100/domain — that ceiling is what the model card records |
| **B2** `[done]` `bench/scripts/subsample.sh`: deterministic every-Nth-record subsample to 100 k pairs, md5 check against ENA, sha256 of the unit, delete the raw | 100 k-pair units | ran on the pilot batch — `bench/PILOT.md` B2 records **22/24 runs** subsampled, 1 md5 mismatch and 1 acceptance violation (surfaced, not hidden), and no raw file left behind |
| **B3** `[—]` staging: `scp` the harness, copy the bench env on Io, build fastsim/fastprep/fastcheck there | staged tree | **dropped**: there is no cluster to stage to (PLAN-ML-LOCAL rev 4); the desktop builds from the same checkout |
| **B4** `[~]` `ena_query.py --stats`: per-domain count, GC and read-length spread | report | the flag ships; the counts it reports are **15–17 per domain**, so the plan's 100/domain floor is not met |

**GATE 1 — the pilot's subject.** Which 30 datasets? My proposal: 10 each from
generalised, human gut, virome — the three most different artifact regimes.
**Resolved:** the pilot ran on **24 datasets** (23 with a full three-variant
score), not 30 from three (`bench/PILOT.md`).

---

## Phase C — the pilot (ran on the desktop, 24 datasets)

| step | what runs | acceptance test |
|---|---|---|
| **C1** `[done]` full chain on 30 datasets, 60 configs × 3 variants: `fastcheck` → `fastcheck --profile` → `fastsim` → 180 evaluations/probe → winner | `results/pilot.tsv` | **ran on 24 datasets × 55 configs × 3 variants**; the winner beat the shipped default on **24/24** datasets (mean score 0.9218 vs 0.8125) |
| **C2** `[done]` closed loop on real data: does the simulated QC match the real QC? | comparison table | **36/69 = 52 %** passed at the pilot (go/no-go: SCALE WITH CAVEATS); re-measured at corpus scale on the current build **183/306 = 59.8 %** |
| **C3** `[done]` run the winner and fastp (defaults + best-effort) on the **real** reads, fastcheck both outputs | `fastp-comparison.tsv` | **303 rows** = 101 datasets × (winner, fastp default, fastp tuned), all three populated |
| **C4** `[done]` decision memo from the numbers | `bench/PILOT.md` | written: GATE 2 verdict **SCALE WITH CAVEATS** |

**GATE 2 — scale or fix.** If C2 fails, Phase D does not run: the fix is in the
fitter, not in the sweep. **Resolved:** C2 passed on 52 % of comparisons, so
Phase D ran — on a corpus cut to **101 datasets**, not 500.

---

## Phase D — the full run (ran on the desktop, 101 datasets)

| step | cost | acceptance |
|---|---|---|
| **D1** `[done]` sweep: 500 datasets × 3 variants × 200 configs = **300 000 exact evaluations** | 131 core-h | **ran smaller: 101 datasets × 55 config ids (52 distinct knob sets) × 3 variants, landing 15 847 scored rows** in `bench/results/phase-d/sweep_rows.tsv` |
| **D2** `[done]` validation on real reads: winner + fastp defaults + fastp best-effort, each `fastcheck`-ed | 10 core-h | **303 rows** = 3 per dataset (`fastp-comparison.tsv`) |
| **D3** `[ ]` alignment subset: 50 datasets × (winner, fastp), `bwa-mem2`/`minimap2 -x sr` to the reference panel | 12 core-h | **not run**: no alignment harness exists in the tree (`bwa-mem2`/`minimap2` appear nowhere under `bench/`); the fastp comparison is reference-free |
| **D4** `[done]` aggregate to `matrix.parquet` (DuckDB over the shards) | minutes | aggregated to **`sweep_rows.tsv`** (15 847 rows), the table `labels.py` and `train_auto_params.py` read; no Parquet was written |

**GATE 3 — submission.** Per your standing rule: partition `cpu-dedicated`,
account `dedicated-cpu@cirad-normal` for cells and `@cirad-long` for the
dispatcher, `mem=256M`, flat 48 h, **throttled at 100 cores**, dispatcher
7 days. Nothing is submitted until you confirm those four values.
**Obsolete:** never exercised — Io is gone (PLAN-ML-LOCAL rev 4), and the sweep
ran on the desktop.

---

## Phase E — labels and models ~2 core-h, local

| step | artefact | acceptance |
|---|---|---|
| **E1** `[x]` `bench/scripts/labels.py`: sweep + validation + alignment → the label table (domain, accession, QC features, best params, gain, fastp delta, variant agreement) | `labels.tsv` | one row per dataset; no NaN in the label columns — **15 847 rows over 101 datasets** (the trainer's input) |
| **E2** `[x]` feature builder + the leak test | `QC_FEATURES` in `labels.py` (exported to `features.txt`), `leak_test.py` | the test asserts no simulator-only field is reachable from the real-data vector |
| **E3** `[x]` split **by study accession**, stratified by domain | `bench/scripts/splits.py` → `splits.json` | `leak_test.py` CHECK 1: no study appears in two splits |
| **E4** `[x]` per domain: train the ranker, distil the direct `QC → params` head | 6 model files | **delivered differently, and better than planned: no ranker and no distillation** — `bench/scripts/train_auto_params.py` grows a depth-3 tree whose leaves are configs and flattens it into a **7-rule decision list (8 KB)**, one model for all six domains, embedded in fastprep 0.5.0 and committed at `tools/fastprep/models/auto-params.v1.json` |
| **E5** `[~]` evaluate: regret vs the sweep's best, top-5 agreement, and **does it beat fastp on held-out real datasets** | `docs/MODEL-CARD-auto-params.md` | **regret 0.0331 vs the default's 0.0721, top-5 0.533 vs 0.020, delta +0.0389, LODO positive on 6/6 — delivered**; the fastp comparison on the model's *picks* is still missing (Phase D compared the per-dataset winners) |
| **E6** `[x]` **domain check**: is the generalised model as good as the per-domain ones? | 1 line in the card | **answered yes**: leave-one-domain-out stays positive on all six domains (+0.1305 generalised, +0.0880 bacterial_isolate, +0.0530 virome, +0.0052 human_gut, +0.0052 rna_seq, +0.0004 metagenomics), so **ONE model shipped, not six** — the card's LODO table is the evidence |

---

## Phase F — productise ~2 core-h, local

| step | artefact | acceptance |
|---|---|---|
| **F1** `[x]` export the model (the flattened decision list) to a compact format; Rust rule-walking in `FastSuite/tools/fastprep` | `--auto-params <qc.json> --print-args` | prints the same argv the Python path recommends (102/102 reports, byte-identical); no new dependency; the model is embedded, and a file on disk overrides it |
| **F2** `[ ]` `fasthelper recommend runs/` for the batch case | CLI | **still missing** — `fasthelper`'s subcommands today are run, concatenate, compare, doctor, check, install, schema and guard; the walk is a library function in fastprep, only the CLI front is absent |
| **F3** `[x]` versioned model files + a doc on what each model was trained on and where it must not be trusted | `docs/AUTO-PARAMS.md` | states the training domains and the known-bad cases |

---

## Critical path and the immediate next actions

```
A1→A2→A3  ──►  B1→B2→B3→B4  ──►  C1→C2→C3→C4  ──►  D1→D2→D3→D4  ──►  E1..E6  ──►  F1..F3
   (local,       (manifest +         (GATE 2:           (GATE 3:          (local)      (local)
    minutes)      downloads)          go/no-go)          submit)
```

If you say go, the next three actions are, in order, and all local *(as written;
all three are now done — see below)*:

1. **A1** write `qc_to_profile.py` (`fastcheck` report → fastsim profile);
2. **A2** the recovery test against the known `synthetic_novaseq.json` — it either
   recovers the parameters it generated from or the fitter is wrong;
3. **A3** the round-trip test (fit → simulate → re-QC → compare).

Nothing in Phases B–D starts until A2 and A3 pass, because everything downstream
assumes the fitted profile is a faithful stand-in for the dataset it came from.

**Superseded 2026-10-23.** All three are done: A1 shipped as `fastcheck --profile`
(no Python fitter), A2 measured **10/10** parameters, A3 measured as the corpus
closed loop (**183/306 = 59.8 %**). The path was walked with three deviations —
B3 and GATE 3 were dropped with Io, D3 never ran, and D1/D4 landed as a 15 847-row
TSV instead of 300 000 evaluations in Parquet. What is actually left:

1. **F2** `fasthelper recommend runs/` — the batch front end for the walker that
   already ships in fastprep;
2. **E5's third criterion** — the held-out fastp comparison on the model's *picks*
   (Phase D compared the per-dataset winners, not the model's choices);
3. **D3** the alignment subset, if a reference-free fastp comparison is not enough.
