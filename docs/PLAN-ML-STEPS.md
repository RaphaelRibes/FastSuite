# PLAN-ML — step-by-step execution

Companion to `PLAN-ML.md` (rev 3). This is the runnable breakdown: each step
names its artefact, its acceptance test, and its cost. Steps that need a decision
from you are marked **GATE**.

Legend: `[done]` already built and verified · `[next]` next up · `[ ]` pending.

---

## Phase A — the one new component: QC → fitted profile ~1 core-h, local

The whole rev-3 premise rests on this: a profile fitted from a real QC must
reproduce that dataset's QC closely enough that sweeping it means something.

| step | artefact | acceptance test |
|---|---|---|
| **A1** `[done]` write `bench/scripts/qc_to_profile.py`: `fastcheck --json` + `--positions-tsv` → a fastsim profile + its CLI knobs | `qc_to_profile.py` | runs on the last fastcheck report without error |
| **A2** `[next]` **recovery test, local and free of downloads**: profile the *known* `synthetic_novaseq.json` dataset, fit it back, and compare the fitted profile to the one that generated it | test in `bench/tests/` | head CDF recovered within 3 % per level; quality curve within 0.5 Q; GC within 0.005; tail run lengths within one level |
| **A3** `[next]` round-trip test: simulate from the *fitted* profile, run fastcheck on the simulation, compare to the original QC | test | per-position G fraction within 0.02; mean Q within 0.5; GC within 0.01; duplicate rate within 2× |

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
| **A4** `[ ]` ensemble for the non-identifiable fields: 2–3 variants per dataset (insert shape, adapter-vs-short-insert split), bounded and documented | `--variants N` flag | variants differ only in the non-identifiable fields; all load in fastsim |
| **A5** `[ ]` unit tests for the mapping (G-run extraction from the position matrix, quality curve, GC, length) | tests | green |

**Why A2 before anything else:** I have a profile whose parameters I set by hand,
so I can check the fitter *recovers* them. That is a free, exact test of the new
component — no downloads, no cluster, seconds to run.

---

## Phase B — corpus and staging ~3 h transfer, ~2 core-h

| step | artefact | acceptance test |
|---|---|---|
| **B1** `[ ]` `bench/scripts/ena_query.py`: build `manifests/real.tsv` from the ENA portal API, stratified by domain (generalised, metagenomics, human gut, virome, isolate, RNA-seq), with `fastq_ftp`, `fastq_md5`, `fastq_bytes`, `read_count` | manifest | 100 datasets per domain, every row has an FTP path and an md5 |
| **B2** `[ ]` `bench/scripts/subsample.sh`: deterministic every-Nth-record subsample to 100 k pairs, md5 check against ENA, sha256 of the unit, delete the raw | 100 k-pair units | read counts equal on both mates; md5 matches; **no raw file left on scratch** |
| **B3** `[ ]` staging: `scp` the harness, copy the bench env on Io, build fastsim/fastprep/fastcheck there | staged tree | `run_pilot.sh` prints `PROBE PASS` on Io |
| **B4** `[ ]` `ena_query.py --stats`: per-domain count, GC and read-length spread | report | no domain below its floor, no accidental duplicate accessions |

**GATE 1 — the pilot's subject.** Which 30 datasets? My proposal: 10 each from
generalised, human gut, virome — the three most different artifact regimes.

---

## Phase C — the pilot, ~2.5 core-h on Io + downloads ~GATE

| step | what runs | acceptance test |
|---|---|---|
| **C1** `[ ]` full chain on 30 datasets, 60 configs × 3 variants: `fastcheck` → `qc_to_profile` → `fastsim` → 180 evaluations/probe → winner | `results/pilot.tsv` | ≥ 95 % of evaluations scored; every dataset produces a winner |
| **C2** `[ ]` closed loop on real data: does the simulated QC match the real QC? | comparison table | the same thresholds as A3 — this is the go/no-go |
| **C3** `[ ]` run the winner and fastp (defaults + best-effort) on the **real** reads, fastcheck both outputs | `fastp-comparison.tsv` | all three runs complete; quality/speed/memory columns populated |
| **C4** `[ ]` decision memo from the numbers | `bench/PILOT.md` | states plainly whether to scale, fix, or stop |

**GATE 2 — scale or fix.** If C2 fails, Phase D does not run: the fix is in the
fitter, not in the sweep.

---

## Phase D — the full run, ~153 core-h on Io ~GATE

| step | what runs | cost | acceptance |
|---|---|---|---|
| **D1** `[ ]` sweep: 500 datasets × 3 variants × 200 configs = **300 000 exact evaluations** | 131 core-h | every shard file present; per-dataset ≥ 190/200 configs scored |
| **D2** `[ ]` validation on real reads: winner + fastp defaults + fastp best-effort, each `fastcheck`-ed | 10 core-h | 3 rows per dataset, no refusals unexplained |
| **D3** `[ ]` alignment subset: 50 datasets × (winner, fastp), `bwa-mem2`/`minimap2 -x sr` to the reference panel | 12 core-h | `mapped_frac`, `proper_pair_frac`, insert IQR for both tools |
| **D4** `[ ]` aggregate to `matrix.parquet` (DuckDB over the shards) | minutes | row count matches the shard sum; one git hash per run |

**GATE 3 — submission.** Per your standing rule: partition `cpu-dedicated`,
account `dedicated-cpu@cirad-normal` for cells and `@cirad-long` for the
dispatcher, `mem=256M`, flat 48 h, **throttled at 100 cores**, dispatcher
7 days. Nothing is submitted until you confirm those four values.

---

## Phase E — labels and models ~2 core-h, local

| step | artefact | acceptance |
|---|---|---|
| **E1** `[ ]` `bench/scripts/build_labels.py`: sweep + validation + alignment → the label table (domain, accession, QC features, best params, gain, fastp delta, variant agreement) | `labels.parquet` | one row per dataset; no NaN in the label columns |
| **E2** `[ ]` feature builder + the leak test | `features.py`, test | a test asserting no simulator-only field is reachable from the real-data vector |
| **E3** `[ ]` split **by study accession**, stratified by domain and instrument | `splits.json` | test: no study appears in two splits |
| **E4** `[ ]` per domain: train the ranker, distil the direct `QC → params` head | 6 model files | each trains in minutes and is < 1 MB |
| **E5** `[ ]` evaluate: regret vs the sweep's best, top-5 agreement, and **does it beat fastp on held-out real datasets** | `MODEL-CARD.md` per domain | the fastp comparison on held-out data, including the strata where fastp wins |
| **E6** `[ ]` **domain check**: is the generalised model as good as the per-domain ones? | 1 line in the card | if yes, ship one model and say so |

---

## Phase F — productise ~2 core-h, local

| step | artefact | acceptance |
|---|---|---|
| **F1** `[ ]` export the tree ensemble to a compact format; Rust tree-walking in `FastSuite/tools/fastprep` | `--auto-params <qc.json> --print-args` | prints the same argv the Python path recommends; no new dependency; refuses to run without a model file |
| **F2** `[ ]` `fasthelper recommend runs/` for the batch case | CLI | works on a directory of QC reports |
| **F3** `[ ]` versioned model files + a doc on what each model was trained on and where it must not be trusted | `docs/AUTO-PARAMS.md` | states the training domains and the known-bad cases |

---

## Critical path and the immediate next actions

```
A1→A2→A3  ──►  B1→B2→B3→B4  ──►  C1→C2→C3→C4  ──►  D1→D2→D3→D4  ──►  E1..E6  ──►  F1..F3
   (local,       (manifest +         (GATE 2:           (GATE 3:          (local)      (local)
    minutes)      downloads)          go/no-go)          submit)
```

If you say go, the next three actions are, in order, and all local:

1. **A1** write `qc_to_profile.py` (`fastcheck` report → fastsim profile);
2. **A2** the recovery test against the known `synthetic_novaseq.json` — it either
   recovers the parameters it generated from or the fitter is wrong;
3. **A3** the round-trip test (fit → simulate → re-QC → compare).

Nothing in Phases B–D starts until A2 and A3 pass, because everything downstream
assumes the fitted profile is a faithful stand-in for the dataset it came from.
