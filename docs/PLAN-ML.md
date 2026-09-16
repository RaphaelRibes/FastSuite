# PLAN-ML — from downloaded datasets to a parameter-recommending model

> **For Hermes:** phased plan; each phase has an acceptance criterion and a cost.
> Do not start a phase before the previous one's criterion is met.

**Goal:** given only the QC profile of a short-read dataset, recommend the
`fastprep` parameter set that maximises preprocessing quality — where "quality"
is defined exactly on simulated reads (ground truth) and verified by alignment
on real reads.

**One sentence of architecture:** a *contextual surrogate* — `f(QC features of
the dataset, parameter vector) → predicted quality` — trained on a large sweep
of configurations, then used as an argmax over a candidate set. Not a direct
mapping from QC to parameters: many configurations are near-optimal for the
same dataset, so regressing the argmax is unstable, while scoring a candidate
is well-posed and reusable off-policy (every sweep row is a training row).

**Why this shape:** the sweep is where the compute goes, so it must produce
*reusable* labels. A ranking model consumes every (dataset, config, score) row;
a direct QC→params regressor would throw away all but one row per dataset.

---

## 1. Measured cost model (this machine, release builds, 150 bp, gz→gz)

Not estimates — numbers produced while writing this plan, on
`/tmp/pairent_r{1,2}.fastq.gz` (2 000 000 pairs each) and 200k-pair output.

| operation | input size | time | rate | RSS |
|---|---|---|---|---|
| `fastprep -j1` (trim, gz in → gz out) | 2 M pairs | **15.48 s** | 258 420 reads/s | bounded by chunk |
| `fastprep -j8` | 2 M pairs | **5.45 s** | 734 272 reads/s | " |
| `fastcheck` (score pass, pair, 16 cores) | 2 M pairs | **3.71 s** | 1.08 M reads/s | 15.6 MB |
| `fastsim` (PE, gz + gt) | 200 k pairs | **19.9 s** | 10 050 pairs/s | small |

**Work unit chosen for the sweep: one dataset = 200 000 read pairs (~37 MB gz,
~4.5 MB gt).** Scaling the table linearly (fixed startup costs are ~0.1 s):

| per (dataset, config) | core-seconds |
|---|---|
| `fastprep -j1` on 200k pairs | ≈ 1.55 |
| tier-1 score (PXMT vs GT, streaming) | ≈ 0.10 |
| **tier-1 evaluation total** | **≈ 1.65** |
| tier-3 score (fastcheck on output) instead | ≈ 0.40 |
| tier-2 score (alignment, validation subset only) | 30–60 |

So: **one exact evaluation ≈ 1.7 core-seconds.** Everything below is priced
from that one number.

---

## 2. The score: three tiers, one plug-in interface

The score is the heart of the project and the easiest place to lie to
ourselves. It is a versioned module with an explicit weights file, and every
result row records the score version that produced it.

### Tier 1 — exact, simulated data only (the training signal)

`fastsim` writes `PREFIX.gt.tsv` with one row per pair:
`pair_id  insert  r1_h  r1_t  r2_h  r2_t  full_flag  is_dup`, where `r{n}_h` /
`r{n}_t` are the head/tail bases a **correct** trimmer must remove from each
mate. `fastprep --metadata` writes a PXMT record for every read it trims or
drops. Comparing the two gives an exact, per-read, reference-free verdict:

| counter | definition |
|---|---|
| `keep_clean_ok` | pair untouched by GT (h=t=0), untouched by fastprep |
| `trim_exact` | pair trimmed, and `(h', t') == (h, t)` for both mates |
| `trim_over` | fastprep removed more true bases than GT: `h' > h` or `t' > t` |
| `trim_under` | fastprep left artefact: `h' < h` or `t' < t` |
| `drop_wrong` | pair dropped by fastprep but GT keeps it |
| `drop_ok` | pair dropped by fastprep and GT's within-pair overlap makes the pair untrimmable (`full_flag` / must-drop class) |

```
score_t1 = w_keep·keep_clean_ok
         + w_exact·trim_exact
         − w_over·trim_over        # deleting real sequence is the sin that matters
         − w_under·trim_under      # leaving artefact is the milder one
         − w_drop·drop_wrong       # losing a pair costs a downstream read pair
           (all normalised by pairs_in; weights in score_weights.toml, sum of the
            negative weights fixed at 1.0 of the pair budget so scores are comparable
            across datasets)
```

Default weights (explicit, arguable, versioned — `score_weights.toml`):
`w_keep = 1.0`, `w_exact = 1.0`, `w_over = 0.5`, `w_under = 0.25`,
`w_drop = 0.5`. A near-tie tolerance `eps = 0.005` is applied when reporting
"best config", so noise-level differences are not treated as signal.

**Why this matters:** it is exact (no proxy), cheap (~0.1 s), needs no
aligner, no reference and no network, and it is the only tier where the
*right answer is known*. It is what makes 10⁶ evaluations affordable.

### Tier 2 — alignment-validated, real data (the transfer check)

No ground truth exists for real reads, so the real tier uses the closest thing
to a downstream truth: **alignment**. On a subset only (50 datasets × 20
configs), align both mates to a small reference panel with `bwa-mem2 mem` (or
`minimap2 -x sr`), and record `mapped_frac`, `proper_pair_frac`,
`mapq30_frac` and the insert-size median/IQR. Reference panel (download once,
~300 MB): *E. coli* K-12, *S. cerevisiae*, *A. thaliana*, human chr20 — one
small genome per organism stratum, plus the organism of each dataset when it is
in the panel.

This tier is not for training; it answers one question: **does the ranking the
model learned on simulated data transfer to real data?** Reported as Spearman ρ
between tier-1-predicted ranking and tier-2-measured ranking, per dataset.

### Tier 3 — cheap real-data proxy (all 700 real datasets)

`fastcheck` on the processed output: `read_retention`, output mean quality and
its tail behaviour (positions 1, 2, last 3), and `fastmaps` corner fractions
(`polyG_frac`, `polyT_frac`, `off_diagonal_frac`) on the output. Cost 0.4 s.
Used for *description* and as a smoke test, not as the training label — its
components are not independent of the trimmer being evaluated.

**Known dependency:** none of the three tiers independently measures *residual
adapter content* on real data. That needs a k-mer over-representation module
(fastcheck phase 3). Until it exists, real-data correctness leans on tier 2.
Listed as an open dependency, not hidden.

---

## 3. Data model

One row per evaluated configuration. Parquet, ~40 columns, one file per
dataset, then concatenated. Nothing in this table is data; it is measurements.

```sql
-- results/sweep/<dataset_id>.parquet, one row per (dataset_id, config_hash)
dataset_id        VARCHAR   -- sha1 of the manifest entry (source, accession, hash)
dataset_kind      VARCHAR   -- 'simulated' | 'real'
study_accession   VARCHAR   -- split unit; NEVER split inside a study
organism          VARCHAR   -- reference panel key, or 'none'
instrument        VARCHAR   -- NovaSeq6000 / NextSeq2000 / ...
stratum           TEXT[]    -- GC bin, read length, adapter prior, coverage bin
config_hash       VARCHAR   -- sha1 of the canonical config json
config_json       VARCHAR   -- canonical, sorted keys, exactly the argv that ran
argv              VARCHAR   -- the literal command line
score_version     VARCHAR   -- e.g. 't1-v1'
score             DOUBLE
score_components  VARCHAR   -- json: every counter above, unaggregated
pairs_in          BIGINT
pairs_out         BIGINT
bases_trimmed     BIGINT
elapsed_s         DOUBLE
cores             INT
tool_versions     VARCHAR   -- json: fastprep/fastcheck/fastsim + git hash each
host              VARCHAR
ran_at            TIMESTAMP
```

Directories (on cluster scratch, **never** in the repo):

```
$SCRATCH/fastsuite-bench/
  raw/          <study>/<run>_R{1,2}.fastq.gz     downloaded, rotating (deleted after subsample+hashing)
  work/         <dataset_id>/{in_r1,in_r2}.fastq.gz, gt.tsv, out/, pxmt/
  results/
    sweep/<dataset_id>.parquet        tier-1 rows
    real/<dataset_id>.parquet         tier-3 rows
    align/<dataset_id>.parquet        tier-2 rows
    matrix.parquet                    the concatenation = the training set
  models/       ranker_v1.txt (LightGBM text), export_v1.json (Rust-inference form), card.md
  logs/
```

`raw/` peak is bounded by a rolling window: download → subsample to 200k pairs
→ sha256 both the raw and the subsample → keep the raw only until the subsample
is verified → delete. Peak ≈ 60 GB, not 1.4 TB.

---

## 4. Phases

### Phase 0 — pilot harness (~2 core-hours, and the only phase to run before approval)

**Objective:** prove the plumbing end-to-end on 10 datasets (5 simulated, 5
real) × 30 configs, and produce the first "predicted vs achieved" scatter.

| task | file | done when |
|---|---|---|
| 0.1 Create the harness repo `~/hermes_scratch/FastSuite-bench` (sibling of FastSuite: the suite repo must stay shippable as bioconda packages, and this one carries heavy Python deps) | `pixi.toml` | `pixi run python -c "import polars, lightgbm, duckdb, snakemake"` |
| 0.2 Declarative sweep space, one entry per dimension with its type and range | `config/config_space.toml` | `pixi run sample-configs --n 30` prints 30 canonical argv sets |
| 0.3 Score module + its weights, unit-tested on hand-built fixtures (a GT with h=t=0, a correct trim, an over-trim, an under-trim, a wrong drop) | `scripts/score_gt.py`, `config/score_weights.toml`, `tests/test_score_gt.py` | `pytest` green; fixture scores match hand-computed values to 1e-9 |
| 0.4 One (dataset, config) runner, idempotent, writing exactly one result row | `scripts/run_config.py` | running it twice leaves one row; `--force` re-runs |
| 0.5 Five real datasets downloaded, subsampled, hashed, raw deleted | `scripts/ena_query.py`, `scripts/subsample.sh` | `manifests/real.tsv` has 5 rows with sha256, no data in the repo |
| 0.6 Snakefile: `acquire → simulate → sweep → score → aggregate` | `Snakefile`, `workflow/rules/*.smk` | `snakemake -n` plans the whole DAG; a full run finishes < 30 min on 8 cores |
| 0.7 Baseline rows for the same datasets: fastp defaults, fastprep `--adapter-preset nextera` alone, and "no preprocessing" | `scripts/run_baseline.py` | baselines appear in the same parquet as ordinary configs |

**Acceptance:** `results/matrix.parquet` with 300 tier-1 rows, a scatter of
predicted (in-sample, 10 points) and measured score, and a written list of the
top-3 configs per dataset. If the plumbing surprises us, it surprises us here,
for 2 core-hours.

### Phase 1 — real dataset corpus (~30 core-hours + network)

**Objective:** 700 real runs across strata, each reduced to a 200k-pair
subsample plus its QC profile.

* Query ENA (not SRA: ENA serves FASTQ over HTTPS, no `sra-tools`, no
  conversion step) via the portal API:
  `result=read_run`, `platform=ILLUMINA`, `library_strategy in (WGS, RNA-Seq,
  WXS)`, `instrument_model in (NovaSeq 6000, NextSeq 2000, HiSeq 2500,
  NovaSeq X)`, `fastq_ftp` non-empty, paired, `read_count` between 5 M and
  50 M.
* **Strata** (each with a floor of ~40 datasets): instrument model × GC bin
  (35/45/55 %) × read length (100/150) × organism class (bacterial, fungal,
  plant, metagenome, human) × library prep (Nextera vs TruSeq, inferred from
  the study description where it can be, `unknown` otherwise — recorded, never
  guessed silently).
* `study_accession` is the **split unit** (train/val/test 70/15/15, stratified
  by instrument). Splitting runs from the same study across train and test is
  the single most likely way to get a beautiful, meaningless result.
* Deterministic subsample: `seqtk sample -s 42` to 200k pairs (seqtk and the
  aligner are the only non-FastSuite dependencies in the whole pipeline).
* Parallel downloads: `aria2c -x8 -j4` or `curl --parallel`, 4 concurrent,
  retry-on-503, checksum against ENA's `fastq_md5` when present.

**Acceptance:** `manifests/real.tsv` (700 rows: accessions, strata, sha256 of
subsample, QC profile path), zero raw FASTQ left on scratch, all 700 tier-3 QC
reports present, and a per-stratum count table with no stratum below its floor.

### Phase 2 — simulated corpus (~30 core-hours)

**Objective:** the training corpus, with exact truth and wide parameter
coverage.

* 5 000 simulated datasets. Each draws *artefact parameters*, not just
  sequences: polyG/polyX rate and run length, adapter readthrough rate (0–15 %),
  insert distribution (150–600 mean), duplicates (0–10 %), per-base quality
  decay curve, GC (30–65 %), read length (50/75/100/150), seed.
* **Randomising the artefact parameters is what stops the model from learning
  `fastsim`'s quirks** instead of the general behaviour of a trimmer. Wide
  coverage is the mitigation for the model-realism risk in §6, and it also
  makes the corpus a genuine benchmark for fastprep rather than a fixture.
* Where a measured profile is available (the polyG scanner + `build_profile.py`
  lineage), include it as *one* point in the space, not as the default.
* `--read-len 100`, `--insert-min` and the drop class are the levers that make
  datasets genuinely different in how hard they are to trim; sample them, and
  keep the GT.

**Acceptance:** 5 000 `work/<dataset_id>/{in_r1,in_r2}.fastq.gz,gt.tsv`
generated and hashed in ≤ 30 core-hours, with a coverage report showing each
artefact dimension spanning its range.

### Phase 3 — the sweep (~620 core-hours, the bulk)

**Objective:** 1.06 M exact evaluations.

| tier | datasets | configs | evaluations | core-hours |
|---|---|---|---|---|
| simulated | 5 000 | 200 | 1 000 000 | 458 |
| real (tier-1 unavailable → tier-3) | 700 | 200 | 140 000 | 76 (fastprep + fastcheck = 1.95 core-s each) |
| **total** | | | **1.14 M** | **~534** |

* Configs per dataset: 150 Latin-hypercube samples from `config_space.toml`
  (covers the 12–14 dimensions without a grid blow-up) + 30 hand-picked
  anchors (fastp defaults, "trim nothing", "trim everything", preset
  nextera/TruSeq, the extremes of each axis) + 20 local-refinement samples
  around the best-so-far for that dataset (a second, cheap stage that fixes
  the "the model can only recommend inside the grid" problem).
* Dimensions (`config_space.toml`, generated from the real `fastprep --help`,
  never hand-invented): adapter mode (auto / preset illumina / nextera /
  small-rna / none), `--adapter-error-rate` (0.03–0.20),
  `--min-adapter-overlap` (4–14), `--min-run` (4–12), `--min-len` (30–0.9L),
  a 5′/3′/interior scan switch as a categorical, `--r2-scan-reverse` on/off,
  `--mean-q` (0 / 10 / 15 / 20 / 25), `--phred-offset` (auto / 33 / 64),
  `--trim-5`/`--trim-3` (0 / 1–5), `--dedup` on/off.
* Parallelism: **two levels, both used** — many samples in flight (the
  scheduler's job) and `fastprep -jN` inside a task (measured 2.8× at `-j8` on
  2 M pairs, 15.48 s → 5.45 s, thanks to the gz output thread pool). One SLURM array task = one dataset × all its
  configs (~200 × 0.6 s ≈ 2 min at `-j4`), so 5 700 tasks, walltime limit
  30 min, 4 cores, 2 GB each.
* Idempotency: a task writes `results/sweep/<dataset_id>.parquet.tmp`, then
  renames. Existence of the final file = done; nothing recomputes, nothing
  half-writes.
* Every row records `tool_versions` with the git hash of each binary: a
  result is only meaningful next to the binary that produced it.

**Acceptance:** `results/matrix.parquet` with ≥ 1.1 M rows, < 0.1 % failed
configs, every row carrying its argv and tool versions, and a per-dataset
check that at least 190 of the 200 configs returned a score.

### Phase 4 — features (~1 core-hour)

**Objective:** turn a QC report into the model's input vector.

* From `fastcheck --json`: reads/bases, mean length and its histogram moments,
  GC, mean quality, the quality-byte census (`phred_offset_certain`, offset),
  N rate, per-position summary of the quality plane (pos 1, 2, mean of the
  first 10, slope over 5→100, last-3 mean, min), per-position base-composition
  summary (max deviation from 0.25 in the first 10 positions and the last 10 —
  this is where polyG and adapter readthrough live), and the over-`max-read-len`
  fraction.
* From `fastmaps`: `polyG_frac`, `polyT_frac`, `polyC_frac`, `polyA_frac`,
  `off_diagonal_frac`, per R1/R2.
* From `fastsim` (simulated only): the true artefact rates. **These are
  training-only features and must be excluded from the real-data feature
  vector** — otherwise the model silently depends on information real data
  never has. The exclusion list is a test, not a convention.
* Config features: the parameter vector, one-hot for the categoricals.

**Acceptance:** a feature builder with a test asserting that every feature is
computable from a real dataset's artifacts alone, and that no `gt`-derived
field appears in the real vector.

### Phase 5 — the model (~2 core-hours, CPU only)

**Objective:** a ranker that beats the baselines on held-out studies.

* **Primary: LightGBM `lambdarank`**, groups = datasets, label = score
  (binned), ~300 trees, 63 leaves, `learning_rate 0.05`, on
  features × config. Trains in minutes on 1.1 M rows, ~1 MB on disk, and its
  text model is trivial to reimplement in Rust (~150 lines, no C dependency,
  keeping the suite's static-binary and no-shared-crate rules). This is the
  choice I would confirm before Phase 5.
* **Alternatives, only if the primary underperforms:** a small MLP
  (2 × 64, ~10 k params) exported to ONNX and run with `tract` (pure Rust —
  still no C dependency); or a direct QC→params regressor distilled from the
  ranker for a one-step `--auto-params`, accepting its instability.
* Loss/eval protocol: group split by study, never by read. Report:
  * **regret** `1 − score(best predicted) / score(true best in sweep)` — the
    business metric, targeted at ≤ 5 % median on the test split;
  * Spearman ρ between predicted and achieved score per dataset;
  * **top-5 overlap**: does the predicted-best fall inside the true top-5?
    (the honest way to handle near-ties);
  * per-stratum breakdown (instrument, GC bin, read length, organism) — a
    model that only works on NovaSeq is a finding, not a success;
  * tier-2 transfer: the same regret computed on the 50 alignment-validated
    real datasets.
* Baselines the model must beat, all present in the sweep as ordinary rows:
  `fastp` defaults, `fastprep` with `--adapter-preset nextera` and defaults,
  "no preprocessing", and the **oracle** (true best config in the sweep) as
  the ceiling.

**Acceptance:** model card (`models/card.md`) with regret, ρ, top-5 overlap,
per-stratum table, baselines, and the transfer gap — including the honest
cases where a baseline wins.

### Phase 6 — productisation (~2 core-hours)

* `fastprep --auto-params <qc.json> --print-args`: load the exported tree
  ensemble, evaluate the candidate set, print the argv. No new runtime
  dependency, works offline, and the model file is versioned next to the tool.
* `fasthelper recommend --qc runs/` for the batch case.
* `docs/` entry: what the model is, what it was trained on, what it must not
  be trusted on (out-of-distribution instruments, long reads — D-7).

**Acceptance:** `fastprep --auto-params` on a held-out QC report prints an
argv whose measured score is within the regret bound, and the flag is off by
default and refuses to run without a model file.

---

## 5. Cluster execution

* **Never compute on the login node** (submission only: `sbatch`, `squeue`,
  `sacct`). Everything below is submitted.
* Snakemake + the SLURM executor plugin, or, if the plugin is unavailable,
  `sbatch --array=1-5700%200` driving `run_dataset.sh` with an index file.
  Arrays are the right primitive: thousands of tiny, independent, 2-minute
  tasks.
* Resources: 1 task = 1 dataset, 4 cores, 2 GB, 30 min walltime. The sweep is
  embarrassingly parallel; the only shared resource is scratch I/O, so keep
  concurrent readers bounded (`%200`) and stage per-dataset work inside
  `$SCRATCH` (never `$HOME`).
* Aggregation after the sweep, on one node: DuckDB over
  `results/sweep/*.parquet` → `matrix.parquet` (a 1.1 M-row, ~200 MB table).
* Budget guard: every phase has a core-hour cap; the driver stops submitting
  when the cap is reached. Check `sacct -S <date> --format=JobID,Elapsed,AllocCPUS`
  at the end of each phase and write the actual cost into the phase log.
* Cleanup: `raw/` emptied after subsampling, `work/` kept only while a
  dataset is in flight, `results/` is small and kept (it is the deliverable).

---

## 6. Risks, and what each one costs if we ignore it

| # | risk | mitigation |
|---|---|---|
| 1 | **Simulator realism is the ceiling.** `fastsim`'s artefact model (hardcoded Nextera adapter, insert distribution as a CLI choice) may not cover real libraries → the model learns fastsim, not trimming | randomise artefact parameters widely (§Phase 2), include real datasets in training, report the tier-2 transfer gap explicitly, prefer a measured profile as one point in the space |
| 2 | **Objective mismatch.** Tier-1 exactness is only exact *within the simulator's model of what is wrong* | tier 2 (alignment) on a real subset, plus the tier-3 description; never claim real-data gains on tier-1 evidence alone |
| 3 | **Study leakage** → inflated scores | group split by `study_accession`, stratified; a test asserting no study appears in two splits |
| 4 | **Argmax instability** (many near-optimal configs) | ranking loss, near-tie tolerance `eps`, top-5-overlap reporting, regret instead of "did we pick the exact best" |
| 5 | **Grid-bound** recommendations | two-stage sweep: LHS + anchors, then local refinement around the best |
| 6 | **Compute overrun** | caps per phase, throttled arrays, pilot first (Phase 0 = 2 core-h), kill switch on `sacct` |
| 7 | **Reproducibility drift** (binary changes mid-sweep → mixed labels) | `tool_versions` + git hash in every row; the aggregation refuses to mix git hashes without `--allow-mixed` |
| 8 | **Data hygiene** | data never enters a repo; manifests carry accessions + sha256 only; `raw/` on scratch is deleted after subsampling; `.gitignore` in the harness blocks `*.fastq*`, `raw/`, `work/` |
| 9 | **fastprep lineage question** (which source tree is the shipped binary — still open in `DECISIONS.md`) | resolve before Phase 3: a sweep against a binary that is not the shipped one is worse than no sweep |

---

## 7. Budget summary

| phase | deliverable | core-hours | wall on 200 cores | wall on 32 cores |
|---|---|---|---|---|
| 0 pilot | plumbing proven, 300 rows | 2 | 1 min | 4 min |
| 1 acquire | 700 real datasets, manifests | 30 + network | ~1 h (bandwidth-bound) | ~1 h |
| 2 simulate | 5 000 datasets + GT | 30 | 10 min | 1 h |
| 3 sweep | 1.1 M scored configs | 536 | 2.7 h | 17 h |
| 4 features | feature builder | 1 | seconds | seconds |
| 5 train | ranker + model card | 2 | seconds | seconds |
| 6 productise | `--auto-params` | 2 | seconds | seconds |
| **total** | | **~600 core-hours** | **~4 h** | **~1.5 days** |

600 core-hours is the whole project. Storage peak: ~60 GB scratch. No GPU
needed at any point — this is the reason the plan prefers a 1 MB tree ensemble
over a neural model.

---

## 8. Open decisions (need a human, not a tool)

1. **Score weights** (`w_over = 0.5`, `w_under = 0.25`): is over-trimming
   really twice as bad as under-trimming *for your downstream use*? If the
   downstream is an assembler, yes; if it is a k-mer counter, less so.
2. **Cluster budget**: is ~600 core-hours available in one run, or should the
   sweep be staged over several days? (Staged is safe: the sweep is
   checkpointed per dataset, so it can be stopped and resumed.)
3. **Strata mix**: 700 real datasets weighted toward Illumina NovaSeq/NextSeq
   (the artefact model in fastsim is NovaSeq-like). Add long-read or 10x
   strata? D-7 says short-read only.
4. **Reference panel** for tier 2: the four genomes proposed, or the real
   organisms of the 50 validation datasets (more faithful, more download and
   more aligner index building).
5. **Where the harness lives**: sibling repo `FastSuite-bench` (my
   recommendation, keeps the suite shippable) vs `FastSuite/bench/`.
6. **LightGBM ranker as the model** — confirm before Phase 5.

---

## 9. Definition of done

A `fastprep --auto-params qc.json` that, on datasets it has never seen
(grouped by study), lands within 5 % median regret of the swept optimum,
beats `fastp` defaults on ≥ 70 % of held-out datasets, and carries a model card
stating its training distribution, its measured transfer gap to real data, and
the cases where it is known to be worse than the baseline.
