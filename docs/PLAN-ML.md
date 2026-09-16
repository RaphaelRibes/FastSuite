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

Every number below is a 3-run mean with < 1 % variance, measured while writing
this plan. Two ways of measuring turned out to matter, so both are given:

| operation | work | wall | CPU actually used | verdict |
|---|---|---|---|---|
| `fastprep -j1` **pinned to one CPU** | 200 k pairs | **3.03 s** | 1.0 core-s | the honest cost of a 1-core task |
| `fastprep -j1`, unpinned | 200 k pairs | 1.27 s | 3.4 core-s | the 2.6× extra is **Rayon spin-waiting**, not work |
| `fastprep -j4` on 4 CPUs | 200 k pairs | **1.27 s** | 4.6 core-s | 2.4× less wall, +50 % core-seconds |
| `fastprep -j8` on 8 CPUs | 200 k pairs | 0.54 s | 5.6 core-s | flat past 4 — the win is the gz output pool |
| `fastcheck` score pass, 4 CPUs | 200 k-pair output | **0.38 s** | 0.86 core-s | |
| `fastsim --pairs` | 200 k pairs + GT | 19.9 s | 1.0 core-s | 10 050 pairs/s |
| fixed cost | — | ~0 | ~0 | 20 k-pair run: 0.135 s wall / 0.337 s core → cost is linear in pairs |

Three operational consequences, each of which changes the plan:

1. **`fastprep -j1` unpinned does not fit in one core** — 3.4 core-seconds of CPU
   for 1.27 s of wall, and `--gz-threads 1` does not change it. The extra CPU is
   Rayon workers spinning instead of parking: burned without doing work. So
   budget by *allocated cores × wall*, and pin tasks (cgroup / taskset) so the
   spin cannot steal from neighbours.
2. **More threads cost more core-hours, not fewer.** `-j4` is 2.4× faster in
   wall time but consumes 4.6 core-s instead of 3.0. On a cluster that bills
   core-hours, a 1-core task is the efficient shape; `-jN` is a latency lever to
   pull only when the queue is the bottleneck and the budget is not.
3. Per-pair cost is ~15–17 µcore-s under either accounting, so **the unit size
   is the budget dial**:

| sweep unit (pairs) | tier-1 evaluation | 1.14 M evaluations | wall on 200 slots | wall on 32 cores |
|---|---|---|---|---|
| 200 k | 3.1 core-s | 1013 core-h | 5 h | 32 h |
| **100 k (chosen)** | **1.57 core-s** | **507 core-h** | **2.5 h** | **16 h** |
| 50 k | 0.79 core-s | 254 core-h | 1.3 h | 8 h |

**Work unit: one dataset = 100 000 read pairs (~18 MB gz, ~2 MB gt).** That is a
0.15–0.22 % standard error on a per-read rate — below the `eps = 0.005` tie
tolerance already used when comparing configurations — and it halves the budget
for no loss of ranking signal. Task 0.8 verifies this empirically instead of
trusting the argument, and the unit size is one constant in the harness if the
measurement disagrees.

| per (dataset, config), 100 k-pair unit | core-seconds |
|---|---|
| `fastprep -j1`, 1 allocated core | 1.52 |
| tier-1 score (PXMT vs GT, streaming) | ≈ 0.05 (measured in task 0.8) |
| **tier-1 evaluation total** | **≈ 1.57** |
| tier-3 instead (fastcheck on the output) | 1.52 + 0.43 = 1.95 |
| tier-2 (alignment, validation subset only) | 30–60 |

---

## 2. The score: three tiers, one plug-in interface

The score is the heart of the project and the easiest place to lie to
ourselves. It is a versioned module with an explicit weights file, and every
result row records the score version that produced it.

### Tier 1 — exact, simulated data only (the training signal)

**As built:** `bench/scripts/score_gt.py`, 18 unit tests, verified against the
producers' own code rather than assumed — the definitions below were corrected
after the first pilot run contradicted an earlier draft.

`fastsim` writes `PREFIX.gt.tsv`, one row per pair:
`pair_id insert r1_h r1_t r2_h r2_t full_flag is_dup`. Per `fastsim/src/truth.rs`:

* `r{n}_h` / `r{n}_t` are the **5′/3′ G-run lengths** of that mate — the artefact
  spans, *not* "the trim a correct trimmer must make".
* The correct kept window is `[h, min(I, L−t))`, so with `L` the read length the
  ideal removals are: 5′ `= h`, 3′ `= L − min(I, L−t)`. When `I < L − t` the 3′
  removal includes adapter readthrough (e.g. `I=58, t=0, L=150` → remove 92).
* `full_flag = 1` → the window is empty and the correct action is **DROP**.

`fastprep --metadata` writes PXMT v3 (`fastprep/src/metadata.rs`) recording
**every** read — kept, trimmed or dropped — with per-end removals split into
adapter and polyX parts (`t5 = t5_adapter + t5_poly`). The join to GT is
asserted, not assumed, two independent ways: mates must alternate R1,R2,R1,R2,
and the read ID must resolve to the same pair index.

> **Flowcell geometry (found by that assertion, on the first real run):**
> fastsim's IDs are `instr:run:flowcell:lane:tile:x:y` with
> **`pair_index = y × TILE_WIDTH + x`**, `TILE_WIDTH = 10 000`. `x` resets to 0
> every 10 000 pairs, so a scorer assuming a monotonic coordinate refuses real
> files. A silent join error here would have mislabelled every row of the sweep.

**The score is config-aware**, because the same read is a correct drop under one
parameter set and a wrong drop under another:

| counter | definition |
|---|---|
| `exact` | kept, and both mates match the ideal removals exactly |
| `over_trim` | kept, real sequence removed beyond the artefact window |
| `under_trim` | kept, artefact left behind |
| `mixed_trim` | kept, one mate over- and the other under-trimmed |
| `drop_correct` | dropped, and dropping was right: `full_flag=1`, **or** the ideal window is shorter than the config's `--min-len` (default `0.6L`), **or** `is_dup=1` with `--dedup` on |
| `drop_wrong` | dropped a pair that should have been kept |
| `keep_wrong` | kept a pair whose window was empty or below `--min-len` |
| `partial` | one mate kept, one dropped — the pair is broken |

```
score = (exact + drop_correct − Σ wᵢ·errorᵢ) / pairs        # 1.0 = perfect
w: over 0.5, under 0.25, mixed 0.4, drop_wrong 0.5, keep_wrong 0.5, partial 0.5
```

Weights live in `config/score_weights.toml` and the score version in every
result row, so a weight change is visible in the data rather than silent.

**First real measurement (100 k pairs, nextera preset, defaults, this profile):**

| score | exact | over | under | mixed | drop_ok | drop_wrong |
|---|---|---|---|---|---|---|
| **0.277** | 0.450 | 0.274 | 0.090 | 0.140 | 0.044 | 0.0015 |

The loss is dominated by over-trim, which is a calibration decision, not a bug
report: fastprep's quality-weighted scan (mercy 3.0, `--min-run 8`) routinely
takes a few bases beyond the injected G-run. Either that is real over-trimming
(the current, strict reading — it is what a quality-agnostic ground truth says)
or the score needs an explicit tolerance of *k* bases at the artefact boundary.
**Decide before the sweep:** the choice sets the optimum the model will learn
toward. Task 0.9 measures how sensitive the ranking is to it.

**Why this matters:** it is exact (no proxy), cheap (~0.05 core-s), needs no
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
| 0.8 **Label-noise calibration** (decides the unit size): same dataset, same config, 10 subsample seeds, at 50 k / 100 k / 200 k pairs | `scripts/label_noise.py` | SD of the score per unit size; keep the smallest whose SD stays under `eps`/2 |
| 0.9 Score-tolerance calibration: re-score the pilot with a k-base tolerance at the artefact boundary (k = 0,1,2,3) | `scripts/score_gt.py --boundary-tolerance k` | ranking of the pilot's configs is reported for each k; the plan's tier-1 text records the decision |

**Acceptance:** `results/matrix.parquet` with 300 tier-1 rows, a scatter of
predicted (in-sample, 10 points) and measured score, and a written list of the
top-3 configs per dataset. If the plumbing surprises us, it surprises us here,
for 2 core-hours.

### Phase 1 — real dataset corpus (~15 core-hours + network)

**Objective:** 700 real runs across strata, each reduced to a 100k-pair
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
* Deterministic subsample: `seqkit sample -s 42` to 100k pairs (seqkit and the
  aligner are the only non-FastSuite dependencies in the whole pipeline; seqkit
  is already present under `~/.pixi/bin`, and it is the pixi-native choice).
* Parallel downloads: `aria2c -x8 -j4` or `curl --parallel`, 4 concurrent,
  retry-on-503, checksum against ENA's `fastq_md5` when present.

### Storage and transfer (from real ENA metadata, not estimates)

Queried the ENA portal API for Illumina WGS runs and read the published
`fastq_bytes`: 199 runs, of which 66 fall in this plan's strata (5 M – 50 M
reads). Compressed bytes per **base** is the stable metric — median **0.645**
(p25 0.532, p75 0.722, max 0.845); bytes per *read* swings 33–212 because read
length varies, so do not budget from it.

| quantity | median case | worst case |
|---|---|---|
| one run, R1+R2 gzipped | 2.00 GB | 4.57 GB |
| **700 runs, full download (transfer)** | **1.40 TB** | 3.20 TB |
| at 100 MB/s | 3.9 h | 8.9 h |
| **peak on disk** (download → subsample → hash → delete raw, 8 in flight) | **16 GB** | 36.6 GB |

Transfer is the cost, not capacity: `/scratch` has ~45 TB, and the *total*
footprint of the whole sweep peaks at **~43 GB** worst case (36.6 GB of raw
reads in flight, 6 GB of sweep working set at 100 concurrent shards — a shard of
10 datasets is 60 MB, measured — and ~0.4 GB of results Parquet).

**Do not try to be clever with byte ranges.** A `urllib`/`curl -r 0-9000000`
prefix is 200× smaller (a 100 k-pair unit needs only 9.2 MB, so ~11 GB for the
whole corpus), but only a *contiguous prefix* of a gzip member is decodable —
deflate state is sequential, so a range at 50 % cannot be decompressed without
everything before it. A prefix unit therefore samples only the start of a run,
where quality and artefact profiles genuinely drift. **Download whole runs
(the 1.4 TB is cheap at cluster bandwidth) and keep the prefix trick as a
fallback** if the link turns out slow — and then validate 50 runs against the
full file before trusting a prefix-derived recommendation.

Sanity check on the cost model: this bench file compresses to 0.624 bytes/base,
within 3 % of the real-data median of 0.645 — the synthetic units cost the same
to move as real ones.

**Acceptance:** `manifests/real.tsv` (700 rows: accessions, strata, sha256 of
subsample, QC profile path), zero raw FASTQ left on scratch, all 700 tier-3 QC
reports present, and a per-stratum count table with no stratum below its floor.

### Phase 2 — simulated corpus (~8 core-hours)

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
generated and hashed in ≤ 8 core-hours (100 k pairs ≈ 10 s each), with a coverage report showing each
artefact dimension spanning its range.

### Phase 3 — the sweep (~620 core-hours, the bulk)

**Objective:** 1.06 M exact evaluations.

| tier | datasets | configs | evaluations | core-s each | core-hours |
|---|---|---|---|---|---|
| simulated | 5 000 | 200 | 1 000 000 | 1.57 | 436 |
| real (tier-1 unavailable → tier-3) | 700 | 200 | 140 000 | 1.95 | 76 |
| **total** | | | **1.14 M** | | **~512** |

* Configs per dataset: **170 shared + 30 per-dataset refinement**. The shared
  core (Latin hypercube + anchors: fastp defaults, "trim nothing", "trim
  everything", the extremes of each axis) is evaluated on *every* dataset, so
  each configuration is observed ~5 700 times and its effect is measured rather
  than inferred; the 30 refinement samples are drawn around that dataset's
  best-so-far. An all-per-dataset draw would leave ~1.14 M mostly-unique rows —
  a ranker fitting noise. Same evaluation count either way.
* Dimensions (declared in `bench/config/config.yaml`, values taken from the real
  `fastprep --help`, never hand-invented): `--adapter-error-rate` (0.03–0.20),
  `--min-adapter-overlap` (4–14), `--min-run` (4–12), `--min-len` (0.3L–0.9L),
  `--mean-q` (0–25), `--phred-offset` (auto/33/64), interior scan on/off,
  `--r1-scan-reverse`, `--r2-scan-reverse`, `--trim-5`/`--trim-3` (0–5),
  `--dedup`.
* **Adapter mode is fixed to `auto`** (decided): `--adapter-preset`,
  `--adapters` and `--no-adapter` are not swept, because auto-detection is what
  a user with no prior knowledge of the library kit actually gets. The
  consequence is handled explicitly: auto-detection **refuses** (exit 1) when
  the sampled reads carry no adapter, so the simulated corpus floors
  `adapter_rate` at 2 %, and a refusal is recorded as its own outcome — never
  silently retried with a preset, which would change what is measured.
* **Only `--no-interior` is swept among the scan flags.** The two end-scans stay
  on (the simulator injects G-runs at both ends of every read, `[0,h)` and
  `[L−t,L)`, so disabling one leaves whole artefact runs → guaranteed
  under-trim), and the scan directions stay forward. Dropping the directions
  takes the swept space from 15 to 10 dimensions (coarse grid 34 M per dataset
  instead of 569 M).
  *Caveat recorded:* `fastsim/src/artefact.rs` models a **per-mate head CDF** and
  calls it "the 5′-R2 regime" — the simulator does inject 5′ R2 artefacts, which
  is exactly what `--r2-scan-reverse` exists for ("Recommended for
  Nextera/NextSeq/NovaSeq libraries" per `fastprep --help`). Leaving it out of
  the sweep means the model can never recommend it, on data that contains the
  artefact it addresses. Cheap to revisit: one boolean, 2 levels.
* Parallelism: **two levels, both used** — many samples in flight (the
  scheduler's job) and `fastprep -jN` inside a task (measured 2.8× at `-j8` on
  2 M pairs, 15.48 s → 5.45 s, thanks to the gz output thread pool). One SLURM array task = one dataset × all its
  configs (200 × 1.57 ≈ 5 min on one core), so 5 700 tasks, walltime limit
  30 min, **1 core**, 2 GB each. Raise to `-j4` only if the queue, not the
  budget, is the bottleneck.
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

**Io-specific mapping — partitions, accounts, env strategy, staging commands,
measured resource requests and the submit-for-confirmation block — is in
`docs/DEPLOY-IO.md`. It was written against the `novaseq-bench` deployment that
already works on Io, not from assumptions.**

* **Never compute on the login node** (submission only: `sbatch`, `squeue`,
  `sacct`). Everything below is submitted.
* Snakemake + the SLURM executor plugin, or, if the plugin is unavailable,
  `sbatch --array=1-5700%200` driving `run_dataset.sh` with an index file.
  Arrays are the right primitive: thousands of tiny, independent, 2-minute
  tasks.
* Resources per cell: **1 core, 256 MB** (measured peaks: fastprep 29.5 MB,
  fastcheck 10.7 MB, fastsim 46.5 MB — a large `--mem` trips the group quota and
  parks jobs in PENDING), flat 48 h walltime per the standing over-request
  policy. Packing: 10 datasets per job → 570 jobs of ~52 min; `-j 100` ≈ 5.5 h
  wall, `-j 200` ≈ 3 h. The sweep is
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
| 2 simulate | 5 000 datasets + GT | 8 | 3 min | 15 min |
| 3 sweep | 1.14 M scored configs | 512 | 2.6 h | 16 h |
| 4 features | feature builder | 1 | seconds | seconds |
| 5 train | ranker + model card | 2 | seconds | seconds |
| 6 productise | `--auto-params` | 2 | seconds | seconds |
| **total** | | **~550 core-hours** | **~3 h** | **~17 h** |

~550 core-hours is the whole project. Storage peak: ~60 GB scratch. No GPU
needed at any point — this is the reason the plan prefers a 1 MB tree ensemble
over a neural model.

---

## 8. Open decisions (need a human, not a tool)

1. **Score weights** (`w_over = 0.5`, `w_under = 0.25`): is over-trimming
   really twice as bad as under-trimming *for your downstream use*? If the
   downstream is an assembler, yes; if it is a k-mer counter, less so.
2. **Cluster budget**: is ~550 core-hours available in one run, or should the
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
