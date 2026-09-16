# PLAN-ML — parameter recommendation from QC (final)

**Revision 2 (final).** Supersedes the v1 narrative: several of its statements were
contradicted by measurement and are corrected here. Every line below is either a
**decision record** or a **measurement** — nothing is an estimate unless it says
so.

**Goal.** `fastprep --auto-params <qc.json>` recommends the preprocessing
parameters for a dataset from its QC profile alone — the manual tuning loop
fastp makes you do by hand.

**Shape.** Learn `f(QC features, config) → score`, then take the argmax over a
candidate set. Not `QC → params`: dozens of configurations are near-optimal for
the same data, so the argmax is a multimodal target, while scoring a candidate
is well-posed and every sweep row becomes a training row.

---

## 1. What is built and verified (workstation)

| component | state | evidence |
|---|---|---|
| `bench/scripts/score_gt.py` | exact tier-1 scorer: PXMT v3 vs ground truth, config-aware | **19 unit tests green**, refuses a misaligned join |
| `bench/scripts/sample_configs.py` | Latin hypercube + anchors, repairs fastprep's own feasibility couplings | 40 sampled points, 0 refusable |
| `bench/scripts/run_shard.py` | N datasets × M configs → one TSV row each; crash-safe, idempotent | 400/400 rows scored in 244 s |
| `bench/scripts/gen_profiles.py` | randomised synthetic libraries | 12 profiles, 0 % – 44 % adapter readthrough, all load |
| `bench/scripts/run_pilot.sh` | cluster probe: build, egress, env, one dataset E2E | `PROBE PASS` |
| `deploy/io/*.sbatch`, `slurm/config.yaml` | the SLURM headers (pilot, dispatcher, acquire) + cell profile | ready, not submitted |

Verified end to end: profile generation → 8-dimension sampling → fastprep →
PXMT parsing → exact scoring → results table.

**Findings that changed the plan** (full detail in `bench/FINDINGS.md`):
the PXMT→GT join is flowcell coordinates (`pair = y×10000 + x`, `x` resets every
10 000) — a monotonic assumption would have mislabelled **100 %** of rows; the
ground truth's `r*_h`/`r*_t` are G-run lengths and the correct window is
`[h, min(I, L−t))`; the strict score is only usable once the corpus contains the
artefact it scores; and a corpus whose artefact model does not vary trains a
model to output a **constant**.

---

## 2. Decisions (frozen)

| decision | value | why |
|---|---|---|
| adapter mode | **`auto` only** | the deployment case: no prior knowledge of the kit. Its refusal on an adapter-free dataset is handled by the corpus floor (≥ 2 % readthrough) and recorded as its own outcome — never retried with a preset, which would change what is measured |
| Phred offset | **`auto`** | fastprep decides from the quality-byte census and reports its certainty |
| interior scan | **always `--no-interior`** | the corpus injects artefact only at the ends and the GT window carries no interior annotation, so nothing measurable is given up |
| dedup | **always `--dedup`** | measured cost 31.5 MB vs 27.5 MB peak, 0.54 s vs 0.59 s — and the scorer then runs with `--dedup`, because dropping an `is_dup` pair is always correct |
| scan directions | **not swept** | *caveat on record:* fastsim models a per-mate head CDF ("the 5′-R2 regime"), so it does inject 5′ R2 artefacts and `--r2-scan-reverse` is what fastprep's help recommends for them. The model can therefore never recommend that flag. Re-adding is one boolean |
| swept dimensions | **8**: `adapter_error_rate` (0.03–0.20), `min_adapter_overlap` (4–14), `min_run` (4–12), `min_len` (0.3L–0.9L), `mean_q` (0/25), **`head_dropout` (0–4)**, `trim_5`, `trim_3` (`{0,2,5,7,10,12,15,17,20}`) | `--head-dropout` is the flag built for the NovaSeq first-cycle G-dropout; a hard trim was a blunt stand-in for it, and 1–3 base runs are too short for the polyX scan to find |
| sampling design | **170 shared + 30 per-dataset refinement** | each shared configuration is then observed ~5 700 times instead of ~1.14 M mostly-unique rows |
| boundary tolerance | **K = 0, strict** | with the corpus now carrying the artefact, exactness is reachable (0.468 on the winner), so no tolerance is needed to make the score informative |
| corpus | **per-dataset randomised synthetic profiles** + 700 real runs | a single profile gave the *same* winner on 10/10 datasets |
| split unit | `study_accession` | splitting runs of one study across train/test is the fastest way to a meaningless result |

Grid if enumerated: **5 056 128 combinations per dataset** (28.8 billion
corpus-wide, ~15 years at 100 cores) — hence sampling, at 200 per dataset.

---

## 3. The remaining work, in order, with cost

| # | step | where | cost | note |
|---|---|---|---|---|
| 1 | Snakefile wiring the shard runner + corpus index | local | ~1 core-h | the last code before the cluster |
| 2 | generate the corpus: 5 000 randomised profiles + datasets | Io | 8 core-h | 100 k-pair units |
| 3 | **the sweep: 1.14 M exact evaluations** | Io | **512 core-h** | 5 700 datasets × 200 configs; ~5 h wall at 100 cores |
| 4 | features + train the ranker | **local** | 3 core-h | LightGBM is not installable on Io and does not need to be |
| 5 | tier-2 alignment validation on a 50-dataset subset | Io | 12 core-h | the transfer check |
| 6 | productise `--auto-params` + model card | local | 2 core-h | ~150 lines of Rust tree-walking, no C dependency |
|  | **total outstanding** |  | **~537 core-h** | peak disk ~43 GB against ~45 TB free |

Cluster shape (partitions/accounts/time are in `deploy/io/DEPLOY-IO.md`, with the
`scp` + `ssh Io` + `srun --time=01:00:00` workflow): one dispatcher under
`@cirad-long`, 570 one-core cells under `@cirad-normal`, 256 MB each, **hard
ceiling 100 cores** enforced by the profile's summed resource cap.

---

## 4. Acceptance criteria

1. **regret** `1 − score(recommended)/score(best in sweep)` ≤ 5 % median on
   held-out studies;
2. beats **fastp defaults** on ≥ 70 % of held-out datasets, with the oracle as
   the ceiling;
3. a **model card** stating the training distribution, the per-stratum
   breakdown, and the measured tier-2 transfer gap — including the strata where
   a baseline wins;
4. the achievable ceiling is reported too: on some libraries no configuration
   does well (measured: −0.08 on the hardest of 10), and saying so is part of
   the deliverable.

---


---

# Appendices — the measured detail behind the plan

## Appendix A — cost model

Every number is a 3-run mean measured on the workstation; the pinned-CPU
method matters (see the note on Rayon spin-waiting).

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

## Appendix B — the three scoring tiers



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

## Appendix C — result schema and layout



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
## Appendix D — corpus storage and transfer (measured)

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


---

## 5. Live risks

| risk | mitigation |
|---|---|
| simulator realism caps everything | artefact parameters randomised per dataset; hand-authored distributions describe mechanisms, never measured data; tier-2 alignment measures the gap |
| the optimum is dataset-independent after all | the 10-dataset test showed 3 distinct winners of 10; if 200 configurations shows one winner, the honest output is a good default, not a model — that is a result, not a failure |
| 12/36 → 0 refused, but more couplings may lurk | the sampler repairs what fastprep refuses; the runner lists any dataset that yields no rows |
| data hygiene | reads never enter a repo; the corpus is generated, the real runs stay on scratch, manifests carry accessions + hashes only |

---

## 6. Falsifiers — what would change this plan

* a corpus where the winner is the same for every dataset ⇒ ship a default, drop the model;
* a tier-2 transfer gap so large the simulated ranking does not predict real ranking ⇒ the training signal must come from alignment-scored real data, at ~30× the cost;
* `--head-dropout` never chosen by the optimizer on profiles that contain dropout ⇒ the flag is not doing the work the data says it should, and the corpus is wrong.
