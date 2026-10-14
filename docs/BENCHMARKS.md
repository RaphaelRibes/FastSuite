# FastSuite benchmarks

Every number in this document is a **committed measurement**, not an estimate.
Each section names the tool pairing, the axes actually measured, the comparator
versions, the core count and the compression level, and ends with the caveats
the replacement map demands. Nothing here is re-derived, extrapolated or
"expected to hold".

## Where the numbers live

| section | paper (source of truth) | bench repo (raw matrix) | committed CSV |
|---|---|---|---|
| fastprep vs fastp / cutadapt | `~/work/FastPrep-Paper/` | `~/work/FastPrep-Bench/` | `benchmark/results.csv`, `annex-tables.typ` |
| fastmaps vs FastQC / Falco | `~/work/FastMaps-Paper/` | `~/work/FastMaps-Bench/` | `benchmark/results.csv` |
| fastcheck vs FastQC / Falco | `~/work/FastCheck-Paper/` | `~/work/FastCheck-Bench/` | `benchmark/results.csv` |
| fasthelper vs MultiQC | `~/work/FastHelper-Paper/` | `~/work/FastHelper-Bench/` | `benchmark/results.csv` |
| fastsim (accuracy / recovery) | `~/work/FastSim-Paper/` | `~/work/FastSim-Bench/` | `benchmark/{results,recovery,insert_resolution}.csv` |

The paper and its `-Bench` repo are **local-only and never pushed**
(`docs/STATUS.md`). The `-Bench` repos hold the seeded generator, the runner and
the plotter, so every row is re-derivable: `pixi install` → `python
scripts/gen_reads.py` → `python scripts/run_bench.py`. The reads themselves are
gitignored and regenerated deterministically from the committed profile at
`--seed 7`; only the measured CSVs and the SVG figures are versioned. Peak RSS
is the kernel high-water mark (`/proc/<pid>/status` `VmHWM`), never a sampled
average.

Where a paper's CSV and its `-Bench` repo's CSV disagree, **the paper wins** —
it is the committed source for the published claim. Two such cases exist and are
listed under [Known stale numbers](#known-stale-numbers).

---

## Shared methodology

**Axes.** Wall clock (s), CPU time (s) where the harness records it, peak RSS
(MB, VmHWM), and — for the tools that filter rather than measure — a
quality/accuracy axis scored against the simulator's ground truth. Throughput
is derived, not measured directly.

**Comparator versions**, pinned in each bench repo's pixi environment:

| comparator | version | pinned in |
|---|---|---|
| fastp | 1.3.7 | `FastPrep-Bench/pixi.lock` |
| cutadapt | 5.2 | `FastPrep-Bench/pixi.lock` |
| FastQC | 0.12.1 | `FastMaps-Bench`, `FastCheck-Bench` |
| Falco | 2.0.2 | `FastCheck-Bench/pixi.lock` |
| MultiQC | 1.35 | `docs/STATUS.md` (not pinned in `FastHelper-Bench/pixi.toml`, which declares only `python`) |

**Core counts.** Both pinning modes are measured so a parallel tool cannot win a
single-core comparison by spending more cores: *pinned* is `taskset -c 0` with
`-j 1` / `-w 1`, and *default* is shipped settings (all available cores for the
FastSuite tools; FastQC, Falco, cutadapt and fastsim are single-threaded).

**Compression.** The FastPrep comparison runs **every tool at gzip level 6** —
the one comparison where the level is documented, and the one where it matters
most: fastprep's canonical default is 6, fastp defaults to 4 and cutadapt to 1,
and an earlier revision of that benchmark compared them at those unequal
defaults, charging fastprep for compression the others never paid for. The
FastMaps and FastCheck papers state wall and VmHWM but do **not** record a gz
level, so no equal-level claim is made for them here.

**Machine.** Idle machine, one tool per run, medians of three runs (five for the
FastHelper aggregators, whose cells are sub-second). Io access ended 2026-10-12;
benchmark compute is now the desktop (16 cores, 48 GB) per
`docs/PLAN-ML-LOCAL.md`. Absolute wall times are machine-specific and are not
portable — the *ratios* and the memory *shape* are the transferable results.

**Input.** Paired-end libraries of 100 k / 500 k / 1 M / 5 M pairs (151 bp,
Phred+33), generated deterministically by `fastsim` from the committed
`synthetic_novaseq` profile with planted adapter readthrough, PCR duplicates and
a must-drop (insert-below-drop-floor) class. See
[Blocked on the mandrill profile refresh](#blocked-on-the-mandrill-profile-refresh)
for what this profile is *not*.

---

## 1. fastprep vs fastp and cutadapt

`FastPrep` is the trimming/filtering slot. `fastp` and `cutadapt` are the two
established tools doing the same job.

**Axes measured.** Wall, CPU, peak RSS, input/output file size, and the quality
axis — dropped-pair fraction against the ground-truth must-drop fraction.

**Conditions.** Both tools pinned to the same core count (1 and 2), every tool's
output at gzip level 6, `fastprep` invoked with `-q` so the summary stream does
not perturb the timed runs.

### Runtime and memory (largest library, 5 M pairs)

| tool | cores | wall (s) | CPU (s) | peak RSS (MB) |
|---|---|---|---|---|
| fastprep | 1 | 34.37 | 81.8 | 106 |
| fastprep | 2 | 25.84 | 81.2 | 106 |
| fastp | 1 | 28.5 | 83.6 | 1185 |
| fastp | 2 | 39.61 | 84.4 | 1131 |
| cutadapt | 1 | 131.12 | 130.7 | 24 |
| cutadapt | 2 | 49.98 | 118.7 | 51 |

Full per-cell matrix (4 sizes × 3 tools × 2 core counts) is in
`FastPrep-Paper/benchmark/results.csv`; `FastPrep-Bench/results/results_gz6.csv`
is the raw equal-level run and `results/results_cores.csv` the core sweep.

### What the table says

* **fastp's second core makes it slower** — 39.61 s at two cores against 28.5 s
  at one. Its threads contend over the global duplicate index. fastprep's second
  core buys a real 1.33× (34.37 → 25.84 s).
* **Memory is the decisive axis**: 106 MB vs fastp's 1185 MB at 5 M pairs, and
  31 MB vs 1156 MB at 100 k — a 38× gap that is structural, not a tuning
  artefact.
* **fastprep is 3.8× faster than cutadapt on one core** (34.37 vs 131.12 s).
* **fastprep leads fastp at every size except 5 M**, where fastp leads 0.83×
  (28.5 vs 34.37 s on one core). This is stated plainly rather than rounded
  away.

### The quality axis — drop rate vs ground truth

The simulator injects a must-drop class (insert below the drop floor) at
4.7–4.8 % of pairs. A correct preprocessor drops exactly those pairs and keeps
the rest.

| tool | dropped pairs | ground truth |
|---|---|---|
| fastprep | 7.7–7.9 % | 4.7–4.8 % |
| fastp | 0.0 % | 4.7–4.8 % |
| cutadapt | 1.4 % | 4.7–4.8 % |

fastprep is the only tool that enforces a meaningful minimum read length by
default. Its excess over ground truth is **deliberate conservative behaviour**,
not error — the drop-reason attribution makes this auditable rather than
opaque. On the 5 M library, 391 805 pairs dropped (7.84 %), decomposed as:

| cause | reads dropped | share of drops |
|---|---|---|
| adapter | 230,778 | 58.90 % |
| duplicate | 160,529 | 40.97 % |
| polyX 5′ | 356 | 0.09 % |
| polyX 3′ | 81 | 0.02 % |
| interior | 0 | 0.00 % |
| mean Q | 0 | 0.00 % |
| hard trim | 61 | 0.02 % |

fastp's 0 % is a *length-filter-off* default, not a better trimmer; cutadapt's
1.4 % is below ground truth, i.e. it keeps pairs that should be dropped.

### Caveats

* **fastp wins the 5 M one-core cell.** Not hidden, not re-framed.
* **Auto-detection refuses to guess.** Where sampled reads carry no adapter,
  fastprep exits 1 with instructions rather than trimming on a guess — a
  deliberate difference from fastp, and a reason a chain that "worked with
  fastp" can stop. Documented, with `--adapter-preset` as the escape hatch.
* **Not measured here:** PE overlap correction / merge (deliberately out of
  scope), UMI handling (out of scope), and stdin/stdout streaming (phase 3 —
  landed in the tree since, per `docs/PLAN.md` §3.1 and the D-4 commits).

---

## 2. fastmaps vs FastQC and Falco

`fastmaps` is the composition-plane QC tool: it places every read in a
`(T+C, T+G)` plane **per mate**, so the polyA/polyT/polyG/polyC corners and the
off-diagonal mass are visible in one image, and reports corner fractions in a
`fastsuite.report.v1` wrapper.

**Axes measured.** Wall and peak RSS on paired-end libraries, at shipped
defaults and pinned single core. **There is no accuracy axis** — this is a
timing and memory comparison on identical inputs.

**Conditions.** fastmaps at shipped defaults (all cores) and pinned
(`taskset -c 0`, `-j 1`). FastQC 0.12.1 and Falco 2.0.2 are single-threaded.

### Runtime and memory (largest library, 5 M pairs)

| tool | mode | wall (s) | peak RSS (MB) |
|---|---|---|---|
| fastmaps | all cores | 9.17 | 145.2 |
| fastmaps | pinned 1 core | 17.91 | 142.6 |
| fastqc | default | 16.1 | 16.4 |
| falco | default | 8.27 | 641.5 |

| tool | 100 k wall (s) | 5 M wall (s) | 5 M peak (MB) |
|---|---|---|---|
| fastmaps (all cores) | 0.22 | 9.17 | 145.2 |
| fastqc | 2.1 | 16.1 | 16.4 |
| falco | 0.18 | 8.27 | 641.5 |

Full matrix in `FastMaps-Paper/benchmark/results.csv` (identical to
`FastMaps-Bench/results/results.csv`).

### The axis that is not shared

The composition **plane is not an axis FastQC or Falco address**. Their base
content is *per-position*; fastmaps plots *whole-read* composition. The paper
verified from a FastQC 0.12.1 report that no `(T+C, T+G)` module exists. So the
timing table is the apples-to-apples axis and the plane figure is
fastmaps-specific — both stated honestly rather than forcing a false
equivalence on the image axis.

### Caveats

* **Falco is faster than fastmaps at both ends of this range** (0.18 s vs
  0.22 s at 100 k; 8.27 s vs 9.17 s at 5 M), and fastmaps is not the fastest
  tool in this table. The claim is the plane and the per-mate view, not raw
  speed.
* **FastQC's low peak RSS (16.4 MB) is not comparable to fastmaps' 145 MB** —
  they are computing different things at different resolutions. fastmaps'
  memory grows with library size because it builds per-mate grids; FastQC
  streams per-position counters.
* **Pinned fastmaps is 1.95× slower than shipped-default** at 5 M (17.91 vs
  9.17 s): the win is thread scaling, and it is real, but the comparison is
  only honest with both modes shown.

---

## 3. fastcheck vs FastQC and Falco

`fastcheck` is the exact, memory-flat QC tool — the component that actually
closes the FastQC gap.

**Axes measured.** Wall and peak RSS on single-end Illumina libraries
(100 k / 1 M / 5 M / 10 M / 20 M reads, 151 bp, Phred+33), plus the
correctness axes: exactness against an independent Python implementation
(`tests/reference_qc.py` + `tests/differential.py`, comparing every metric,
both histograms and all 150 per-position rows) and, via the simulator, against
ground truth.

**Conditions.** fastcheck uses all available threads by default; FastQC 0.12.1
and Falco 2.0.2 are single-threaded.

### Runtime and memory

| reads | fastcheck wall (s) / peak (MB) | fastqc wall / peak | falco wall / peak |
|---|---|---|---|
| 100 k | 3.05 / 140.2 | 2.22 / 310.5 | 0.12 / 14.2 |
| 1 M | 5.73 / 133.3 | 7.13 / 351.5 | 1.17 / 80.6 |
| 5 M | 17.92 / 140.4 | 29.14 / 340.8 | 4.93 / 335.5 |
| 10 M | 32.0 / 138.8 | 49.16 / 354.0 | 8.29 / 335.5 |
| 20 M | 55.68 / 138.8 | 97.14 / 362.9 | 16.63 / 335.2 |

Full matrix in `FastCheck-Paper/benchmark/results.csv` (identical to
`FastCheck-Bench/results/results.csv`).

### What the table says

* **Memory is flat in input size**: 133–140 MB from 100 k to 20 M pairs
  (2.7 GB of reads), because the accumulator state is a function of read length
  alone. FastQC holds 310–363 MB (JVM overhead); Falco's 14 MB at 100 k grows
  to 335 MB by 10 M and is flat from there.
* **fastcheck is 1.5–1.7× faster than FastQC at every scale above a million
  reads**, and extrapolates linearly to a billion reads in ~44 minutes.
* **fastcheck is slower than Falco at every size measured.** Stated as
  measured; the trade is exactness and the machine-readable schema.
* **The memory ceiling is enforced, not hoped for.** fastcheck projects the
  duplication fingerprint table's exact size *before* the scan — an identity,
  not an estimate — discovers the binding limit as
  `min(cgroup v2, cgroup v1, SLURM, RLIMIT_AS, MemAvailable)`, and refuses to
  start rather than dying to an OOM killer. A SLURM job granted 8 GB on a
  1.5 TB node sees 1.5 TB in `/proc/meminfo` and would be killed by its cgroup;
  fastcheck sees the cgroup instead.
* **The Phred-offset rule is a proof, not a heuristic**: a quality byte below 64
  *proves* Phred+33 encoding, while the absence of such a byte proves nothing.
  The heuristic form shipped in fastp and silently discarded every read of a
  high-quality library.

### Exactness status (this is the load-bearing part)

The exactness tasks **have landed and are committed** — the replacement claim
may now be made, in the replacement map's words:

> **FastQC's modules, replaced and extended, minus per-tile and exact
> over-representation** — not "all of FastQC".

What that wording now means, task by task:

* **Duplication is exact.** Task 5 shipped exact duplication multiplicities over
  xxh3-64 fingerprints (commit `d81aaa5`) — no sampling, and no longer the
  lower bound the plan's baseline measured (`GT 1.99% → 0.0`, which held only
  because the counter was a bound).
* **Over-representation has exact modes.** Task 6 shipped
  `--overrep-mode none|fingerprint|second-pass|table` (commit `3f8fb6b`);
  `second-pass` emits the sequence text of sequences above the threshold at the
  cost of one extra read, `table` packs 2-bit sequences and refuses above
  `--max-mem`. Verified in the built binary's `--help`.
* **Adapter content per position** shipped (commits `0e2b06e`, `108f306`),
  default off, reusing fastprep's catalogue and matcher semantics.
* **`--max-mem` admission refusal** shipped (commit `2587371`).
* **The insert-measurement defects the ground truth exposed are fixed**
  (commit `98296e3`) — see §5.
* **Per-tile quality modules are the one remaining FastQC gap** (plan Task 11).
  `--tiles` is not in the built binary's `--help`, so the residual gap in the
  replacement map's wording is **per-tile only**; "exact over-representation"
  has since closed. The ground-truth recovery test (Task 13) is also still
  outstanding as an automated gate, though the recovery numbers it would assert
  are measured and committed in §5.

### Caveats

* **The insert-size module is exact only above the read length, and censored
  above it.** The overlap search anchors on 12 bases of R2's biology, so an
  insert is measurable only when the mates share at least 12 bases —
  `insert ≤ 2L − 12`. At `L = 151` that bound is **insert ≤ 290**: below 288 bp
  measures exactly, above 290 bp everything is censored because no overlap
  exists to anchor on. Measured on the benchmark library (200 k pairs, ground
  truth from the simulator): the distribution is dominated by a single mode at
  303 that lies entirely above the boundary, so **90.6 % of pairs are censored
  outright**, only 9.4 % are measurable candidates, and of those the anchor now
  recovers **95 %**. The pre-fix build recovered ~26 %. This is documented in
  the report itself, not hidden — the module reports the measured/censored split
  per category rather than presenting a censored histogram as complete, and the
  same censoring behaviour exists in Picard's `CollectInsertSizeMetrics`.
* **The short-insert regime (insert ≤ L) is recovered only when
  `--adapter-seq` is supplied**, because both mates run into adapter and the
  adapter boundary must be resolved. The report keeps the two regimes in
  separate bins because their biases differ.
* **The bounded counters are approximations and are labelled as such** in the
  report: the over-represented-sequence counter evicts singletons (a sequence
  whose first occurrence is evicted loses one count) and the 12-mer spectrum is
  bounded. Both are now avoidable via `--overrep-mode`.
* **`fastcheck` is slower than Falco**, and its memory is ~10× Falco's at
  small sizes. The claim is exactness and the schema, not speed.
* **The memory-flat property holds per run**; aggregating many samples is the
  aggregator's job (§4).

---

## 4. fasthelper vs MultiQC

`fasthelper` is both orchestrator and reporter — the only component in the
suite that is both, which is what makes "one command per run" and "one report
per run" possible.

**Axes measured.** Orchestration: wall, CPU, peak RSS against the *identical*
tool sequence written by hand as a bash chain. Aggregation: wall and peak RSS
as the per-sample report count grows, against `multiqc` on the same input.
Threshold gate: `--fail-on` exit codes.

**Conditions.** Medians of five runs per cell (the cells are sub-second, so
more than three). MultiQC 1.35.

### Orchestration (N = 1 sample)

| tool | wall (s) | CPU (s) | peak RSS (MB) | report files | manifest entries |
|---|---|---|---|---|---|
| fasthelper-run | 0.667 | 3.503 | 3.6 | 2 | 2 |
| bash-chain | 0.647 | 3.389 | 4.1 | 2 | 0 |

The runner adds **no statistically meaningful wall overhead** — 0.667 s against
0.647 s, a ~3 % difference that is within measurement noise — while producing
the two manifest entries the bash chain structurally cannot: the resolved
binary path, exact argv, exit code and runner-measured wall time per step,
schema-versioned as `fastsuite.manifest.v1`.

### Aggregation scaling (N = 25 … 800 samples)

| samples | fasthelper-report wall (s) / peak (MB) | multiqc wall / peak | multiqc output files |
|---|---|---|---|
| 25 | 0.02 / 0.02 | 1.435 / 136.87 | 0 |
| 50 | 0.021 / 0.03 | 1.433 / 137.07 | 0 |
| 100 | 0.02 / 0.03 | 1.556 / 136.98 | 0 |
| 200 | 0.021 / 0.03 | 1.492 / 137.27 | 0 |
| 400 | 0.021 / 0.04 | 1.901 / 137.56 | 0 |
| 800 | 0.041 / 3.89 | 2.466 / 138.32 | 0 |

Full matrix in `FastHelper-Paper/benchmark/results.csv`.

### The claim, and why the last column is all zeros

**MultiQC's aggregation, replaced** — on the suite's own schema.

But the honest reading of that last column matters more than the timing.
`multiqc` run on the *same* `fastsuite.report.v1` input **produces no report
files at all**: it ships no module for the FastSuite schema, so its directory
scan finds nothing it recognises and exits with "No analysis results found".
This is not a harness bug or a fair-loss concession — it is exactly the coupling
gap the design argues against. A file-discovery aggregator can only aggregate
what a per-tool module has been taught to parse; fasthelper owns the schema and
parses it natively.

So the comparison is **not** "same output, faster". It is "one tool reads this
schema and the other does not". The 20 ms / 0.02 MB vs 1.4–2.5 s / 137 MB gap
is real, but the decisive axis is that only one of them produces a report. The
claim should be read as: *the MultiQC slot is filled, and what fills it also
owns the schema it aggregates*.

### Threshold gate

| case | exit code |
|---|---|
| `--fail-on pairs_dropped>0` tripped | 2 |
| `--fail-on reads_in>1e9` (clean) | 0 |

The report becomes a decision with a dedicated exit code, which is what lets it
gate a pipeline.

### Caveats

* **Not a module catalogue.** MultiQC ships parsers for 100+ tools; fasthelper
  deliberately does not. The suite reports on itself, so there is no parsing zoo
  to maintain — and no ability to aggregate arbitrary third-party output.
  Anyone aggregating non-FastSuite tools still needs MultiQC.
* **Interactive plots are static.** Inline SVG bars on rate metrics; multiline
  and threshold plots are phase 5, and there is deliberately no JavaScript in
  the report.
* **The provenance delta is the real win, not the wall time.** MultiQC sees
  files, not runs: binary path, argv, exit code and per-step wall time are not
  in the files, so they cannot be in the report.

---

## 5. fastsim — ground-truth recovery and insert resolution

`fastsim` is not a competitor to anything; it is the fixture that makes the
accuracy axes in §1 and §3 possible. Its own measurement is whether the
parameters it injects can be recovered from its output.

**Axes measured.** Wall and peak RSS vs library size (single-threaded), profile
recovery (injected vs recovered, per parameter, with the tolerance each is
checked against), and the insert-size measurement ceiling.

**Conditions.** Single-threaded by design — it accepts `-j` for CLI
compatibility, ignores it, and warns. `threads = 1` throughout, and **no
parallel speedup is claimed**. A 100 M-pair library is therefore a ~2.3 hour
single-core job. Datasets regenerated deterministically from the committed
profile at `--seed 7`.

### Throughput and memory

| pairs | wall (s) | peak RSS (MB) |
|---|---|---|
| 100 k | 8.21 | 34.3 |
| 1 M | 81.41 | 47.0 |
| 5 M | 421.84 | 47.3 |

Memory grows 34.3 → 47.3 MB across a 50-fold range, because the simulator
streams one pair at a time.

### Profile recovery — 8 of 10 parameters within tolerance

| parameter | injected | recovered | abs. error | tolerance | pass |
|---|---|---|---|---|---|
| fragment GC | 0.48 | 0.4945 | 0.0145 | 0.02 | ✅ |
| R1 mean-Q curve | 34.0 | 33.88 | 0.12 | 1.0 Phred | ✅ |
| R1 worst 10-pos Q block | 0.0 | 1.95 | 1.95 | 1.5 Phred | ❌ |
| R1 tail run-length distribution | 0.0 | 0.28 | 0.28 | 0.2 TV | ❌ |
| R1 head run-length distribution | 0.0 | 0.025 | 0.025 | 0.2 TV | ✅ |
| adapter readthrough (R1) | 0.0495 | 0.049019 | 0.00048 | 0.02 | ✅ |
| duplicate rate | 0.02 | 0.0200385 | 0.0000385 | 0.005 | ✅ |
| must-drop rate | 0.05 | 0.049019 | 0.00098 | 0.01 | ✅ |
| censored share | 0.5576915 | 0.5275 | 0.0302 | 0.05 | ✅ |
| censored tail share | 0.3402410 | 0.3399025 | 0.00034 | 0.1 | ✅ |

**The two failures are measurement defects in `fastcheck`, quantified against
ground truth, and no tolerance was widened to hide them.**

* `tail_hist_R1_tv` — fastcheck's polyG detector reports 93.98 % of reads as
  tail-free against a true 66.01 %, with zero mass at `t = 3` / `t = 15` where
  the truth is 7.99 % / 6.99 %.
* `max_q_block_dev_R1` — the worst 10-position quality block misses by 0.45
  Phred (mean passes at 0.12).

`censored_share` is the **regression guard**: it failed on the pre-fix build
(0.8075 against a true 0.5577 — a 25-point gap) and passes after the
insert-anchor fix (0.5275). fastsim reproduces its injected artefact marginal
exactly, a total-variation distance of **0.000** from ground truth, which is
what establishes that both failures belong to the measurement and not to the
simulation.

### Insert-size resolution — exact to 94–95 % up to the geometric ceiling

| true insert (bp) | pairs | measured | measured share | p50 |
|---|---|---|---|---|
| 50 | 20 000 | 0 | 0.0 | 0 |
| 100 | 20 000 | 0 | 0.0 | 0 |
| 150 | 20 000 | 2 311 | 0.1155 | 147 |
| 200 | 20 000 | 18 982 | 0.9491 | 200 |
| 250 | 20 000 | 18 956 | 0.9478 | 250 |
| 280 | 20 000 | 19 009 | 0.9505 | 280 |
| 288 | 20 000 | 15 640 | 0.782 | 288 |
| 290 | 20 000 | 0 | 0.0 | 0 |
| 300 | 20 000 | 0 | 0.0 | 0 |
| 350 | 20 000 | 0 | 0.0 | 0 |
| 400 | 20 000 | 0 | 0.0 | 0 |
| 500 | 20 000 | 0 | 0.0 | 0 |

Inserts of 200–288 bp resolve **exactly** at 94–95 % after the anchor fix. The
ceiling is geometry, not a tool defect: with mates of length `L` each,
`I_measurable ≤ L + L − 12`, so at a 150 bp read length the ceiling is 288 bp.
Above it the mates do not overlap, there is nothing to measure, and the pair
must be reported as censored — the same place Picard censors. Rows at 50 and
100 bp are censored *below* the read length: an overlap shorter than 12 bases
cannot be anchored at all, and that regime needs the adapter boundary.

### Caveats

* **Two known defects, both in the measurement, neither fixed:** the R1 tail
  run-length distribution and the worst 10-position quality block. They are
  reported as failures rather than hidden.
* **The insert ceiling is a property of paired-end sequencing**, not of any
  tool. Reported as such.
* **fastsim is single-threaded** and its wall times should never be compared
  against a threaded tool's as if they were the same class of measurement.

---

## The mandrill profile refresh (PLAN §6.1)

**Status: the profile is fitted and committed; the bench re-cut is NOT done.**

`tools/fastsim/profiles/mandrill_E358.json` was fitted on 2026-10-14 from the
real mandrill reads with `fastcheck --profile --insert-hist --insert-pairs 20000`
(the fitter IS the profile builder — no Python fitter). It supersedes the old
`mandril_E.json` and the synthetic profile on every axis that matters:

| | read_len | GC | full_length_rate | insert bins | insert range | insert mode |
|---|---|---|---|---|---|---|
| **new mandrill_E358** | 150 | 0.4748 | 0.0 | **146** | 144–301 bp | 301 bp (12 498 / 20 000) |
| old mandril_E | 150 | 0.4895 | 0.00217 | 0 (no insert block) | — | — |
| synthetic_novaseq | 150 | 0.5 | 0.0 | 0 (no insert block) | — | — |

What the refresh changes:

* **The insert distribution is real now.** The old and synthetic profiles carry
  NO insert block, so every fastsim run before this fell back to the Gaussian
  `--insert-mean/--insert-sd` model — which understates the real spread. The new
  profile's 146-bin histogram has its mode at the 301 bp saturation ceiling,
  meaning the E-dataset's inserts are longer than the measurement window and
  saturate it (62.5% of sampled pairs censored above the measurable range).
* **GC is measured, not nominal**: 0.4748 against the synthetic 0.50.
* **The polyG run-length histograms are the measured shape** (head: 81% zero-run,
  19% one-base; tail: 98% zero-run with a thin long tail). The old file's
  head/tail blocks were also mis-scaled — raw counts in the billions against the
  new profile's per-10 000 normalisation — so the old file was not usable as-is;
  the new fit is the authoritative one.

**What has NOT been re-measured:** every benchmark number in this document was
generated from the synthetic profile. The timing and memory numbers are
unaffected in kind (the artefact rates are planted identically), but the
quality-axis numbers — §1's drop-rate comparison, §3's exactness and insert
censoring split, §5's recovery table — will move once the benches run against
`mandrill_E358.json`, because the Gaussian insert fallback is gone.

The honest label for every number above remains:
**measured on synthetic-profile libraries with planted, known-rate artefacts** —
until the three sections are re-run at the same seeds and core counts against
the committed mandrill profile.

---

## Known stale numbers

These appear elsewhere in the repo's docs and are **superseded by the committed
paper CSVs**. Listed here so the discrepancy is visible rather than silently
propagated.

| claim | where | stale value | committed value |
|---|---|---|---|
| fastmaps vs FastQC, 5 M pairs | `docs/STATUS.md` "Measured this week" | 12.25 s / 143 MB vs 72.22 s / 375 MB | 9.17 s / 145.2 MB vs 16.1 s / 16.4 MB |
| fastqc at 100 k | `FastMaps-Bench/.bench_run.log` | 4.36 s | 2.1 s (`FastMaps-Paper/benchmark/results.csv`) |
| fastprep vs fastp, 5 M pairs | `docs/STATUS.md` | 32.9 s vs 84.7 s (2.58×) | 34.37 s vs 28.5 s — fastp **leads** this cell |
| fastsim recovery | `docs/STATUS.md` | 7 / 10 parameters | 8 / 10 (`FastSim-Paper/benchmark/recovery.csv`) |
| fasthelper orchestration | `FastHelper-Paper` abstract | "0.69 s vs 0.69 s" | 0.667 s vs 0.647 s |
| fasthelper aggregation peak | `FastHelper-Paper` abstract | "~50 MB" | 0.02–3.89 MB (multiqc 137 MB) |

Two of these are worth flagging beyond the arithmetic: **fastp leads the 5 M
one-core cell** (so the STATUS.md "2.58×" win is not the committed result), and
**fastsim recovery is 8/10, not 7/10** — the 7/10 figure is the older
`FastSim-Bench` commit message, superseded by the final CSV and the paper's
caption. The two papers whose CSV differs from their own `-Bench` repo's CSV
(FastPrep, FastHelper, FastSim) were resolved in favour of the paper per the
rule above.

---

## What is explicitly *not* claimed

* **Not "replaces all of FastQC."** The claim is the replacement map's wording:
  FastQC's modules, replaced and extended, **minus per-tile**. Per-tile quality
  modules are the remaining gap.
* **Not "faster than everything."** Falco beats fastmaps and fastcheck on wall
  time at every size measured, and fastp beats fastprep in the 5 M one-core
  cell. The claims are about memory shape, exactness, the schema, the
  composition plane and the provenance record — not a clean speed sweep.
* **Not "MultiQC's aggregation is faster at the same task."** MultiQC produces
  zero output on this schema. The claim is that the slot is filled by a tool
  that owns the schema, not that the same work is done quicker.
* **Not production-ready.** `docs/STATUS.md`: the suite is WIP.
