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

| tool | cores | wall (s) | CPU (s) | CPU/wall | peak RSS (MB) |
|---|---|---|---|---|---|
| fastprep | 1 | 73.99 | 73.83 | 1.00 | 103 |
| fastprep | 2 | 37.55 | 74.84 | 1.99 | 106 |
| fastp | 1 | 80.03 | 79.83 | 1.00 | 1171 |
| fastp | 2 | 40.55 | 80.32 | 1.98 | 1132 |
| cutadapt | 1 | 130.60 | 130.27 | 1.00 | 24 |
| cutadapt | 2 | 68.51 | 114.47 | 1.67 | 49 |

Full per-cell matrix (4 sizes × 3 tools × 2 core counts) is in
`FastPrep-Paper/benchmark/results.csv`. **That file is the source of truth, and
its `cpu_s/wall_s` ratio is the audit: a pinned run sits at 1.00, an escaped one
above it.**

`FastPrep-Bench/results/results_gz6.csv` is the earlier raw equal-level run and
**must not be quoted for a pinned cell**: its rows marked `pinned` carry
`cpu_s/wall_s` = 2.43–2.44 for fastprep (33.0 s wall, 80.5 s CPU) while fastp and
cutadapt sit at 0.98–0.99, so only the fastprep arm escaped `taskset`. Every
number derived from those rows — including the one-core win they attributed to
fastp — is superseded below.

### What the table says

* **fastprep leads fastp at every size, on one core and on two.** One core:
  1.51 vs 1.83 s at 100 k (1.21×), 7.41 vs 8.21 at 500 k, 14.79 vs 16.09 at 1 M,
  73.99 vs 80.03 at 5 M (1.08×). Two cores: 0.79 vs 1.01, 3.87 vs 4.19, 7.49 vs
  8.20, 37.55 vs 40.55 (1.08–1.28×). The previous revision of this table had
  fastp *leading* the 5 M one-core cell; that was the unpinned-harness artefact
  described above, not a fastprep deficit, and it is withdrawn.
* **Both tools scale ~1.98× on the second core when pinned.** The earlier claim
  that fastp's second core makes it *slower* came from the same unpinned rows
  (39.61 s at "two cores" against 28.5 s at "one") and is likewise withdrawn.
* **Memory is the decisive axis**: 103 MB vs fastp's 1171 MB at 5 M pairs, and
  32 MB vs 1148 MB at 100 k — a 36× gap that is structural, not a tuning
  artefact.
* **fastprep is 1.77× faster than cutadapt on one core** (73.99 vs 130.60 s),
  and cutadapt's second core is the one that genuinely helps it (130.60 →
  68.51 s).

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
* **PE overlap correction is now measured** (opt-in `--pe-correct`, off by
  default): posterior consensus over the mate overlap leaves 5.2% fewer
  residual fragment base errors than `fastp` 1.4.0 `--correction` on the
  20 000-pair truth corpus (22 190 vs 23 405; 13 vs 10 corrupted bases),
  per `bench/results/consensus-headroom*.tsv` and `docs/PE-CORRECTION.md`.
  **Not measured here:** read merging (`fastp`'s `-m/--merge` exists;
  merging is a different output object, not a comparison axis), UMI
  handling (out of scope), and stdin/stdout streaming (phase 3 --- landed
  in the tree since, per `docs/PLAN.md` §3.1 and the D-4 commits).

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
**Every fastmaps cell must name its `--buffer-mb`**: the 512 MB figure was the
default until 2026-10-18 and is 4× the peak of the 8 MB default, so a row that
omits it is unreadable and looks like a 9× memory gap against FastQC that the
shipped tool does not have.

### Runtime and memory (process-tree peak RSS, sampled every 20 ms)

| reads | tool | wall (s) | peak RSS (MB) |
|---|---|---|---|
| 100 k | fastmaps (8 MB buffer, all cores) | 0.18 | 10.6 |
| 100 k | fastqc 0.13.0 | 3.22 | 343.5 |
| 100 k | falco 2.0.2 | 0.27 | 18.5 |
| 5 M | fastmaps (8 MB buffer, all cores) | 5.78 | 22.2 |
| 5 M | fastmaps (512 MB buffer, all cores) | 5.83 | 123.8 |
| 5 M | fastqc 0.13.0 | 20.31 | 479.1 |
| 5 M | falco 2.0.2 | 8.48 | 641.2 |

**Read the RSS column as the tree, not the root.** `fastqc` is a Python
launcher that forks the JVM: its root process reads 16.4 MB and can never
produce the 479.1 MB cell, so a root-only poll reports a launcher, not a
footprint. Every number here is the summed `VmHWM` of every process in the
run's tree, which is also what `FastMaps-Paper/benchmark/results.csv` records.
Pinned one-core rows are in that CSV.

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
* **fastmaps' 145 MB peak RSS is decomposed and bounded, and it is not the
  grid.** Measured on a real 5 M-pair library (50 corpus datasets concatenated,
  10 M reads, 22 064 grid pixels; peak RSS read from `wait4(rusage)` of the
  child):

  | run | peak RSS before | peak RSS after |
  |---|---|---|
  | shipped buffer (512 MB until 2026-10-18; 8 MB since) | 202.6 MB | 141.2 MB |
  | `--buffer-mb 8` (now the default) | 100.5 MB | 36.1 MB |
  | `--preview-reads 0` | 126.3 MB | 126.3 MB |
  | `--buffer-mb 1 --preview-reads 0` | 16.5 MB | 16.0 MB |

  The terms are ~8 MB of process base, the worker index buffer, and the preview
  payloads. The index buffer holds one **8-byte** `(x, y, file_id, read_num)`
  entry per read (2+2+4; the budget charges 16) until `--buffer-mb` (default
  **8 MB**, divided across workers) is
  reached. At the old 512 MB default with two input files each worker's share was
  256 MB, so an entire
  5 M-read file (~60 MB of entries, ~110 MB with `Vec` growth) was retained and
  nothing spilled — the "budget" never bound. The preview store used to keep
  `preview_n` full read payloads per (pixel, mate) for the whole run (~76 MB,
  ≈15 MB per `--preview-reads` step). FastQC keeps neither a per-read record
  nor a read payload, yet it still peaks at 479 MB at 5 M reads: the 8 MB
  index buffer is not what makes fastmaps smaller. The comparison is 22.2 MB
  against 479.1 MB, not the reverse; the old "9× gap" came from reading only
  FastQC's launcher process.
* **The preview payloads are now streamed, not resident** (fixed 2026-10-18).
  The spiller writes each preview read to its row's temp file as it arrives
  (the worker flushes previews on their own 4 MB threshold, independent of the
  index buffer), and the writer reads one row file back at a time. The preview
  term fell from ~76 MB to ~3 MB per step, flat in library size, with
  **byte-identical** `heatmap.data` / `.preview` / `.json` / `.png` on the same
  inputs (only `heatmap_3d.html` varies, and it already varied run-to-run
  because it serialises a `HashMap`). The remaining term is the index buffer,
  a `--buffer-mb` policy, and the policy was pulled on 2026-10-18: at the old
  512 MB default a 5 M-pair run was ~141 MB, while the shipped 8 MB default
  gives ~36 MB with wall time unchanged (5.3–5.9 s both, i.e. within
  run-to-run noise). Measured independently on a 2 M-read library: 48.6 MB at
  512 MB against **22.7 MB at the shipped 8 MB**, 1.0 s both, outputs
  byte-identical.
* **What stays structural.** One 8-byte index record per read is written to a
  per-row temp file and re-read to build `heatmap.data` — that is the tool's
  contract (every read is addressable in the grid), so fastmaps cannot go below
  roughly one 8-byte record per read plus the 8 MB buffer — 22.2 MB at 5 M
  reads. That is already 22x below FastQC's measured 479 MB; the PLAN-QC
  rewrite is only needed for a sub-10 MB target, not for parity. Costed in
  `docs/PLAN-QC.md` §1.
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

**Conditions.** fastcheck *accepts* `-j` and defaults to all cores, but a single
input file is scanned on **one thread** whatever the flag says (`main.rs` takes
the serial branch when `files.len() == 1`, and the decompression pool is only
used for BGZF). So the single-file cells below are single-core numbers on both
sides, which is the fair comparison; the flag only parallelises across files.
FastQC 0.12.1 and Falco 2.0.2 are single-threaded.

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
| fasthelper-run (adapter guard on) | 4.437 | 5.300 | 110.97 | 2 | 2 |
| fasthelper-run (guard off) | 1.125 | 2.023 | 50.21 | 2 | 2 |
| bash-chain | 1.124 | 1.994 | 50.24 | 2 | 0 |

The runner adds **no statistically meaningful wall overhead** — 1.125 s against
1.124 s with the guard off, a difference inside measurement noise; the post-trim
adapter guard is an opt-in extra pass that costs 3.3 s and 61 MB — while producing
the two manifest entries the bash chain structurally cannot: the resolved
binary path, exact argv, exit code and runner-measured wall time per step,
schema-versioned as `fastsuite.manifest.v1`.

### Aggregation scaling (N = 25 … 800 samples)

| samples | fasthelper-report wall (s) / peak (MB) | multiqc wall / peak | multiqc output files |
|---|---|---|---|
| 25 | 0.003 / 17.18 | 1.304 / 137.28 | 0 |
| 50 | 0.005 / 17.31 | 1.270 / 137.96 | 0 |
| 100 | 0.010 / 17.31 | 1.318 / 137.74 | 0 |
| 200 | 0.019 / 17.43 | 1.473 / 138.07 | 0 |
| 400 | 0.037 / 17.81 | 1.840 / 138.43 | 0 |
| 800 | 0.075 / 18.43 | 2.439 / 139.45 | 0 |

The 17.2 MB floor is the harness process itself (`fasthelper-startup` with no
input); fasthelper's own share grows by 1.25 MB (+7 %) from 25 to 800 samples.

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
schema and the other does not". The 3–75 ms / 17.2–18.4 MB vs 1.3–2.4 s / 137–139 MB gap
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

### Profile recovery — 10 of 10 parameters within tolerance

| parameter | injected | recovered | abs. error | tolerance | pass |
|---|---|---|---|---|---|
| fragment GC | 0.48 | 0.4952 | 0.0152 | 0.02 | ✅ |
| R1 mean-Q curve | 34.0 | 34.0 | 0.0 | 1.0 Phred | ✅ |
| R1 worst 10-pos Q block | 0.0 | 0.05 | 0.05 | 1.5 Phred | ✅ |
| R1 tail run-length distribution | 0.0 | 0.013 | 0.013 | 0.2 TV | ✅ |
| R1 head run-length distribution | 0.0 | 0.005 | 0.005 | 0.2 TV | ✅ |
| adapter readthrough (R1) | 0.0495 | 0.04895 | 0.00055 | 0.02 | ✅ |
| duplicate rate | 0.02 | 0.01998 | 0.0000185 | 0.005 | ✅ |
| must-drop rate | 0.05 | 0.04895 | 0.00105 | 0.01 | ✅ |
| censored share | 0.5582 | 0.5318 | 0.0264 | 0.05 | ✅ |
| censored tail share | 0.34024 | 0.34027 | 0.00013 | 0.1 | ✅ |

**The three rows that ever failed were measurement defects in `fastcheck`,
quantified against ground truth, and no tolerance was widened to hide them.
All three are now fixed and the two detector fixes are recorded with their
mechanisms:**

* `tail_hist_R1_tv` (**0.28 → 0.013**) — fastcheck's profile fitter estimated
  the head/tail run-length distributions from the per-position G fraction,
  monotone-clamped and capped at 12 bases, so it reported 93.98 % of reads as
  tail-free against a true 66.01 % and had no bin at all at `t = 3` / `t = 15`
  where the truth is 7.99 % / 6.99 %. It now reads the exact per-base G-run
  histograms and subtracts the biological boundary run, giving
  {0: 66.6 %, 3: 7.98 %, 5: 13.6 %, 10: 4.75 %, 15: 6.38 %} against the true
  {66.0 %, 7.99 %, 14.0 %, 5.0 %, 7.0 %}.
* `max_q_block_dev_R1` (**1.95 → 0.05 Phred**) — the fitter stored the measured
  per-position mean as `pos_q`, but the simulator consumes `pos_q` as the
  out-of-run curve and applies the flat in-run mean itself, so the tail
  collapse was counted twice. It now emits the deconvolved out-of-run curve
  plus the measured `run_q_mean` / `deep_head_q` / `pos_q_uncond`.

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

* **The two detector defects are fixed** (R1 tail run-length distribution
  0.28 → 0.013; worst 10-position quality block 1.95 → 0.05 Phred). They were
  reported as failures while they were open; both now pass with the mechanism
  recorded in §5.
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
| fastmaps vs FastQC, 5 M pairs | `docs/STATUS.md` "Measured this week" | 12.25 s / 143 MB vs 72.22 s / 375 MB | **5.78 s / 22.2 MB vs 20.31 s / 479.1 MB**, tree-summed RSS (`FastMaps-Paper/benchmark/results.csv`); the old comparison paired a 512 MB-buffer cell with FastQC's launcher-only 16.4 MB |
| fastqc at 100 k | `FastMaps-Bench/.bench_run.log` | 4.36 s | 2.1 s (`FastMaps-Paper/benchmark/results.csv`) |
| fastprep vs fastp, 5 M pairs | `docs/STATUS.md` | 32.9 s vs 84.7 s (2.58×) | **73.99 s vs 80.03 s on one properly pinned core** — fastprep leads, 1.08× (`FastPrep-Paper/benchmark/results.csv`; the `pinned` rows of `results_gz6.csv` are NOT pinned, they burn 2.4 cores) |
| fastsim recovery | `docs/STATUS.md` | 7 / 10 parameters | 9 / 10 — the GC row is a **false pass** on the campaign build (`profile.gc` was inert); the field is honoured since 2026-10-21 but the row is unverified until the recovery is re-run (`FastSim-Paper/benchmark/recovery.csv`) |
| fasthelper orchestration | `FastHelper-Paper` abstract | "0.69 s vs 0.69 s" | 1.125 s vs 1.124 s (guard off); 4.437 s with the adapter guard |
| fasthelper aggregation peak | `FastHelper-Paper` abstract | "~50 MB" | 17.18–18.43 MB, of which 17.2 MB is the harness floor (multiqc 137–139 MB) |

Two of these are worth flagging beyond the arithmetic. **The fastprep one-core
cell has been re-measured twice for two different reasons**: the old harness was
not pinned (its `pinned` rows burn 2.4 cores), and the first correction of this
table carried the *unpinned* verdict forward as "fastp leads". With `taskset`
actually applied and `cpu_s/wall_s` ≈ 1.00 on every row, fastprep leads fastp at
every size on one and two cores. Second, **fastsim recovery is 10/10, not 7/10**
— the 7/10 figure is an older `FastSim-Bench` commit message and the 8/10 one the
pre-fix fitter, superseded by the final CSV and the paper's caption. The two papers whose CSV
differs from their own `-Bench` repo's CSV
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
