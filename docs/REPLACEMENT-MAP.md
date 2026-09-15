# Replacing FastQC/Falco, MultiQC and fastp

What the suite covers today, what it does not, and where each gap is closed.
Anything claimed here is either implemented and exercised in this repo, or
explicitly marked as a gap with the phase that closes it.

## fastp → fastprep (+ fastdedup)

| fastp capability | FastSuite | status |
|---|---|---|
| Illumina adapter trimming (auto-detect + presets) | `fastprep --adapter-preset illumina\|nextera\|small-rna`, `--adapters FASTA`, catalogue auto-detection | ✅ implemented; **auto-detection refuses to guess** (exit 1 with instructions) when the sampled reads carry no adapter — a deliberate difference from fastp, which trims on a guess |
| polyX / polyG trimming at 3′ | `fastprep` (quality-weighted, 5′ **and** 3′ **and** interior) | ✅ superset: fastp cannot touch 5′ polyG, the NovaSeq R2 artefact this suite was written for |
| polyX with error tolerance | `--match-leniency μ`, `--mismatch-strength ν`, `--mercy Δ_m`, `--gate-shape`, `--phred-sat` | ✅ |
| quality filtering | `--mean-q` | ✅ |
| hard trimming | `--trim-5`, `--trim-3` | ✅ |
| length filtering | `--min-len N|0.6L` | ✅ |
| PCR deduplication | `fastdedup` (separate binary) | ✅ and better: deterministic fingerprints, ~4.6 GB at 3e8 reads, no Bloom-filter false positives (`prinseq++` loses >90% of reads at that depth) |
| per-read metadata/traceability | `fastprep --metadata` (binary PXMT) + `fastprep-meta` reader | ✅ fastp has nothing equivalent |
| PE overlap correction / merge | **not implemented** | ⚠ gap — deliberately removed from the polyx lineage (it doubled memory for <1% gain). Sequence-merging users are not served yet; phase 4 decides re-add vs. declare out of scope |
| `-`/pipe in and out | not yet | ⚠ gap, phase 3 |
| HTML report per run | `fastprep --json` + `fasthelper` | ✅ different shape (one aggregate report for N samples) |
| UMI handling | no | ⚠ out of scope for now |

## FastQC / Falco → fastmaps (+ fasthelper)

| FastQC/Falco module | FastSuite | status |
|---|---|---|
| per-base sequence composition | `fastmaps` — reads placed in a `(T+C, T+G)` plane | ✅ but **not the same view**: FastQC plots per-position base fractions; fastmaps plots whole-read composition. The library-level questions (polyX contamination, R1/R2 asymmetry, off-composition mass) are answered in one image, the position-level questions are not |
| polyX/polyG contamination detection | `fastmaps` `polyG_frac`, `r2_polyG_frac` | ✅ arguably clearer than FastQC's per-base plot for this specific failure |
| per-base quality | not implemented | ⚠ gap, phase 4 — needs per-position accumulation; `fastsim` already models the quality curves it would be compared against |
| per-sequence quality / GC distribution | not implemented (GC implicit in the composition plane) | ⚠ gap, phase 4 |
| adapter content module | implicit: `fastprep`'s adapter stats | partial |
| duplication levels | `fastdedup` exact rate (not FastQC's estimate) | ✅ better: exact, no sampling |
| overrepresented sequences | not implemented | ⚠ gap, phase 4 (cheap once fastdedup's fingerprint table is reused) |
| k-mer content | not implemented | ⚠ gap, low priority |
| per-sample HTML/zip report | `fasthelper report/report.html` | ✅ different shape (multi-sample, one file, no zip) |
| multi-file / batch QC | `fasthelper run` | ✅ one config, N samples |

Honest summary: **fastmaps is not a drop-in FastQC replacement yet.** It replaces
the parts that matter for this pipeline (composition, polyX, duplication) and
beats FastQC on them, but per-base quality and k-mer modules are not there. The
phase-4 tasks in `PLAN.md` are what close it; until then the replacement claim
is "FastQC's composition modules, replaced and extended", not "all of FastQC".

## MultiQC → fasthelper

| MultiQC behaviour | FastSuite | status |
|---|---|---|
| discover per-tool outputs in a run directory | `fasthelper report <dir>` walks for `*.json` and auto-detects the producing tool | ✅ |
| aggregate heterogeneous tools into one report | normalises into one metric space; three adapters (v1 wrapper, fastprep raw stats, fastmaps grid) | ✅ |
| per-sample + per-module tables | `report.html` (per-tool sections), `report.tsv`, `report.md`, `report.json` | ✅ |
| interactive plots | static inline SVG bars on rate metrics | ⚠ multiline/threshold plots are phase 5 — deliberately no JavaScript in the report |
| module catalogue for 100+ tools | **not the goal** | ✅ by design: the suite reports on itself, so no parsing zoo to maintain |
| provenance (which tool version, which command) | `manifest.json` records resolved binary path, argv, exit code, wall time, log path for every step | ✅ stronger than MultiQC, which sees files, not runs |

## fasthelper beyond MultiQC

MultiQC only aggregates. Fasthelper also **runs** the pipeline
(`fasthelper run -c fastsuite.toml`), records the manifest, deletes each
intermediate once its last consumer has read it, and fails loudly with the log
path when a step exits non-zero. It is the only component here that is both an
orchestrator and a reporter — which is what makes "one command per run" and
"one report per run" possible.

## What is *not* in the suite's scope

* Assembly/binning/annotation (`megahit`, `kraken2`, …) — the suite is QC and
  preprocessing only.
* Alignment-dependent QC (`qualimap`, `mosdepth`) — no reference involved here.
* Workflow engine semantics (SLURM arrays, conda envs, resumability). Fasthelper
  deliberately stays a thin runner; a Snakemake/Nextflow wrapper around
  `fasthelper run` is the composable answer, not a re-implementation.
