# Replacing FastQC/Falco, MultiQC and fastp

What the suite covers today, what it does not, and where each gap is closed.
Anything claimed here is either implemented and exercised in this repo, or
explicitly marked as a gap with the phase that closes it.

## fastp → fastprep (+ fastdedup)

| fastp capability | FastSuite | status |
|---|---|---|
| Illumina adapter trimming (auto-detect + presets) | `fastprep --adapter-preset illumina\|nextera\|small-rna`, `--adapters FASTA`, catalogue auto-detection | ✅ implemented; when the sampled reads carry no adapter, auto-detection **takes the safe action instead of guessing** — no adapter trimming, the assumption recorded in the report and on stderr, exit 0 (D-6; `--strict` restores the refusal). Still a deliberate difference from fastp, which trims on a guess |
| polyX / polyG trimming at 3′ | `fastprep` (quality-weighted, 5′ **and** 3′ **and** interior — interior is opt-in and off by default) | ✅ superset: fastp cannot touch 5′ polyG, the NovaSeq R2 artefact this suite was written for |
| polyX with error tolerance | `--match-leniency μ`, `--mismatch-strength ν`, `--mercy Δ_m`, `--gate-shape`, `--phred-sat` | ✅ |
| quality filtering | `--mean-q` | ✅ |
| hard trimming | `--trim-5`, `--trim-3` | ✅ |
| length filtering | `--min-len N|0.67L` (the floor is also the trim budget: 0.67L keeps 100 bp of a 150 bp read and drops a read that would need more trim) | ✅ |
| PCR deduplication | `fastprep --dedup` (native, shipped ON) **or** the standalone `fastdedup` — same fingerprint model (64/128-bit, deterministic, ~4.6 GB at 3e8 reads, no Bloom-filter false positives — `prinseq++` loses >90% of reads at that depth). Use `fastprep --dedup` for the one-pass chain; use `fastdedup` when you need dedup as a separate, pipeable step | ✅ **two entry points, one implementation model** — fastprep's is the integrated one; fastdedup (v1.3.0, published) is the standalone one |
| per-read metadata/traceability | `fastprep --metadata` (binary PXMT) + `fastprep-meta` reader | ✅ fastp has nothing equivalent |
| PE overlap correction / merge | `--pe-correct [GATE]` (opt-in) corrects base errors in the doubly-retained mate overlap by posterior consensus at a measured margin over `fastp --correction` (5.2% fewer residual base errors at gate 2.0 on the 20k-pair truth corpus); `--pe-overlap` infers the adapter boundary from the mate pair's own overlap. **No sequence MERGER** (no `--merge`, no merged-output mode) | ⚠ partial — correction shipped, merging still a gap (deliberately removed from the polyx lineage: it doubled memory for <1% gain). Sequence-merging users are not served |
| `-`/pipe in and out | not yet | ⚠ gap, phase 3 |
| HTML report per run | `fastprep --json` + `fasthelper` | ✅ different shape (one aggregate report for N samples) |
| UMI handling | no | ⚠ out of scope for now |

## FastQC / Falco → fastmaps (+ fasthelper)

| FastQC/Falco module | FastSuite | status |
|---|---|---|
| per-base sequence composition | `fastmaps` — reads placed in a `(T+C, T+G)` plane | ✅ but **not the same view**: FastQC plots per-position base fractions; fastmaps plots whole-read composition. The library-level questions (polyX contamination, R1/R2 asymmetry, off-composition mass) are answered in one image, the position-level questions are not |
| polyX/polyG contamination detection | `fastmaps` `polyG_frac`, `r2_polyG_frac` | ✅ arguably clearer than FastQC's per-base plot for this specific failure |
| per-base quality | `fastcheck` — exact per-position quality matrices + quantiles from the histogram | ✅ implemented (`--positions-tsv`, `mean_q/q10/q50/q90` columns) |
| per-sequence quality / GC distribution | `fastcheck` — `read_mean_q_hist`, `read_length_hist`, per-read GC (`gc_fraction`, `gc_hist`) | ✅ implemented |
| adapter content module | `fastcheck --adapter-content` (opt-in) plus `fastprep`'s adapter stats | ✅ implemented; fastcheck's block is opt-in, fastprep's stats are on every run |
| duplication levels | `fastdedup` exact rate (not FastQC's estimate) | ✅ better: exact, no sampling |
| overrepresented sequences | `fastcheck` bounded counters (singleton eviction, top-20 above `max(10, 0.1%)`) | ⚠ **approximate by default**: a sequence whose first occurrence is evicted loses 1 count (documented). Exact modes ship — `--overrep-mode none\|fingerprint\|second-pass\|table`, default `none` — see `docs/plans/2026-10-12_fastcheck-exactness.md` |
| k-mer content | `fastcheck` bounded 12-mer spectrum (top-20) | ⚠ **partial**: bounded approximate for the same reason; per-position exact k-mer (`kmer_pos[p][kmer]`) is planned in the same fastcheck plan |
| per-sample HTML/zip report | `fasthelper concatenate` → `report.html` | ✅ different shape (multi-sample, one file, no zip) |
| multi-file / batch QC | `fasthelper run` | ✅ one config, N samples |

Honest summary: **fastmaps alone is not a drop-in FastQC replacement — but
`fastcheck` is what closes the gap.** fastmaps replaces the composition
modules (and beats FastQC on the polyX/R1-R2 view); fastcheck adds exact
per-base quality, per-sequence quality, GC, length distributions and
duplication levels. What is still missing from the FastQC slot: per-position
k-mer content (the 12-mer spectrum is whole-read, top-20). Per-tile modules
(`--tiles`, off by default) and exact over-representation modes
(`--overrep-mode`) are no longer gaps. The replacement claim is therefore
"FastQC's modules, replaced and extended, minus per-position k-mer content —
see `docs/plans/2026-10-12_fastcheck-exactness.md`", not "all of FastQC".

## MultiQC → fasthelper

| MultiQC behaviour | FastSuite | status |
|---|---|---|
| discover per-tool outputs in a run directory | `fasthelper concatenate <dir>` walks for `*.json` and auto-detects the producing tool | ✅ |
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
