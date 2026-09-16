# FastSuite

Rust tools for FASTQ quality control and preprocessing, built as one suite:
**one flag grammar, one report format, one config, one command.**

The goal is to replace the `FastQC`/`Falco` + `MultiQC` + `fastp` trio with
faster, single-binary, native-Rust components that share an interface instead
of merely coexisting in a pipeline.

| tool | replaces | what it does | status |
|---|---|---|---|
| **fastdedup** | `fastp --dedup` (and standalone `seqkit rmdup`, `prinseq++`) | exact PCR-duplicate removal by 64/128-bit fingerprint; O(N) memory, no Bloom-filter false positives | v1.3.0, shipped |
| **fastprep** | `fastp` (trimming + adapter + filters) | asymmetric quality-weighted polyX trimming (5′ **and** 3′ **and** interior), adapter trimming with catalogue auto-detection, hard trim, mean-Q filter, optional dedup | working (was `rpolyx`) |
| **fastsim** | *no equivalent* (benchmark scaffolding) | NovaSeq-like PE read simulator with full ground truth: artefact CDFs, quality curves, adapters, duplicates, must-drop class — driven by a profile measured from real data | working (was `novaseq-sim`) |
| **fastmaps** | FastQC/Falco composition plots | 2-D heatmap of reads by base composition `(x = T+C, y = T+G)`; natively exposes polyA/polyT/polyG/polyC corners and off-diagonal mass, per mate | working (was `fastmaps` in `heatmaps/`) |
| **fastcheck** | FastQC / Falco | exact one-pass QC: basic statistics, per-position base/N/quality matrices, per-read length/GC/mean-Q histograms, duplication levels, adapter content, k-mer content — no sampling, no verdicts, memory printed before the run | v0.1 working: statistics + per-position modules; duplication/adapters/k-mers next (`docs/PLAN-QC.md`) |
| **fasthelper** | MultiQC (+ workflow glue) | a toolbox: `concatenate` (per-sample reports → one table/report), `run` (the chain from one TOML config), `doctor`, `schema` (json/tsv/md/html), checks the install, owns the report schema | new, this repo |

Still outside the suite on purpose: `fastx-diff` (visual before/after
inspection) and `gzippy` (compression engine) are natural candidates for
`fastview` / a shared compression backend — see `docs/PLAN.md` phase 5.

## Layout

```
FastSuite/
├── README.md
├── docs/
│   ├── CLI-SPEC.md          the unified option grammar (the contract)
│   ├── REPORT-SCHEMA.md     fastsuite.report.v1
│   ├── RENAME-MAP.md        every old name → new name, and what changed in the code
│   ├── REPLACEMENT-MAP.md   FastQC/Falco/MultiQC/fastp → FastSuite, feature by feature
│   └── PLAN.md              phased plan to full compatibility, with estimates
├── examples/
│   ├── fastsuite.toml       a complete run config (synthetic control sample)
│   └── demo.sh              end-to-end demo on generated data
└── tools/
    ├── fastcheck/           src, tests (incl. the Python differential harness)
    ├── fastdedup/           src, tests, benchmarks, paper/
    ├── fastprep/            src, paper/
    ├── fastsim/             src, tests, profiles/, scripts/, validation/
    ├── fastmaps/            src/
    └── fasthelper/          src (new)
```

## Build

```fish
for t in fastcheck fastdedup fastprep fastsim fastmaps fasthelper
    cd tools/$t; and cargo build --release --offline; and cd ../..
end
```

Then check the install:

```fish
./tools/fasthelper/target/release/fasthelper doctor --tools-dir tools
```

## Use

One command for a whole preprocessing run:

```fish
./tools/fasthelper/target/release/fasthelper run -c examples/fastsuite.toml --tools-dir tools
```

That writes `report.{json,tsv,md,html}`, `manifest.json` (every command, exit
code and wall time) and per-step logs; intermediates are deleted once their
last consumer has read them.

Aggregate reports you already have (including from the pre-suite binaries):

```fish
# per-sample JSON reports -> one table (one row per sample, one column per metric)
fasthelper concatenate runs/ -o report/ --format html,tsv

# per-position tables (from `fastcheck --positions-tsv`) -> one positions x samples matrix
fasthelper concatenate --positions runs/ -o matrix/                  # metric mean_q by default
fasthelper concatenate --positions runs/ --metric depth --stdout     # any column of the table
fasthelper concatenate --positions runs/ --all-metrics -o matrix/     # tidy: sample,pos,metric,value
```

Every tool is still a normal Unix filter — the suite adds an interface, it does
not hide one:

```fish
fastdedup -1 R1.fq.gz -2 R2.fq.gz -o d1.fq.gz -p d2.fq.gz -j 16 --json dedup.json
fastprep  -1 d1.fq.gz -2 d2.fq.gz -o c1.fq.gz -p c2.fq.gz -j 16 --gz-level 6 --json prep.json
fastmaps  -o maps/ c1.fq.gz c2.fq.gz
```

## Conventions

* `-1/-2` inputs, `-o/-p` outputs, `-j/--threads`, `--gz-level`, `--json <PATH>`,
  `--force`, `--dry-run`, `-q/--quiet`, `-v/--verbose`, `-V/--version`.
  Full spec: `docs/CLI-SPEC.md`.
* Every tool writes the same report schema with `--json`; `fasthelper` reads it.
  Full spec: `docs/REPORT-SCHEMA.md` (also `fasthelper schema`).
* Long names win; the pre-rename spellings stay as aliases for one release, so
  existing pipelines keep working.
