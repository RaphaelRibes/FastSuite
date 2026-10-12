# FastSuite — Work In Progress

> **This is work in progress: nothing works yet.** The tools build, the tests
> pass and the demo chain runs, but the suite is pre-release, the CLI is still
> converging (see `docs/CLI-SPEC.md` §4) and every tool is WIP. Do not use it
> in production. The goal and the design are settled; the implementations are
> not.

## What this is

Rust tools for FASTQ quality control and preprocessing, built as one suite:
**one flag grammar, one report format, one config, one command.**

The goal is to replace the `FastQC`/`Falco` + `MultiQC` + `fastp` trio with
faster, single-binary, native-Rust components that share an interface instead
of merely coexisting in a pipeline. The suite is a **naming scheme, not an
architecture** (decision D-1): the tools are independent binaries with no
shared crate, kept consistent by tests rather than by a compiler.

| tool | replaces | what it does | maturity |
|---|---|---|---|
| **fastdedup** | `fastp --dedup` as a separate step (and `seqkit rmdup`, `prinseq++`) | the same exact fingerprint dedup as `fastprep --dedup`, as a standalone pipeable binary — use it when dedup is its own pipeline step | v1.3.0, WIP (published: bioRxiv 2026.04.29.721745) |
| **fastprep** | `fastp` (trimming + adapter + filters) | asymmetric quality-weighted polyX trimming (5′ **and** 3′ **and** interior), adapter trimming with catalogue auto-detection, hard trim, mean-Q filter, **native deduplication** (`--dedup`, on by default — the same exact fingerprint model as fastdedup) | v0.1.0, WIP |
| **fastsim** | *no equivalent* (benchmark scaffolding) | NovaSeq-like PE read simulator with full ground truth: artefact CDFs, quality curves, adapters, duplicates, must-drop class — driven by a profile measured from real data or fitted by `fastcheck` | v0.1.0, WIP |
| **fastmaps** | FastQC/Falco composition plots | 2-D heatmap of reads by base composition `(x = T+C, y = T+G)`; natively exposes polyA/polyT/polyG/polyC corners and off-diagonal mass, per mate | v0.1.0, WIP |
| **fastcheck** | FastQC / Falco | exact one-pass QC: basic statistics, per-position base/N/quality matrices, per-read length/GC/mean-Q histograms, duplication levels, polyX runs, bounded repeated-sequence and k-mer counters, insert-size histogram — no sampling, no verdicts, memory printed before the run. `--profile` writes a fastsim profile directly | v0.1.0, WIP |
| **fasthelper** | MultiQC (+ workflow glue) | a toolbox: `run` (the chain from one TOML config, with a provenance manifest), `concatenate` (per-sample reports → one table/report), `compare` (raw vs processed side by side), `doctor`, `schema` — owns the report schema | v0.1.0, WIP |

Release history and per-tool maturity: `docs/RELEASES.md`.
What the suite replaces today, feature by feature, honestly: `docs/REPLACEMENT-MAP.md`.

## Quick start

```fish
# build every tool (the registry cache has every dep, offline works)
for t in fastcheck fastdedup fastprep fastsim fastmaps fasthelper
    cd tools/$t; and cargo build --release --offline; and cd ../..
end

# check the install
./tools/fasthelper/target/release/fasthelper doctor --tools-dir tools

# one command, the whole preprocessing chain (synthetic control sample)
bash examples/demo.sh 20000
```

Then use it:

```fish
# a whole run from one config: fastsim -> fastdedup -> fastprep -> fastmaps
# (fastprep dedups natively too: drop the fastdedup step for the one-pass chain)
./tools/fasthelper/target/release/fasthelper run -c examples/fastsuite.toml --tools-dir tools

# the closed loop: measure a dataset, then simulate it back
fastcheck --insert-hist --profile profile.json R1.fq.gz R2.fq.gz
fastsim --profile profile.json --pairs 20000 -o sim/

# aggregate reports you already have (the MultiQC slot)
fasthelper concatenate runs/ -o report/ --format html,tsv
fasthelper concatenate --positions runs/ -o matrix/            # positions x samples matrix
fasthelper compare --before raw/ --after processed/ -o cmp/    # before/after side by side
```

Every tool is still a normal Unix filter:

```fish
# dedup as its own pipeable step (same fingerprint model as fastprep's)
fastdedup -1 R1.fq.gz -2 R2.fq.gz -o d1.fq.gz -p d2.fq.gz -j 16 --json dedup.json
# ...or let fastprep do it in the same pass (--dedup is ON by default)
fastprep  -1 R1.fq.gz -2 R2.fq.gz -o c1.fq.gz -p c2.fq.gz -j 16 --gz-level 6 --json prep.json
fastcheck R1.fq.gz R2.fq.gz --json qc.json --html qc.html --positions-tsv pos.tsv
fastmaps  -o maps/ c1.fq.gz c2.fq.gz
```

## Conventions

* `-1/-2` inputs, `-o/-p` outputs, `-j/--threads`, `--gz-level`, `--json <PATH>`,
  `--force`, `--dry-run`, `-q/--quiet`, `-v/--verbose`, `-V/--version`.
  Full spec: `docs/CLI-SPEC.md`; the remaining gaps and the plans that close
  them: `docs/CLI-SPEC.md` §4 and `docs/plans/`.
* Every tool writes the same report schema with `--json`; `fasthelper` reads
  it. Spec: `docs/REPORT-SCHEMA.md` (also `fasthelper schema`).
* Long names win; the pre-rename spellings stay as aliases for one release.
* **Exact numbers, no verdicts.** Counts are integers read off the whole
  file. No sampling presented as measurement, no PASS/FAIL from magic
  constants — the decision is yours (`fasthelper` supports `--fail-on`).
* **Never ship any data** (D-11): every public demo uses generated data
  (`fastsim`) or hand-written fixtures. Real profiles live outside the repo.

## Layout

```
FastSuite/
├── README.md
├── LICENSE                MIT — Copyright 2026 Raphaël Ribes
├── docs/
│   ├── CLI-SPEC.md        the unified option grammar (the contract)
│   ├── REPORT-SCHEMA.md   fastsuite.report.v1
│   ├── RENAME-MAP.md      every old name → new name, and what changed
│   ├── DECISIONS.md       the decision log (D-1 … D-12)
│   ├── REPLACEMENT-MAP.md FastQC/Falco/MultiQC/fastp → FastSuite, feature by feature
│   ├── PLAN.md            phased plan to full compatibility
│   ├── PLAN-QC.md         fastcheck module design
│   ├── PLAN-ML.md         the parameter-recommendation research plan (design)
│   ├── PLAN-ML-STEPS.md   the runnable breakdown, with acceptance tests (historical — Io-era)
│   ├── PLAN-ML-LOCAL.md   the current ML execution plan: desktop-only benchmarks
│   ├── DEPLOY-IO.md       running the campaign on the Io cluster (historical — access ended 2026-10-12)
│   ├── RELEASES.md        release train + per-tool maturity
│   └── plans/             per-tool implementation plans (bite-sized TDD)
├── examples/              fastsuite.toml + demo.sh (synthetic control)
├── deploy/io/             historical SLURM scripts (Io access ended 2026-10-12)
└── tools/
    ├── fastcheck/         src, tests, CITATION.cff (incl. the Python differential harness)
    ├── fastdedup/         src, tests, benchmarks, paper/ (published)
    ├── fastprep/          src, tests, paper/
    ├── fastsim/           src, tests, profiles/, scripts/
    ├── fastmaps/          src/
    └── fasthelper/        src, tests
```

## Testing

```fish
for t in fastcheck fastdedup fastprep fastsim fastmaps fasthelper
    cd tools/$t; and cargo test --offline; and cd ../..
end
# 345 tests, 0 failures (as of the 2026-10-12 status)
```

## Citing

Each tool prints a one-line citation request at the end of every run and is
independently citable via its `CITATION.cff` (GitHub renders these
automatically). **fastdedup** and **fastprep** belong to the Mandrillus
Project and carry the ISDM/ISEM and Céline Mandier / Alice Baniel credits;
fastdedup is published as [Ribes, Mandier & Baniel, bioRxiv 2026](https://www.biorxiv.org/content/10.64898/2026.04.29.721745v1).
Everything else is [Raphaël Ribes](https://www.raphaelrib.es)'s own work.

The license stays **MIT** (bioconda requires an OSI-approved license) with a
citation request appended as etiquette, not a condition — see `LICENSE`.

## Out of scope (deliberately)

PE overlap/merge, UMI handling, alignment-based QC, workflow-engine
semantics (arrays, resume — that is Snakemake/Nextflow's job), long reads
for now. Each with its reason: `docs/DECISIONS.md` D-12.

## License

MIT — see [LICENSE](LICENSE). `tools/fastprep/LICENSE` carries the
Mandrillus Project attribution. Every LICENSE carries the same citation
request: asked for, not enforced — that is what keeps the tools OSI-compatible
for bioconda (D-3) while making the ask visible on every run.
