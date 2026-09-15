# Rename map and what changed in the code

Every rename made by the FastSuite reorganisation. **Nothing was moved out of
its original location** — the source trees in `~/hermes_scratch/` (and
`~/FastDedup`, `~/novaseq-sim`, …) are untouched. `FastSuite/tools/*` are
copies with the new names applied.

## 1. Names

| old | new | kind |
|---|---|---|
| `rpolyx` / `fastrpolyx` (crate `polyx`) | **fastprep** | directory, crate, binary |
| `rpolyx --metadata`'s reader binary `pxmt` | `fastprep-meta` | binary |
| `FastDedup` (binary `fdedup`) | **fastdedup** | directory's binary renamed (crate was already `fastdedup`) |
| `novaseq-sim` (crate `novaseq_sim`) | **fastsim** | directory, crate, lib, binary |
| `heatmaps/` (crate/bin `fastmaps`) | **fastmaps** | directory moved into the suite; crate already named `fastmaps` |
| `fastq_heatmap` (clap display name) | `fastmaps` | help/usage string |
| — | **fasthelper** | new crate + binary |

Name rationale: `fast<purpose>` for the binaries (short, unambiguous on
`$PATH`), project names in prose (`FastPreprocessing` for fastprep). The suite
directory is `FastSuite`; it prints `fastsuite.report.v1` as the schema.

## 2. Source of truth used for each copy

| suite tool | copied from | why that tree |
|---|---|---|
| fastprep | `~/hermes_scratch/fastrpolyx/` | newest lineage (2026-08-11), git-tracked, and a **superset**: it carries the parallel-BGZF / zero-copy-chunk / composition-prefilter work plus the catalogue + pigeonhole + dedup stages that the older `rpolyx-code/` snapshot (2026-08-06) has, 6,986 vs 6,037 src lines |
| fastdedup | `~/FastDedup/` | the only tree with `parallel.rs`, `writer.rs` and abundance/pre-count support (matches the `fdedup-dev` binary's flags); `~/work/FastDedup`, `~/test/FastDedup` and the `Downloads/` copies are older |
| fastsim | `~/hermes_scratch/novaseq-sim/` | sole tree |
| fastmaps | `~/hermes_scratch/heatmaps/` | sole tree |

⚠ **One decision to confirm.** `rpolyx-bench` builds its binary from a
`rpolyx-code` checkout (`RPOLYX_REPO`), while the newer development lineage is
`fastrpolyx` — the same class of two-lineage divergence the earlier three-repo
split dealt with. `FastSuite/tools/fastprep` was taken from the **newer**
lineage. If the bench must keep building `rpolyx-code`, either pin the bench to
`tools/fastprep` and re-verify the accuracy numbers, or keep `rpolyx-code` as
`legacy/` inside the suite. Not decided here; it is phase 0 of `PLAN.md`.

Also left alone deliberately: `fastrpolyx/data/` (6.4 GB of reads),
`novaseq-sim/local_work/` (15 GB), `heatmaps/*.data/.preview/.png` outputs,
every `target/`, `.pixi/` and nested `.git/`.

## 3. Code changes applied to the copies

### fastprep (from fastrpolyx)

| file | change |
|---|---|
| `Cargo.toml` | `[package] name = "polyx"` → `"fastprep"`; `[[bin]] rpolyx` → `fastprep`; `[[bin]] pxmt` → `fastprep-meta` |
| `src/main.rs` | `use polyx::…` → `use fastprep::…` (3 imports) |
| `src/cli.rs` | clap `name = "polyx"` → `"fastprep"`; runtime progress prefix `"polyx: "` → `"fastprep: "`; `--gz-level` gained `visible_alias = "compression"` |
| `src/pipeline.rs`, `src/bin/*` | warning prefixes, doc comments, `about =` string |
| `pixi.toml` | workspace name `polyx-trim` → `fastprep` |

Verified: `cargo check --all-targets`, `cargo build`, `cargo test`.

### fastdedup (from ~/FastDedup)

| file | change |
|---|---|
| `Cargo.toml` | `[[bin]] name = "fdedup"` → `"fastdedup"`; added `serde_json = "1"` |
| `src/cli.rs` | long aliases `--in1/--in2/--out1/--out2` added to `--input/--input-r2/--output/--output-r2`; `--compression` renamed to `--gz-level` (old name kept as a visible alias, `-c` unchanged); **new** `--json <PATH>` and `-q/--quiet` |
| `src/main.rs` | new end-of-run summary line `<tool>: N reads in, M out, D duplicates (x%), T threads, time`; new report writer emitting `fastsuite.report.v1` |
| `tests/*.rs` | `cargo_bin_cmd!("fdedup")` → `"fastdedup"` |
| `src/processor.rs` | temp file prefix `fdedup_abundance_` → `fastdedup_abundance_` |

### fastsim (from novaseq-sim)

| file | change |
|---|---|
| `Cargo.toml` | `name = "novaseq-sim"` → `"fastsim"`, description prefixed `FastSim:` |
| `src/*.rs`, `tests/*.rs` | `use novaseq_sim::` → `use fastsim::`; `env!("CARGO_BIN_EXE_novaseq-sim")` → `…_fastsim`; temp-file prefixes; doc/help strings |
| `src/cli.rs` | **new** `-o` (alias of `--out`, plus `--out1`/`--out-prefix`), **new** `-j/--threads` (accepted, warns that simulation is single-threaded), **new** `--json <PATH>`, **new** `-V/--version`; help text updated |
| `src/main.rs` | new v1 report writer for the synthetic control (pairs, artefact/duplicate/must-drop rates, seed) |
| `pixi.toml`, `README.md` | workspace name, title, binary path |

Left as-is on purpose: `scripts/*.sh`, `*.sbatch` still reference
`/scratch/users/ribesr/novaseq-sim` and `~/hermes_scratch/novaseq-sim` — those
paths are the *deployed* HPC checkout, not this copy. Renaming them in place
would silently break a working deployment; update them the day the checkout is
renamed.

### fastmaps (from heatmaps/)

| file | change |
|---|---|
| `src/cli.rs` | clap `name = "fastq_heatmap"` → `"fastmaps"`; `-o` short added to `--output-dir` |
| `pixi.toml` | workspace name `wta` → `fastmaps` |

### fasthelper (new)

Four subcommands: `run` (chain from `fastsuite.toml`), `report` (aggregate
existing JSON), `doctor` (install check), `schema` (print the report contract).
Zero non-`cargo` dependencies beyond `anyhow`, `clap`, `serde_json`, `toml`.

## 4. Compatibility status

| surface | break? |
|---|---|
| existing command lines (`rpolyx -1 …`, `fdedup -1 …`, `novaseq-sim --out …`) | **no** — every old spelling is still accepted, either unchanged or as an alias |
| `--compression` | accepted (alias of `--gz-level`) |
| binary names on `$PATH` | yes — `rpolyx`/`fdedup`/`novaseq-sim` no longer exist under those names; the suite ships `fastprep`/`fastdedup`/`fastsim`. Symlinks are the migration path |
| `--json` payloads | fastdedup/fastsim now emit the v1 wrapper; fastprep/fastmaps still emit their legacy shapes (fasthelper adapts both) |
| pipelines/cron calling `rpolyx` | need the symlink or the flag rename; both are one-liners |
