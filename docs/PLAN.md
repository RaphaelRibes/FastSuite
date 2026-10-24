# FastSuite plan — from four tools to one suite

Goal: replace `FastQC`/`Falco` + `MultiQC` + `fastp` with native-Rust
components that share one CLI grammar, one report schema and one runner, and
that are demonstrably faster on the same machine and the same data.

Estimates are for one developer working with an agent, inclusive of
verification. Phase 0 must be answered before phase 1 starts; phases 1–4 are
sequential, 5–6 partially parallel.

---

## Phase 0 — Decisions to lock (30 min, no code)

| # | decision | why it blocks |
|---|---|---|
| 0.1 | **fastprep lineage.** Take `fastrpolyx` (2026-08-11, superset, +950 src lines) as the suite's fastprep, or freeze on `rpolyx-code` (2026-08-06) that `rpolyx-bench` builds via `RPOLYX_REPO`? | decides which tree every later phase edits, and whether the bench numbers must be re-run |
| 0.2 | **Names.** `fastprep` (binary) for FastPreprocessing, or the long name? `FastSuite` vs `FastSuit`? | every path, doc and report title |
| 0.3 | **Migration shims.** Ship symlinks (`rpolyx → fastprep`, `fdedup → fastdedup`), a `fasthelper install` command, or nothing? | decides whether existing cron/SLURM jobs break |
| 0.4 | **fastsim scope.** Is the simulator a shipped tool in the suite, or an internal benchmark fixture? | whether it gets the full CLI treatment or stays a fixture |

Deliverable: the four answers written into `docs/RENAME-MAP.md`, plus a
`legacy/` decision if 0.1 goes the other way.

---

## Phase 1 — Interface convergence (1–1.5 days)

Make the CLI spec (`docs/CLI-SPEC.md`) true. Everything here is small, local and
independently verifiable.

**1.1 fastprep: v1 report wrapper** (~2 h)

* `src/stats.rs`: add a `to_report_v1()` that emits the documented shape
  (`schema/tool/tool_version/command/mode/inputs/outputs/elapsed_s/metrics`),
  reusing the existing `Stats` struct; keep the raw object under
  `"legacy": {...}` for one release.
* `src/main.rs`: write the wrapper from both the single and paired paths.
* Test: `tests/report_v1.rs` — run on a 4-read fixture, assert `schema`,
  `tool`, `metrics.reads_in`, and that `metrics` has no nested value.
* Verify: `fasthelper concatenate <out.json>` no longer needs the legacy adapter
  (assert it stops adding the `legacy raw stats` note).

**1.2 fastprep: `--force`, `--dry-run`, `--quiet`** (~3 h)

* `--force`: refuse to overwrite an existing `-o/-p` unless set (today it
  clobbers silently). Errors name the file and suggest `--force`.
* `--dry-run`: read, compute every counter, write no FASTQ and no metadata; the
  JSON report is still emitted, flagged `dry_run: true` (delivered — see
  `docs/CLI-SPEC.md` §1).
* `--quiet`: suppress the progress/`eprintln!` stream, keep warnings and errors.
* Tests: output path exists → exit 1 without `--force`, exit 0 with it;
  `--dry-run` leaves the directory byte-identical (`find -newer` check).

**1.3 fastmaps: `--json`, `-j/--threads`, `--force`, `-q`** (~3 h)

* `--json <PATH>`: write `heatmap.json` where told (default stays
  `<output-dir>/heatmap.json`) and additionally emit the v1 wrapper around the
  metrics fasthelper already derives — so the corner maths stops being a
  consumer-side convention.
* `-j/--threads`: wire to the worker pool size; `0` = all cores.
* `--force` for the output directory, `-q` for the summary line.
* Tests: existing 27 tests stay green; new test asserts the v1 file's keys.

**1.4 Shared env defaults + `--threads 0`** (~1 h)

* `FASTSUITE_THREADS` read by all four binaries when `-j` is absent.
* `0` = all cores in fastprep and fastmaps (fastdedup already does it).
* Test: `FASTSUITE_THREADS=4 fastprep --json r.json` → `metrics.threads == 4`.

**1.5 fastsim: clap adoption** ~~(~2 h, optional/cosmetic)~~ **[done 2026-10-12 —
see docs/plans/2026-10-12_fastsim-cli-convergence.md]**

* clap 4.5 derive replaced the hand-rolled parser; `-h`, `--help`, `-V` and the
  error messages now match the rest of the suite. Every existing flag name is
  kept (`-o` plus the `--out1`/`--out-prefix`/legacy `--out` aliases).
* `--force`/`-f` (refuse to clobber the output prefix) and `-q`/`--quiet`
  (suppress the summary line) were added at the same time.
* `FASTSUITE_THREADS` is read when `-j` is absent (CLI-SPEC gap #6); the value
  is reported-and-ignored — fastsim is single-threaded.
* Verification: 58 tests green (45 unit + 13 integration), the 19-form parity
  gate passes on the release build, and the fastcheck→fastsim closed loop runs
  end to end on a fastcheck-built profile.
* The deployment scripts were migrated off the pre-rename `novaseq-sim`
  checkout to a `FASTSUITE_ROOT` indirection with a
  `/scratch/users/$USER/FastSuite` fallback (Task 6/7 of the same plan).

Verification for the whole phase: `cargo test --offline` in all five crates
(378 tests at the time; 852 today, see below) plus `fasthelper run -c examples/fastsuite.toml`
on the synthetic control, plus `fasthelper doctor` showing every tool at the
new version.

---

## Phase 2 — FastHelper hardening (1 day)

**2.1 Parallel samples** (~3 h). Run samples concurrently with a worker pool
(`-P/--jobs`, default 1). Per-sample step order stays strict; threads divide
across the pool. Real win: a 86-sample mandrill run is embarrassingly parallel.

**2.2 Resume and `--from <step>`** (~3 h). Skip steps whose output exists and
whose recorded argv+tool version match the manifest. Makes a 34 TB run
restartable without a workflow engine.

**2.3 Config validation** (~2 h). `fasthelper check -c config.toml`: verify input
files exist, that every requested tool is present and supports the flags the
config passes (`--help` scrape), and that `steps` names are known. Fail before
starting, not in the middle of hour three.

**2.4 `fasthelper install --prefix DIR`** (~2 h). Symlink the five binaries into
one `bin/`, optionally with legacy names (`rpolyx`, `fdedup`, `novaseq-sim`)
per phase 0.3.

**2.5 Report polish** (~2 h). Sample-level ordering, a "cross-sample" section
(worst retention, highest polyG, highest duplicate rate), and threshold flags
(`--fail-on polyG_frac>0.05`) so the report can gate a pipeline.

---

## Phase 3 — Streaming and composition (2–3 days)

**3.1 stdin/stdout (`-`)** for fastprep and fastmaps (~1 day). Needs a
reader that does not seek: fastprep's chunked pipeline already reads
sequentially, so the work is in `io.rs` (format detection from magic bytes
instead of extension) and in the output writer (no filename-derived format).
This is the single biggest blocker to replacing `fastp` in existing pipes.

**3.2 Chain-equivalence check** (~0.5 day). `fastdedup + fastprep` vs
`fastprep --dedup` must be byte-identical on the bench fixture; a test in
fasthelper (`--selftest`) that fails loudly if the two paths drift.

**3.3 Memory envelope** (~0.5 day). `--max-mem` honoured by fastprep and
fastmaps (fastdedup has `--dedup-max-mem` already); refuses to start rather
than OOM-killing a 300 M-read job.

**3.4 Compression backend** (~1 day, optional). Adopt the in-process `gzippy`
engine (already used by fastdedup) in fastprep/fastmaps so no tool shells out
to `pigz` and the `--gz-level` semantics are identical everywhere.

---

## Phase 4 — FastQC/Falco module parity (3–4 days)

> **Superseded by `docs/PLAN-QC.md`** (a full QC binary, `fastcheck`, rather
> than modules bolted onto fastmaps). The table below is the module checklist
> that plan inherits; read PLAN-QC.md for the design, memory model and phases.

Add the missing QC modules — as `fastmaps` submodules writing into the same
report schema, so `fasthelper` needs no change:

| module | metric keys | estimate |
|---|---|---|
| per-base quality | `mean_q_pos.<n>` / a compact percentile curve | 1 day |
| per-sequence quality distribution | `mean_q_hist.*`, `q20_frac`, `q30_frac` | 0.5 day |
| GC distribution | `gc_hist.*`, `gc_mean` | 0.5 day |
| adapter content per position | `adapter_pos.<n>` (reuses fastprep's catalogue) | 1 day |
| overrepresented sequences | `top_seq.<i>.count` (reuses fastdedup's fingerprint table) | 1 day |
| k-mer content (k=5) | `kmer.<kmer>.frac` (top-N only) | 0.5 day |

End of phase: `fastmaps` answers every FastQC question that matters for
short-read QC, in one pass, in one file. Position-level accumulation is the only
structural change (per-position arrays instead of a 2-D grid); the spill-to-disk
machinery already exists and is reused.

---

## Phase 5 — Suite completion (1–2 days, parallelisable)

* **fastview** (from `fastx-diff`): before/after visual inspection as a suite
  member — the natural companion to fastprep, currently a separate tool with an
  incompatible CLI.
* **Report plots** (~0.5 day): composition heatmap thumbnail, retention bars,
  per-sample QC traffic lights in `report.html`.
* **Docs/CI** (~0.5 day): a `pixi.toml` at the suite root with
  `build`/`test`/`demo` tasks; a one-command `cargo test` across all tools.

---

## Phase 6 — Benchmarks and the public claim (2 days)

1. Regenerate/refresh the fastsim profile from the mandrill data.
2. Run the head-to-head on Io (SLURM): `fastdedup`/`fastprep` vs `fastp`,
   `fastmaps` vs `FastQC`/`Falco`, `fasthelper` vs `MultiQC` — same data,
   same threads, wall clock and peak RSS, three reps.
3. Accuracy on the simulator's ground truth: retention, retention of reads that
   must be dropped, adapter false-positive rate.
4. Publish a `docs/BENCHMARKS.md` with the raw matrix and the caveats; store the
   run under `bench/` so the numbers are reproducible, not folklore.

Do not claim "replaces FastQC" before phase 4 and "replaces MultiQC" before
phase 2.5 — the replacement map is written to survive scrutiny.

---

## Current state (what is already true, verified)

* The six binaries build in this tree; `cargo test` green: fastcheck 143,
  fastdedup 135, fastprep 306, fastsim 93, fastmaps 96, fasthelper 78, plus the
  pathology bundle — 852 tests, 0 failures (2026-10-23; the per-tool manifest
  `scripts/expected-test-counts.txt` is authoritative and `pixi run ci` asserts
  it). The 436 this bullet used to quote was the 2026-10 figure.
* `fasthelper run -c examples/fastsuite.toml` executes
  `fastsim → fastdedup → fastprep → fastmaps`, records a manifest with argv,
  exit codes and wall times, and aggregates four reports into
  `report.{json,tsv,md,html}`.
* The fastmaps corner adapter is verified against a hand-built ground truth
  (3 polyT + 2 polyG + 5 random reads → `polyT_frac 0.30`, `polyG_frac 0.20`).
* Score against the CLI spec: fastdedup ~95%, fastsim ~80%, fastprep ~70%,
  fastmaps ~65% (gap list in `docs/CLI-SPEC.md` §4). The CLI-convergence
  plans under `docs/plans/2026-10-12_*.md` close those gaps task by task.

## Risks

| risk | impact | mitigation |
|---|---|---|
| lineage divergence (phase 0.1) resurfaces | bench numbers disagree with the shipped binary | pick one tree, delete the other from the suite, re-run the accuracy sweep once |
| adapter auto-detection refuses on thin data | a chain that "worked with fastp" now stops | documented + `--adapter-preset`; add a `--adapter=auto:lenient` (guess and warn) only if users ask |
| streaming (`-`) touches the reader hot path | performance regression on the main use case | byte-identical output gate + throughput matrix before/after (the existing bench harness) |
| position-level QC modules | memory blow-up on 300 M reads | reuse the existing spill-to-disk path; add a hard `--max-mem` test |
| report schema drift across tools | aggregation silently loses metrics | `fasthelper check` asserts each tool's `--json` keys against the spec |
