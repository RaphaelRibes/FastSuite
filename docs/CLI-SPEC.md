# FastSuite CLI specification (v1)

The contract every FastSuite binary follows. It exists so that one runner
(`fasthelper`) can drive every tool with the same arguments, and so a user
learns the interface once.

Guiding rule: **the canonical name is the long option; every pre-suite spelling
survives as an alias for one release.** Nothing breaks, everything converges.

## 1. The common core

| purpose | flag | value | notes |
|---|---|---|---|
| R1 / single-end input | `-1, --in1 <PATH>` | path | `-` = stdin (target, see §5) |
| R2 input | `-2, --in2 <PATH>` | path | presence selects paired mode |
| R1 / single-end output | `-o, --out1 <PATH>` | path | `-` = stdout (target) |
| R2 output | `-p, --out2 <PATH>` | path | required when `-2` is given |
| threads | `-j, --threads <N>` | int | `0` = all cores (target); `-t` is **not** threads |
| gzip level | `--gz-level <1..9>` | int | canonical; `--compression`, `-c` accepted |
| output format on a stream | `--out-format fastq\|fasta` | enum | `-o -` only; named paths use the extension (D-4) |
| output compression on a stream | `--out-compression none\|gz` | enum | `-o -` only; default `none` (D-4) |
| overwrite | `--force` / `-f` | flag | without it, refuse to clobber |
| dry run | `--dry-run` | flag | report only, write nothing |
| report | `--json <PATH>` | path | `fastsuite.report.v1` |
| quiet | `-q, --quiet` | flag | suppress the summary line |
| verbose | `-v, --verbose` | flag | per-step detail on stderr |

### `--gz-level`: the measured cost of the default

The default stays at **6** in every tool. The numbers below are the measured
trade for `--gz-level 1`, so the choice is visible rather than folklore. All at
5 M pairs, `-j 8`, one machine, three reps per arm (2026-10-21).

| tool | `--gz-level 1` vs `6` | output bytes | peak RSS |
|---|---|---|---:|
| fastprep | **2.07× faster** (3.81 s vs 7.90 s) | 1.08× (441 vs 407 MB) | unchanged |
| fastdedup | 1.16× faster (17.5 s vs 20.3 s) | **1.87×** (729 vs 390 MB) | **25× less** (152 MB vs 3 780 MB) |
| fastsim | not applicable — its bytes are the ground-truth corpus | — | — |

Two facts that are not obvious:

* **fastprep's parallel writer has no level floor**, so level 1 keeps full
  parallel compression (busy cores 11.79 → 10.42). **fastdedup's does**: its
  parallel gzip engages only at level ≥6, so level 1 falls back to
  single-threaded flate2 (busy cores 6.16 → 3.70) — and is still faster,
  because the level-1 path is cheap enough to win outright.
* **At level 1, `-j 1` beats `-j 8` for fastdedup** (12.9 s vs 17.5 s). Once
  compression is cheap the hashing is I/O-bound and extra threads only contend.

Both level-1 outputs pass `gzip -t`: the non-standard-format hazard applies
only to a `-j >1` gzippy build *below* level 6, which the code does not do.
| temp location | `--temp-dir <PATH>` | path | where spill/intermediate files go |
| strict mode | `--strict` | flag | an assumption becomes a fatal error (D-6) |
| help / version | `-h, --help`, `-V, --version` | flag | `-V` prints `<binary> <semver>` |

Rules:

0. **`-` means the stream** (D-4). FASTQ input from stdin and FASTQ/report
   output to stdout, in a pipe or on a file, with no intermediate on disk:
   `zcat in.fq.gz | fastdedup -1 - -o - | fastprep -1 - -o clean.fq.gz`.
   Compression is decided from **magic bytes**, never from a filename, so a
   stream and a mislabelled file behave identically. When stdout carries a
   report (`--json -`), the summary line moves to stderr.
   **Exception: `fastmaps` is file-only**, by design — taking *many* files and
   classifying R1/R2 by name is its job, and a stream has no name. It is the
   one tool where a pipe would remove information rather than save I/O.
1. **No positional inputs** except tools whose job is many files (`fastmaps`,
   and `fastcheck`, whose input is a whole-library QC scan).
2. Short flags are a scarce resource: `-1 -2 -o -p -j -q -v -f -V -h` are the
   suite's; anything else a tool uses for its own options (`fastdedup -t`
   threshold, `-H` hash, `-P` read length) stays tool-local.
3. `--threads 0` means "all available cores" everywhere (fastdedup already does
   this; fastprep/fastmaps treat 0 as 1 today — gap listed below).
4. Success is exit 0; a correct refusal (e.g. fastprep declining to guess an
   adapter) is a non-zero exit **with a message on stderr saying what to pass**.
5. Every tool prints its end-of-run summary as a single line starting with
   `<tool>: `, and `--quiet` suppresses it.

Environment variables (read when the flag is absent):

| variable | meaning |
|---|---|
| `FASTSUITE_TOOLS_DIR` | where the tool checkouts/binaries live (`fasthelper`) |
| `FASTSUITE_THREADS` | default thread count (target, see §4) |

## 2. Per-tool conformance

Legend: **✓** already conforms · **+** added in this reorganisation ·
**→** alias added (old spelling still works) · **gap** still to do (§4).

### fastdedup

| canonical | state | old spelling |
|---|---|---|
| `-1 --in1` | → | `--input` |
| `-2 --in2` | → | `--input-r2` |
| `-o --out1` | → | `--output` |
| `-p --out2` | → | `--output-r2` |
| `-j --threads` | ✓ | — |
| `--gz-level` | → | `--compression` / `-c` |
| `--force -f` | ✓ | — |
| `--dry-run -s` | ✓ | — |
| `--json` | + | — |
| `-q --quiet` | + | — |
| `-v --verbose` | ✓ | — |
| `--temp-dir` | ✓ | — |
| tool-local | `-t --threshold`, `-H --hash`, `-P --read-length`, `--abundance`, `--abundance-sep`, `--pre-count` | — |

### fastprep (ex-rpolyx)

| canonical | state | note |
|---|---|---|
| `-1 -2 -o -p` | ✓ | already `--in1/--in2/--out1/--out2` |
| `-j --threads` | ✓ | `1` = serial; `0` not yet accepted |
| `--gz-level` | → | added `--compression` alias |
| `--json` | ✓ | writes the legacy stats object; wrapper pending (§4) |
| `--temp-dir` | gap | uses its own chunking, no spill flag |
| `--force -f` | ✓ | refuses to clobber an existing `-o`/`-p` (and `--metadata`/`--json`/`--singletons*`), naming the file |
| `--dry-run` | ✓ | every counter computed, nothing written; the report is still emitted |
| `-q --quiet` | ✓ | suppresses the summary stream; warnings and errors still print |
| `--pe-overlap` | + | paired-end only, **default `auto`**: infer the adapter boundary from the mate pair's own overlap (R1's 3′ suffix vs the reverse-complement R2 prefix) and trim from there. Takes a mode — `never` (the mode never runs), `auto` (**the default**; runs only when an explicit adapter set is configured, i.e. `--adapter-preset` or `--adapters`; off with auto-detection or `--no-adapter`, because the safety guard needs a set to check the boundary against), `always` (runs whenever an overlap is detected, even with no adapter set — the caller explicitly opts in). The bare spelling `--pe-overlap` (no value) is accepted as an alias for `always`, so existing invocations that asked for the mode keep working. Byte-identical to the pre-flag behaviour on every pair where detection returns nothing — no overlap, no anchor, a chance anchor hit, or an ambiguous (more than one admissible) alignment — and on any run whose mode resolved to `never`. The detector anchors on a 12-base match and verifies the whole shared window. A boundary is applied **only if** the readthrough past it is adapter sequence from the configured set, starting exactly at the boundary; with no set configured the mode refuses a boundary outright, structurally, so even `--pe-overlap=always` cannot trim biology. Refuses in single-end mode. The value is unmeasured on real libraries: the 102-dataset corpus has no adapter readthrough at all (`adapter.mode: none`; 0.01–0.26% of mates carry adapter), so it is a guard, not a value test (`docs/plans/2026-10-19_004500-pe-overlap-with-without.md`). |
| tool-local | `--metadata`, `--singletons1/2`, `--adapter-preset`, `--adapters`, `--no-adapter`, `--min-len`, `--mercy`, `--phred-sat`, `--gate-shape`, `--match-leniency`, `--mismatch-strength`, `--min-run`, `--head-dropout`, `--min-adapter-overlap`, `--adapter-error-rate`, `--mean-q`, `--phred-offset`, `--auto-phred`, `--trim-5/3`, `--chunk`, `--zst-level`, `--gz-threads`, `--dedup*` | |

### fastsim (ex-novaseq-sim)

| canonical | state | note |
|---|---|---|
| `-o` | ✓ | `--out1`, `--out-prefix`, legacy `--out` all accepted (clap visible aliases) |
| `-j --threads` | + | accepted; simulator is single-threaded, warns. Defaults to `FASTSUITE_THREADS` when absent (§4 gap #6, fastsim side done) |
| `--json` | ✓ | writes the suite report for the synthetic control |
| `-V --version` | ✓ | clap-generated |
| `-h --help` | ✓ | clap-generated, 29 lines; `-f/--force` and `-q/--quiet` present |
| `--force -f` | ✓ | refuses to clobber an existing `-o` prefix (probes `.gt.tsv`, `_R1/_R2.fastq[.gz]`), naming the file |
| `-q --quiet` | ✓ | suppresses the end-of-run summary line; warnings and errors still print |
| tool-local | `--pairs`, `--profile`, `--seed`, `--read-len`, `--insert-*`, `--gc`, `--reference`, `--adapter-seq`, `--no-adapter`, `--no-gzip`, `--dup-rate`, `--dup-pool`, `--drop-rate`, `--drop-min-len` | |
| gap | `-1/-2` meaningless for a simulator (it *generates* inputs); `--dry-run` deliberately absent (a generator always writes — the report records `"dry_run": false`) | |

### fastmaps

| canonical | state | note |
|---|---|---|
| `-o` | + | short for `--output-dir` |
| positional `FILE...` | ✓ | deliberate: many inputs is the tool's job |
| `--temp-dir` | ✓ | per-row spill directory |
| `--json` | ✓ | writes the `fastsuite.report.v1` report at PATH; the grid stays in `-o/heatmap.json` |
| `-j --threads` | ✓ | worker count; `0` = all cores |
| `--force` / `-q` | ✓ | clobber guard and quiet mode |
| tool-local | `--corners`, `--preview-reads`, `--buffer-mb`, `--svg`, `--r1-pattern`, `--r2-pattern` | |

## 3. Composition rules

The chain is `fastdedup → fastprep → fastmaps`, each step a separate process,
each writing its own report. This is intentional:

* no tool embeds another (a bug in the trimmer cannot corrupt dedup);
* every step is restartable and its log/exit code is recorded by `fasthelper`;
* the shared flags make the chain trivially scriptable without a runner.

`fastprep --dedup` exists for the single-binary case; **when fasthelper runs the
chain it uses `fastdedup` + `fastprep --no-dedup` semantics instead**, so the
deduplication rate is reported by the tool that owns the fingerprint table.

## 4. Known gaps (conformance backlog)

Ordered by user-visible value:

1. ~~`--json` wrapper (`fastsuite.report.v1`) for **fastprep**.~~ **Done** — it
   emits v1 with the raw stats under `legacy`.
2. ~~`--force` + `--dry-run` in **fastprep**.~~ **Done 2026-10-12** (per-tool
   commits; `--force`/`-f`, `--dry-run`).
3. ~~`--json <PATH>` and `-j/--threads` in **fastmaps**.~~ **Done** (2026-10-21:
   `--json` writes the v1 report at PATH; the grid stays in `-o`).
4. `--threads 0 = all cores` in **fastprep** — fastprep documents `>= 1` and
   refuses 0; `fasthelper` resolves the sentinel for it.
5. ~~`--quiet` in **fastprep**/**fastmaps**.~~ **Done** (fastmaps gained `-q`).
6. `FASTSUITE_THREADS` honoured as a default by all four binaries. **fastsim
   done 2026-10-12** (read when `-j` is absent; the value is reported-and-ignored
   because the simulator is single-threaded); fastprep/fastmaps/fastdedup to do.
7. `-`/stdin-stdout support (streaming mode) — **fastdedup done** (2026-10-12:
   `-` on `-1/-2/-o/-p`, plus `--out-format`/`--out-compression`); fastprep and
   fastmaps remain, needed to compete with `fastp`'s pipe usage.
8. ~~`fastsim` adopting clap for generated help (cosmetic; hand-parser
   works).~~ **Done 2026-10-12** — clap 4.5 derive, every flag name kept (the
   `-o` aliases are clap visible aliases), errors now in suite style
   (`error: unexpected argument '--bogus' found`, exit 2). See
   `docs/plans/2026-10-12_fastsim-cli-convergence.md`.

### fastcheck-specific gaps (closed by docs/plans/2026-10-12_fastcheck-exactness.md)

* `--max-mem` admission refusal — ported from fastprep's proven model.
* `--overrep-mode none|fingerprint|second-pass|table` — the four honest
  trade-offs; default `none` (PLAN-QC §5.9).
* adapter content per position — reusing fastprep's catalogue/matcher
  semantics, copied per D-1.
* per-tile quality modules — default off, ~5 MB per 100 tiles.
* `--only <module>` gating — the parser materialises only what enabled
  modules need.
* ground-truth recovery test (PLAN-QC §7.3) — the decisive correctness gate.

## 5. Compatibility promises

* Adding a flag: minor. Removing or renaming one: major, and the old spelling
  must remain as an alias for at least one release.
* `--json` payload: see `REPORT-SCHEMA.md` for the versioning rule.
* Exit codes: `0` success — a warning (e.g. a residual-adapter trip *without*
  `--fail-on-adapter-residual`) does not change it; `1` input/run error (a
  missing or unreadable config, a config `check` rejects, or a step that
  failed); `2` "ran fine but a threshold tripped" (a `--fail-on` violation in
  `concatenate`, or a residual-adapter trip under `--fail-on-adapter-residual`
  in `run` or in the standalone `guard`; the report and manifest are still
  written first), and equally the usage errors clap reports before the run
  (unknown flag, missing argument — see §4, item 8); `>2` reserved. All of `0`,
  `1` and `2` are produced today.

### `fasthelper guard` — the adapter-residual contract

One command, one exit code, one sentence of what it proves:

```
fasthelper guard -1 trimmed_R1.fastq.gz [-2 trimmed_R2.fastq.gz] \
    --adapter-preset illumina --fail-on-adapter-residual
```

It measures adapter content in a file the user already has — the same
measurement `run --fail-on-adapter-residual` applies to a step — so "adapter
survived trimming" is a machine-checkable failure rather than a judgement call.
Exit `2` when the residual exceeds `--threshold` (default 1 %), `0` otherwise.
The one-line report names the measured percentage and the read count, e.g.
`fasthelper guard: 100.00% of 20 reads still carry adapter (threshold 1.00%) -- RESIDUAL TOO HIGH`.
Without `--fail-on-adapter-residual` a trip is a warning and the exit code is
`0`. A file that cannot be measured is an error (`1`), never a clean bill of
health (D-6).

## fastcheck

`fastcheck [FILES...] [--json PATH] [--html PATH] [--positions-tsv PATH]
[--profile PATH] [--insert-hist] [--tiles] [--processed-1 FQ] [--processed-2 FQ]`

* `--json PATH` — the `fastsuite.report.v1` QC report (metrics + rich blocks).
* `--html PATH` — a self-contained FastQC-style report (inline SVG, no JS):
  per-base content/quality/N, distributions, insert sizes, polyX run lengths
  (4 base tabs), repeated sequences + 12-mer spectrum, and the heatmap.
* `--positions-tsv PATH` — per-position matrix (base counts, quality, `g_frac`,
  `polyx3_frac`, `polyx5_frac`) for `fasthelper concatenate --positions`.
* `--profile PATH` — a fastsim profile fitted directly from this run (no Python
  fit step); `--insert-hist` makes its insert distribution measured.
* `--processed-1/--processed-2 FQ` — optional processed files; the report then
  carries a `processed` block and the HTML gains a Before/After section.
* `--insert-hist` — the paired lockstep pass that measures insert sizes.
* `--tiles` — per-tile quality (FastQC's per-tile modules), off by default: quality by
  tile x position x quality bin, parsed from 7-field Illumina headers
  (`@instrument:run:flowcell:lane:tile:x:y`), emitted as a `per_tile` block. Reads whose
  header has no tile field are counted in `reads_without_coordinates`, never dropped.

Compare two runs side by side with `fasthelper compare --before <set>
--after <set> -o out/`.
