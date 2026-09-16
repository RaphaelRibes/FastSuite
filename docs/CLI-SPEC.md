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
| overwrite | `--force` / `-f` | flag | without it, refuse to clobber |
| dry run | `--dry-run` | flag | report only, write nothing |
| report | `--json <PATH>` | path | `fastsuite.report.v1` |
| quiet | `-q, --quiet` | flag | suppress the summary line |
| verbose | `-v, --verbose` | flag | per-step detail on stderr |
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
1. **No positional inputs** except tools whose job is many files (`fastmaps`).
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
| `--force` | gap | currently overwrites silently |
| `--dry-run` | gap | no stats-only mode |
| `-q --quiet` | gap | everything already goes to stderr as progress |
| tool-local | `--metadata`, `--singletons1/2`, `--adapter-preset`, `--adapters`, `--no-adapter`, `--min-len`, `--mercy`, `--phred-sat`, `--gate-shape`, `--match-leniency`, `--mismatch-strength`, `--min-run`, `--head-dropout`, `--min-adapter-overlap`, `--adapter-error-rate`, `--mean-q`, `--phred-offset`, `--auto-phred`, `--trim-5/3`, `--chunk`, `--zst-level`, `--gz-threads`, `--dedup*` | |

### fastsim (ex-novaseq-sim)

| canonical | state | note |
|---|---|---|
| `-o` | + | `--out1`, `--out-prefix`, legacy `--out` all accepted |
| `-j --threads` | + | accepted; simulator is single-threaded, warns |
| `--json` | + | writes the suite report for the synthetic control |
| `-V --version` | + | was missing entirely |
| tool-local | `--pairs`, `--profile`, `--seed`, `--read-len`, `--insert-*`, `--gc`, `--reference`, `--adapter-seq`, `--no-adapter`, `--no-gzip`, `--dup-rate`, `--dup-pool`, `--drop-rate`, `--drop-min-len` | |
| gap | help text is hand-rolled; `-1/-2` meaningless for a simulator (it *generates* inputs) | |

### fastmaps

| canonical | state | note |
|---|---|---|
| `-o` | + | short for `--output-dir` |
| positional `FILE...` | ✓ | deliberate: many inputs is the tool's job |
| `--temp-dir` | ✓ | per-row spill directory |
| `--json` | gap | writes `heatmap.json` into `-o`; no path control |
| `-j --threads` | gap | worker count not exposed |
| `--force` / `-q` | gap | — |
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

1. `--json` wrapper (`fastsuite.report.v1`) for **fastprep** — today it emits the
   legacy stats object; `fasthelper` adapts it, but new consumers should not
   have to.
2. `--force` + `--dry-run` in **fastprep**.
3. `--json <PATH>` and `-j/--threads` in **fastmaps**.
4. `--threads 0 = all cores` in **fastprep** and **fastmaps**.
5. `--quiet` in **fastprep**/**fastmaps**.
6. `FASTSUITE_THREADS` honoured as a default by all four binaries.
7. `-`/stdin-stdout support (streaming mode) — needed to compete with
   `fastp`'s pipe usage; blocked on the reader architecture in fastprep.
8. `fastsim` adopting clap for generated help (cosmetic; hand-parser works).

## 5. Compatibility promises

* Adding a flag: minor. Removing or renaming one: major, and the old spelling
  must remain as an alias for at least one release.
* `--json` payload: see `REPORT-SCHEMA.md` for the versioning rule.
* Exit codes: `0` success, `1` input/usage error, `2` (target) "ran fine but a
  threshold tripped", `>2` reserved. Today only 0/1 are used.
