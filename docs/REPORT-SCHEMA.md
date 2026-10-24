# fastsuite.report.v1

The single machine-readable output format every FastSuite tool writes with
`--json <PATH>`. `fasthelper` reads exactly this (plus two legacy shapes) and
turns it into a multi-sample report.

Print the normative copy with:

```fish
fasthelper schema
```

## Payload

```json
{
  "schema":       "fastsuite.report.v1",
  "tool":         "fastdedup",
  "tool_version": "1.3.0",
  "command":      "fastdedup -1 R1.fq -2 R2.fq -o d1.fq -p d2.fq -j 16",
  "mode":         "paired",
  "dry_run":      false,
  "inputs":       { "r1": "R1.fq", "r2": "R2.fq" },
  "outputs":      { "r1": "d1.fq", "r2": "d2.fq" },
  "elapsed_s":    12.53,
  "metrics": {
    "reads_in":        1000000,
    "reads_out":       998000,
    "reads_duplicate": 2000,
    "duplicate_rate":  0.002,
    "hash_bits":       64
  }
}
```

## Rules

1. `metrics` is a **flat numeric namespace**. Booleans as `0`/`1`. Never nest.
   Anything structured (per-adapter counts, per-position curves) goes in a
   tool-specific block *outside* `metrics`.
2. Metric names are `snake_case`, and prefixed only when a name would be
   ambiguous across tools (`r2_polyG_frac`, not `fastmaps_r2_polyG_frac`).
3. Recommended shared keys — tools emit whatever they can compute, none are
   mandatory: `reads_in`, `reads_out`, `reads_dropped`, `bases_in`, `bases_out`,
   `bases_dropped`, `bases_trimmed`. `elapsed_s` is top-level, not a metric.
   `fastprep` emits all of these except the drop-cause split's siblings, and
   satisfies `bases_in == bases_out + bases_trimmed + bases_dropped` (a
   dropped read can still have been trimmed first). Its `reads_out` counts
   RECORDS written, not trim decisions. No other tool emits `bases_out` yet. No
   key is mandatory, and `fasthelper check --report <json>` enforces only what
   the schema states: the `schema` tag, the identifying keys (`tool`,
   `tool_version`, `command`, `metrics`), that every metric is a number, and the
   base identity above whenever `bases_out` and `bases_dropped` are both
   present. A missing recommended key is a warning, never an error.
4. Wall-clock `elapsed_s` is measured by the tool itself; `fasthelper` keeps its
   own measurement in `manifest.json` for cross-checking.
5. `command` is the exact argv, space-joined: enough to reproduce the run
   without guessing which flags produced the numbers.
6. **Record ORDER within a FASTQ/FASTA output is not part of the contract.**
   `fastdedup` keeps a stable SET across thread counts but does not promise an
   order; anything gating on byte-identity must compare sorted records. See
   `tools/fastdedup/README.md` for the measurement.

## Versioning

* **Minor** (adding a metric, adding a top-level key) → keep `v1`.
* **Major** (renaming/removing a metric, changing a type) → new schema string;
  emit the old key as an alias for one release and note the rename in the tool's
  CHANGELOG.
* Consumers must ignore unknown keys and must not require any metric to be
  present — the adapter layer in `fasthelper` is written that way.

## Producers, today

| tool | status | notes |
|---|---|---|
| `fastdedup` | native v1 | `--json report.json` |
| `fastsim` | native v1 | synthetic control: pairs, artefact/duplicate/must-drop rates |
| `fastprep` | native v1 | `--json report.json`; the raw stats object is preserved under `legacy` for one release |
| `fastmaps` | native v1 | `--json report.json` writes the v1 wrapper at PATH; the grid stays at `<output-dir>/heatmap.json`, where the legacy adapter still finds it |
| `fastcheck` | native v1 | `--json report.json`; carries `assumptions` and the profile/fit blocks |

## Legacy adapter (what fasthelper accepts and derives)

### fastprep raw stats

Drop causes are per-read and partition `reads_dropped`. `reads_dropped_length`
(added 2026-10-20) is the post-trim length floor: a read that survived trimming
but ended below `--min-len` (an absolute floor above the read length is the usual
cause). It is distinct from the trim-budget causes — the budget caps how much
trim is *allowed*, this one rejects a read that was too short to begin with.

Only reached by report files written before fastprep emitted the v1
wrapper (PLAN.md §1.1). A current `fastprep --json` file is a v1 report
and keeps its raw object under `legacy`.

| derived metric | source |
|---|---|
| `reads_in`, `reads_out`, `reads_dropped` | sum over `r1`/`r2` |
| `reads_dropped_hard`, `_mean_q`, `_adapter`, `_duplicate`, `_polyx` | per-cause counters, **summed across mates** (halving them by taking only the last mate is a real bug the adapter avoids) |
| `bases_trimmed`, `bases_trimmed_frac` | `bases_trimmed_5 + bases_trimmed_3` over both mates |
| `reads_adapter_trimmed` | `r1.adapter.reads_trimmed + r2.adapter.reads_trimmed` |
| `read_retention` | `reads_out / reads_in` |
| `pairs_out`, `pairs_dropped`, `singletons_r1/r2` | top-level, paired runs |
| `bases_corrected` | paired runs with `--pe-correct`: retained bases rewritten by the posterior consensus correction, both mates summed (0 when the flag is off). Mirrored in `legacy` with the rest of the raw `PairedStats` object |
| `margin_bases_5p`, `margin_bases_3p`, `reads_with_margin` | aggregate-only `--overtrim-margin` counters: bases the deliberate overtrim extension added at each end, and reads that received at least one (a read with no margin is not counted). Emitted in **both** `metrics` and `legacy`, mirrored, so a v1 consumer and the pre-v1 adapter read the same numbers |

There is deliberately **no per-read margin record and no PXMT field**. PXMT v3 is
unpacked at fixed offsets by the bench scripts, so inserting a field would
silently misalign every one of them; the margin is attributed from the run
configuration plus these counters (D-15, and the plan's Tier 3a guardrail).

### fastmaps `heatmap.json`

Grid is row-major with `index = y*size + x`, where `x = #(T+C)` and
`y = #(T+G)` (verified against `fastmaps/src/format.rs`). The adapter computes,
for `all`, `r1` and `r2` separately, with a corner window of `max(2, size/50)`
bases:

| metric | definition |
|---|---|
| `polyA_frac` | mass in the low-x/low-y corner (0,0) |
| `polyT_frac` | mass in the high-x/high-y corner |
| `polyG_frac` | mass in the low-x/high-y corner — **the NovaSeq 5′ polyG signature** |
| `polyC_frac` | mass in the high-x/low-y corner |
| `off_diagonal_frac` | mass with `|x − y| > size/10` — reads whose composition is nothing like the library's |
| `total_reads`, `grid_size` | scalars from the file |
| `min_seq_len`, `max_seq_len` | the length span the grid was built over (raw reads, per mate) |
| `corner_window` | the window actually applied, `max(2, size/50)`; published so the bias below is visible from the report alone |

**Length bias, stated.** The grid is sized for the longest read (`size =
max_seq_len + 1`), so a read shorter than `size − 1 − corner_window` can never
reach a corner. On a mixed-length run the corner fractions therefore describe
only the longest reads, and a single untrimmed read can move them by more than
any biological effect. Since 2026-10-20 the tool writes one stderr warning when
`max_seq_len − min_seq_len > corner_window`, and both outputs carry the three
scalars above so a consumer can detect the case without stderr. The corner
definition itself is unchanged; a length-relative plane is the open fix.

Ground-truth check of the adapter (10 reads: 5 random, 3 polyT, 2 polyG) gives
`polyT_frac = 0.30`, `polyG_frac = 0.20`, `off_diagonal_frac = 0.20` — the two
polyG reads are the only off-diagonal mass, the three polyT reads sit on the
diagonal at (150,150). Reproduce with `examples/demo.sh`.

## Assumptions (D-6)

Every run states what it had to assume instead of measuring, as data rather
than a stderr aside:

```json
"assumptions": [{"what": "quality offset", "value": "phred+33 assumed",
                 "basis": "no explicit --phred-offset and no census that excludes +64",
                 "consequence": "a +64 library would read ~31 Q high"}]
```

The array is empty when nothing was assumed. `fastcheck` populates it from its
base/quality census; `fastprep` populates it for the Phred offset (and records
`phred_offset` in `metrics`). `fasthelper` carries the block into the aggregate
as notes on the report row, so a published table cannot silently drop the
caveat.

## Aggregate report (fasthelper output)

`fasthelper concatenate` (alias `report`) / `fasthelper run` write:

| file | content |
|---|---|
| `report.json` | `{schema, generator, title, reports:[{sample, tool, tool_version, source, metrics, notes}]}` |
| `report.tsv` | wide matrix: one row per `(sample, tool)`, one column per metric, empty cell = not applicable |
| `report.md` | one table per tool |
| `report.html` | self-contained page (no scripts, no external assets), inline bars on rate metrics |
| `manifest.json` | only for `run`: run config, resolved binaries, every argv, exit code, wall time, log path, failures |

Sample naming: `fasthelper run` records each report together with the sample it
belongs to (the run already knows the name), so a report can never be
misattributed. `fasthelper concatenate` on loose files still derives the name
from the path: the file's stem, unless it is a generic report name
(`heatmap`, `fastcheck`, `fastprep`, `report`, …) — then the nearest non-generic
parent directory is used. So `runs/demo/E94/maps/heatmap.json` is attributed to
`E94`, and a QC-only run attributes `alpha/fastcheck.json` to `alpha`.

## fastcheck rich report (v1, additive)

fastcheck emits the QC metrics plus the tool-specific blocks below and three
extra per-position columns. All are additive (minor): consumers that ignore
unknown keys are unaffected, and every numeric metric flows through `fasthelper
concatenate report` with no adapter change. Module-gated blocks follow the same
absent-not-null rule: `adapter_content` exists only when `--adapter-content` ran,
`per_tile` only with `--tiles`, and a `modules` array is emitted only when
`--only`/`--no-*` gating was used.

### New `metrics` (all numeric)

`polyx3_max`, `polyx5_max`, `polyx3_frac_gt_5`, `polyx5_frac_gt_1`,
`head_g_frac` (position-1 G fraction), `head_dropout` (library G baseline minus
position-1 G fraction; >0 = 2-colour first-cycle dropout), `repeated_seq_count`,
`repeated_seq_frac`, `distinct_seqs`, `duplication_level` (1 - distinct/reads),
`kmer_top_count`, `kmer_top_frac`, `heatmap_peak`, `heatmap_peak_x`,
`heatmap_peak_y`.

### New tool-specific blocks (outside `metrics`)

| block | shape |
|---|---|
| `polyx` | `{3p: {A|C|G|T: {len: count}}, 5p: {A|C|G|T: {len: count}}}` — homopolymer run lengths |
| `repeated_sequences` | `[{seq, count, frac}]` — over-represented full sequences (>= max(10, 0.1% of reads)) |
| `kmer_spectrum` | `[{kmer, count, frac}]` — top 12-mers |
| `heatmap` | `{size, counts}` — the 256x256 grid, row-major `y*size+x`, `x = #(T+C)`, `y = #(T+G)` (same encoding as fastmaps) |
| `processed` | `null`, or a mirror of the above for the `--processed-*` reads (before/after) |

### New positions TSV columns

`g_frac` (per-position G fraction), `polyx3_frac` (fraction of reads whose
trailing run reaches this position), `polyx5_frac` (fraction whose leading run
reaches this position). `fasthelper concatenate --positions` reads the header
at runtime, so these become matrix metrics automatically (`--metric g_frac`, or
`--all-metrics`).
