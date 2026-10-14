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
   `bases_trimmed`. `elapsed_s` is top-level, not a metric.
4. Wall-clock `elapsed_s` is measured by the tool itself; `fasthelper` keeps its
   own measurement in `manifest.json` for cross-checking.
5. `command` is the exact argv, space-joined: enough to reproduce the run
   without guessing which flags produced the numbers.

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
| `fastmaps` | **legacy** | emits `heatmap.json` (grid + `counts[_r1/_r2]`); adapter derives the corner fractions |

## Legacy adapter (what fasthelper accepts and derives)

### fastprep raw stats

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

Ground-truth check of the adapter (10 reads: 5 random, 3 polyT, 2 polyG) gives
`polyT_frac = 0.30`, `polyG_frac = 0.20`, `off_diagonal_frac = 0.20` — the two
polyG reads are the only off-diagonal mass, the three polyT reads sit on the
diagonal at (150,150). Reproduce with `examples/demo.sh`.

## Aggregate report (fasthelper output)

`fasthelper concatenate` (alias `report`) / `fasthelper run` write:

| file | content |
|---|---|
| `report.json` | `{schema, generator, title, reports:[{sample, tool, tool_version, source, metrics, notes}]}` |
| `report.tsv` | wide matrix: one row per `(sample, tool)`, one column per metric, empty cell = not applicable |
| `report.md` | one table per tool |
| `report.html` | self-contained page (no scripts, no external assets), inline bars on rate metrics |
| `manifest.json` | only for `run`: run config, resolved binaries, every argv, exit code, wall time, log path, failures |

Sample naming: a report file's stem is the sample name, unless the stem is a
generic report name (`heatmap`, `fastprep`, `report`, …) — then the nearest
non-generic parent directory is used. So
`runs/demo/E94/maps/heatmap.json` is correctly attributed to `E94`.

## fastcheck rich report (v1, additive)

fastcheck emits the QC metrics plus four tool-specific blocks and three extra
per-position columns. All are additive (minor): consumers that ignore unknown
keys are unaffected, and every numeric metric flows through `fasthelper
concatenate report` with no adapter change.

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
