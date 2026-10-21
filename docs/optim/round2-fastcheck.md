# fastcheck — round-2 speed propositions (2026-10-21)

Evidence: `docs/optim/raw/round2/_salvage-fastcheck.md` (52 probe results),
`bench/profile/pipeline_5m.csv` (5 M pairs, `-j 8`: wall 10.304 s, cpu 21.826 s,
**peak 785.0 MB** — vs fastmaps' 22.7 MB on the same input).

### FC2-1 — Parallelise the single-input path (re-price FC-4)
work-partitioning · [V] (`main.rs` takes the serial branch when
`files.len() == 1 || threads <= 1`; `reader::open` consults `gz_threads` only in
its BGZF branch, so plain FASTQ/ordinary .gz get no inflate pool — optim-dossier
item 3, unchanged) · behaviour-change only if float reductions reorder (they
must not: deterministic merge order)
The 2x500k pair at `-j 8` runs **1.93 busy cores** (8.159 cpu / 4.219 wall).
Chunked single-file parallelism (split on record boundaries, merge in input
order) is the missing axis for the most common invocation. A/B: same 1M-read
file, `-j 1` vs `-j 8`; cpu/M-read must stay flat and every report field
byte-identical.

### FC2-2 — Shrink the 785 MB peak
allocation · [V] (accum.rs `BASE_POSITIONS` 512-wide after the fitter fix; RSS
probe: a single 32 MB allocation dominates VmHWM — 35.3 MB vs 4.5 MB for 1 MB) ·
byte-identical reports
Compact the duplication table (fixed-width key + 32-bit counter, open addressing
instead of HashMap) and fuse the per-base censuses into ONE scan (composition +
quality + polyX + k-mer rolling register — the fusion pattern already used for
`--only kmer`). A/B: peak RSS at 5 M pairs and `--only` arm timings.

### FC2-3 — Profile fitter width from observed read length
algorithmic · [E] · behaviour-change (flag loudly: fitted values shift within
tolerance)
512-wide rows on short libraries are mostly padding; derive the width from the
observed read length (clamped), keep 512 only where the data reaches it.

### FC2-4 — Fixed startup cost on tiny inputs
I/O · [V] (attributable rows only — the seven-variant table's labels were
truncated in salvage and is NOT quoted; the `--only` arms are: `--only dup` on
1M reads 0.437 s / 103.9 MB, on 100k 0.039 s / 9.9 MB) · byte-identical
The full pass costs far more than any single census on the same input — the
evidence that independent censuses should fuse (FC2-2). Startup itself is not
the bottleneck; keep this as a one-look item.

| ID | title | tier | gate | target | expected effect |
|---|---|---|---|---|---|
| FC2-1 | single-file parallelism | V | byte-identical (deterministic merge) | main.rs serial branch | up to ~4-6x wall single-file |
| FC2-2 | 785 MB peak: compact dup table + fused scan | V | byte-identical | accum/dup tables | large RSS cut |
| FC2-3 | fitter width from data | E | behaviour-change | profile fitter | memory + minor wall |
| FC2-4 | tiny-input cost | V | byte-identical | startup | one-look |
