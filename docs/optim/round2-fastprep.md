# fastprep — round-2 speed propositions (2026-10-21)

Evidence: `docs/optim/raw/round2/_salvage-fastprep-speed.md` (62 probe results),
`~/work/FastPrep-Bench/results/results_gz6.csv`, `bench/profile/pipeline_5m.csv`.
Tree: `30e8fb3` + `af268d4`; round 1 (FP-2..FP-8, FP-14) landed in `154a39f`.
Tiers: [M] measured (artifact named) / [E] estimated (pricing A/B named) /
[V] verified on today's source or binary. Proposals only — nothing implemented.

## Where the wall time goes

Probe (200k PE pairs, generated input, built binary; task-0 salvage):

| run | wall s | cpu s | busy cores | peak MB |
|---|---|---|---|---|
| gz6 output, base | 0.722 | 5.145 | 7.13 | 38.5 |
| gz6 + `--metadata` | 0.830 | 5.349 | 6.45 | 44.0 |
| `-j 2` | 1.340 | 3.726 | 2.78 | 34.0 |
| `-j 4` | 0.754 | 4.229 | 5.61 | 39.9 |
| `-j 8` | 0.649 | 5.004 | 7.71 | 43.0 |
| `-j 12` | 0.657 | 5.317 | 8.10 | 48.3 |
| `-j 16` | 0.651 | 5.631 | 8.66 | 44.6 |
| `-j 8 --chunk 4096` | 0.796 | 4.793 | 6.02 | 47.1 |

Saturation is at ~8 workers on this box: CPU per pair rose 12.6 % from `-j 8`
to `-j 16` (5.004 -> 5.631 cpu-s) while wall did not move. The remaining levers
are per-read, not more concurrency. Metadata costs +0.11 s wall / +0.20 cpu
at 38.4 B per read.

Against fastp (`results_gz6.csv`, 5 M pairs, gz6): fastprep **16.36 CPU-s per M
pairs vs fastp 17.64 (7 % less work per pair)** and **101 MB peak vs 1141 MB
(11x lighter)**. The wall gap is utilisation, not per-pair work. **FP-0 caveat:
the harness's `pinned` rows are broken** (fastprep rows burned 2.44 cores inside
a supposedly pinned run — optim-dossier item 1); quote CPU-seconds and busy
cores, never that wall column.

## Propositions

### F2-1 — Three-stage pipeline: parse | trim | ordered commit (re-price FP-1)
work-partitioning · [E] ~26–28 s at 5 M vs 33.5 s · behaviour-change (byte-identical FASTQ intended; order of emit preserved by the commit thread)
`pipeline.rs:445-501`: parse of chunk N+1 shares the main thread with the write
of chunk N, `bounded(2)` in front; the dedup probe, stats update and staged emit
sit on the same thread. Split into parse thread -> `ChunkBuf` ring (depth 4) ->
trim pool -> ordered commit thread (dedup probe + stats + 256 KiB staging,
FP-4). The 200k probe already reaches 7.7 busy cores at `-j 8`; the point is
decoupling parse from the ordered write so a slow write no longer stalls parse.
A/B: 5 M pairs, pinned `-j 8`, busy cores + CPU/M-pair vs `results_gz6.csv`.

### F2-2 — Output gzip level: default to a fast level and document the trade
default-behaviour · [M] (fastsim probe gz6 13.7–16.4 s vs gz1 3.1–3.7 s at 200k;
fastdedup `gz_levels.csv` gz1 18.6 s / 154 MB vs gz6 22.79 s / 4372 MB at 5 M)
· behaviour-change (flag loudly: .gz bytes change; decompressed FASTQ identical)
`--gz-level` is in the canonical grammar; the proposition is the **default** and
the documented trade (fastp-class speed at level 1/4, size +18–88 %). On every
tool that writes .gz this is the single largest wall lever.

### F2-3 — Clamp `-j` past saturation; report busy cores
default-behaviour · [V] (probe table: j8->j16 = +12.6 % CPU, +0.2 % wall) ·
byte-identical
Document the ~8-worker plateau and stop spawning past it (or say so in
`--help`); optionally print `cpu/wall` in `--json` so users see utilisation.

### F2-4 — Per-pair allocation trims in the PE path
allocation · [V] (`pe_overlap.rs:177` allocates a `Vec<u8>` via `revcomp(r2)`
per pair; `pe_overlap.rs:196` `admissible: Vec<usize>` per pair) · byte-identical
Reuse a thread-local buffer, or compare in complement space without materialising
the reverse complement. Small per pair, whole-chunk in aggregate.

### F2-5 — Metadata writer staging
I/O · [V] (38.4 B/read record; +0.11 s wall at 200k) · byte-identical to the
PXMT format
`metadata.rs` writes variable-length records; give it the same 256 KiB staging
FP-4 gave the FASTQ writer.

## Deliberately NOT proposed
* Re-checking prefilter widths (FP-3 landed; `ms_safe` is sound per `af268d4`).
  Any loosening is an undertrim/overtrim decision — see `docs/UNDERTRIM-REDUCTION.md`.
* SIMD quality scan (FP-11): the gate LUT already indexes raw bytes; a SIMD
  variant needs its own equivalence proof for marginal gain on a 151-base scan.
* FP-12 hash-set sharding: conflicts with the ordered writer (round-1 reason
  still holds).

| ID | title | tier | gate | target | expected effect |
|---|---|---|---|---|---|
| F2-1 | three-stage pipeline | E | behaviour-change (bytes intended identical) | pipeline.rs:445-501 | +0.7 busy cores, ~20 % wall at 5 M |
| F2-2 | fast default gzip level | M | behaviour-change (.gz bytes) | writer path | up to 4x wall on the output stage |
| F2-3 | clamp -j past saturation | V | byte-identical | cli/main | no wall change; honest UX |
| F2-4 | PE-path alloc trims | V | byte-identical | pe_overlap.rs:177,196 | small, free |
| F2-5 | metadata staging | V | byte-identical | metadata.rs | -0.1 s / 200k |
