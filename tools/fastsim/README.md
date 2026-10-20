# FastSim (formerly novaseq-sim)

NovaSeq-like paired-end read simulator with **complete ground truth**, driven
by a measured profile (built from the 82-billion-read Mandrillus E-dataset
polyG profile).

The goal: simulated reads that carry maximal information — every base is
knowable as *biology* or *artefact* — so we can (a) understand what is
biology and what is not in two-colour chemistry data, and (b) score trimming
tools (rpolyx, fastp, cutadapt, ...) **exactly** against the correct answer,
not relative to another tool.

## Output

```
PREFIX_R1.fastq[.gz]   simulated R1 reads (150 bp, NovaSeq-style headers)
PREFIX_R2.fastq[.gz]   simulated R2 reads
PREFIX.gt.tsv          ground truth:  pair_id  insert  r1_h  r1_t  r2_h  r2_t  full_flag  is_dup
```

`is_dup` marks PCR duplicates: a re-sequencing of a previously emitted
fragment with fresh quality, errors and artefacts (per-cluster dropout),
matching the measured 1–4% duplication rate in real NovaSeq data.

**Reads are always exactly `read_len` bases.** The simulator models
sequencer output, not fragments: a short fragment still produces 150 bp
reads (adapter readthrough fills the read). The `--drop-rate` class is a
fragment too short to survive trimming — its reads are still 150 bp, but
the GT keep window `[h, min(I, L−t))` is below `min_length` on both mates,
so a correct trimmer must drop the pair.

The correct trimmed window per mate is `[h, min(I, L−t))`:
`[0, h)` = 5′ G-run (artefact), `[h, min(I, L−t))` = biology (keep),
`[min(I, L−t), L−t)` = adapter readthrough (artefact, only when `I < L−t`),
`[L−t, L)` = 3′ G-run (artefact). `full_flag=1` means the pair is entirely
artefact — the correct answer is to **drop** it.

## What it reproduces (measured from 82.3B reads)

| property | R1 | R2 |
|---|---|---|
| 3′ polyG tail rate (≥1 G) | 35.4% | 34.5% |
| mean tail length | 4.01 bp | 3.57 bp |
| tail p90 / p99 | 16 / 44 bp | 13 / 40 bp |
| 5′ head rate (≥1 G) | 45.0% | 45.2% |
| mean head length | 0.57 bp | 1.02 bp |
| quality inside tail run | 26.5 | 37.1 |
| quality outside tail run | 38.2 | 38.8 |
| full-length spike (read = all G) | ~0% | 0.22% |

Key effects reproduced:
- **Correlated R1×R2 tail lengths** (joint distribution, not independent)
- **R1-specific quality collapse** (−11.4 Phred inside tails vs −1.6 for R2)
- **Per-position quality curves** (R1/R2 separate) + **per-read variance**
- **Adapter readthrough** when insert < read length (Nextera by default)
- **Ambiguous bases (`N`)** from the profile's measured no-call / padding
  model (`n_run`: read rate, run-length histogram, edge shares). Absent from a
  profile → the simulator emits none (the pre-N behaviour, byte for byte).

## The profile's N model

```
"n_run": { "read_rate": 0.024, "len": {"1": 9000, "150": 400},
           "leading": 0.05, "trailing": 0.60 }
```

`read_rate` is the measured share of reads carrying an `N`, `len` the measured
run-length histogram, `leading`/`trailing` the share of runs touching the
read's first / last base (the rest interior). All four are measured by
fastcheck's fitter (merged over both mates, like `dup_rate`); the standalone
measurement lives in `bench/scripts/n_stats.py`. `r1_n`/`r2_n` in the GT say
how many `N` the overlay put in each mate.

## Usage

```bash
fastsim --pairs 20_000_000 --profile profiles/mandril_E.json \
            --out /scratch/.../sim --seed 7
```

| flag | default | meaning |
|---|---|---|
| `--pairs N` | required | number of read pairs |
| `--profile FILE` | required | measured profile JSON |
| `--out PREFIX` | required | output prefix |
| `--seed S` | 42 | RNG seed (deterministic output) |
| `--read-len L` | 150 | read length (must match profile) |
| `--insert-mean M` | 300 | insert length mean (bp) |
| `--insert-sd S` | 50 | insert length sd |
| `--insert-min M` | min(150, read_len) | insert length floor |
| `--gc G` | 0.49 | GC fraction for random DNA |
| `--reference FILE` | — | sample inserts from a FASTA instead of random DNA |
| `--adapter-seq SEQ` | profile | override adapter sequence |
| `--no-adapter` | — | disable adapter readthrough |
| `--no-gzip` | — | write uncompressed FASTQ |
| `--gz-level L` | 6 | DEFLATE level for the gzip outputs, 1..=9 (see below) |
| `--dup-rate R` | 0.02 | PCR duplicate fraction (mandril measures 1–4%) |
| `--dup-pool N` | 131072 | ring-buffer capacity for the duplicate pool |
| `--drop-rate R` | 0.05 | fraction of pairs whose fragment is too short to survive trimming (must-drop class) |
| `--drop-min-len L` | 100 | insert ceiling for the forced-drop class (bench `min_length`) |
| `--stats` | — | add the instrument counters (RNG words drawn, substitutions, mean Q, gzip bytes) to `--json` |

### `--gz-level`: the compression trade

The level changes only the DEFLATE encoding — the decompressed FASTQ is
identical at every level, so changing it cannot invalidate a campaign. Deflate
is the bulk of fastsim's wall time at level 6, and FASTQ quality strings are
near-incompressible noise, so the search mostly fails. Measured on this tool's
own output, 200 000 pairs, `taskset -c 0` — in-process backend first, then the
review's system-`gzip` CLI sweep on the same output:

| level | backend | input rate | size (both mates) | |
|---|---|---|---|---|
| 1 | flate2/miniz_oxide | ~200 MB/s | 52.63 MB (37.2 % of plain) | fastest, **+35.8 % size** |
| 6 | flate2/miniz_oxide | ~11 MB/s | 38.75 MB (27.4 %) | **default** — 0.1.0's bytes |
| 1 | system `gzip` proxy | 109 MB/s | 23.04 MB/mate (32.6 %) | the review's sweep |
| 4 | system `gzip` proxy | 59.7 MB/s | 20.87 MB/mate (29.5 %) | balanced |
| 6 | system `gzip` proxy | 9.2 MB/s | 19.54 MB/mate (27.7 %) | 0.1.0's level |

The two backends disagree and the difference matters: miniz's fast mode packs
worse than C zlib's, so **size a campaign with the +35.8 % figure**, not the
proxy's +17.9 %. End to end, 200 000 pairs go from ~14.8 s to ~2.3 s at level 1
(~6x), with the generator's ~1.8 s as the floor.

The default stays 6 so existing corpora keep their exact `.gz`; `--gz-level 1`
is the end-to-end win, `--gz-level 4` the balanced choice. `--json` records the
level it used (`metrics.gz_level`) next to the files it shaped.

## Build & test

```bash
pixi install -e build
pixi run -e build test      # 89 tests (65 unit + 24 integration)
pixi run -e build build     # release binary at target/release/fastsim
```

## Closed-loop validation

```bash
pixi run -e analysis build-profile   # 700 profiler JSONs -> profiles/mandril_E.json
# simulate, then re-profile with polyg-profiler and compare:
python3 scripts/validate_sim.py profiles/mandril_E.json \
    sim_R1.fastq.gz sim_R2.fastq.gz sim.gt.tsv validation/closed_loop
```

Acceptance: tail rate within 0.5%, tail-length KS D small, per-position Q
within ±1 Phred, collapse within ±1.

## Project layout

```
src/profile.rs    profile JSON -> typed CDFs (heads, tails, joint, quality)
src/randutil.rs   CDF + joint sampling (seeded ChaCha8)
src/biology.rs    random DNA (GC%) or reference-FASTA inserts
src/artefact.rs   heads, correlated tails, full-length class, adapters
src/quality.rs    position curve + per-read baseline + collapse + errors
src/fastq.rs      FASTQ writer (plain/gzip, NovaSeq headers)
src/truth.rs      ground-truth TSV writer
scripts/          build_profile.py, validate_sim.py, plot_heads_tails.py
profiles/         mandril_E.json (measured), test_small.json
```

See `docs/ARCHITECTURE.md`, `docs/ALGORITHMS.md`, `docs/CLI.md`.
