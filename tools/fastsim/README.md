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
| `--dup-rate R` | 0.02 | PCR duplicate fraction (mandril measures 1–4%) |
| `--dup-pool N` | 131072 | ring-buffer capacity for the duplicate pool |
| `--drop-rate R` | 0.05 | fraction of pairs whose fragment is too short to survive trimming (must-drop class) |
| `--drop-min-len L` | 100 | insert ceiling for the forced-drop class (bench `min_length`) |

## Build & test

```bash
pixi install -e build
pixi run -e build test      # 49 tests (43 unit + 6 integration)
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
