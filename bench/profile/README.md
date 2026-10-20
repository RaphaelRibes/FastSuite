# Whole-suite pipeline profile

The measurements behind "FastSuite — performance across the six tools"
(`~/.hermes/cache/scratch/fsreview/performance.pdf`, 20 Oct 2026).

Unlike the per-tool campaigns in `bench/`, these numbers are **same machine, same
input, same I/O mode**, so the stages are comparable to each other rather than to
their own papers.

## What was measured

`profile_pipeline.py` generates one library with `fastsim` and carries it through
`fastdedup → fastprep → fastmaps → fastcheck`, each stage reading the previous
stage's output. `-j 8` everywhere, gzip level 6 in and out.

`dedup_threads.py` re-runs `fastdedup` on the same 5 M-pair input with only `-j`
changed, which is where the memory surprise is.

Peak RSS is the reaped child's `ru_maxrss` from `wait4` — not a `/proc` poll, which
under-reports anything that finishes between polls (this suite has been bitten by
that once already: see the FastHelper benchmark).

## Results

| file | what it holds |
|---|---|
| `pipeline_1m.csv` | five stages at 1 M pairs |
| `pipeline_5m.csv` | five stages at 5 M pairs |
| `dedup_threads.csv` | `fastdedup` at 5 M pairs, `-j` 1/2/4/8 |

Machine: AMD Ryzen 7 7730U, 16 threads, 15 GB RAM, NVMe. Build `595ef4f`.

## Reproducing

```bash
/usr/bin/python3 bench/profile/profile_pipeline.py 1000000 8 /tmp/pipe1m
/usr/bin/python3 bench/profile/profile_pipeline.py 5000000 8 /tmp/pipe5m
/usr/bin/python3 bench/profile/dedup_threads.py 1 2 4 8
```

Each needs the release binaries built (`pixi run build`) and roughly 0.5 GB of input
per million pairs; `dedup_threads.py` reuses `pipe5m/sim`, so run the 5 M profile
first.

## The two findings

1. **`fastdedup`'s memory is a function of input size, not thread count, once
   `-j > 1`.** At `-j 1` it is 15.7 B/pair (the paper's model, reproduced). Every
   `-j > 1` carries a ~4 GB floor here (78.7 MB → 4 025 MB from one thread to two),
   and 2 and 4 threads are *slower* than 1 while using 51x the memory.
2. **`fastsim` is ~80% of the chain's wall time** and single-threaded by
   construction, so `-j` does nothing for it. The analysis stages total 15 s per
   million pairs; the simulator alone is 63 s.
