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
| `dedup_threads.csv` | `fastdedup` at 5 M pairs, `-j` 1/2/4/8 (post-fix) |
| `gz_levels.csv` | the same run at `--gz-level` 1, 6 and 9 (compression vs memory) |

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

## The findings

1. **`fastsim` is ~80% of the chain's wall time** and single-threaded by
   construction, so `-j` does nothing for it. The analysis stages total 15 s per
   million pairs; the simulator alone is 63 s.
2. **`fastdedup`'s hashing is ~10 s and ~100 MB at 5 M pairs, at any thread count.**
   With plain (uncompressed) output every `-j` in 1/2/4/8 lands in 9.8–10.8 s and
   77–126 MB, so the hashing pipeline itself does not scale with threads because it
   is I/O-bound on the read side.
3. **The gzip writer is what costs memory, and it is a trade, not a leak.**
   `--gz-level 6` (gzippy, multi-block) peaks at 4.37 GB for 888 MB of output;
   `--gz-level 1` (single-threaded flate2) peaks at 153 MB for 1 681 MB of output and
   is *faster* on this box (18.3 s against 29.5 s). The level-6 footprint scales with
   input, so a 100 M-pair gzipped run needs ~80 GB.
4. **A second, independent contributor was an unbounded hand-off queue** between the
   hashers and the writer: when compression is the bottleneck the writer blocks, and
   that queue absorbed the whole stream (4.0 GB at `-j 2` and `-j 4`, and those two
   settings were slower than `-j 1`). It is now bounded (`53d6686`), and the same runs
   measure 101 MB and 105 MB, with `-j 4` 22% *faster* than before.
