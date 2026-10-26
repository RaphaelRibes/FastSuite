# What the second-wave fixes cost, and what the harness found on the way

Measured on the desktop (`cachyos-x8664`) on 2026-10-26, after the wave, because
two of the fixes sit in code paths the closed loop does not time at all:

* **fastdedup** now re-orders worker results by input ordinal before deciding
  which copy of a duplicated record survives (`82c46af`). That is a `BTreeMap`
  insert/remove per record in the writer.
* **fastprep** now measures a budget-capped run's extent past the cap
  (`d6a866e`). That is extra scan work inside the same read.

Baseline = `4f65573` (the published 1.3.0 tree), built in a detached worktree at
`/mnt/Fast/cost-bench/base`. Every number below is the tool's **own** elapsed
time from its summary line, not a shell clock.

## The reorder costs nothing measurable

Matched `-j`, base vs head, five interleaved reps on two corpora
(`scale3.tsv`, elapsed from the summary line):

| corpus | j=1 | j=2 | j=4 | j=8 | j=16 |
|---|---|---|---|---|---|
| 2.4 M reads (`corpus6`) | +4.6 % | +7.5 % | +18.4 % | −30.8 % | −1.1 % |
| 4.8 M reads (`corpus24`) | +1.9 % | +8.9 % | +18.0 % | −44.3 % | +1.5 % |

The spread between binaries at fixed `-j` (−44 % to +18 %) is larger than any
plausible cost of a map insert, and the sign is not stable: the runs at
`j = 8` are *faster* for the fixed binary, which cannot be a property of the
change. The only cells where the comparison is stable are `j = 1` (no scheduler
noise, and the reorder receives records already in order) — +4.6 %, +1.9 %, and
−2.2 % on a third corpus. **Conclusion: the ordering buffer is inside the noise
floor of this machine at `j = 1`, and unmeasurable above it.**

## The harness found a much larger defect: `-j` makes fastdedup slower

The same measurements, speedup against `-j 1` on one binary:

| corpus | binary | j=1 | j=2 | j=4 | j=8 | j=16 |
|---|---|---|---|---|---|---|
| 1.2 M reads, 2 mates | head | 1.00× | 0.44× | 0.44× | 0.40× | **0.09×** |
| 4.8 M reads, 2 mates | head | 1.00× | 0.41× | 0.44× | 0.42× | **0.09×** |
| 200 k reads, 1 dataset | head | 1.00× | 0.41× | 0.37× | 0.31× | **0.12×** |
| 1.6 M reads, gz input | head | 1.00× | 0.53× | 0.55× | 0.57× | **0.59×** |

The published 1.3.0 binary shows the same curve (`scale3.tsv`, tag `base`: 1.00×
/ 0.43× / 0.52× / 0.23× / 0.09×), so **this is not caused by the wave** — but it
contradicts the tool's purpose and the qualitative advice in its own paper
("one worker thread per physical core", §5.4 of
`tools/fastdedup/paper/paper.typ`).

Measured details (`scale3.tsv`, `scale4.tsv`, `scale5.tsv`; 3–5 reps, min and
median within 2 % of each other, so this is systematic, not contention noise):

* `-j 1` is the fastest setting on every corpus tried, gz or not.
* Degradation is monotone in `j` up to 8 and then falls off a cliff at 16
  (`0.09×`), i.e. 11 s where `-j 1` takes 1.03 s on 4.8 M reads.
* It is not duplication load: the 4.8 M-read corpus is 0.33 % duplicates, the
  gz corpus is ~99 % duplicates (eight copies of one unit), and both degrade.
* It is not corpus construction: it reproduces on a single unmodified unit
  (200 k reads, 0.01 % duplicates).

**Not fixed here.** A thread-scaling defect of this shape is a change to the
dedup engine's concurrency design (sharded verifier, per-thread tables, or a
serialised decision stage), not a fix-pass edit, and it needs its own before/after
campaign. It is recorded in `docs/STATUS.md` under known defects with the
reproduction below.

## Reproduction

```bash
# base worktree (the published 1.3.0 tree)
git -C ~/FastSuite worktree add --detach /mnt/Fast/cost-bench/base 4f65573
cd /mnt/Fast/cost-bench/base && ~/.pixi/bin/pixi run --as-is \
    cargo build --release --manifest-path tools/fastdedup/Cargo.toml

# one unmodified unit, uncompressed and gz
zcat units/human_gut/DRR425489.R1.fastq.gz > corpus-one/R1.fq   # + R2
for j in 1 2 4 8 16; do
  tools/fastdedup/target/release/fastdedup -1 corpus-one/R1.fq -2 corpus-one/R2.fq \
      -o /tmp/o.fq -p /tmp/o2.fq -j $j --force 2>&1 | grep 'reads in'
done
```

`bench/scripts/thread-scaling.sh` runs exactly this, for both binaries, and
prints the table above.

## Harness traps this file exists to prevent repeating

Recorded because each one produced a *plausible* number that was pure garbage:

* **`bc` is not installed** on the desktop, and the locale uses a comma decimal:
  `$(echo "$t1 - $t0" | bc -l)` yields an empty string, and Python then fails on
  `0,000`. Time in integer nanoseconds, reduce in one Python pass with
  `LC_ALL=C`.
* **`/usr/bin/time` is not installed.** `-f '%M'` silently made every run fail in
  0.7 ms, and the table read "0.74 ms, 0 MB RSS" for a 584 MB input.
* **fastdedup's summary goes to STDOUT** for file outputs (only a `-` stream moves
  it to stderr), so `2>&1 >/dev/null` captured nothing.
* **`-q` suppresses the summary entirely**, and `elapsed_s` lives at the report's
  TOP level, not under `metrics` — both cost a run.
