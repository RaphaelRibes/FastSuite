# Optimisation bench (Snakemake, desktop)

Measures the six tools **before vs after** the optimisation round described in
`docs/optim/`. Written for the desktop (16 cores, 48 GB), because that is the
only machine the project benchmarks on; the laptop stays interactive-only.

## Why this workflow and not the paper benches

The `~/work/Fast*-Bench` repos measure one tool against a comparator on fixed
sizes, one script per tool. This workflow answers a different question — *did
this change help, and did it move memory or only wall time?* — for all six tools
in one DAG, with both arms in the same run so the machine state is shared.

The rules it enforces, each a documented failure of an earlier bench:

* one **row per file**, aggregated with a **count assertion**: a shared CSV that
  a driver opens for a header and never appends to silently empties a campaign,
  and an appending driver silently doubles it on a re-run;
* `-j N` always paired with `taskset -c 0..N-1`, and a `cpu_over_wall` column:
  an escaped child shows up as > 1.05 on a pinned cell. That column is how the
  fastprep rows in `FastPrep-Bench/results/results_gz6.csv` were caught claiming
  to be pinned while burning 2.44 cores;
* the child's **own** `ru_maxrss` from `wait4` plus a live `VmHWM` poll, never
  `getrusage(RUSAGE_CHILDREN)` — that is a never-resetting high-water mark and
  once made three tools report the same 1209.8 MB;
* **no row unless** exit 0, non-empty output, and no unexpected stderr. The
  citation line and tool progress banners are allow-listed; anything matching
  `error|panicked|failed to|cannot|no such file` fails the cell;
* the read corpus is **generated once** by the frozen baseline `fastsim` at
  `--seed 42` from the committed profile, so a change to fastsim cannot move the
  input under its own comparison.

## Layout

```
bench/optim/Snakefile        DAG: reads -> cell (one per tool x size x threads x arm) -> aggregate
bench/optim/config.yaml      desktop paths and the cell grid
bench/optim/scripts/measure.py    one cell, refuses to emit a fake row
bench/optim/scripts/helper_tree.py    synthesises the fasthelper report tree
```

## Running it on the desktop

```bash
# 0. corpus + scratch on the data partition, never $HOME
ssh desktop 'bash -c "mkdir -p /mnt/Fast/fastsuite-optim/{results,reads,baseline-bin}"'

# 1. baseline binaries: built from the pre-change commit, frozen
#    (on the laptop: git worktree add /tmp/fs-baseline <sha>; build the six;
#     scp them to /mnt/Fast/fastsuite-optim/baseline-bin)

# 2. the tree under test, built in place
ssh desktop 'bash -c "cd ~/FastSuite && ~/.pixi/bin/pixi run build"'
ssh desktop 'bash -c "mkdir -p /mnt/Fast/fastsuite-optim/new-bin && for t in fastcheck fastdedup fastprep fastsim fastmaps fasthelper; do cp ~/FastSuite/tools/\$t/target/release/\$t /mnt/Fast/fastsuite-optim/new-bin/; done"'

# 3. the run, detached so an ssh drop cannot kill it
ssh desktop 'bash -c "cd ~/FastSuite/bench/optim && setsid nohup ~/.pixi/bin/pixi run snakemake -c 6 --rerun-triggers=mtime all > /mnt/Fast/fastsuite-optim/bench.log 2>&1 & "'
```

`config.yaml` points at `/mnt/Fast/fastsuite-optim/repo` for the profile; on a
checkout elsewhere, edit `repo` and `profile`.

## What a cell is

`tools x sizes x threads x arms`, where the arms are `baseline` (frozen binary)
and `new` (tree under test). `fastsim` is single-threaded by construction, so it
carries the gzip-level axis instead of threads; the baseline binary predates
`--gz-level`, so its arm only has the default-level cell (measuring a flag the
binary does not have would be a fabricated comparison, and `measure.py` refuses
it rather than silently running the default).

`fasthelper` does not read FASTQ: its cells aggregate a synthesised tree of real
reports (`helper_tree.py` clones one genuine fastcheck v1 report and attaches a
real fastmaps grid to every tenth sample), at 2 000 and 5 000 samples.

## Output

`/mnt/Fast/fastsuite-optim/results/results.tsv`, one row per cell:

| column | meaning |
|---|---|
| tool, kind | tool and arm (`baseline` / `new`) |
| reads, threads, gz_level | the cell |
| wall_s, cpu_s, cpu_over_wall | wall clock, the child's own CPU time, and their ratio |
| peak_rss_kb | the child's `ru_maxrss`, cross-checked by the live `VmHWM` poll |
| out_bytes, files, rc | output size, file count, exit code |

Peak RSS on a pinned cell above ~1.05 in `cpu_over_wall` means the child kept
threads of its own, which is a finding about the tool, not about the bench.

## Why the campaign runs with `--cores 1`

Each cell is pinned (`taskset -c 0..N-1` for `-j N`) so that `cpu_s / wall_s` near
1.00 proves no work escaped the pin. That only means something if nothing else is
using those cores. With `--cores 8` Snakemake happily ran several cells at once,
and because every `-j 1` cell is pinned to core 0 they queued on one core: the
wall times came out inflated by contention, `cpu_over_wall` fell below 1, and the
numbers were not comparable between cells.

Two rules follow.

* Run the campaign serially (`--cores 1`). The whole matrix costs ~30-45 min
  serial, which is cheaper than a contaminated table.
* Before launching, check for orphans from an earlier run: a killed Snakemake
  leaves its `taskset`-pinned grandchildren alive, still holding core 0, and the
  next campaign measures against them. `pkill -f 'fastsi[m] --pairs'` (the
  bracket keeps the pattern from matching your own shell) then confirm
  `pgrep -c` reports 0.
* Snakemake also refuses to rebuild after a killed run until `--unlock`, and
  marks the interrupted cells incomplete; `--rerun-incomplete` is required.

## Force the cells, not the corpus

`-F` forces every job in the DAG, including the three `reads` rules, and
regenerating the 5 M-pair corpus costs ~10 min for files that are already there
and byte-identical (fastsim is deterministic at a fixed seed). Target the cells
instead: `--forcerun cell` plus the aggregate, or list the targets explicitly.
`-F all` is only right on a cold workspace.

## Which config produces which table

`config.yaml` is the full matrix — three sizes, three thread counts, two comparators — and serially it
is hours of compute. `config-papers.yaml` is the targeted campaign: one size, one thread count, both
arms, the gzip-level axis, and it is what produced the fifteen rows quoted in `docs/optim/report/` and
in the six papers (archived at `results/2026-10-20-targeted.tsv`). It takes six minutes serially.

`run_bench.sh` defaults to `config-papers.yaml` and `BENCH_CORES=1` for exactly that reason: the first
campaign of this round ran the full matrix at `--cores 8`, which is the one combination that cannot
produce a sound row, and 90 of its cells had to be discarded. Override with `BENCH_CONFIG` and
`BENCH_CORES` only when you know why.

The systemd unit is in `systemd/fs-bench.service` with install instructions at the top.

## Three traps that make a row unreadable rather than wrong

* *The two arms can report the same version.* `baseline-bin/fastdedup` and
  `new-bin/fastdedup` both print `fastdedup 1.3.0` — the version was never bumped
  between the frozen commit and the tree under test. A cell is therefore
  identified by its path, never by `--version`, and a report should quote the
  commit, not the version string.
* *A reboot in the middle of a campaign is silent.* One run was interrupted by a
  host reboot; Snakemake left `IncompleteFilesException` markers on three cells
  and the aggregate never ran. Check `uptime` against the campaign's start time
  before blaming the harness, and re-run only the incomplete cells.
* *The synthetic helper tree is not a sample-count axis if the reports have no
  sample name.* `helper_tree.py` rewrites the `sample` field to name each
  directory, but the fastcheck report on this desktop carries no `sample` key, so
  the rewrite is a no-op and every row is labelled `qc`. That makes
  `helper__2000` and `helper__5000` the same tree, not two sizes. Fix the tree
  builder before reading that column as a scaling axis.
