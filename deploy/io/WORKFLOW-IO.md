# Working on Io: scp for transfer, `ssh Io` + `srun` for anything that runs

**The rule (corrected 2026-09):** no `rsync` onto the server — use `scp`. Then
`ssh Io`, take a slice of a partition with `srun` (1 hour maximum), and run
things in there. Nothing runs on the login node.

## 1. Transfer with scp

```fish
# single files (the common case)
scp FastSuite/tools/fastcheck/src/accum.rs Io:/scratch/users/ribesr/fastsuite-bench/FastSuite/tools/fastcheck/src/

# a whole subdirectory that is safe to overwrite
scp -r bench/scripts Io:/scratch/users/ribesr/fastsuite-bench/bench/

# verify content, not just presence (a file can land but be stale)
scp Io:/scratch/users/ribesr/fastsuite-bench/bench/Snakefile /tmp/chk; diff /tmp/chk bench/Snakefile
```

Never `scp -r` a Rust checkout as-is: `target/` is gigabytes and the binary is
rebuilt on the cluster anyway. Push `src/`, `Cargo.toml`, `Cargo.lock` and the
docs you need — not the build directory.

## 2. Get a shell inside a partition (1 h cap)

```fish
ssh Io
# then, on the login node:
srun --partition=cpu-dedicated --account=dedicated-cpu@cirad-normal \
     --cpus-per-task=1 --mem=2G --time=01:00:00 --pty bash
```

Inside that shell the env is used **by directory** (`pixi install` deadlocks on
this cluster):

```bash
cd /scratch/users/ribesr/fastsuite-bench
export PATH=$PWD/.pixi/envs/bench/bin:$PATH
snakemake --version
```

Leave with `exit`, and cancel anything forgotten with `scancel -u ribesr`.

## 3. What the 1 h cap means for this project

| work | fits in one `srun` (≤1 h)? | how to run it |
|---|---|---|
| pilot job (build + smoke test) | **yes**, typically minutes | `bash slurm/pilot.sbatch` style steps, or `sbatch` |
| **one sweep shard** (10 datasets × 200 configs) | **yes** — ~52 min by design | the shard size was chosen to fit this cap |
| full sweep (570 shards, ~5.5 h) | no | the `sbatch` dispatcher holds the controller; it submits the cells |
| corpus download (1.4 TB) | no | `sbatch slurm/acquire.sbatch`, resumable |
| training the ranker | — | runs locally, not on Io |

The shard size is the hinge: 10 datasets ≈ 52 min keeps every unit of real work
inside your 1 h `srun`, while the 570-shard matrix is driven by the dispatcher.

## 4. Debugging a shard by hand

```bash
# inside the srun shell
cd /scratch/users/ribesr/fastsuite-bench/bench
export PATH=$PWD/../.pixi/envs/bench/bin:$PATH
export FASTPREP_BIN=$PWD/../FastSuite/tools/fastprep/target/release/fastprep
python3 scripts/run_shard.py --shard 0007 --config config/config_space.toml --out results/sweep/
```

`run_shard.py` is idempotent: it writes `results/sweep/0007.parquet.tmp` and
renames on success, so an interrupted hand-run costs only the elapsed time.
