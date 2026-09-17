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
docs you need — not the build directory. The campaign wrapper `stage-to-io.sh`
does exactly this (dry-run first, `--run` to actually scp).

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
export PATH=$PWD/.pixi/envs/build/bin:$PATH
cargo --version
```

Leave with `exit`, and cancel anything forgotten with `scancel -u ribesr`.

## 3. What the 1 h cap means for this project

| work | fits in one `srun` (≤1 h)? | how to run it |
|---|---|---|
| pilot job (build + smoke test) | **yes**, typically minutes | `srun ... --pty bash` then cell steps, or `sbatch` |
| **one campaign cell** (1 dataset × 600 configs, vs sim gt) | **yes** — ~10 min | the cell is sized to fit this cap |
| full campaign (18 datasets) | no | `sbatch --array=0-17 deploy/io/sweep-cell.sbatch` (18 cores, ≤100 ceiling) |
| training the ranker | — | runs locally, not on Io |

The cell is the hinge: 1 dataset × 600 configs ≈ 10 min keeps every unit of real
work inside your 1 h `srun`, while the 18-cell matrix is driven by the array.
The campaign is **light** — it ships code + measured profiles + the config grid,
not reads (each dataset is *simulated* from its profile on the cluster), so there
is no TB-scale transfer.

## 4. The campaign (steps 1+2 of the validation)

Stage once (dry-run first):
```fish
bash deploy/io/stage-to-io.sh          # print the plan
bash deploy/io/stage-to-io.sh --run    # actually scp (needs your OK)
```
Then on Io (in an `srun` shell or via `sbatch`):
```bash
cd /scratch/users/ribesr/fastsuite-bench
export PATH=$PWD/.pixi/envs/build/bin:$PATH      # env by directory (pixi deadlocks)
sbatch --array=0-17 deploy/io/sweep-cell.sbatch  # 18 cells × 600 configs
```
Per-dataset results land in `results/sweep/<RUN>.tsv`; pull them back with
`scp` and correlate best-config-per-dataset against shape features locally.

## 5. Debugging a cell by hand
```bash
# inside the srun shell
cd /scratch/users/ribesr/fastsuite-bench
export PATH=$PWD/.pixi/envs/build/bin:$PATH
python3 bench/scripts/resweep_reference.py --run ERR10123282 --out results/sweep/ERR10123282.tsv --pairs 100000
```
`resweep_reference.py` is idempotent: it skips fastprep for any config whose
`.pxmt` already exists, so an interrupted hand-run costs only the elapsed time.
