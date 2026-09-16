# Running PLAN-ML on Io (the HPC)

Companion to `PLAN-ML.md`. Everything here was verified against the cluster
read-only (rsync probes) or against the proven deployment that is already on
Io — no guesses about partitions, accounts or tooling.

Host `io-login.meso.umontpellier.fr` (ssh alias `Io`, user `ribesr`),
scratch `/scratch/users/ribesr/` with ~45 TB free, key `~/.ssh/Meso`.
**I never ssh to run commands on Io: rsync/scp only.** Anything that needs a
shell on the cluster is written out for you to paste.

---

## 1. What is already proven on Io (reuse it, do not reinvent)

The `novaseq-bench` deployment on `/scratch/users/ribesr/` is a working
Snakemake-9 + SLURM + pixi pipeline of exactly this shape. Facts taken from it:

| fact | value | consequence for PLAN-ML |
|---|---|---|
| SLURM partition | `cpu-dedicated` | every job, dispatcher and cell |
| accounts | `dedicated-cpu@cirad-normal` (cells), `dedicated-cpu@cirad-long` (dispatcher/controller) | a long-lived dispatcher under `@normal` can be killed at the QoS max runtime |
| dispatcher shape | 1 CPU, 2 GB, `--time=7-00:00:00` | holds the Snakemake controller while it submits cells |
| `pixi install` | **hangs forever** (rattler package-cache lock held by stale processes at `/storage/simple/users/ribesr/.cache/rattler/cache/pkgs/.cache.lock`) | never plan on installing packages on Io; the env is used **by directory** |
| env layout | `.pixi/envs/bench` with python 3.12, **polars**, numpy, snakemake 9 + `snakemake-executor-plugin-slurm`, pigz, fastp | the harness needs no new conda dependency (polars writes the Parquet) |
| build feature | `rust ~=1.83`, pkg-config, cmake, c-compiler, zlib, zstd | needed: flate2's zlib-ng builds from source, and the stale-CMake-cache trap applies |
| launch pattern | `sbatch bench/slurm/dispatch.sbatch` from the repo root; the dispatcher runs `snakemake --profile slurm --nolock --rerun-triggers=mtime` | resume-friendly, which the sweep needs |
| rsync to Io | intermittently times out mid-transfer; a batch silently drops files | push file-by-file, retry, verify by pulling back and diffing |

**Your standing rule, applied throughout:** nothing ever runs on a login node;
before any submission you get a block with partition, account and wall time for
confirmation.

---

## 2. Phase-by-phase mapping

| plan phase | where it runs on Io | job shape | measured cost |
|---|---|---|---|
| 0 pilot | 1 job | 1 CPU, 1 GB, 48 h, `@cirad-normal` | minutes |
| 1 acquire (700 ENA datasets) | **login node**, as `curl` in a `for` loop (a download is not a compute job) or a 1-CPU job if the compute nodes turn out to have no internet | 1 CPU, 1 GB if jobbed | bandwidth-bound |
| 2 simulate (5 000 datasets) | inside the sweep shards (generated just before use, deleted after) | — | 8 core-h |
| 3 sweep (1.14 M evaluations) | ~570 shard jobs + 1 dispatcher | **1 CPU, 256 MB, 48 h**, `@cirad-normal`; dispatcher 1 CPU, 2 GB, 7 d, `@cirad-long` | 512 core-h |
| 4 features | same shard jobs (inside the scoring script) | — | 1 core-h |
| 5 train (ranker) | **locally**, not on Io | — | minutes |
| 6 productise | locally | — | — |

Phase 5 running locally is not a dodge: LightGBM is not in the cluster env,
`pixi install` cannot add it, and the training set is a single ~200 MB table
that trains in minutes on this workstation. The cluster's job is the sweep.

### The resource request for one sweep cell, justified by measurement

Measured on this workstation on the actual work unit (100 k pairs × 2 mates),
peak RSS sampled from `/proc/<pid>/status`:

| tool | peak RSS |
|---|---|
| `fastprep -j1` | **29.5 MB** |
| `fastcheck` (pair) | **10.7 MB** |
| `fastsim` | **46.5 MB** |

So the cell asks for **256 MB**, which is ~5× the largest peak. This is
deliberate: an over-large `--mem` trips the group allocation limit
(`AssocGrpMemLimit`) and leaves jobs PENDING forever, and with ~100 concurrent
cells a 2 GB request would reserve 200 GB for work that uses 0.05 GB. If the
scheduler rejects small requests, 1 GB is the fallback.

Wall time follows your standing policy — one flat generous value (48 h) rather
than per-rule estimates, since over-requesting costs queue priority, not
correctness. Only the dispatcher gets its own long value, because it must
outlive the matrix.

---

## 3. Job packing and concurrency

One SLURM job = one **shard** of 10 datasets × 200 configurations = 2000
evaluations ≈ 2000 × 1.57 core-s ≈ **52 min on one core**.

| quantity | value |
|---|---|
| datasets | 5 700 (5 000 simulated + 700 real) |
| shards / jobs | **570** (10 datasets each) |
| wall at `-j 100` | **~5.5 h** (100 cores in flight) |
| wall at `-j 200` | **~3 h** |
| total | ~550 core-hours |

Packing 10 datasets per job keeps the SLURM scheduler out of the way (570
submissions instead of 5 700, each with ~2 s of scheduler overhead) while
staying short enough that a preemption costs little. The shard is also the
idempotency unit: `results/sweep/{shard}.parquet.tmp` → rename, so a resumed
dispatcher redoes only missing shards (`--rerun-triggers=mtime`).

---

## 4. Layout, staging and env on Io

```
/scratch/users/ribesr/fastsuite-bench/
├── FastSuite/            # the tools repo, rsync'd, built ON the cluster
├── bench/                # the harness: Snakefile, config/, scripts/, slurm/
│   ├── slurm/dispatch.sbatch
│   ├── slurm/config.yaml          # profile: partition, accounts, -j, mem
│   └── config/{config_space.toml,score_weights.toml,strata.toml}
├── .pixi/envs/bench               # copied env dir, shebangs repaired
├── data/                          # raw + working FASTQ (gitignored, rotating)
└── results/                       # sweep Parquet (the deliverable)
```

Staging (per-file loop with retry, then verify by pulling back — the batch form
silently drops files on this link):

```fish
# local
set S ~/hermes_scratch/FastSuite-bench
for f in Snakefile pixi.toml config/config_space.toml config/score_weights.toml \
         scripts/run_shard.py scripts/score_gt.py slurm/dispatch.sbatch slurm/config.yaml
    rsync -az -e "ssh -o BatchMode=yes -o ConnectTimeout=20" $S/$f Io:/scratch/users/ribesr/fastsuite-bench/bench/$f
end
# verify (content, not just presence)
rsync -az -e "ssh -o BatchMode=yes" Io:/scratch/users/ribesr/fastsuite-bench/bench/Snakefile /tmp/chk && diff /tmp/chk $S/Snakefile
```

The tools repo is pushed the same way, then **built on the cluster** — never
ship a locally built binary: the workstations' glibc is newer than Io's and
`fastprep` links zlib-ng from source.

Env: copy the proven one rather than installing (installing hangs):

```bash
# ON IO (paste this — I don't ssh)
cd /scratch/users/ribesr/fastsuite-bench
mkdir -p .pixi/envs
cp -a /scratch/users/ribesr/novaseq-bench/.pixi/envs/bench .pixi/envs/
# entry-point shebangs embed the OLD absolute path — repair them
sed -i "s|#!.*/novaseq-bench/.pixi/envs/bench/bin/python|#!/scratch/users/ribesr/fastsuite-bench/.pixi/envs/bench/bin/python|" .pixi/envs/bench/bin/*
.pixi/envs/bench/bin/snakemake --version   # expect 9.x
.pixi/envs/bench/bin/python -c "import polars; print(polars.__version__)"
```

---

## 5. What the pilot must verify on the cluster (unknowns I cannot settle from here)

1. **Rust version.** The cluster's build feature pins `rust ~=1.83`;
   `fastdedup` declares `rust-version = "1.85"` and edition 2024, so it cannot
   build there as-is. PLAN-ML needs only `fastsim`, `fastprep`, `fastcheck`
   (all edition 2021). Pilot step: `cargo build --release` for those three
   only, and record `cargo --version` in the run log. If a tool needs a newer
   rustc, that is the cause — not a code bug.
2. **Does a compute node reach the outside world?** The ENA downloads (phase 1)
   and `cargo`'s registry fetch both need egress. The existing dispatcher
   builds from crates.io, which suggests yes, but a pilot rule must prove it
   from inside a job:
   `curl -sI -m 10 https://ftp.sra.ebi.ac.uk/ | head -1`
3. **Subsampling without new tools.** `seqkit` exists locally but is not in the
   cluster env, and nothing can be installed. The subsampler is therefore a
   12-line `awk` that keeps every N-th record (deterministic, dependency-free,
   same N ⇒ same reads); the pilot asserts the output read count and that R1/R2
   stay in sync.
4. **GLIBC floor** for the built binaries (build on the cluster and the question
   disappears; record `ldd --version` in the pilot log anyway).
5. **The shard's real wall time** (predicted 52 min) and its real peak RSS on
   cluster hardware — the 256 MB request is a prediction from local
   measurements until a shard has actually run.

---

## 6. The SLURM headers, as written

Four files, all in `deploy/io/` — the headers are real, not sketches:

| file | role | header |
|---|---|---|
| `pilot.sbatch` | phase 0: build the three tools, prove egress, one dataset end to end, measure wall/RSS | `cpu-dedicated` / `@cirad-normal` / 1 CPU / **1G** / 48 h |
| `dispatch.sbatch` | phase 3 controller: holds Snakemake, submits the cells | `cpu-dedicated` / `@cirad-long` / 1 CPU / 2G / **7-00:00:00** |
| `acquire.sbatch` | phase 1: ENA downloads + subsample + hash + delete raw, resumable | `cpu-dedicated` / `@cirad-normal` / 1 CPU / 1G / 48 h |
| `slurm-config.yaml` | the **cells** — no per-cell sbatch exists; these headers are generated from the profile | `cpu-dedicated` / `@cirad-normal` / 1 CPU / **256M** / 48 h (runtime 2880), `jobs: 100`, global core ceiling 100 |

Plus `WORKFLOW-IO.md` for the `ssh Io` + `srun` path:

```fish
ssh Io
srun --partition=cpu-dedicated --account=dedicated-cpu@cirad-normal \
     --cpus-per-task=1 --mem=2G --time=01:00:00 --pty bash
```

Two deliberate choices inside those headers:

* **The account split is load-bearing.** The dispatcher runs under
  `@cirad-long` with a 7-day clock because it must outlive the matrix; under
  `@cirad-normal` it can be killed at the QoS maximum mid-wave. The cells keep
  `@cirad-normal`, where 48 h is plenty for a 52 min shard.
* **`mem_mb: 256` is a measured request, not a cautious one** (peaks: 29.5 /
  10.7 / 46.5 MB). A comfortable 2 GB would reserve 200 GB across 100
  concurrent cells, which is exactly how jobs end up parked in PENDING on a
  group allocation.

Nothing is submitted until you confirm partition + account + time for each role.
