# PLAN-ML — execution on the desktop

> **Revision 4, 2026-10-12: Io is gone, and all benchmark compute runs on the
> desktop — nothing else.** The ISDM-MESO cluster access that PLAN-ML-STEPS
> (rev 3) assumed is no longer available, and the laptop stays out of the
> benchmark path entirely (it is the interactive machine: manifests, fitting,
> model training, and the free A2/A3 tests only).
>
> | machine | cores | RAM | role |
> |---|---|---|---|
> | **desktop** | 16 | 48 GB | **the only benchmark machine** — pilot, sweep, fastp comparison, alignment subset |
> | laptop | 16 | 15 GB | interactive only: A2/A3 tests, ENA manifests, label/model training |
>
> **What changed:** no `sbatch`, no `srun`, no `cpu-dedicated`, no
> `@cirad-*` accounts, no `scp` to a shared scratch. The desktop runs cells
> from a queue directory with `flock`; the laptop prepares the queue and
> collects results over SSH. **What did not change:** the science (A2/A3
> first, GATE 2 before the full run), the cost model (~8 core-min per
> 600-evaluation cell), or the rule that nothing runs unsupervised on a
> machine the user is using — the desktop workers are started by hand.
>
> The historical SLURM scripts and `docs/DEPLOY-IO.md` are kept for
> reference under `deploy/io/` but are **not** the execution path.

---

## The compute budget, rescaled (desktop-only)

Measured cost of one config evaluation (from the bench harness): ~0.8 core-s
per 50 k pairs (`fastprep + fastsim gt scoring`). A cell = 1 dataset × 200
configs × 3 variants ≈ 600 evaluations ≈ **8 core-min**.

| phase | work | wall on the desktop (16 cores) |
|---|---|---|
| A2/A3 (recovery + round-trip) | seconds × N tests | minutes — **on the laptop, free** |
| B (corpus manifests) | ENA API queries, no reads | ~1 h — on the laptop |
| C (pilot) | 30 datasets × 60 configs × 3 variants | **~3 h** overnight |
| D (full run) | 500 datasets × 200 configs × 3 variants | **~5–7 days** of continuous running |
| D3 (alignment subset) | 50 datasets × (winner, fastp), bwa-mem2 | ~1 day — RAM-rich, the desktop handles it |
| E (labels + models) | local, minutes | ~2 h — on the laptop |
| F (productise) | local, minutes | ~1 h — on the laptop |

**Honest framing:** 300 000 evaluations at 8 core-min per 600 → 4 000 core-h
total. The desktop gives 16 cores ≈ 384 core-h per 24 h, so the full sweep is
**~5–7 days of continuous running** — spread over a week of evenings, not a
weekend. The 48 GB RAM is the second gift: `fastprep`'s dedup admission
(4.6 GB at 3e8 reads) never refuses, and `bwa-mem2` indexes a human
reference panel comfortably.

**The dial if that is too slow:** the sweep does not need 500 datasets to
answer "does a fitted profile reproduce a dataset's QC". Phase C (30
datasets) already answers the go/no-go. Phase D can be cut to the top 100
datasets by domain (~20 000 core-min ≈ 21 core-h ≈ **one weekend**) and the
model retrained when the corpus grows. **Decide the corpus size at GATE 2
with the pilot numbers in hand** — do not pre-commit to 500.

## Coordination: no scheduler, just a queue file

The laptop prepares; the desktop consumes:

```fish
# on the laptop (the coordinator): the queue is a directory of pending cell files
mkdir -p campaign/queue campaign/done
for run in (cat corpus.tsv); touch campaign/queue/$run.cell; end

# on the desktop (the ONLY benchmark machine), one command, repeatable:
rsync -a laptop:campaign/queue/ campaign/queue/          # see pending cells
for cell in campaign/queue/*.cell
    if mkdir campaign/locks/(basename $cell).lock 2>/dev/null   # flock-free claim
        pixi run cell $cell
        rsync -a results/ laptop:campaign/done/
        rm campaign/queue/$cell campaign/locks/$cell.lock
    end
end
```

No job scheduler, no shared filesystem, no cluster account. A worker that
dies mid-cell leaves a `.partial` file; the next run redoes that cell (cells
are idempotent — the output is content-addressed by dataset + config hash).
`campaign/README.md` (to be written at GATE 2) documents the exact rsync
layout; this document is the plan, not the runbook. The desktop is the only
machine that ever runs `fastprep`/`fastsim` at benchmark scale.

## Phase-by-phase (the steps, with what replaced what)

### Phase A — unchanged, local, free, on the laptop

A2 (recovery test) and A3 (round-trip) run on the laptop in seconds. They are
now written up as TDD tasks in
`docs/plans/2026-10-12_fastcheck-exactness.md` (Task 13) and the suite plan's
Phase A. **Nothing in Phases B–D starts until they pass** — the fitter is the
hinge of the whole design.

### Phase B — corpus (laptop; replaces the Io staging steps)

* **B1** `bench/scripts/ena_query.py`: build `manifests/real.tsv` from the ENA
  portal API (100 datasets per domain, every row with `fastq_ftp` +
  `fastq_md5`). Runs on the laptop; the ENA API needs no cluster.
* **B2** `bench/scripts/subsample.sh`: deterministic every-Nth-record
  subsample to 100 k pairs, md5 check against ENA, delete the raw. **D-11
  still applies**: no raw file may persist on either machine beyond the
  subsample's verification.
* **B3 (was "stage to Io")** → **distribute the harness**: the desktop gets
  the same checkout via `git pull` (the campaign code is in the repo) and
  runs its own `pixi run build`. The corpus manifest and the queue directory
  are rsynced from the laptop. Nothing is "built on the cluster".
* **B4** `ena_query.py --stats` per-domain sanity check (unchanged).

**GATE 1 — the pilot subject:** 30 datasets, 10 each from generalised, human
gut, virome (unchanged proposal).

### Phase C — the pilot (desktop only)

* **C1** full chain on 30 datasets, 60 configs × 3 variants: the same
  `qc_to_profile → fastsim → fastprep` evaluation the old plan ran, driven by
  the queue-file worker on the desktop instead of `sbatch --array`.
* **C2** closed loop on real data (unchanged — the go/no-go).
* **C3** winner vs fastp on the real reads (unchanged; fastp runs on the
  desktop, it is the bench reference).
* **C4** decision memo `bench/PILOT.md` (unchanged).

**GATE 2 — scale or fix.** With the pilot's measured cell cost in hand, pick
the Phase D corpus size: 500 (5–7 days) or 100 (a weekend) — the model is
retrainable, so the smaller corpus is a legitimate choice, not a compromise
to hide.

### Phase D — the full run (desktop only)

* **D1** the sweep, driven by the queue on the desktop. Acceptance
  unchanged: every cell file present, per-dataset ≥ 190/200 configs scored.
* **D2** validation on real reads (unchanged).
* **D3** alignment subset: 50 datasets × (winner, fastp) with `bwa-mem2` /
  `minimap2 -x sr` — runs on the desktop (48 GB; a human reference panel
  indexes comfortably). If the desktop is unavailable for a stretch, cut the
  reference panel or run it last.
* **D4** aggregate to `matrix.parquet` with DuckDB/polars — **on the
  laptop** (it is interactive work, minutes, and needs no benchmark RAM).

### Phase E and F — unchanged, laptop

Labels, leak test, study-level splits, per-domain rankers, model cards,
`fastprep --auto-params` (Rust tree-walking) and `fasthelper recommend`.
All laptop work, all minutes.

## What the old SLURM material is now

* `deploy/io/*.sbatch`, `deploy/io/stage-to-io.sh`, `deploy/io/WORKFLOW-IO.md`,
  `docs/DEPLOY-IO.md`: **historical**. They document a working deployment
  that no longer has a target. Kept in the tree (they cost nothing and the
  design — cells, manifests, idempotent reruns — is reusable if another
  cluster appears), but the README no longer advertises Io and no plan
  depends on them.
* The campaign's coordination model (cells, queue, idempotency) is borrowed
  from them deliberately — the design outlived the cluster.

## Risks

| risk | mitigation |
|---|---|
| the 5–7 day full sweep never finishes | GATE 2 picks the corpus size with pilot numbers; 100 datasets is a legitimate target |
| the desktop is in use (it is your gaming machine) | workers are started by hand; the queue tolerates the machine joining late; cells are resumable |
| rsync drops files (the old Io lesson) | every cell verifies its output hash before being marked done |
| ENA downloads are slow or fail | manifests carry md5s; failures are per-dataset, not per-run |
