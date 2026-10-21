# fasthelper — round-2 speed propositions (2026-10-21)

Evidence: `docs/optim/raw/round2/_salvage-fasthelper.md` (52 results).

### FH2-1 — Sample-level parallelism in `concatenate` / render
work-partitioning · [M] (shard probe: 16 shards x 625 grids, serial sum 8.90 s
-> 16 concurrent processes 1.79 s = **5.0x on 16 cores**) · byte-identical
outputs
`-P/--jobs` exists for steps; extend it to report/concatenate aggregation over
samples. This is the measured shape of the win.

### FH2-2 — Adapter-guard cost: reuse the step's own fastcheck pass
I/O · [E] (the guard measures fastcheck adapter content on reads entering AND
leaving every trimming step — an extra fastcheck pass per step per sample;
sample-size probe: 20k-read sample 0.216 s vs full-pass order of magnitude) ·
proven-equivalent (same numbers, fewer processes)
Where the chain already runs fastcheck on that output, reuse it; otherwise keep
the sampling knob.

### FH2-3 — Crash-resume cleanliness
correctness · [V] (probe: SIGKILL/timeout leaves partial sample dirs —
`a_clean_R1.fastq` etc. with no terminal manifest record; `--from` resume
exists) · behaviour-change (flag loudly)
Write step outputs to `*.tmp` then rename; a killed run leaves the step's full
output or nothing. Makes the resume contract honest.

### FH2-4 — Parse each report once
allocation · [E] · byte-identical
check -> aggregate -> concatenate -> html each re-read reports; parse once into
the in-memory model and render from it.

### FH2-5 — The two reproduced output defects (verify before fixing)
correctness · [V] (round-1 dossier) · behaviour-change
(1) `concatenate` counts the `heatmap.report` wrapper as a phantom sample
(stem not in the GENERIC name set -> "4 report(s) from 3 sample(s)"). (2)
`check -c` accepts duplicate `[[samples]]` names (shared work dir/manifest/
resume key) — must refuse. The 2026-10-21 wave may have fixed (1) via
sample-name propagation — verify on the built binary before proposing a patch.

| ID | title | tier | gate | target | expected effect |
|---|---|---|---|---|---|
| FH2-1 | sample-level parallelism | M | byte-identical | concatenate/render | 5.0x measured |
| FH2-2 | reuse fastcheck pass in guard | E | proven-equivalent | adapter_guard.rs | one process per step saved |
| FH2-3 | tmp+rename step outputs | V | behaviour-change | runner.rs | honest resume |
| FH2-4 | parse reports once | E | byte-identical | report path | modest |
| FH2-5 | two output defects | V | behaviour-change | concatenate/check | correctness |
