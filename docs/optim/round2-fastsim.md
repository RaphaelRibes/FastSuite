# fastsim — round-2 speed propositions (2026-10-21)

Evidence: `docs/optim/raw/round2/_salvage-fastsim.md` (34 results),
`bench/profile/pipeline_5m.csv` (**fastsim is 283.5 s wall / 283.0 s cpu of the
~307 s chain — 80 % of total wall**, single-threaded by construction: `-j` does
nothing for it).

### FS2-1 — Runner-level default `--gz-level 1` (or 4) for simulated corpora
default-behaviour · [M] (probe, 200k pairs, mandrill_E358, built binary: gz6
13.67/16.40 s, gz4 4.18/4.78 s, gz1 3.11/3.68 s, `--no-gzip` 1.85 s — deflate is
~86 % of the gz6 run) · behaviour-change (flag loudly)
The flag landed (FS-1) and defaults to 6 "so existing campaigns keep their exact
bytes" (`cli.rs`). Keep the tool default for byte-compat; make **fasthelper's
runner** select level 1 for intermediates and document the trade. Deflate is the
single largest wall lever in the whole suite.

### FS2-2 — Parallel chunked simulation with per-pair deterministic RNG
work-partitioning · [E] (~8x on 8 cores; pairs are independent given a per-pair
stream) · behaviour-change (the emitted stream changes even at `-j 1` if the RNG
call pattern changes — needs a versioned GT format and a closed-loop gate)
Derive each pair's stream from `(seed, pair_index)`: reproducible AND parallel.
Gate: the closed loop (fastcheck --profile -> fastsim -> fastcheck) must still
pass at the same tolerance, and the GT schema must carry the RNG contract.

### FS2-3 — Fuse the per-base generation loops
algorithmic · [E] · byte-identical if the RNG draw order is kept
`quality.rs` + `biology.rs` each walk the read; one walk with hoisted
distribution objects (no per-base RNG construction, LUT the quality CDF).
A/B: wall on the 200k probe.

### FS2-4 — `--stats` counters zero-cost when absent
allocation · [V] (FS-11 landed) · byte-identical
Verify the counters sit behind the flag in the hot loop, not collected then
discarded.

| ID | title | tier | gate | target | expected effect |
|---|---|---|---|---|---|
| FS2-1 | runner gz-level 1 | M | behaviour-change (.gz bytes) | runner/profile default | up to 4x wall |
| FS2-2 | parallel per-pair RNG | E | behaviour-change (versioned GT) | simulation core | up to ~8x wall |
| FS2-3 | fused generation loops | E | byte-identical (draw order kept) | quality/biology | modest wall |
| FS2-4 | stats no-op path | V | byte-identical | stats | verify only |
