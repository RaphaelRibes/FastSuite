# fastmaps — round-2 speed propositions (2026-10-21)

Evidence: `docs/optim/raw/round2/_salvage-fastmaps.md` (42 results incl. a
working cheap-hasher prototype), `bench/profile/pipeline_5m.csv` (5 M pairs:
3.856 s wall / 7.552 s cpu / **22.7 MB peak**), `docs/BENCHMARKS.md:217-219`.

### FM2-1 — Cheap hasher for the per-pixel maps
allocation · [M] (hasher prototype: std SipHash-1-3 9.4 ns/read -> identity+mul
5.9 ns/read -> dense vec accumulator 1.0 ns/read) · byte-identical
`worker.rs` hashes `HashMap<(u16,u16), PixelAcc>` (`counts`, `counts_r1`,
`counts_r2`) and `HashMap<u16,u64>` (`per_file_counts`, `worker.rs:64,158,184,396`)
per read. SipHash's DoS resistance is pointless for derived (u16,u16) keys; and
`per_file_counts` should be a `Vec` indexed by `file_id`. A/B: wall at 5 M pairs.

### FM2-2 — Release profile: lto = "fat", codegen-units = 1
build · [V] (fastmaps' own `Cargo.toml` `[profile.release]` carries ONLY
`opt-level = 3` — verified this session; the scratch prototype measured the
other flags separately) · byte-identical
Add the two flags; price the build-time cost once against the runtime win.

### FM2-3 — Fused per-read passes
algorithmic · [E] · byte-identical
One pass per read with three accumulator updates beats three scans in fewer
cache lines. A/B: wall at 5 M pairs.

### FM2-4 — Index record 12 B padded -> 8 B packed
buffer-width · [V] (the doc comment claims 8 B; the struct is (u16,u16,u16,u32)
= 10 B padded to 12; the spilled record already drops y) · proven-equivalent
Pack two u16 + u32 or drop y from the in-memory record too. Halves index memory.

| ID | title | tier | gate | target | expected effect |
|---|---|---|---|---|---|
| FM2-1 | cheap hasher / Vec for per-file | M | byte-identical | worker.rs | ~3-9x on map ops |
| FM2-2 | lto + codegen-units 1 | V | byte-identical | Cargo.toml | free runtime win |
| FM2-3 | fused passes | E | byte-identical | scan path | modest wall |
| FM2-4 | 8 B index record | V | proven-equivalent | index struct | -33 % index memory |
