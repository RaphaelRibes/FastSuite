# fastdedup — round-2 speed propositions (2026-10-21)

Evidence: `docs/optim/raw/round2/_salvage-fastdedup.md` (55 results),
`bench/profile/{dedup_threads,gz_levels,pipeline_5m}.csv`.

### FD2-1 — Record ORDER is not deterministic in parallel mode
correctness · [V] (probe: a.fq vs c.fq sorted-seq md5 identical
`f07fb28e2037567da0170f9ecf766620`, 204 675 distinct each, but first differing
record at line 826 568 / record 206 641) · behaviour-change (flag loudly) OR
documented contract
The KEEP SET is identical (same sequences and counts) but OUTPUT ORDER differs
across runs/modes. Anything gating on byte-identity (including our own gates!)
must compare sorted records — which the probes had to do. Either serialise emit
order (the writer is already the single ordered stage; sequence the chunk order)
or write "order is not part of the contract" into REPORT-SCHEMA.md.

### FD2-2 — Default `--gz-level` to 1 for intermediates
I/O · [M] (`gz_levels.csv`, 5 M pairs: gz1 18.6 s / 154 MB / 1679 MB out vs
gz6 22.79 s / 4372 MB / 888 MB out; gz9 121 s) · behaviour-change (.gz bytes)
The writer is the memory hog (4.37 GB peak at gz6 — 28x gz1) and gz1 is *faster*
on this box. The flag exists; the default is the proposition.

### FD2-3 — Deflate backend swap
I/O · [E] (typically 2-3x faster deflate at equal-or-better ratio for
zlib-rs/libdeflate-class backends) · behaviour-change (.gz bytes)
flate2's default backend -> a faster one; pairs with FD2-2. A/B on the
`gz_levels.csv` protocol.

### FD2-4 — Explain the -j 1..4 flatline, then fix or document
work-partitioning · [V] (`dedup_threads.csv`: j1 94.84 s, j2 92.13 s, j4 92.2 s
all ~92 s wall / ~100 MB vs j8 22.73 s / 4368 MB) · byte-identical
With uncompressed output the hash pipeline is ~10 s at any `-j` (README); the
92 s flatline is read-side I/O plus writer stall. One A/B with uncompressed
output separates parse from compression cost; then either the queue bound
(`53d6686`) needs tuning at low `-j` or the flatline is documented.

### FD2-5 — 128-bit hash visibility
default-behaviour · [V] (review: the 128-bit path is unreachable without
`--pre-count`) · byte-identical (collision analysis)
Expose it or document the 64-bit collision bound in the report schema.

| ID | title | tier | gate | target | expected effect |
|---|---|---|---|---|---|
| FD2-1 | emit-order determinism | V | behaviour-change or documented | writer | honest byte-gates |
| FD2-2 | default gz-level 1 | M | behaviour-change (.gz bytes) | writer default | 18 % faster wall, 28x less RAM |
| FD2-3 | faster deflate backend | E | behaviour-change (.gz bytes) | writer | 2-3x on the output stage |
| FD2-4 | -j flatline | V | byte-identical | pipeline | explanation or fix |
| FD2-5 | 128-bit hash | V | byte-identical | hasher | documented bound |
