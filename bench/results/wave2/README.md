# The 2026-10-26 second-wave re-run (desktop, `cachyos-x8664`)

Re-derived after the second fix wave (commits `82c46af`…`d6a866e`, 13:29–13:30),
because those changes sit on the path that produces the reads, the fitted
profiles and the trims:

* fastcheck (`402429c`) — `head_dropout` and the profile's `g_mean` divided a
  `--max-read-len`-truncated G sum by every base in the file (the same file read
  0.05 / 0.10 / 0.30 at caps 20 / 40 / 1000); they now divide by the bases the
  per-position plane holds. The same commit publishes `busy_cores` next to
  `threads`, so a report can no longer claim `-j 8` bought eight cores.
* fastsim (`7f26d71`) — the GT TSV's adapter columns stop publishing template
  `N` padding as adapter bases.
* fastprep (`d6a866e`) — a 3′ artefact is no longer kept whole when the 5′ one
  already spent the trim budget (the two-sided boundary case).
* fastdedup (`82c46af`) — the kept copy of a duplicated sequence is the first in
  input order; the winner used to depend on thread scheduling, so a `-j 8` run
  could pick a different header per run.

The binaries under `tools/` were built at 13:31, i.e. from the tree at `6aee6dc`
(13:30:31); the next commit, `4649a81` (13:40, fastdedup creating a missing
parent directory), landed during the run and cannot touch fastcheck / fastsim /
fastprep output. `work-wave2/` and `results-wave2/` were fresh; the previous
waves' trees were left in place as the baseline.

## Commands (desktop)

```bash
# stage 1: closed loop, 102 datasets x 3 variants (fastcheck -> fastsim -> fastcheck)
cd /mnt/Fast/fastsuite-scratch
WORKROOT=/mnt/Fast/fastsuite-scratch/work-wave2 \
RESULTS=/mnt/Fast/fastsuite-scratch/results-wave2 \
BIN=/home/raphael/FastSuite/tools JOBS=8 THREADS=2 SIM_PAIRS=100000 \
  bash rerun_closed_loop_gt2026.sh          # 306 sim reports, 13:37:17 -> 13:42:20

# stage 2: trim accuracy against the same GT, shipped default
/usr/bin/python3 /home/raphael/FastSuite/bench/scripts/trim_score.py \
  --work /mnt/Fast/fastsuite-scratch/work-wave2 \
  --fastprep /home/raphael/FastSuite/tools/fastprep/target/release/fastprep \
  --out /mnt/Fast/fastsuite-scratch/results-wave2/trim-accuracy-shipped.tsv
```

Both stages ran as `campaign-wave2.service` (systemd `--user`, oneshot):
`Result=success`, `ExecMainStatus=0`, closed loop 13:37:17→13:42:20, trim
13:42:20→13:44:40.

## Closed loop — `closed_loop.tsv`

**183/306 comparisons pass (59.8%)**, **60/102 datasets pass every variant**
(62/102 pass at least one), median per-position G-fraction difference
**0.0092** against the 0.02 bound (raw median
0.00925). 123 rows carry a failing gate.

| domain | PASS | share |
|---|---|---|
| human gut | 51/51 | 100.0% |
| generalised | 47/51 | 92.2% |
| bacterial isolates | 40/51 | 78.4% |
| RNA-seq | 21/51 | 41.2% |
| viromes | 21/51 | 41.2% |
| metagenomes | 3/51 | 5.9% |

Gates over the 123 failing rows: G fraction 121, mean quality
8, GC 4.

**Nothing moved.** `closed_loop.tsv` is byte-identical (sha256
`2076d6a075e11d22265943e8adfab5d9e79f46898321b01dde55d3ae0717ff5b`) to `bench/results/fix2026/closed_loop.tsv` *and* to
`bench/results/gt2026/closed_loop.tsv`: 0/306 rows changed verdict, and not just
over the 17 non-scored columns — all 25 columns match, `verdict`, `failing`,
`verdict_lm` and `failing_lm` included.

That invariance is expected, not disappointing: every scored quantity here is
produced by the **same** fastcheck invocation on both sides — once on the real
reads, once on the reads fastsim generated from the first run's profile — and
closed_loop.py compares the two reports. A systematic report-side estimator
error therefore cancels between numerator and denominator of the comparison; it
moves the absolute values both reports publish, not the difference between them.
The fix wave changes what the reports *say*, not the agreement they measure, so
this stage is a regression check and its reproduction of 183/306 is the result.

## Trim accuracy — `trim-accuracy-shipped.tsv`

Shipped default (`0.67L`, `--min-run 12`, `--mercy 3`, `--head-dropout 4`,
boundary-anchored, dedup on), 102 datasets, tolerance 0:

| | mates | of scored | of all |
|---|---|---|---|
| exact | 13,914,963 | **77.88%** | 68.21% |
| overtrim | 3,574,241 | 20.00% | 17.52% |
| undertrim | 378,184 | 2.12% | 1.85% |
| dropped (dedup + filters) | 2,532,612 | — | 12.41% |
| scored | 17,867,388 | | |

Total mates 20,400,000, unaligned 0, tolerance 0.

Here the table does move — by exactly one mechanism. 39 of 102 rows differ from
`fix2026`, and in every one of them the movement is mates leaving the undertrim
bucket: 2,143 mates net, of which 1,690 land in exact, 433 in overtrim and 20 in
dropped (1,690 + 433 + 20 = 2,143, the identity holds to the mate). No row moves
undertrim the other way. That is the two-sided polyX fix (`d6a866e`): a read whose
5′ cut had already spent the trim budget kept its 3′ artefact whole, and is now
trimmed. Largest single row: 560 mates of its 200,000 (`rna_seq/DRR002364`).

Net effect at the headline: exact/scored 77.87% → 77.88%,
**+0.0095 pp** unrounded (prints as +0.01 pp at two decimals); undertrim
−0.0120 pp, overtrim +0.0024 pp, dropped +0.0001 pp — each delta taken on its own
row's denominator (scored for the three buckets, all 20,400,000 mates for dropped).

## Comparison against `fix2026`

| metric | fix2026 | wave2 | delta |
|---|---|---|---|
| closed-loop PASS | 183/306 (59.8%) | 183/306 (59.8%) | 0 rows changed verdict |
| datasets all variants pass | 60/102 | 60/102 | 0 |
| median g_frac_max_diff | 0.0092 | 0.0092 | 0 |
| fail gates (g_frac / mean_q / gc) | 121 / 8 / 4 | 121 / 8 / 4 | 0 |
| closed_loop.tsv sha256 | `2076d6a075e1…` | `2076d6a075e1…` | identical file |
| trim exact of scored | 77.87% | 77.88% | +0.0095 pp |
| trim overtrim of scored | 20.00% | 20.00% | +0.0024 pp |
| trim undertrim of scored | 2.13% | 2.12% | -0.0120 pp |
| trim dropped of all | 12.41% | 12.41% | 0.0001 pp |
| scored mates | 17,867,408 | 17,867,388 | -20 |

## What did and did not move

- Closed loop: **nothing** — byte-identical file across all 25 columns, both
  against `fix2026` and against `gt2026`. Regression check, not a new result.
- Trim: one movement, −2,143 undertrim mates redistributed as +1,690 exact / +433
  overtrim / +20 dropped, i.e. +0.0095 pp exact-of-scored. Real, attributable,
  and small: it cannot be quoted as an accuracy *gain* beyond 0.01 pp.
- Dropped rate is unchanged at 12.41% of all mates; the 20 extra dropped mates are
  part of the same reclassification, not a change in filter behaviour.

## Note on the driver's summary line

The driver prints `closed_loop: 0/306 comparisons pass (0.0%)` and
`0/102 datasets pass all arms`. Those lines are wrong: the inline Python compares
`verdict == "pass"` (lowercase) while `closed_loop.py` writes `PASS`. The TSV is
correct; every count here is computed from `closed_loop.tsv`. The bug is still
present in `rerun_closed_loop_gt2026.sh` (same as the `fix2026` wave) and should
be fixed there, not worked around in prose.

Raw counts live in `summary.json` next to this file.
