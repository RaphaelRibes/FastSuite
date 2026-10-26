# The 2026-10-26 post-fix re-run (desktop, `cachyos-x8664`)

Re-derived after the fix wave that closed the 2026-10-26 review's reproduced
defects (`3c9cc28`), because two of them change the data behind published
numbers:

* fastsim's GT `r1_n`/`r2_n` now carry a template `N` (it was netted to zero and
  booked as adapter padding), and `full_flag` counts both mates;
* fastcheck's profile fitter no longer emits `gc: 0.0` for an all-readthrough
  R1, and its insert sampler reports what it actually scanned;
* fastprep trims a 3′ polyX run that lands at the trim-budget boundary instead
  of emitting it whole.

## Commands (desktop)

```bash
# closed loop: 102 datasets x 3 variants, fastcheck -> fastsim -> fastcheck
cd /mnt/Fast/fastsuite-scratch
WORKROOT=/mnt/Fast/fastsuite-scratch/work-fix2026 \
RESULTS=/mnt/Fast/fastsuite-scratch/results-fix2026 \
BIN=/home/raphael/FastSuite/tools JOBS=8 THREADS=2 SIM_PAIRS=100000 \
  bash rerun_closed_loop_gt2026.sh          # 306 sim reports in 294 s

# trim accuracy against the same GT, shipped default, tolerance 0
/usr/bin/python3 trim_score_fix2026.py \
  --work /mnt/Fast/fastsuite-scratch/work-fix2026 \
  --fastprep /home/raphael/FastSuite/tools/fastprep/target/release/fastprep \
  --out results-fix2026/trim-accuracy-shipped.tsv
```

## Closed loop — `closed_loop.tsv`

**183/306 comparisons pass (59.8%)**, 60/102 datasets pass every variant,
median per-position G-fraction difference **0.0092** against the 0.02 bound.
Per domain: human gut 51/51, generalised 47/51, bacterial isolates 40/51,
RNA-seq 21/51, viromes 21/51, metagenomes 3/51. Gates over the 123 failing rows:
G fraction 121, mean quality 8, GC 4.

The first 17 columns are **byte-identical** to `bench/results/gt2026/` (diff of
the two files over those columns is empty): the fix wave repairs report and GT
integrity, and does not move any of the four scored gates on this corpus. That
is the honest reading — a re-run that reproduced the old numbers is a
regression check, not a new result.

The driver's own summary line prints `0/306` because its inline Python compares
`verdict == "pass"` while `closed_loop.py` writes `PASS`; the TSV is correct and
the counts above are computed from it.

## Trim accuracy — `trim-accuracy-shipped.tsv`

Shipped default (`0.67L`, `--min-run 12`, `--mercy 3`, `--head-dropout 4`,
boundary-anchored, dedup on), 102 datasets, tolerance 0:

| | mates | of scored | of all |
|---|---|---|---|
| exact | 13,913,273 | **77.87%** | 68.20% |
| overtrim | 3,573,808 | 20.00% | 17.52% |
| undertrim | 380,327 | 2.13% | 1.86% |
| dropped (dedup + filters) | 2,532,592 | — | 12.41% |
| scored | 17,867,408 | | |

The previous campaign records 75.5 / 22.5 / 2.0 for the same configuration
(`docs/TRIM-ACCURACY.md`): the budget-boundary fix moves ~2.4 pp of mates from
overtrim to exact, which is the defect it closes. `docs/TRIM-ACCURACY.md` cites
`bench/scripts/overtrim_vs_fastp.py` as the producer of that table; **that file
does not exist in the repository**, so the fastp columns were not re-derived
here — only the fastprep row, with the scorer committed as
`bench/scripts/trim_score.py`.
