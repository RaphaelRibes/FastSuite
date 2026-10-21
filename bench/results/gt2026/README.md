# Closed loop, re-derived on the post-2026-10-21 build

`closed_loop.tsv` is the same 306 comparisons (102 real ENA datasets x 3
variants) and the same four thresholds as `../phase-d/closed_loop.tsv`, but
every report behind it was regenerated:

    run   2026-10-21, desktop (16 cores), build 202990d
    driver  /mnt/Fast/fastsuite-scratch/rerun_closed_loop_gt2026.sh
            (fastcheck --json --profile -> fastsim -> fastcheck --json per
             dataset/variant; 306 jobs, 8 in parallel, 300 s)
    scorer  bench/scripts/closed_loop.py --force
    work    /mnt/Fast/fastsuite-scratch/work-gt2026 (not committed)

Why: the phase-d table came from binaries where `fastsim` ignored
`profile.gc` (emitted the `--gc` default 0.49), `fastcheck` measured insert
length with a leftmost-anchor bias, and its profile fitter estimated the polyG
and quality curves from the per-position G fraction instead of the exact
histograms it held.

Result: **183/306 (59.8%)**, 62/102 datasets on all three variants; median
per-position G-fraction difference 0.0093 against the 0.02 bound (phase-d:
142/306, 48/102, 0.0218). Gate counts: G fraction 121, mean quality 8, GC 4,
duplicate ratio 0.

The metric definitions live in `bench/scripts/closed_loop.py`; the tolerance
constants are embedded in the scorer because the aggregated file's published
thresholds come from there.
