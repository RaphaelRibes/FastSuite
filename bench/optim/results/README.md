# Campaign outputs

* `2026-10-20-targeted.tsv` — the campaign the dossier and the six papers quote. Fifteen rows,
  strictly serial, `cpu_s / wall_s` 0.984–0.996 on every row. Reproduce with
  `bash run_bench.sh` (defaults: `config-papers.yaml`, `--cores 1`).
* `2026-10-20-discarded-cores8.tar.gz` — the 91 cell files of the first campaign, kept as evidence
  rather than as data. It ran the full matrix at `--cores 8` while every cell is pinned with
  `taskset`, so cells queued on the same core: `-j 1` rows show `cpu_s / wall_s` between 0.13 and
  0.28 and `-j 4` rows up to 2.96. Those rows are unreadable and are cited nowhere except as the
  example of a measurement that looks complete and is not.

The 20 GB of tool output those cells produced has been deleted; the TSVs are what carries the
evidence, and the outputs are regenerable.
