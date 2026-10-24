# CLI reference

```
fastsim — NovaSeq-like paired-end read simulator with ground truth

USAGE:
    fastsim --pairs N --profile FILE --out PREFIX [OPTIONS]

REQUIRED:
    --pairs N              number of read pairs to simulate
    --profile FILE         profile JSON (artefact CDFs, quality curves, adapter)
    -o, --out PREFIX       output prefix (writes PREFIX_R1.fastq.gz,
                           PREFIX_R2.fastq.gz, PREFIX.gt.tsv)

OPTIONS:
    --seed S               RNG seed (default 42)
    --read-len L           read length (default: the profile's read_len; an
                           explicit value must match it)
    --insert-mean M        insert length mean (default 300; ignored when the
                           profile carries a measured insert CDF)
    --insert-sd S          insert length sd (default 50; ignored likewise)
    --insert-min M         insert length floor (default min(150, read_len);
                           ignored likewise)
    --gc G                 GC fraction for random DNA (default: the profile's
                           gc, falling back to 0.49)
    --reference FILE       sample inserts from a reference FASTA instead of
                           random DNA
    --adapter-seq SEQ      adapter sequence for readthrough (default: from the
                           profile, Nextera CTGTCTCTTATACACATCT)
    --no-adapter           disable adapter readthrough
    --no-gzip              write uncompressed FASTQ
    --gz-level L           DEFLATE level for the gzip outputs, 1..=9
                           (default 6)
    --dup-rate R           PCR duplicate fraction (default 0.02; a profile's
                           measured dup_rate overrides it)
    --dup-pool N           ring-buffer capacity for the duplicate pool
                           (default 131072)
    --drop-rate R          fraction of pairs whose fragment is too short to
                           survive trimming (default 0.05)
    --drop-min-len L       insert ceiling for the forced-drop class
                           (default 100)
    --single               one read per fragment (no R2)
    --json PATH            write the FastSuite report (schema v1)
    --stats                add the instrument counters to the --json report
    -j, --threads N        accepted for FastSuite CLI compatibility; the
                           simulator is single-threaded and warns (default 1,
                           or FASTSUITE_THREADS)
    -f, --force            overwrite existing output files (default: refuse)
    -q, --quiet            suppress the end-of-run summary line
    -V, --version          print the version
    -h, --help             this help
```

## Notes

- `--read-len` must equal the profile's `read_len` (the profile's per-position
  quality curves are that length). With no flag the profile's value is used.
- `--insert-mean`/`--insert-sd`/`--insert-min` are used only when the profile
  has no measured `insert` CDF; `--gc` and `--dup-rate` likewise defer to the
  profile's measured `gc`/`dup_rate` when present.
- `--reference` loads a plain FASTA (contigs concatenated per header); inserts
  are sampled from random contigs/positions, skipping positions that would run
  past a contig end.
- `--adapter-seq ""` (via `--no-adapter`) disables readthrough; with no flag,
  the profile's adapter is used.
- `-o/--out` also accepts `--out1` and `--out-prefix` (visible aliases).
- `--single` writes no R2 and a single-end GT header (`pair_id insert r1_h
  r1_t r1_adapter full_flag is_dup r1_n r1_full`).
- Output is deterministic for a fixed (seed, profile, args).
