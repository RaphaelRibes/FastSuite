# CLI reference

```
fastsim — NovaSeq-like paired-end read simulator with ground truth

USAGE:
    fastsim --pairs N --profile FILE --out PREFIX [OPTIONS]

REQUIRED:
    --pairs N              number of read pairs to simulate
    --profile FILE         profile JSON (artefact CDFs, quality curves, adapter)
    --out PREFIX           output prefix (writes PREFIX_R1.fastq.gz,
                           PREFIX_R2.fastq.gz, PREFIX.gt.tsv)

OPTIONS:
    --seed S               RNG seed (default 42)
    --read-len L           read length (default 150)
    --insert-mean M        insert length mean (default 300)
    --insert-sd S          insert length sd (default 50)
    --insert-min M         insert length floor (default 150)
    --gc G                 GC fraction for random DNA (default 0.49)
    --reference FILE       sample inserts from a reference FASTA instead of
                           random DNA
    --adapter-seq SEQ      adapter sequence for readthrough (default: Nextera
                           CTGTCTCTTATACACATCT from profile)
    --no-adapter           disable adapter readthrough
    --no-gzip              write uncompressed FASTQ
    -h, --help             this help
```

## Notes

- `--read-len` must equal the profile's `read_len` (the profile's per-position
  quality curves are that length).
- `--reference` loads a plain FASTA (contigs concatenated per header); inserts
  are sampled from random contigs/positions, skipping positions that would run
  past a contig end.
- `--adapter-seq ""` (via `--no-adapter`) disables readthrough; with no flag,
  the profile's adapter is used.
- Output is deterministic for a fixed (seed, profile, args).
