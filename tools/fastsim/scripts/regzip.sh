#!/bin/bash
# Re-gzip the plain sim FASTQs (overwrite truncated gz), keep originals.
cd /home/raphael/hermes_scratch/novaseq-sim/local_work/sim
gzip -kf sim_R1.fastq
gzip -kf sim_R2.fastq
echo "GZIP_DONE"
ls -la sim_R1.fastq.gz sim_R2.fastq.gz
