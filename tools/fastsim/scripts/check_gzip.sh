#!/bin/bash
# Check gzip state after the timeout.
cd /home/raphael/hermes_scratch/novaseq-sim/local_work/sim
ls -la sim_R1.fastq.gz sim_R2.fastq.gz
echo "--- testing ---"
gzip -t sim_R1.fastq.gz 2>&1 && echo "R1 OK" || echo "R1 BAD"
gzip -t sim_R2.fastq.gz 2>&1 && echo "R2 OK" || echo "R2 BAD"
