#!/bin/bash
#SBATCH --job-name=novaseq-validate2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-normal
#SBATCH --output=logs/validate2_%j.out
#SBATCH --error=logs/validate2_%j.err

cd /scratch/users/ribesr/novaseq-sim
python3 scripts/validate_sim.py \
  profiles/mandril_E.json \
  validation/sim_mandril_R1.fastq.gz \
  validation/sim_mandril_R2.fastq.gz \
  validation/sim_mandril.gt.tsv \
  validation/closed_loop
echo "EXIT=$?"
ls -la validation/
