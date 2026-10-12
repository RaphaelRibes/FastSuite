#!/bin/bash
#SBATCH --job-name=fastsim-validate2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-normal
#SBATCH --output=logs/validate2_%j.out
#SBATCH --error=logs/validate2_%j.err

# Cluster-side: no $HOME here, so the root is the Io scratch checkout
# (override with FASTSUITE_ROOT if the checkout lives elsewhere).
FASTSUITE_ROOT="${FASTSUITE_ROOT:-/scratch/users/$USER/FastSuite}"
cd "$FASTSUITE_ROOT/tools/fastsim"
python3 scripts/validate_sim.py \
  profiles/mandril_E.json \
  validation/sim_mandril_R1.fastq.gz \
  validation/sim_mandril_R2.fastq.gz \
  validation/sim_mandril.gt.tsv \
  validation/closed_loop
echo "EXIT=$?"
ls -la validation/
