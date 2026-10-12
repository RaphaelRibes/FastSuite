#!/bin/bash
#SBATCH --job-name=build-profile
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-normal
#SBATCH --output=logs/build_profile_%j.out
#SBATCH --error=logs/build_profile_%j.err

# Cluster-side: no $HOME here, so the root is the Io scratch checkout
# (override with FASTSUITE_ROOT if the checkout lives elsewhere).
FASTSUITE_ROOT="${FASTSUITE_ROOT:-/scratch/users/$USER/FastSuite}"
cd "$FASTSUITE_ROOT/tools/fastsim"
mkdir -p logs profiles
python3 scripts/build_profile.py /scratch/users/$USER/mandril-polyg-2026/stats profiles/mandril_E.json
echo "EXIT=$?"
ls -la profiles/
