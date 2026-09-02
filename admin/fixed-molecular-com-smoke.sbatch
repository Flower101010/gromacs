#!/usr/bin/env bash
# Minimal V100 smoke test for the experimental GMX_FIXED_MOLECULAR_COM path.
# This script intentionally lives outside CGWorkflow and accepts an existing TPR.
# Usage: sbatch fixed-molecular-com-smoke.sbatch /path/to/input.tpr /path/to/output-dir [nsteps]
#SBATCH --job-name=fixed-com-smoke
#SBATCH --partition=gpu2-32GB
#SBATCH --account=sait
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
#SBATCH --time=00:10:00

set -euo pipefail

input_tpr="$1"
output_dir="$2"
nsteps="${3:-1000}"
source "$HOME/MeOH-H2O/cluster/site_env.sh"
gmx="$HOME/opt/src/gromacs-2026.3-fixed-com/build/bin/gmx_mpi"

mkdir -p "$output_dir"
export GMX_FIXED_MOLECULAR_COM=1
"$gmx" mdrun -s "$input_tpr" -deffnm "$output_dir/run" -nsteps "$nsteps" \
    -ntomp "$SLURM_CPUS_PER_TASK" -pin on -update cpu -nb gpu -pme gpu
