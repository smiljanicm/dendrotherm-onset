#!/bin/bash
#SBATCH --job-name=compile_boot
#SBATCH --output=hpc/logs/compile_%A_%a.out
#SBATCH --error=hpc/logs/compile_%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

module purge
export PATH=$PATH:$HOME/bin
eval "$(micromamba shell hook -s bash)"

micromamba activate r-env

cd $SLURM_SUBMIT_DIR

Rscript hpc/R/compile_boots.R "$1"

