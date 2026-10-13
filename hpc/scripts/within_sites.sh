#!/bin/bash
#SBATCH --job-name=WS
#SBATCH --array=1-100
#SBATCH --output=hpc/logs/WS_%A_%a.out
#SBATCH --error=hpc/logs/WS_%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

module purge
export PATH=$PATH:$HOME/bin
eval "$(micromamba shell hook -s bash)"

micromamba activate r-env

cd $SLURM_SUBMIT_DIR

BOOT_START=$(( ($SLURM_ARRAY_TASK_ID - 1) * 20 + 1 ))
BOOT_END=$(( $SLURM_ARRAY_TASK_ID * 20 ))

echo "Running bootstrap $BOOT_START to $BOOT_END"

Rscript hpc/R/block_boots_within_site_hpc.R Ladebow $BOOT_START $BOOT_END ./data/gs_within_site.csv Ladebow
Rscript hpc/R/block_boots_within_site_hpc.R Demmin $BOOT_START $BOOT_END ./data/gs_within_site.csv Demmin
