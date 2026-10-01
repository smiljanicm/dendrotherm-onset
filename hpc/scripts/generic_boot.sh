#!/bin/bash
#SBATCH --array=1-100
#SBATCH --output=hpc/logs/%x_%A_%a.out
#SBATCH --error=hpc/logs/%x_%A_%a.err
#SBATCH --time=16:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

module purge
export PATH=$PATH:$HOME/bin
eval "$(micromamba shell hook -s bash)"

micromamba activate r-env

cd $SLURM_SUBMIT_DIR

# Standardized for 100 array tasks yielding 2000 total bootstraps
BOOT_START=$(( ($SLURM_ARRAY_TASK_ID - 1) * 20 + 1 ))
BOOT_END=$(( $SLURM_ARRAY_TASK_ID * 20 ))

echo "Running bootstrap $BOOT_START to $BOOT_END for $SPECIES"

Rscript hpc/R/block_boots_hpc.R "$SPECIES" $BOOT_START $BOOT_END "$GS_DATA" "$TYPE"
Rscript hpc/R/block_boots_hpc.R "$SPECIES" $BOOT_START $BOOT_END "$WIDE_DATA" "$TYPE"
Rscript hpc/R/block_boots_hpc.R "$SPECIES" $BOOT_START $BOOT_END "$NARROW_DATA" "$TYPE"