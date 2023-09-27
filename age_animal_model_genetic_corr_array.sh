#!/bin/bash
#SBATCH --job-name=genetic_corr
#SBATCH --output=slurm_outs/genetic_corr_slurm_out_%A_%a
#SBATCH --error=slurm_outs/genetic_corr_slurm_err_%A_%a
#SBATCH --array=1-6
#SBATCH --mem-per-cpu=9G
#SBATCH --time=39-00:00:00

# Note: do this before running: mkdir slurm_outs stdouts
module load R/4.0.3

Rscript --vanilla \
Oct2022_genetic_corr.R \
${SLURM_ARRAY_TASK_ID} \
> stdouts/genetic_corr_stdout_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} \
2>&1
