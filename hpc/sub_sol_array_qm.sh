#!/bin/bash
#SBATCH --job-name=cbrdb
##SBATCH --output=output_%A_%a.out   # Output file for each array task
##SBATCH --error=error_%A_%a.err     # Error file for each array task
#SBATCH --array=0-4                  # Array of 10 jobs (adjust to your range) 15200 17023
#SBATCH --time=0-04:00:00            # Time limit
#SBATCH --ntasks=1                   # Number of tasks (1 per array task)
#SBATCH --cpus-per-task=8            # Number of CPU cores per task
#SBATCH --mem=8G                     # Memory per task 4 or 16
#SBATCH --partition=htc              # general lightwork highmem
#SBATCH --qos=public                 # private public

cd $SLURM_SUBMIT_DIR

module load mamba/latest
source activate cbrdb

# Run the Python script, passing the SLURM_ARRAY_TASK_ID to the script
srun python3 calc_qm_spectrum.py "${SLURM_ARRAY_TASK_ID}" CBRdb_C.csv
