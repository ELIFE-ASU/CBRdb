#!/bin/bash --login
#SBATCH --job-name=trace_atom
#SBATCH --array=0-100
#SBATCH --time=0-4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=htc
#SBATCH --qos=public
#SBATCH -o /home/mshahjah/AtomTracking/Logs/out/slurm.%A_%a.out
#SBATCH -e /home/mshahjah/AtomTracking/Logs/err/slurm.%A_%a.err
#SBATCH --mem=100G
module conda activate rxn_env
srun $HOME/.conda/envs/rxn_env/bin/python3.8 /home/mshahjah/AtomTracking/AtomTracking.py --row_index ${SLURM_ARRAY_TASK_ID} --dynamic_cofactors 