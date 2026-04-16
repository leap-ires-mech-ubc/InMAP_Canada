#!/bin/bash
#SBATCH --job-name=gemmach_lcc_fix
#SBATCH --account=def-rscholes
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-10
#SBATCH -o /home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/slurm/logs/lcc_fix_%A_%a.out
#SBATCH -e /home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/slurm/logs/lcc_fix_%A_%a.err

module load python/3.12 scipy-stack
source ~/inmapenv/bin/activate #mv home/inmapenv/bin/activate ~/inmapenv

mkdir -p /home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/slurm/logs

python /home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/gemmach_reproject_lcc.py $SLURM_ARRAY_TASK_ID
