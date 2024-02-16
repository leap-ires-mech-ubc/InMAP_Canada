#!/bin/bash
#SBATCH --time=23:59:00
#SBATCH --account=def-rscholes #def-agiang01
#SBATCH --job-name=plotemiss
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=4G
#SBATCH --nodes=1
#SBATCH --array=2-5,9
python plot_emissions.py $SLURM_ARRAY_TASK_ID
#python emiss_to_shp.py