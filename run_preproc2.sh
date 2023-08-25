#!/bin/bash
#SBATCH --time=167:59:0
#SBATCH --account=def-agiang01
#SBATCH --ntasks-per-node=64
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name='Preprocs'
#SBATCH --array=1-6
./inmap preproc --config=cmd/inmap/configGEMMACH_working_$SLURM_ARRAY_TASK_ID.toml
#./inmap preproc --config=cmd/inmap/config_test_geos.toml