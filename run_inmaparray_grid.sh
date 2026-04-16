#!/bin/bash
#SBATCH --time=71:59:0
#SBATCH --account=def-agiang01 #def-rscholes #def-agiang01 #def-rscholes #
#SBATCH --ntasks-per-node=64
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=3984M
#SBATCH --job-name='Staticgrid_20251002'
#SBATCH --array=6
./inmap run steady --config=cmd/inmap/configGEMMACH_scenariogrid$SLURM_ARRAY_TASK_ID.toml
#./inmap preproc --config=cmd/inmap/config_test_geos.toml./inmap run steady --config=cmd/inmap/configGEMMACH_scenario1.toml