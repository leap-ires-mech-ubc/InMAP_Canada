#!/bin/bash
#SBATCH --time=71:59:59
#SBATCH --account=def-rscholes #def-agiang01 #def-rscholes #def-agiang01 #
#SBATCH --ntasks-per-node=64
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=3984M
#SBATCH --job-name='20260310_run'
#SBATCH --array=2
./inmap run steady --config=cmd/inmap/configGEMMACH_scenario$SLURM_ARRAY_TASK_ID.toml
#./inmap preproc --config=cmd/inmap/config_test_geos.toml ./inmap run steady --config=cmd/inmap/20260226_configGEMMACH_test.toml