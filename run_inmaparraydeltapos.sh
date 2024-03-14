#!/bin/bash
#SBATCH --time=71:59:0
#SBATCH --account=def-rscholes #def-agiang01 #def-agiang01 #def-rscholes #
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4600M
#SBATCH --job-name='pos'
#SBATCH --array=2-9
./inmap run steady --config=cmd/inmap/configGEMMACH_deltaposscenario$SLURM_ARRAY_TASK_ID.toml
#./inmap preproc --config=cmd/inmap/config_test_geos.toml./inmap run steady --config=cmd/inmap/configGEMMACH_scenario1.toml