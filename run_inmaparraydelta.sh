#!/bin/bash
#SBATCH --time=71:59:0
#SBATCH --account=def-agiang01 #def-rscholes #def-agiang01 #def-rscholes #
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4600M
#SBATCH --job-name='delarea_40k'
#SBATCH --array=6-8
./inmap run steady --config=cmd/inmap/configGEMMACH_deltascenario$SLURM_ARRAY_TASK_ID.toml
#./inmap preproc --config=cmd/inmap/config_test_geos.toml./inmap run steady --config=cmd/inmap/configGEMMACH_scenario1.toml