#!/bin/bash
#SBATCH --time=167:59:0
#SBATCH --account=def-agiang01
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4500M
#SBATCH --job-name='run_inmap_fine'
./inmap run steady --config=cmd/inmap/configGEMMACH_testcoarser2.toml
#./inmap preproc --config=cmd/inmap/config_test_geos.toml