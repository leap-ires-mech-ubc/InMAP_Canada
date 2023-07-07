#!/bin/bash
#SBATCH --time=71:59:0
#SBATCH --account=def-agiang01
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4500M
#SBATCH --job-name='run_vargrid_NoMay'
./inmap grid --config=cmd/inmap/configGEMMACH_testcoarser.toml
#./inmap preproc --config=cmd/inmap/config_test_geos.toml