#!/bin/bash
#SBATCH --time=23:0:0
#SBATCH --account=def-agiang01
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name='run_preproc'
./inmap preproc --config=cmd/inmap/configGEMMACH_Feb.toml
#./inmap preproc --config=cmd/inmap/config_test_geos.toml