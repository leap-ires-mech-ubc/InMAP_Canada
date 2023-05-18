#!/bin/bash
#SBATCH --time=1:0:0
#SBATCH --account=def-agiang01
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1024
#SBATCH --ntasks=1
#SBATCH --job-name='test_preproc'
./inmap preproc --config=cmd/inmap/configExampleGEMMACH.toml
#./inmap preproc --config=cmd/inmap/config_test_geos.toml