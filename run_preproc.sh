#!/bin/bash
#SBATCH --time=0:05:0
#SBATCH --account=def-agiang01
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=512
#SBATCH --ntasks=1
#SBATCH --job-name='test_preproc'
./inmap preproc --config=cmd/inmap/config_test_GEMMACH.toml
#./inmap preproc --config=cmd/inmap/config_test_geos.toml