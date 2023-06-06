#!/bin/bash
#SBATCH --time=4:0:0
#SBATCH --account=def-agiang01
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name='run_vargrid'
./inmap grid steady --config=cmd/inmap/configExampleGEMMACH.toml
#./inmap preproc --config=cmd/inmap/config_test_geos.toml