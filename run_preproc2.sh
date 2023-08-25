#!/bin/bash
#SBATCH --time=23:59:0
#SBATCH --account=def-agiang01
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name='may1_preproc'
./inmap preproc --config=cmd/inmap/configGEMMACH_working2.toml
#./inmap preproc --config=cmd/inmap/config_test_geos.toml