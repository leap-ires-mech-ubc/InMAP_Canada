#!/bin/bash
#SBATCH --time=167:59:0
#SBATCH --account=def-rscholes #def-agiang01 #def-rscholes #def-agiang01
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4600M
#SBATCH --job-name='20km'
./inmap run steady --config=cmd/inmap/configGEMMACH_testcoarser3.toml
#./inmap preproc --config=cmd/inmap/config_test_geos.toml