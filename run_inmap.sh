#!/bin/bash
#SBATCH --time=71:59:0
#SBATCH --account=def-agiang01 #def-rscholes #def-agiang01
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4600M
#SBATCH --job-name='NewDataNewMe'
./inmap run steady --config=cmd/inmap/configGEMMACH_testsavegrid.toml
#./inmap preproc --config=cmd/inmap/config_test_geos.toml
#/home/tfmrodge/projects/def-agiang01/tfmrodge/InMAP_Canada/cmd/inmap/configGEMMACH