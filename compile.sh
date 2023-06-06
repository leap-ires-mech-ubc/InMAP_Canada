#!/bin/bash
#SBATCH --time=0:30:0
#SBATCH --account=def-agiang01
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1024
#SBATCH --job-name='compile_inmap'
#./inmap preproc --config=cmd/inmap/configExampleGEMMACH.toml
module load go/1.18.3
GO111MODULE=on go build ./cmd/inmap