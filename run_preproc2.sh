#!/bin/bash
#SBATCH --time=167:59:0
#SBATCH --account=def-agiang01
#SBATCH --ntasks-per-node=64
#SBATCH --nodes=1
#SBATCH --mem=500G
#SBATCH --job-name='Preprocs'
#SBATCH --array=5-6
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    ID=1
else
    ID=$SLURM_ARRAY_TASK_ID
fi
echo $ID
./inmap preproc --config=cmd/inmap/configGEMMACH_preproc$ID.toml
#./inmap preproc --config=cmd/inmap/config_test_geos.toml