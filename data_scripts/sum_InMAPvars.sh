#!/bin/bash
#SBATCH --time=01:30:00
#SBATCH --account=def-agiang01
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=3
#SBATCH --array=1-9
#SBATCH --job-name='sum_conc'

#module load nco cdo
#prefix=BASEGM_2015_017
#Sum all files in a folder
SLURM_ARRAY_TASK_ID=2
prefixes=/home/tfmrodge/projects/def-agiang01/tfmrodge/InMAP_Canada/data_scripts/config_gemmachouts.txt
prefix=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $prefixes)
if [ "$prefix" = "BASEGM_2015_017" ]; then
    prefix=BASEGM_surface
fi
echo "Array task ${SLURM_ARRAY_TASK_ID}, averaging ${prefix}."
#Set the path to the folder
fpath=/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/$prefix/
#Limit number of files
MAXFILES=""
if [ -n "$MAXFILES" ]; then
    echo "Limiting to first $MAXFILES files for testing"
    files=$(ls "${fpath}"*.nc | sort | head -n "$MAXFILES")
    #echo $files
    outfile="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/scenario_sums/20250704_${prefix}_test.nc"
else
    echo "Using all files"
    files=$(ls "${fpath}"*.nc | sort)
    outfile="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/scenario_sums/20251017_${prefix}.nc"
fi

#Combine variables
if [ "$prefix" = "BASEGM_surface" ]; then
    cdo -L -sellevel,0.99875
fi

ncap2 -O -s BASEPM25=AF
ncap2 -O -s BASEPNO3=(TNI1)*RHO
ncap2 -O -s BASEPNH4=(TAM1)*RHO
ncap2 -O -s BASEPSO4=(TSU1)*RHO ${fpath}*.nc ${outfile}.nc
ncap2 -O -s BASESOA=(TOC1)*RHO ${fpath}*.nc ${outfile}.nc





#ncra ${fpath}*.nc ${outfile}.nc
