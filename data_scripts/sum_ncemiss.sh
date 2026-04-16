#!/bin/bash
#SBATCH --time=23:59:00
#SBATCH --account=def-agiang01
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --job-name='Emiss_Process2'
#SBATCH --array=1-9
export IGNORE_ATT_COORDINATES=1
#https://code.mpimet.mpg.de/boards/1/topics/3274
#https://code.mpimet.mpg.de/boards/2/topics/3487
module load nco cdo
#fpath=/home/tfmrodge/scratch/GEMMACH_data/data/BASEGM_2015_E010/
#prefix=BASEGM_2015_E010/home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/sum_ncemiss.sh
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
SLURM_ARRAY_TASK_ID="${SLURM_ARRAY_TASK_ID:-${1:-1}}"
set -euo pipefail
prefixes=/home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/config_gemmachemiss.txt
prefix=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $prefixes)
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the emissions year is ${prefix}."
fpath=/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/$prefix/
gridpth="/home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/gridfile_lcc.txt"
#Check for year - these two were based in 2020, which is a leap year
yrs=(6 8)
# Check if $SLURM_ARRAY_TASK_ID is in the list of COVID 2020 scenarios
if [[ " ${yrs[*]} " == *" $SLURM_ARRAY_TASK_ID "* ]]; then
    year=2020
    hrs=8784
else
    year=2017
    hrs=8760
fi

#Function to count the occurrences of a specific day in a given month and year, made with help of chat GPT
count_day_occurrences() {
    local year="$1"
    local month="$2"
    local day_idx="$3" # your files: 1=Sun,2=Mon,...,7=Sat

    # map to ISO 8601 Mon=1..Sun=7 used by `date +%u`
    local target=$((day_idx - 1))
    if (( target == 0 )); then target=7; fi

    # days in month
    local ndays
    ndays=$(cal "$month" "$year" | awk 'NF {d=$NF} END {print d}')

    local count=0
    for (( d=1; d<=ndays; d++ )); do
        if (( $(date -d "$year-$month-$d" +%u) == target )); then
            ((count++))
        fi
    done

    echo "$count"
}

for month in "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12"
do
	for i in 1 2 3 4 5 6 7
    do
        num_day=$(count_day_occurrences "$year" "$month" "$i")
        #printf $num_day
        mulcx=$(bc -l <<< "$num_day/$hrs")
        cdo -mulc,$mulcx ${fpath}${prefix}_major_${month}_${i}.nc ${fpath}${prefix}_major_${month}_w${i}.nc &
		cdo -mulc,$mulcx ${fpath}${prefix}_area_${month}_${i}.nc ${fpath}${prefix}_area_${month}_w${i}.nc
    done
    wait
    #Sum all 7*24 hours into a single file.
    ncra -O -y ttl ${fpath}${prefix}_area_${month}_w[1234567].nc ${fpath}${prefix}_area_${month}.nc
    rm ${fpath}${prefix}_area_${month}_w[1234567].nc
    ncra -O -y ttl ${fpath}${prefix}_major_${month}_w[1234567].nc ${fpath}${prefix}_major_${month}.nc
    rm ${fpath}${prefix}_major_${month}_w[1234567].nc
done
wait 
#Sum everything to one time step. All values are the weighted average, so units come back to g/s
ncra -O -y ttl -n 12,2,1 ${fpath}${prefix}_area_01.nc ${fpath}${prefix}_areasum.nc &
ncra -O -y ttl -n 12,2,1 ${fpath}${prefix}_major_01.nc ${fpath}${prefix}_majorsum.nc
wait
#Then, sum into the required values for InMap
#VOCs
ncap2 -O -s "VOC=EA2+EA3+ECRE+EISO+EARO+EMEK+ETOL+EC38+EETH+EALD+EHCH" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc &
ncap2 -O -s "VOC=EA2+EA3+ECRE+EISO+EARO+EMEK+ETOL+EC38+EETH+EALD+EHCH" ${fpath}${prefix}_areasum.nc ${fpath}${prefix}_areasum.nc
wait
#NOx
ncap2 -O -s "NOx=ENO+ENO2+EHON" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc &
ncap2 -O -s "NOx=ENO+ENO2+EHON" ${fpath}${prefix}_areasum.nc ${fpath}${prefix}_areasum.nc
wait
#SOx
ncap2 -O -s "SOx=ESO2+ESO4" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc &
ncap2 -O -s "SOx=ESO2+ESO4" ${fpath}${prefix}_areasum.nc ${fpath}${prefix}_areasum.nc
wait
#PM25 - No fugitive (FXX1) for Major
ncap2 -O -s "PM25=EAM1+ENT1+EEC1+EPC1+ECM1+ESU1" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc &
ncap2 -O -s "PM25=EAM1+ENT1+EEC1+EPC1+ECM1+ESU1+FAM1+FNT1+FEC1+FPC1+FCM1+FSU1" ${fpath}${prefix}_areasum.nc ${fpath}${prefix}_areasum.nc
wait
#NH3
ncap2 -O -s "NH3=ENH3" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc &
ncap2 -O -s "NH3=ENH3" ${fpath}${prefix}_areasum.nc ${fpath}${prefix}_areasum.nc
wait
#Just for major. Need to confirm these values, if they should be set now or later. 
ncap2 -O -s "Height=0*NH3" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc
wait
ncap2 -O -s "Diam=0*NH3" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc
wait
ncap2 -O -s "Temp=0*NH3" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc
wait
ncap2 -O -s "Velocity=0*NH3" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc
wait
#Remap to chosen (LCC) grid, no need to convert to netCDF 3 as we will be converting to shapefile later
cdo -O remapbil,"$gridpth" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_major_lcc.nc &
cdo -O remapbil,"$gridpth" ${fpath}${prefix}_areasum.nc ${fpath}${prefix}_area_lcc.nc
wait
#We are going to remove the height variable as it causes problems downstream. 
ncwa -O -a height ${fpath}${prefix}_major_lcc.nc ${fpath}${prefix}_major_lcc.nc &
ncwa -O -a height ${fpath}${prefix}_area_lcc.nc ${fpath}${prefix}_area_lcc.nc
wait
ncks -O -C -x -v height ${fpath}${prefix}_major_lcc.nc ${fpath}${prefix}_major_lcc.nc &
ncks -O -C -x -v height ${fpath}${prefix}_area_lcc.nc ${fpath}${prefix}_area_lcc.nc
wait
#ncks -O -6 data/BASEGM_2015_E010/BASEGM_2015_E010_major_lcc.nc data/BASEGM_2015_E010/BASEGM_2015_E010_major_lcc.nc
# base=data/BASEGM_2015_E010/BASEGM_2015_E010
# #Next, we are going to find the difference between baseGM and this file (scenario-base)
# ncdiff -O ${fpath}${prefix}_major_lcc.nc ${base}_major_lcc.nc ${fpath}${prefix}_majordiff_lcc.nc
# ncdiff -O ${fpath}${prefix}_area_lcc.nc ${base}_area_lcc.nc ${fpath}${prefix}_areadiff_lcc.nc
#ncdiff -O data/BAU_2020_E108/BAU_2020_E108_area_lcc.nc data/BASEGM_2015_E010/BASEGM_2015_E010_area_lcc.nc /home/tfmrodge/scratch/GEMMACH_data/data/BAU_2020_E108/BAU_2020_E108_areadiff_lcc.nc