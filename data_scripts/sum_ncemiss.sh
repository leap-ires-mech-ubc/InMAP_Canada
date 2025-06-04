#!/bin/bash
#SBATCH --time=01:30:00
#SBATCH --account=def-agiang01
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=3
#SBATCH --job-name='Emiss_Process2'
#SBATCH --array=1-9
export IGNORE_ATT_COORDINATES=1
#https://code.mpimet.mpg.de/boards/1/topics/3274
#https://code.mpimet.mpg.de/boards/2/topics/3487

#fpath=/home/tfmrodge/scratch/GEMMACH_data/data/BASEGM_2015_E010/
#prefix=BASEGM_2015_E010
prefixes=/home/tfmrodge/scratch/GEMMACH_data/config_jobs.txt
prefix=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $prefixes)
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the emissions year is ${prefix}."
fpath=/home/tfmrodge/scratch/GEMMACH_data/data/$prefix/
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
    local day_to_count="$3"
    day_to_count=$(bc -l <<< "$day_to_count-1")
    if [ "$day_to_count" -eq 0 ]; then
        day_to_count=7
    fi
    
    # Get the day of the week for the first day of the month
    local first_day_of_month=$(date -d "$year-$month-01" +%u)
    
    # Initialize a counter for the desired day
    local count=0
    
    # Iterate through the days of the month and count occurrences of the desired day
    for (( day = 1; day <= $(cal $month $year | awk 'NF {DAYS = $NF} END {print DAYS}'); day++ )); do
        local day_of_week=$(date -d "$year-$month-$day" +%u)
        if [ "$day_of_week" -eq "$day_to_count" ]; then
            ((count++))
        fi
    done
    
    return $day_to_count
}

for month in "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12"
do
	for i in 1 2 3 4 5 6 7
    do
        count_day_occurrences $year $month $i
        num_day=$?
        #printf $num_day
        mulcx=$(bc -l <<< $num_day'/$hrs')
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
#NOx
ncap2 -O -s "NOx=ENO+ENO2+EHON" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc &
ncap2 -O -s "NOx=ENO+ENO2+EHON" ${fpath}${prefix}_areasum.nc ${fpath}${prefix}_areasum.nc
#SOx
ncap2 -O -s "SOx=ESO2+ESO4" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc &
ncap2 -O -s "SOx=ESO2+ESO4" ${fpath}${prefix}_areasum.nc ${fpath}${prefix}_areasum.nc
#PM25
ncap2 -O -s "PM25=EPC1+ECM1+ESU1" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc &
ncap2 -O -s "PM25=EAM1+ENT1+EEC1+EPC1+ECM1+ESU1+FAM1+FNT1+FEC1+FPC1+FCM1+FSU1" ${fpath}${prefix}_areasum.nc ${fpath}${prefix}_areasum.nc
#NH3
ncap2 -O -s "NH3=ENH3" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc &
ncap2 -O -s "NH3=ENH3" ${fpath}${prefix}_areasum.nc ${fpath}${prefix}_areasum.nc
#Just for major. Need to confirm these values, if they should be set now or later. 
ncap2 -O -s "Height=0*NH3" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc
ncap2 -O -s "Diam=0*NH3" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc
ncap2 -O -s "Temp=0*NH3" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc
ncap2 -O -s "Velocity=0*NH3" ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_majorsum.nc
#Remap to chosen (LCC) grid, no need to convert to netCDF 3 as we will be converting to shapefile later
cdo -O remapbil,gridfile_lcc.txt ${fpath}${prefix}_majorsum.nc ${fpath}${prefix}_major_lcc.nc &
cdo -O remapbil,gridfile_lcc.txt ${fpath}${prefix}_areasum.nc ${fpath}${prefix}_area_lcc.nc
#We are going to remove the height variable as it causes problems downstream. 
ncwa -O -a height ${fpath}${prefix}_major_lcc.nc ${fpath}${prefix}_major_lcc.nc &
ncwa -O -a height ${fpath}${prefix}_area_lcc.nc ${fpath}${prefix}_area_lcc.nc
wait
ncks -O -x -v height ${fpath}${prefix}_major_lcc.nc ${fpath}${prefix}_major_lcc.nc &
ncks -O -x -v height ${fpath}${prefix}_area_lcc.nc ${fpath}${prefix}_area_lcc.nc

#ncks -O -6 data/BASEGM_2015_E010/BASEGM_2015_E010_major_lcc.nc data/BASEGM_2015_E010/BASEGM_2015_E010_major_lcc.nc
base=data/BASEGM_2015_E010/BASEGM_2015_E010
#Next, we are going to find the difference between baseGM and this file (scenario-base)
ncdiff -O ${fpath}${prefix}_major_lcc.nc ${base}_major_lcc.nc ${fpath}${prefix}_majordiff_lcc.nc
ncdiff -O ${fpath}${prefix}_area_lcc.nc ${base}_area_lcc.nc ${fpath}${prefix}_areadiff_lcc.nc
#ncdiff -O data/BAU_2020_E108/BAU_2020_E108_area_lcc.nc data/BASEGM_2015_E010/BASEGM_2015_E010_area_lcc.nc /home/tfmrodge/scratch/GEMMACH_data/data/BAU_2020_E108/BAU_2020_E108_areadiff_lcc.nc