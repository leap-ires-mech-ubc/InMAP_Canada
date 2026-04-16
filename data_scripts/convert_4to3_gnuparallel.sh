#!/bin/bash
#SBATCH --time=2:59:59
#SBATCH --job-name=convert_4to3
#SBATCH --account=def-agiang01
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=10G
#This script should be run from the same folder as process_date.sh and gridfile.lcc
#Make sure you correctly point to the directories you want to access, too - the default is data/[directory]
#Then, run these commands if you haven't before in the session:
#module load nco cdo
#cd GEMMACH_data (or wherever you have all the folders and files)

#salloc --account=def-agiang01 --time=03:00:00 --job-name=test_convert --ntasks=1 --cpus-per-task=8 --mem-per-cpu=10G

###############################################################
# USER SETTINGS — EDIT THESE LIKE YOU EDIT DATES IN ORIGINAL
###############################################################
module load nco cdo
### Output roots
OUT_GEMMACH_ROOT="/home/tfmrodge/scratch/GEMMACH_data/data/nc3/BASEGM_2015_017"
OUT_RDPS_ROOT="/home/tfmrodge/scratch/GEMMACH_data/data/nc3/RDPS_QC"

### Grid file
GRIDFILE="/home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/gridfile_lcc.txt"

### RDPS remap (1 = yes, 0 = no)
DO_RDPS_REMAP=1

### Date ranges (exact same style you currently use)
year=2019
start_month=05
end_month=05
day_s=01
day_e=31
hr_s=00
hr_e=23

### Input roots (NO ARGS — set here) #/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/BASEGM_2015_017/03_Mar/2019032000_001.nc
IN_GEMMACH_ROOT="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/BASEGM_2015_017/05_May"
IN_RDPS_ROOT="/home/tfmrodge/scratch/GEMMACH_data/data/RDPS_QC"

###############################################################
# ENV EXPORTS — used by process_date.sh
###############################################################
export IN_GEMMACH_ROOT IN_RDPS_ROOT OUT_GEMMACH_ROOT OUT_RDPS_ROOT GRIDFILE DO_RDPS_REMAP

###############################################################
# CPU COUNT (same logic as your original)
###############################################################
if [ -z "$SLURM_CPUS_PER_TASK" ]; then
    CORES=4
else
    CORES="$SLURM_CPUS_PER_TASK"
fi

#Change to directory with process_date.sh and gridfile_lcc.txt
cd /home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/
chmod +x process_date.sh
###############################################################
# Start timer
###############################################################
start=$(date +%s)

echo "Using $CORES cores"
echo "Processing year $year months $start_month–$end_month days $day_s–$day_e hours $hr_s–$hr_e"

###############################################################
# RUN — EXACT SAME GNU PARALLEL STRUCTURE AS YOUR SCRIPT
###############################################################

parallel -j $CORES ./process_date.sh ::: $year ::: $(seq -w $start_month $end_month) ::: $(seq -w "$day_s" "$day_e") ::: $(seq -w "$hr_s" "$hr_e")

echo "Finished processing all dates."

###############################################################
# Timer end
###############################################################
end=$(date +%s)
runtime=$(( end - start ))

echo "Elapsed time: $runtime seconds"