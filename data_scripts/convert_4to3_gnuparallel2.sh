#!/bin/bash
#SBATCH --time=12:30:00
#SBATCH --job-name=convert_4to3
#SBATCH --account=def-agiang01 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem-per-cpu=10G
#This script should be run from the same folder as process_date1.sh,etc. and gridfile.lcc
#Make sure you correctly point to the directories you want to access, too - the default is data/[directory]
#Then, run these commands if you haven't before in the session:
#module load nco cdo
#cd GEMMACH_data (or wherever you have all the folders and files)
#You only need to run these once, I believe
#chmod +x process_date1.sh
#chmod +x process_date2.sh
#chmod +x process_date3.sh
# salloc --account=def-agiang01 --time=01:00:00 --job-name=test_convert --ntasks=1 --cpus-per-task=4 --mem-per-cpu=10G

#9s with 4 cores, 28s with 1 core
#parallel --citation # #SLURM_NTASKS_PER_NODE=4
if [ -z "$SLURM_CPUS_PER_TASK" ]; then
    # Not in Slurm → set default number of cores for local testing
    CORES=4
else
    # In Slurm → use Slurm‑assigned CPU count
    CORES="$SLURM_CPUS_PER_TASK"
fi

# Start the timer
start=$(date +%s)

# Set the range of months and days to process
start_month=03
end_month=03
year=2019
day_s=16
day_e=31
hr_s=00
hr_e=23

#source /home/tfmrodge/scratch/GEMMACH_data/process_date.sh
# Generate the commands using GNU Parallel #bash -c 'source /home/tfmrodge/scratch/GEMMACH_data/process_date.sh; process_date1 "$@"' _
#bash -c '. /home/tfmrodge/scratch/GEMMACH_data/process_date.sh; process_date1 "$@"' _
parallel -j $CORES /home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/process_date1.sh ::: $year ::: $(seq -w $start_month $end_month) ::: $(seq -w "$day_s" "$day_e") ::: $(seq -w "$hr_s" "$hr_e")
printf "done job 1"
parallel -j $CORES /home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/process_date2.sh ::: $year ::: $(seq -w $start_month $end_month) ::: $(seq -w "$day_s" "$day_e") ::: $(seq -w "$hr_s" "$hr_e")
#Remove tmp1 files
rm /home/tfmrodge/scratch/GEMMACH_data/data/nc3/BASEGM_2015_017/2019-*-tmp1
printf "done job 2"
parallel -j $CORES /home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/process_date3.sh ::: $year ::: $(seq -w $start_month $end_month) ::: $(seq -w "$day_s" "$day_e") ::: $(seq -w "$hr_s" "$hr_e")
#Remove tmp2 files
rm /home/tfmrodge/scratch/GEMMACH_data/data/nc3/BASEGM_2015_017/2019-*-tmp2
printf "done job 3"
#Reproject gem_geophy file
#cdo remapbil,gridfile_lcc.txt data/GEOPHY_VF/Gem_geophy.nc data/GEOPHY_VF/gem_geophy_lcc.nc
# Calculate the elapsed time
end=$(date +%s)
runtime=$((end - start))

echo "Elapsed time: $runtime seconds"