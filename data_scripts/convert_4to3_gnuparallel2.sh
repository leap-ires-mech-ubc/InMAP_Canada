#!/bin/bash
#SBATCH --time=12:30:00
#SBATCH --job-name=convert_4to3
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
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

#9s with 4 cores, 28s with 1 core
#parallel --citation # #SLURM_NTASKS_PER_NODE=4
# Start the timer
start=$(date +%s)

# Set the range of months and days to process
start_month=05
end_month=05
year=2019

#source /home/tfmrodge/scratch/GEMMACH_data/process_date.sh
# Generate the commands using GNU Parallel #bash -c 'source /home/tfmrodge/scratch/GEMMACH_data/process_date.sh; process_date1 "$@"' _
#bash -c '. /home/tfmrodge/scratch/GEMMACH_data/process_date.sh; process_date1 "$@"' _
parallel -j $SLURM_NTASKS ./process_date1.sh ::: $year ::: $(seq -w $start_month $end_month) ::: {01..31} ::: {00..23}
printf "done job 1"
parallel -j $SLURM_NTASKS ./process_date2.sh ::: $year ::: $(seq -w $start_month $end_month) ::: {01..31} ::: {00..23}
#Remove tmp1 files
rm /home/tfmrodge/scratch/GEMMACH_data/data/nc3/BASEGM_2015_017/2019-*-tmp1
printf "done job 2"
parallel -j $SLURM_NTASKS ./process_date3.sh ::: $year ::: $(seq -w $start_month $end_month) ::: {01..31} ::: {00..23}
#Remove tmp2 files
rm /home/tfmrodge/scratch/GEMMACH_data/data/nc3/BASEGM_2015_017/2019-*-tmp2
printf "done job 3"
#Reproject gem_geophy file
cdo remapbil,gridfile_lcc.txt data/GEOPHY_VF/Gem_geophy.nc data/GEOPHY_VF/gem_geophy_lcc.nc
# Calculate the elapsed time
end=$(date +%s)
runtime=$((end - start))

echo "Elapsed time: $runtime seconds"