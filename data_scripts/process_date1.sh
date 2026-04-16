# Define the function to process a single date

#read -r year month day hour

# Function to process a single date
year=$1
month=$2
day=$3
hour=$4
suf=_00_00_gemmach.nc
suf2=_00_00_rdpsqc.nc
mid=00_0

# Check if the date is valid
if ! date -d "${year}-${month}-${day}" >/dev/null 2>&1; then
    #echo "Invalid date: ${year}-${month}-${day}"
    return
fi

#Check if output files exist, run where they don't
file1="/home/tfmrodge/scratch/GEMMACH_data/data/nc3/BASEGM_2015_017/${year}-${month}-${day}_${hour}${suf}"
file2="/home/tfmrodge/scratch/GEMMACH_data/data/nc3/RDPS_QC/${year}-${month}-${day}_${hour}${suf2}"
ingem="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/BASEGM_2015_017/$year$month$day$mid$hour.nc"
outgem="/home/tfmrodge/scratch/GEMMACH_data/data/nc3/BASEGM_2015_017/${year}-${month}-${day}_${hour}${suf}-tmp1"
inrdps="/home/tfmrodge/scratch/GEMMACH_data/data/RDPS_QC/$year$month$day$mid$hour.nc"
outrdps="/home/tfmrodge/scratch/GEMMACH_data/data/nc3/RDPS_QC/${year}-${month}-${day}_${hour}${suf2}"
if [[ ! -f "$file1" ]]; then
    #Process GEMMACH datafile
    ncap2 --no_tmp_fl -O -s  'pressure=float(exp(a_1+b_1*ln(level1))/100)' $ingem $outgem
else
    echo "$file1 already exists"
fi
if [[ ! -f "$file2" && -f "$inrdps" ]]; then
    #Process RDPS Datafile
    cdo -b F32 remapbil,gridfile_lcc.txt $inrdps $outrdps
else
    echo "$file2 exists or no input RDPS file."
fi