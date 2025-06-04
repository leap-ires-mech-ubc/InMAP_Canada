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

file1="data/nc3/BASEGM_2015_017/${year}-${month}-${day}_${hour}${suf}"
file2="data/nc3/RDPS_QC/${year}-${month}-${day}_${hour}${suf2}"

if [[ ! -f "$file1" ]]; then
    #Process the GEMMACH output data
    cdo remapbil,gridfile_lcc.txt data/nc3/BASEGM_2015_017/${year}-${month}-${day}_${hour}${suf}-tmp1 data/nc3/BASEGM_2015_017/${year}-${month}-${day}_${hour}${suf}-tmp2
elif [[ ! -f "$file2" ]]; then
    #Process RDPS Datafile
    ncks -6 -O data/nc3/RDPS_QC/${year}-${month}-${day}_${hour}${suf2} data/nc3/RDPS_QC/${year}-${month}-${day}_${hour}${suf2}
else
    echo "Both $file1 and $file2 exist."
fi
# Process the date
#echo "$day"

#RDPS

