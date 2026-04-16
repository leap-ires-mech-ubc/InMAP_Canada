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
    #Define filepaths
    file1="data/nc3/BASEGM_2015_017/${year}-${month}-${day}_${hour}${suf}"
    file2="data/nc3/RDPS_QC/${year}-${month}-${day}_${hour}${suf2}"
    ingem="/home/tfmrodge/scratch/GEMMACH_data/data/nc3/BASEGM_2015_017/${year}-${month}-${day}_${hour}${suf}-tmp1"
    outgem="/home/tfmrodge/scratch/GEMMACH_data/data/nc3/BASEGM_2015_017/${year}-${month}-${day}_${hour}${suf}-tmp2"
    inrdps="/home/tfmrodge/scratch/GEMMACH_data/data/nc3/RDPS_QC/${year}-${month}-${day}_${hour}${suf2}"
    outrdps="/home/tfmrodge/scratch/GEMMACH_data/data/nc3/RDPS_QC/${year}-${month}-${day}_${hour}${suf2}"

    if [[ ! -f "$file1" ]]; then
        #Process the GEMMACH output data
        cdo -b F32 remapbil,gridfile_lcc.txt $ingem $outgem
    else
        echo "$file1 already exists"
    fi
    if [[ ! -f "$file2" && -f "$inrdps" ]]; then
        #Process RDPS Datafile
        ncks -6 -O $inrdps $outrdps
    else
        echo "$file2 exists or no RDPS input file"
    fi
    # Process the date
    #echo "$day"

    #RDPS

