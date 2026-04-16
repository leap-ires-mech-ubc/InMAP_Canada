# Define the function to process a single date

#read -r year month day hour

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
ingem="/home/tfmrodge/scratch/GEMMACH_data/data/nc3/BASEGM_2015_017/${year}-${month}-${day}_${hour}${suf}-tmp2"
outgem="/home/tfmrodge/scratch/GEMMACH_data/data/nc3/BASEGM_2015_017/${year}-${month}-${day}_${hour}${suf}"

if [[ ! -f "$file1" ]]; then
    #Process the GEMMACH output data 
    ncks --no_tmp_fl -6 -O $ingem $outgem
    #Add compression, halves the size but not great
    #ncatted -O -a _FillValue,,o,f,-32767 "$ingem"
    #ncpdq -P pack --ppc default=6 -6 -O "$ingem" "$outgem"
else
    echo "$file1 exists"
fi

# Process the date
#echo "$hour"
#Process the GEMMACH output data
#May need to specify no temp files here - 

#rm data/nc3/BASEGM_2015_017/${year}-${month}-${day}_${hour}${suf}-tmp
#No RDPS - all done!