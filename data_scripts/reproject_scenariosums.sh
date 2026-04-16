#!/bin/bash
input_dir="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/scenario_sums"
output_dir="${input_dir}/georeferenced"
mkdir -p "$output_dir"

# Projection metadata
proj4="+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +a=6378137 +rf=298.257222101 +units=m +no_defs"
ulx=-4184312.05377675
uly=4270177.170000
lrx=3165687.946000
lry=-2029822.82977676

for file in "$input_dir"/*.nc; do
    base_filename=$(basename "$file" .nc)
    output_file="$output_dir/${base_filename}_georeferenced.nc"
    # Add projection metadata as global attributes
    ncatted -O \
        -a proj4,global,o,c,"$proj4" \
        -a geospatial_lon_min,global,o,d,$ulx \
        -a geospatial_lon_max,global,o,d,$lrx \
        -a geospatial_lat_min,global,o,d,$lry \
        -a geospatial_lat_max,global,o,d,$uly \
        "$temp_file"

    #ncks 
    # Move to final output
    mv "$temp_file" "$output_file"
done