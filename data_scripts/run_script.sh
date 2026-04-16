#!/bin/bash
#SBATCH --time=2:59:59
#SBATCH --account=def-agiang01 #def-rscholes #def-agiang01 #def-rscholes #
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3900M #--memory=150G 
#SBATCH --job-name='process_outputs'
#SBATCH --array=2-9
module load python/3.12 scipy-stack proj/9.4.1
source ~/inmapenv/bin/activate
export PROJ_LIB="$EBROOTPROJ/share/proj"
export PROJ_DATA="$PROJ_LIB"
#python /home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/gemmach_on_inmap_base_zonalstats.py $SLURM_ARRAY_TASK_ID
python /home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/process_inmap_outputs.py $SLURM_ARRAY_TASK_ID
#salloc --time=3:0:0 --mem-per-cpu=10G --ntasks=10 --nodes=1 --account=def-rscholes
#python /home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/combine_preprocs.py
#module load gdal
# for f in /home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/20250930*.shp; do
#     ogr2ogr -t_srs EPSG:4326 "${f%_lcc.shp}EPSG4326.shp" "$f"
# done
#ogrinfo -so /home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BAU_2020_E108_majorpts.shp file
# for f in /home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/20250926*__*; do
#     mv "$f" "${f/__/}"
# done

# cdo setmisstoc,0 /home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Preproc/inmapData_GEMMACH_BASEGM_2015_017.nc /home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Preproc/inmapData_GEMMACH_BASEGM_2015_017_complete.nc
# cdo -O fillmiss /home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Preproc/inmapData_GEMMACH_BASEGM_2015_017.nc /home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Preproc/inmapData_GEMMACH_BASEGM_2015_017_fillnn.nc