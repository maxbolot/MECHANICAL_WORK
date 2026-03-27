#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 12:00:00 # time required
#SBATCH --mail-user=mbolot@princeton.edu

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

source_dir=/scratch/gpfs/mbolot/data/20191020.00Z.C3072.L79x2_pire/$code
source_dir_area=/scratch/gpfs/mbolot/results/GLOBALFV3/gridarea

target_dir_a_m30p30=/scratch/gpfs/mbolot/results/GLOBALFV3/column_precip_diagnostics.m30p30.iwp.fix/area
target_dir_w_m30p30=/scratch/gpfs/mbolot/results/GLOBALFV3/column_precip_diagnostics.m30p30.iwp.fix/work
target_dir_l_m30p30=/scratch/gpfs/mbolot/results/GLOBALFV3/column_precip_diagnostics.m30p30.iwp.fix/water_lift


[[ -d $target_dir_a_m30p30 ]] || mkdir -p $target_dir_a_m30p30
[[ -d $target_dir_w_m30p30 ]] || mkdir -p $target_dir_w_m30p30
[[ -d $target_dir_l_m30p30 ]] || mkdir -p $target_dir_l_m30p30


# create temporary directory and cd into it

tmp_dir=`mktemp -d --tmpdir=/tmp`
cd $tmp_dir


##### COMPUTE CONVECTIVE WORK

# create temporary file to pass to the fortran program

module load anaconda3/2021.11
. "$CONDA_PREFIX/etc/profile.d/conda.sh"
conda activate nco_env

ncdump -h $source_dir/ps_coarse_C3072_1440x720.fre.nc | ncgen -k4 -o work.nc

ncks -A -v grid_xt_coarse,grid_xt_coarse_bnds,grid_yt_coarse,grid_yt_coarse_bnds,time \
    $source_dir/ps_coarse_C3072_1440x720.fre.nc work.nc

ncrename -v ps_coarse,work work.nc

conda deactivate 
module unload anaconda3

module load anaconda3/2021.11
. "$CONDA_PREFIX/etc/profile.d/conda.sh"
conda activate cdo_env

cdo splitgrid work.nc work   # work contained in grid 01 (necessary in average_DT, etc. are present in file structure)

conda deactivate 
module unload anaconda3

module load anaconda3/2021.11
. "$CONDA_PREFIX/etc/profile.d/conda.sh"
conda activate nco_env

ncatted -a units,grid_xt_coarse,m,c,"degrees_E" work01.nc
ncatted -a units,grid_yt_coarse,m,c,"degrees_N" work01.nc

function ncdmnsz { ncks --trd -m -M ${2} | grep -E -i ": ${1}, size =" | cut -f 7 -d ' ' | uniq ; }
tsize=$(ncdmnsz time work01.nc)

conda deactivate 
module unload anaconda3


# CREATE OUTPUTS FOR THERMAL AND VAPOR PART OF MECHANICAL WORK

module load anaconda3/2021.11
. "$CONDA_PREFIX/etc/profile.d/conda.sh"
conda activate nco_env

cp work01.nc work01_th.nc
ncrename -v work,work_th work01_th.nc
cp work01.nc work01_qv.nc
ncrename -v work,work_qv work01_qv.nc

conda deactivate 
module unload anaconda3


# interpolate ice flux at EDDY resolution and write result to iceflux_eddy_plev01
# create namelist with control parameters to be read by interpolation program

module load intel-oneapi/2024.2
module load netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2
module load hdf5/oneapi-2024.2/1.14.4

cat << EOF > config.nml
&config
    path_dz         = '$source_dir/DZ_C3072_1440x720.fre.nc',
    path_temp       = '$source_dir/temp_coarse_C3072_1440x720.fre.nc',
    path_omega      = '$source_dir/ptend_coarse_C3072_1440x720.fre.nc',
    path_qv         = '$source_dir/sphum_coarse_C3072_1440x720.fre.nc',
    path_omt        = '$source_dir/omT_coarse_C3072_1440x720.fre.nc',
    path_omqv       = '$source_dir/omqv_coarse_C3072_1440x720.fre.nc',
    path_work_out   = '$tmp_dir/work01.nc',
    path_work_th_out = '$tmp_dir/work01_th.nc',
    path_work_qv_out = '$tmp_dir/work01_qv.nc',
    nx              = 1440,
    ny              = 720,
    nplev           = 79,
    nt              = $tsize,
/
EOF

/home/mbolot/Projects/FV3_PIRE/programs/mechanical_work/work config.nml


module unload intel-oneapi/2024.2
module unload netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2
module unload hdf5/oneapi-2024.2/1.14.4


##### COMPUTE WATER LIFT WORK

# create temporary file to pass to the fortran program

module load anaconda3/2021.11
. "$CONDA_PREFIX/etc/profile.d/conda.sh"
conda activate nco_env

ncdump -h $source_dir/ps_coarse_C3072_1440x720.fre.nc | ncgen -k4 -o water_lift.nc

ncks -A -v grid_xt_coarse,grid_xt_coarse_bnds,grid_yt_coarse,grid_yt_coarse_bnds,time \
    $source_dir/ps_coarse_C3072_1440x720.fre.nc water_lift.nc

ncrename -v ps_coarse,water_lift water_lift.nc

conda deactivate 
module unload anaconda3

module load anaconda3/2021.11
. "$CONDA_PREFIX/etc/profile.d/conda.sh"
conda activate cdo_env

cdo splitgrid water_lift.nc water_lift   # work contained in grid 01 (necessary in average_DT, etc. are present in file structure)

conda deactivate 
module unload anaconda3

module load anaconda3/2021.11
. "$CONDA_PREFIX/etc/profile.d/conda.sh"
conda activate nco_env

ncatted -a units,grid_xt_coarse,m,c,"degrees_E" water_lift01.nc
ncatted -a units,grid_yt_coarse,m,c,"degrees_N" water_lift01.nc

function ncdmnsz { ncks --trd -m -M ${2} | grep -E -i ": ${1}, size =" | cut -f 7 -d ' ' | uniq ; }
tsize=$(ncdmnsz time water_lift01.nc)

conda deactivate 
module unload anaconda3

# CREATE OUTPUTS FOR VAPOR, LIQUID WATER, RAIN, ICE, SNOW AND GRAUPEL LIFT


module load anaconda3/2021.11
. "$CONDA_PREFIX/etc/profile.d/conda.sh"
conda activate nco_env

cp water_lift01.nc vapor_lift01.nc
ncrename -v water_lift,vapor_lift vapor_lift01.nc
cp water_lift01.nc liqwat_lift01.nc
ncrename -v water_lift,liqwat_lift liqwat_lift01.nc
cp water_lift01.nc rain_lift01.nc
ncrename -v water_lift,rain_lift rain_lift01.nc
cp water_lift01.nc icewat_lift01.nc
ncrename -v water_lift,icewat_lift icewat_lift01.nc
cp water_lift01.nc snow_lift01.nc
ncrename -v water_lift,snow_lift snow_lift01.nc
cp water_lift01.nc graupel_lift01.nc
ncrename -v water_lift,graupel_lift graupel_lift01.nc

conda deactivate 
module unload anaconda3


# interpolate ice flux at EDDY resolution and write result to iceflux_eddy_plev01
# create namelist with control parameters to be read by interpolation program

module load intel-oneapi/2024.2
module load netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2
module load hdf5/oneapi-2024.2/1.14.4

cat << EOF > config.nml
&config
    path_dz         = '$source_dir/DZ_C3072_1440x720.fre.nc',
    path_omega      = '$source_dir/ptend_coarse_C3072_1440x720.fre.nc',
    path_qv         = '$source_dir/sphum_coarse_C3072_1440x720.fre.nc',
    path_qw         = '$source_dir/liq_wat_coarse_C3072_1440x720.fre.nc',
    path_qr         = '$source_dir/rainwat_coarse_C3072_1440x720.fre.nc',
    path_qi         = '$source_dir/ice_wat_coarse_C3072_1440x720.fre.nc',
    path_qs         = '$source_dir/snowwat_coarse_C3072_1440x720.fre.nc',
    path_qg         = '$source_dir/graupel_coarse_C3072_1440x720.fre.nc',
    path_omqv       = '$source_dir/omqv_coarse_C3072_1440x720.fre.nc',
    path_omqw       = '$source_dir/omql_coarse_C3072_1440x720.fre.nc',
    path_omqr       = '$source_dir/omqr_coarse_C3072_1440x720.fre.nc',
    path_omqi       = '$source_dir/omqi_coarse_C3072_1440x720.fre.nc',
    path_omqs       = '$source_dir/omqs_coarse_C3072_1440x720.fre.nc',
    path_omqg       = '$source_dir/omqg_coarse_C3072_1440x720.fre.nc',
    path_lift_out   = '$tmp_dir/water_lift01.nc',
    path_lift_qv_out = '$tmp_dir/vapor_lift01.nc',
    path_lift_qw_out = '$tmp_dir/liqwat_lift01.nc',
    path_lift_qr_out = '$tmp_dir/rain_lift01.nc',
    path_lift_qi_out = '$tmp_dir/icewat_lift01.nc',
    path_lift_qs_out = '$tmp_dir/snow_lift01.nc',
    path_lift_qg_out = '$tmp_dir/graupel_lift01.nc',
    nx              = 1440,
    ny              = 720,
    nplev           = 79,
    nt              = $tsize,
/
EOF

/home/mbolot/Projects/FV3_PIRE/programs/water_lift/lift config.nml


module unload intel-oneapi/2024.2
module unload netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2
module unload hdf5/oneapi-2024.2/1.14.4


### BINNING BY PRECIP

lon1=0
lon2=360

module load matlab/R2021a
matlab -nodesktop -nosplash << EOF
addpath('/home/mbolot/Projects/FV3_PIRE/programs')
bin_area_by_precip('$source_dir/PRATEsfc_coarse_C3072_1440x720.fre.nc','$source_dir_area/gridarea.nc','$code','$target_dir_a_m30p30',$lon1,$lon2,-30,30)
joint_distribution_precip_work('$source_dir/PRATEsfc_coarse_C3072_1440x720.fre.nc','$tmp_dir/work01.nc','work','$source_dir_area/gridarea.nc','$code','$target_dir_w_m30p30',0,360,-30,30)
joint_distribution_precip_work('$source_dir/PRATEsfc_coarse_C3072_1440x720.fre.nc','$tmp_dir/water_lift01.nc','water_lift','$source_dir_area/gridarea.nc','$code','$target_dir_l_m30p30',0,360,-30,30)
exit
EOF
module unload matlab/R2021a

