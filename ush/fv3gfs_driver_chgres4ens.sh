#!/bin/bash
set -eux
#-------------------------------------------------------------------------------------------------
# Makes ICs on fv3 globally uniform cubed-sphere grid using operational GFS initial conditions.
# Fanglin Yang, 09/30/2016
#  This script is created based on the C-shell scripts fv3_gfs_preproc/IC_scripts/DRIVER_CHGRES.csh
#  and submit_chgres.csh provided by GFDL.  APRUN and environment variables are added to run on
#  WCOSS CRAY.  Directory and file names are standaridized to follow NCEP global model convention.
#  This script calls fv3gfs_chgres.sh.
# Fanglin Yang and George Gayno, 02/08/2017
#  Modified to use the new CHGRES George Gayno developed.
# Fanglin Yang 03/08/2017
#  Generalized and streamlined the script and enabled to run on multiple platforms.
#-------------------------------------------------------------------------------------------------

if [ ${member} =  mem1 ]; then m_enkf=001; fi
if [ ${member} =  mem2 ]; then m_enkf=002; fi
if [ ${member} =  mem3 ]; then m_enkf=003; fi
if [ ${member} =  mem4 ]; then m_enkf=004; fi
if [ ${member} =  mem5 ]; then m_enkf=005; fi
if [ ${member} =  mem6 ]; then m_enkf=006; fi
if [ ${member} =  mem7 ]; then m_enkf=007; fi
if [ ${member} =  mem8 ]; then m_enkf=008; fi
if [ ${member} =  mem9 ]; then m_enkf=009; fi
if [ ${member} = mem10 ]; then m_enkf=010; fi

ulimit -a
ulimit -s unlimited

export CASE=${CASE:-C768}                    # resolution of tile: 48, 96, 192, 384, 768, 1152, 3072
export CDATE=$PDY$CYC

export machine=${machine:-WCOSS_C}
export NODES=1
export OMP_NUM_THREADS_CH=${OMP_NUM_THREADS_CH:-24}
export OMP_NUM_THREADS=6
export OMP_STACKSIZE=2048m

export BASE_GSM=$HOMEfv3href
export script_dir=$HOMEfv3href/ush
export DATA=$WORK_DIR

export gtype=nest	          # grid type = uniform, stretch, or nest

#---------------------------------------------------------
#---------------------------------------------------------

if [ $gtype = uniform ];  then
  echo "creating uniform ICs"
  name=${CASE}
  export ntiles=6
elif [ $gtype = stretch ]; then
  stetch_fac=1.5       	                 # Stretching factor for the grid
  rn=$( echo "$stetch_fac * 10" | bc | cut -c1-2 )
  name=${CASE}r${rn}       		 # identifier based on refined location (same as grid)
  export ntiles=6
  echo "creating stretched ICs"
elif [ $gtype = nest ]; then
  stetch_fac=1.5	                     # Stretching factor for the grid
  rn=$( echo "$stetch_fac * 10" | bc | cut -c1-2 )
  refine_ratio=3   	                     # Specify the refinement ratio for nest grid
  name=${CASE}r${rn}n${refine_ratio}     # identifier based on nest location (same as grid)
  export ntiles=7
  echo "creating nested ICs"
else
  echo "Error: please specify grid type with 'gtype' as uniform, stretch, or nest"
fi

#---------------------------------------------------------------

echo "processing " $CDATE
ymd=`echo $CDATE |cut -c 1-8`
cyc=`echo $CDATE |cut -c 9-10`

# export gfs_dir=${inidir:-$COMROOTp2/gfs/prod/gfs.$ymd}  #directory that contains GFS input file(s)
if [ $member = ctl ];  then
  gfs_dir=/gpfs/hps/nco/ops/com/gfs/para/gfs.$ymd  #directory that contains GFS input file(s)
  gdas_dir=/gpfs/hps/nco/ops/com/gfs/para/gdas.$ymd  #directory that contains GFS input file(s)
else
  enkf_dir=/gpfs/hps/nco/ops/com/gfs/para/enkf.$ymd/$cyc  #directory that contains GFS input file(s)
fi

#if [ -s ${gfs_dir}/siganl.gfs.$CDATE ]; then
#  export ATMANL=$gfs_dir/siganl.gfs.$CDATE
#  export SFCANL=$gfs_dir/sfcanl.gfs.$CDATE
#else
#  export ATMANL=$gfs_dir/gfs.t${cyc}z.sanl
#  export SFCANL=$gfs_dir/gfs.t${cyc}z.sfcanl
#fi

if [ $member = ctl ];  then

   if [ -s ${gdas_dir}/gdas.t${cyc}z.atmanl.nemsio ]; then
     export ATMANL=$gdas_dir/gdas.t${cyc}z.atmanl.nemsio
     export SFCANL=$gdas_dir/gdas.t${cyc}z.sfcanl.nemsio
   else
     export ATMANL=$gfs_dir/gfs.t${cyc}z.atmanl.nemsio
     export SFCANL=$gfs_dir/gfs.t${cyc}z.sfcanl.nemsio
   fi

else

   if [ -s ${enkf_dir}/gdas.t${cyc}z.atmanl.mem${m_enkf}.nemsio ]; then
     export ATMANL=$enkf_dir/gdas.t${cyc}z.atmanl.mem${m_enkf}.nemsio
     export SFCANL=$enkf_dir/gdas.t${cyc}z.sfcanl.mem${m_enkf}.nemsio
   else
     export ATMANL=$enkf_dir/gdas.t${cyc}z.atmanl.mem${m_enkf}.nemsio
     export SFCANL=$enkf_dir/gdas.t${cyc}z.sfcanl.mem${m_enkf}.nemsio
   fi

fi

# make atmospheric and surface data
$script_dir/fv3gfs_chgres.sh

