#!/bin/bash

#----WCOSS_CRAY JOBCARD
#BSUB -L /bin/sh
#BSUB -P GFS-T2O
#BSUB -oo log.grid.%J
#BSUB -eo log.grid.%J
#BSUB -J grid_fv3
#BSUB -q debug
#BSUB -M 2400
#BSUB -W 00:30
#BSUB -extsched 'CRAYLINUX[]'
#----THEIA JOBCARD
#PBS -N fv3_grid_driver
#PBS -A glbss
#PBS -o log.grid
#PBS -e log.grid
#PBS -l nodes=1:ppn=24
#PBS -q debug
#PBS -l walltime=00:30:00


set -eux

export envir=para
source ../parm/fv3href_para_config_cray

machine=${machine:-WCOSS_C}

ulimit -a
ulimit -s unlimited

#-----------------------------
# Makes FV3 cubed-sphere grid
#-----------------------------

#----------------------------------------------------------------

res=768         # resolution of tile: 48, 96, 192, 384, 768, 1152, 3072
gtype=nest      # grid type: uniform, stretch, or nest

#----------------------------------------------------------------
# The orography code runs with threads.  On Cray, the code is
# optimized for six threads.  Do not change.
#----------------------------------------------------------------

export OMP_NUM_THREADS=6
export OMP_STACKSIZE=2048m

if [ $machine = WCOSS_C ]; then
 set +x
 . $MODULESHOME/init/sh
 module load PrgEnv-intel cfp-intel-sandybridge/1.1.0
 module list
 set -x
 export NODES=1
 export APRUN="aprun -n 1 -N 1 -j 1 -d 1 -cc depth"
 export KMP_AFFINITY=disabled
elif [ $machine = THEIA ]; then
 set +x
 . /apps/lmod/lmod/init/sh
 module purge
 module load intel/16.1.150
 module load impi
 module load hdf5/1.8.14
 module load netcdf/4.3.0
 module list
 set -x
 export APRUN=time
fi
#----------------------------------------------------------------
home_dir=$HOMEfv3href
script_dir=$home_dir/ush
topo=$FIXfv3href/topo
export exec_dir=$home_dir/exec

out_dir=$home_dir/fix/fix_fv3/C${res}

mkdir -p $out_dir $TMPDIR
cd $TMPDIR ||exit 8

#----------------------------------------------------------------
if [ $gtype = uniform ];  then
  echo "creating uniform ICs"
elif [ $gtype = stretch ]; then
  stetch_fac=1.5     # Stretching factor for the grid
  target_lon=-97.5   # center longitude of the highest resolution tile
  target_lat=35.5    # center latitude of the highest resolution tile
  title=c96s         # identifier based on refined location
  echo "creating stretched grid"
elif [ $gtype = nest ]; then
  stetch_fac=1.5     # Stretching factor for the grid
  target_lon=-97.5   # center longitude of the highest resolution tile
  target_lat=35.5    # center latitude of the highest resolution tile
  refine_ratio=3     # Specify the refinement ratio for nest grid
  istart_nest=193    # Specify the starting i-direction index of nest grid in parent tile supergrid(Fortran index)
  jstart_nest=329    # Specify the starting j-direction index of nest grid in parent tile supergrid(Fortran index)
  iend_nest=1344     # Specify the ending i-direction index of nest grid in parent tile supergrid(Fortran index)
  jend_nest=1288     # Specify the ending j-direction index of nest grid in parent tile supergrid(Fortran index)
  halo=3             # halo size to be used in the atmosphere cubic sphere model. It only needs to be specified when --nest_grid is set
  title=c768s        # identifier based on nest location
  echo "creating nested grid"
else
  echo "Error: please specify grid type with 'gtype' as uniform, stretch, or nest"
fi

#----------------------------------------------------------------
#filter_topo parameters. C192->50km, C384->25km, C768->13km, C1152->8.5km, C3072->3.2km
if [ $res -eq 48 ]; then
 cd4=0.12;  max_slope=0.12; n_del2_weak=4;   peak_fac=1.1
elif [ $res -eq 96 ]; then
 cd4=0.12;  max_slope=0.12; n_del2_weak=8;   peak_fac=1.1
elif [ $res -eq 192 ]; then
 cd4=0.15;  max_slope=0.12; n_del2_weak=12;  peak_fac=1.05
elif [ $res -eq 384 ]; then
 cd4=0.15;  max_slope=0.12; n_del2_weak=12;  peak_fac=1.0
elif [ $res -eq 768 ]; then
 cd4=0.15;  max_slope=0.12; n_del2_weak=16;  peak_fac=1.0
elif [ $res -eq 1152 ]; then
 cd4=0.15;  max_slope=0.16; n_del2_weak=20;  peak_fac=1.0
elif [ $res -eq 3072 ]; then
 cd4=0.15;  max_slope=0.30; n_del2_weak=24;  peak_fac=1.0
else
 echo "grid C$res not supported, exit"
 exit
fi

#----------------------------------------------------------------
# Make grid and orography

if [ $gtype = uniform ];  then
  ntiles=6
  name=C${res}
  grid_dir=$TMPDIR/$name/grid
  orog_dir=$TMPDIR/$name/orog
  filter_dir=$TMPDIR/$name/filter_topo
  rm -rf $TMPDIR/$name
  mkdir -p $grid_dir $orog_dir $filter_dir

  echo
  echo "............ execute fv3gfs_make_grid.sh ................."
  $script_dir/fv3gfs_make_grid.sh $res $grid_dir $script_dir
#
  echo "Begin uniform orography generation at `date`"
#
#on WCOSS and WCOSS_C use cfp to run multiple tiles simulatneously for the orography
#
if [ $machine = WCOSS_C ]; then
#
#create input files for cfp in order to run multiple tiles of the orography simulataneously
#
  echo "$script_dir/fv3gfs_make_orog.sh $res 3 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >$TMPDIR/orog.file1
  echo "$script_dir/fv3gfs_make_orog.sh $res 6 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
  echo "$script_dir/fv3gfs_make_orog.sh $res 1 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
  echo "$script_dir/fv3gfs_make_orog.sh $res 2 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
  echo "$script_dir/fv3gfs_make_orog.sh $res 4 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
  echo "$script_dir/fv3gfs_make_orog.sh $res 5 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1

  aprun -j 1 -n 4 -N 4 -d 6 -cc depth cfp $TMPDIR/orog.file1
  rm $TMPDIR/orog.file1
elif [ $machine = THEIA ]; then
  for tile in 1 2 3 4 5 6 ; do
    echo
    echo "............ execute fv3gfs_make_orog.sh for tile $tile .................."
    $script_dir/fv3gfs_make_orog.sh $res $tile $grid_dir $orog_dir $script_dir $topo $TMPDIR
  done
fi
  echo "End uniform orography generation at `date`"
#
  echo
  echo "............ execute fv3gfs_filter_topo.sh .............."
  $script_dir/fv3gfs_filter_topo.sh $res $grid_dir $orog_dir $filter_dir $cd4 $peak_fac $max_slope $n_del2_weak $script_dir
  echo "Grid and orography files are now prepared for uniform grid"

elif [ $gtype = stretch ]; then
  ntiles=6
  rn=$( echo "$stetch_fac * 10" | bc | cut -c1-2 )
  name=C${res}r${rn}_${title}
  grid_dir=$TMPDIR/${name}/grid
  orog_dir=$TMPDIR/$name/orog
  filter_dir=$TMPDIR/${name}/filter_topo
  rm -rf $TMPDIR/$name
  mkdir -p $grid_dir $orog_dir $filter_dir

  echo
  echo "............ execute fv3gfs_make_grid.sh ................."
  $script_dir/fv3gfs_make_grid.sh $res $grid_dir $stetch_fac $target_lon $target_lat $script_dir
#
  echo "Begin stretch orography generation at `date`"
#
#on WCOSS and WCOSS_C use cfp to run multiple tiles simulatneously for the orography
#
if [ $machine = WCOSS_C ]; then
#
#create input files for cfp in order to run multiple tiles of the orography simulataneously
#
  echo "$script_dir/fv3gfs_make_orog.sh $res 1 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >$TMPDIR/orog.file1
  echo "$script_dir/fv3gfs_make_orog.sh $res 3 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
  echo "$script_dir/fv3gfs_make_orog.sh $res 4 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
  echo "$script_dir/fv3gfs_make_orog.sh $res 5 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
  echo "$script_dir/fv3gfs_make_orog.sh $res 2 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
  echo "$script_dir/fv3gfs_make_orog.sh $res 6 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1

  aprun -j 1 -n 4 -N 4 -d 6 -cc depth cfp $TMPDIR/orog.file1
  rm $TMPDIR/orog.file1
elif [ $machine = THEIA ]; then
  for tile in 1 2 3 4 5 6 ; do
    echo
    echo "............ execute fv3gfs_make_orog.sh for tile $tile .................."
    $script_dir/fv3gfs_make_orog.sh $res $tile $grid_dir $orog_dir $script_dir $topo $TMPDIR
  done
fi
  echo "End stretch orography generation at `date`"
#
  echo
  echo "............ execute fv3gfs_filter_topo.sh .............."
  $script_dir/fv3gfs_filter_topo.sh $res $grid_dir $orog_dir $filter_dir $cd4 $peak_fac $max_slope $n_del2_weak $script_dir
  echo "Grid and orography files are now prepared for stretched grid"

elif [ $gtype = nest ]; then
  ntiles=7
  rn=$( echo "$stetch_fac * 10" | bc | cut -c1-2 )
  name=C${res}r${rn}n${refine_ratio}_${title}
  grid_dir=$TMPDIR/${name}/grid
  orog_dir=$TMPDIR/$name/orog
  filter_dir=$orog_dir   # nested grid topography will be filtered online
  rm -rf $TMPDIR/$name
  mkdir -p $grid_dir $orog_dir $filter_dir

  echo
  echo "............ execute fv3gfs_make_grid.sh ................."
  $script_dir/fv3gfs_make_grid.sh $res $grid_dir $stetch_fac $target_lon $target_lat $refine_ratio $istart_nest $jstart_nest $iend_nest $jend_nest $halo $script_dir

  echo "Begin stretch nest orography generation at `date`"
#
#on WCOSS and WCOSS_C use cfp to run multiple tiles simulatneously for the orography
#
if [ $machine = WCOSS_C ]; then
#
#create input files for cfp in order to run multiple tiles of the orography simulataneously
#
  export APRUN=time
  echo "$script_dir/fv3gfs_make_orog.sh $res 1 $grid_dir $orog_dir $script_dir $topo $TMPDIR " > $TMPDIR/orog.file1
  echo "$script_dir/fv3gfs_make_orog.sh $res 3 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
  echo "$script_dir/fv3gfs_make_orog.sh $res 4 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
  echo "$script_dir/fv3gfs_make_orog.sh $res 2 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
  echo "$script_dir/fv3gfs_make_orog.sh $res 5 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
  echo "$script_dir/fv3gfs_make_orog.sh $res 6 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1
  echo "$script_dir/fv3gfs_make_orog.sh $res 7 $grid_dir $orog_dir $script_dir $topo $TMPDIR " >>$TMPDIR/orog.file1

  aprun -j 1 -n 4 -N 4 -d 6 -cc depth cfp $TMPDIR/orog.file1
  rm $TMPDIR/orog.file1
elif [ $machine = THEIA ]; then
 for tile in 1 2 3 4 5 6 7; do
    echo
    echo "............ execute fv3gfs_make_orog.sh for tile $tile .................."
    $script_dir/fv3gfs_make_orog.sh $res $tile $grid_dir $orog_dir $script_dir $topo $TMPDIR
  done
fi
  echo "End stretch nest orography generation at `date`"
  echo "Grid and orography files are now prepared for nested grid"
fi
#----------------------------------------------------------------



echo
echo "--------------------------------------------------"

tile=1
while [ $tile -le $ntiles ]; do
 cp $filter_dir/oro.C${res}.tile${tile}.nc $out_dir/C${res}_oro_data.tile${tile}.nc
 cp $grid_dir/C${res}_grid.tile${tile}.nc  $out_dir/C${res}_grid.tile${tile}.nc
 tile=`expr $tile + 1 `
done
#cp $filter_dir/C${res}_mosaic.nc $out_dir/grid_spec.nc
#cp $filter_dir/C${res}_mosaic.nc $out_dir/C${res}_mosaic.nc
cp $grid_dir/C${res}_mosaic.nc $out_dir/C${res}_mosaic.nc

exit
