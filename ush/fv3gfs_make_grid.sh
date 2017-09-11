#!/bin/bash
set -eux

nargv=$#
if [ $nargv -ne 3 -a $nargv -ne 6 -a $nargv -ne 12 ]; then
   echo "number of arguments must be 3 (regular cubic sphere grid), 6 (stretched grid) or 12 (nested grid)"
   echo "Usage for regular cubic sphere grid: "
   echo "  $0 resolution out_dir script_dir"
   echo "Usage for Stretched grid: "
   echo "  $0 resolution out_dir stetch_fac target_lon target_lat script_dir"
   echo "Usage for Nested grid: "
   echo "  $0 resolution out_dir stetch_fac target_lon target_lat refine_ratio istart_nest jstart_nest iend_nest jend_nest halo script_dir"
   exit 1
fi

res=$1
outdir=$2
if [ ! -s $outdir ]; then
  mkdir -p $outdir
fi

nx=`expr $res \* 2 `
cd $outdir

executable=$exec_dir/make_hgrid
if [ ! -s $executable ]; then
  echo "executable does not exist"
  exit 1
fi

if [ $nargv -eq 3 ]; then
  ntiles=6
  script_dir=$3
  $APRUN $executable --grid_type gnomonic_ed --nlon $nx --grid_name C${res}_grid

elif  [ $nargv -eq 6 ]; then
  stretch_fac=$3
  target_lon=$4
  target_lat=$5
  script_dir=$6
  ntiles=6
  $APRUN $executable --grid_type gnomonic_ed --nlon $nx --grid_name C${res}_grid --do_schmidt --stretch_factor ${stretch_fac} --target_lon ${target_lon} --target_lat ${target_lat}

elif  [ $nargv -eq 12 ]; then
  stretch_fac=$3
  target_lon=$4
  target_lat=$5
  refine_ratio=$6
  istart_nest=$7
  jstart_nest=$8
  iend_nest=$9
  jend_nest=${10}
  halo=${11}
  script_dir=${12}
  ntiles=7
  $APRUN $executable --grid_type gnomonic_ed --nlon $nx --grid_name C${res}_grid --do_schmidt --stretch_factor ${stretch_fac} --target_lon ${target_lon} --target_lat ${target_lat} \
     --nest_grid --parent_tile 6 --refine_ratio $refine_ratio --istart_nest $istart_nest --jstart_nest $jstart_nest --iend_nest $iend_nest --jend_nest $jend_nest --halo $halo --great_circle_algorithm

fi

if [ $? -ne 0 ]; then
  echo "ERROR in running create C$res grid without halo "
  exit 1
fi

#---------------------------------------------------------------------------------------
executable=$exec_dir/make_solo_mosaic
if [ ! -s $executable ]; then
  echo "executable does not exist"
  exit 1
fi

if [ $ntiles -eq 6 ]; then
  $APRUN $executable --num_tiles $ntiles --dir $outdir --mosaic C${res}_mosaic --tile_file C${res}_grid.tile1.nc,C${res}_grid.tile2.nc,C${res}_grid.tile3.nc,C${res}_grid.tile4.nc,C${res}_grid.tile5.nc,C${res}_grid.tile6.nc

elif [ $ntiles -eq 7 ]; then
  $APRUN $executable --num_tiles $ntiles --dir $outdir --mosaic C${res}_mosaic --tile_file C${res}_grid.tile1.nc,C${res}_grid.tile2.nc,C${res}_grid.tile3.nc,C${res}_grid.tile4.nc,C${res}_grid.tile5.nc,C${res}_grid.tile6.nc,C${res}_grid.tile7.nc

  $APRUN $executable --num_tiles 6 --dir $outdir --mosaic C${res}_coarse_mosaic --tile_file C${res}_grid.tile1.nc,C${res}_grid.tile2.nc,C${res}_grid.tile3.nc,C${res}_grid.tile4.nc,C${res}_grid.tile5.nc,C${res}_grid.tile6.nc  $executable --num_tiles 1 --dir $outdir --mosaic C${res}_nested_mosaic --tile_file C${res}_grid.tile7.nc
fi
