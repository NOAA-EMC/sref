#!/bin/bash

set -eux
set -o pipefail

FIXDIR=${FIXfv3href}/fix_fv3
FIX_AM=${FIXfv3href}/fix_am
CO2DIR=$FIX_AM/fix_co2_proj

mkdir INPUT
mkdir RESTART

cp -p $FIX_AM/global_solarconstant_noaa_an.txt  solarconstant_noaa_an.txt
cp -p $FIX_AM/global_o3prdlos.f77               INPUT/global_o3prdlos.f77
cp -p $FIX_AM/global_sfc_emissivity_idx.txt     sfc_emissivity_idx.txt
cp -p $FIX_AM/global_co2historicaldata_glob.txt co2historicaldata_glob.txt
cp -p $FIX_AM/co2monthlycyc.txt                 co2monthlycyc.txt
cp -p $FIX_AM/global_climaeropac_global.txt     aerosol.dat

for file in `ls $CO2DIR/global_co2historicaldata* ` ; do
 cp $file $(echo $(basename $file) |sed -e "s/global_//g")
done
#
#copy tile data and orography
#
res=768
ntiles=7
tile=1
while [ $tile -le $ntiles ]; do
 cp -p $FIXDIR/C${res}/C${res}_oro_data.tile${tile}.nc INPUT/oro_data.tile${tile}.nc
 cp -p $FIXDIR/C${res}/C${res}_grid.tile${tile}.nc     INPUT/C${res}_grid.tile${tile}.nc
 cp -p $COMOUT/gfs_data.tile${tile}.nc                 INPUT/gfs_data.tile${tile}.nc
 cp -p $COMOUT/sfc_data.tile${tile}.nc                 INPUT/sfc_data.tile${tile}.nc
 tile=`expr $tile + 1 `
done
cp -p $FIXDIR/C${res}/C${res}_mosaic.nc INPUT/grid_spec.nc
cp -p $COMOUT/gfs_ctrl.nc               INPUT/gfs_ctrl.nc
#
# the next 4 links are a hack required for running a nest
#
cd INPUT
ln -sf C${res}_grid.tile7.nc C${res}_grid.nest02.tile7.nc
ln -sf oro_data.tile7.nc oro_data.nest02.tile7.nc
ln -sf gfs_data.tile7.nc gfs_data.nest02.tile7.nc
ln -sf sfc_data.tile7.nc sfc_data.nest02.tile7.nc
cd ..
#
# copy configuration files
#
cp -p ${PARMfv3href}/data_table .
cp -p ${PARMfv3href}/field_table .
cp -p ${PARMfv3href}/nems.configure .
cp -p ${PARMfv3href}/input.nml .
cp -p ${PARMfv3href}/input_nest02.nml .

START_YYYY=$( echo $PDY | cut -c 1-4 )
START_MM=$( echo $PDY | cut -c 5-6 )
START_DD=$( echo $PDY | cut -c 7-8 )
START_HH=$cyc

cat ${PARMfv3href}/model_configure.IN | sed s/_START_YYYY_/$START_YYYY/ \
                                      | sed s/_START_MM_/$START_MM/ \
                                      | sed s/_START_DD_/$START_DD/ \
                                      | sed s/_START_HH_/$START_HH/ \
                                      | sed s/_FLENGTH_/$FLENGTH/ > model_configure

cat ${PARMfv3href}/diag_table.IN | sed s/_START_YYYY_/$START_YYYY/ \
                                 | sed s/_START_MM_/$START_MM/ \
                                 | sed s/_START_DD_/$START_DD/ \
                                 | sed s/_START_HH_/$START_HH/ > diag_table

#
# copy and run the executable
#
cp ${EXECfv3href}/fv3_32bit.exe fv3_gfs.x

${APRUN} fv3_gfs.x 1>out 2>err

err=$?
ERRSCRIPT=${ERRSCRIPT:-'eval [[ $err = 0 ]]'}
$ERRSCRIPT || exit $err

file=RESTART/sfc_data.nest02.tile7.nc
if [ -f "$file" ]
then
  mv RESTART/sfc_data.nest02.tile7.nc RESTART/sfc_data.tile7.nc
else
  echo "$file not found."
  exit 1
fi

mv atmos*.nc $COMOUT
mv nggps*.nc $COMOUT
