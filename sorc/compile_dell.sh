#!/bin/sh
set -eux

base=$(dirname $PWD)

set +x
module purge
module use ${base}/modulefiles
module load SREF/v7.0.0_dell
module list
set -x

date

#######################################

SORCsref=${base}/sorc
EXECsref=${base}/exec
mkdir -m 775 -p $EXECsref


if [ 1 == 1 ]; then
###
### NMMB
###
cd ${SORCsref}/sref_nmmb.fd/src
cp ${SORCsref}/sref_nmmb.fd/src/atmos/phys/module_LS_NOAHLSM.F90_driersoil ${SORCsref}/sref_nmmb.fd/src/atmos/phys/module_LS_NOAHLSM.F90
./configure dell
set +x
source conf/modules.nems
set -x
gmake clean
gmake nmm
cd ../exe
cp NEMS.x        ${EXECsref}/sref_wrf_nmb_driersoil

cd ${SORCsref}/sref_nmmb.fd/src
cp ${SORCsref}/sref_nmmb.fd/src/atmos/phys/module_LS_NOAHLSM.F90_normalsoil ${SORCsref}/sref_nmmb.fd/src/atmos/phys/module_LS_NOAHLSM.F90
./configure dell
set +x
source conf/modules.nems
set -x
gmake clean
gmake nmm
cd ../exe
cp NEMS.x        ${EXECsref}/sref_wrf_nmb_normalsoil

cd ../src
gmake clean
fi

if [ 1 == 1 ]; then
###
### ARW WRF
###

cd ${SORCsref}/sref_wrf_v3.5.1.fd
cp ${SORCsref}/sref_wrf_v3.5.1.fd/phys/module_sf_noahlsm.F_driersoil ${SORCsref}/sref_wrf_v3.5.1.fd/phys/module_sf_noahlsm.F
./compile_arw_dell.sh
ls -l ${EXECsref}/sref_wrf_arw
mv ${EXECsref}/sref_wrf_arw ${EXECsref}/sref_wrf_arw_driersoil

cp ${SORCsref}/sref_wrf_v3.5.1.fd/phys/module_sf_noahlsm.F_normalsoil ${SORCsref}/sref_wrf_v3.5.1.fd/phys/module_sf_noahlsm.F
./compile_arw_dell.sh
ls -l ${EXECsref}/sref_wrf_arw
mv ${EXECsref}/sref_wrf_arw ${EXECsref}/sref_wrf_arw_normalsoil

cd ..

fi

if [ 1 == 1 ]; then
###
### WPS
###

export JASPERINC=${JASPER_INC}
export JASPERLIB=${JASPER_LIBDIR}

cd ${SORCsref}/sref_wps_v3.5.1.fd
./clean -a
./configure <<EOF
23
EOF
./compile wps > compile_wps.log 2>&1

cp ungrib.exe  ${EXECsref}/sref_ungrib
cp geogrid.exe ${EXECsref}/sref_geogrid
cp metgrid.exe ${EXECsref}/sref_metgrid

./compile mod_levs
cp util/src/mod_levs.exe  ${EXECsref}/sref_ARW_levs

./clean -a

cd ..

fi

date
