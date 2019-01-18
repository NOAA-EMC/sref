#!/bin/bash
set -eux

BUILD_NMMB=false
BUILD_NPS=true
BUILD_WRF_ARW=false
BUILD_WPS=false

BUILD_PERTURB=false

base=$(dirname $PWD)

set +x
module purge
module use ${base}/modulefiles
module load SREF/v7.0.0_dell
module list
set -x

date

SORCsref=${base}/sorc
EXECsref=${base}/exec
mkdir -m 775 -p $EXECsref


if ${BUILD_NMMB}; then
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

fi

if ${BUILD_NPS}; then
###
### NPS
###

export JASPERINC=${JASPER_INC}
export JASPERLIB=${JASPER_LIBDIR}

cd ${SORCsref}/sref_nps.fd
./conf dell
./clean -a
./configure <<EOF
5
EOF
./compile nps > compile_nps.log 2>&1

cp ungrib.exe     ${EXECsref}/sref_ungrib_nems
cp ungribp.exe    ${EXECsref}/sref_ungribp_nems
cp geogrid.exe    ${EXECsref}/sref_geogrid_nems
cp metgrid.exe    ${EXECsref}/sref_metgrid_nems
cp nemsinterp.exe ${EXECsref}/sref_real_nmb

./compile mod_levs
cp util/src/mod_levs.exe  ${EXECsref}/sref_NMB_levs

fi

if ${BUILD_WRF_ARW}; then
###
### WRF_ARW
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

fi

if ${BUILD_WPS}; then
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

fi

if ${BUILD_PERTURB}; then
###
### PERTURB
###

cd ${SORCsref}/sref_perturb.fd
cd dio
sh compile_dell.sh
cd ..
make -f Makefile clean
make -f Makefile_dell
cp breeding_arw.exe    ${EXECsref}/sref_breeding_arw
#cp breeding_nmm.exe    ${EXECsref}/sref_breeding_nmm
cp breeding_nmb.exe    ${EXECsref}/sref_breeding_nmb
cp lbc_perturb_wrf.exe ${EXECsref}/sref_lbc_perturb_wrf
cp lbc_perturb_nmb.exe ${EXECsref}/sref_lbc_perturb_nmb

fi


date
