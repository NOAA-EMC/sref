#!/bin/bash
set -eux

BUILD_NMMB=true
BUILD_NPS=true
BUILD_WRF_ARW=true
BUILD_WPS=true
BUILD_PERTURB=true
BUILD_COLDSTART=true
BUILD_POST=true
BUILD_PRDGEN=true

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

# transutil.fd
cd ${SORCsref}/transutil.fd/v1.0.0
rm -f libtransutil_?.a
cd src
./makelibtransutil.sh

# sref_cluster_NCEP.fd
cd ${SORCsref}/sref_cluster_NCEP.fd
gmake -f makefile_cluster clean
gmake -f makefile_cluster
mv sref_cluster_NCEP ${EXECsref}/.
gmake -f makefile_cluster clean
gmake -f makefile_weight clean
gmake -f makefile_weight
mv sref_clusterweight ${EXECsref}/.
gmake -f makefile_weight clean

# sref_cluster_OU.fd
cd ${SORCsref}/sref_cluster_OU.fd
gmake clean
gmake
mv sref_cluster_OU ${EXECsref}/.
gmake clean

# sref_calfcsterr.fd
cd ${SORCsref}/sref_calfcsterr.fd
gmake clean
gmake
mv sref_calfcsterr ${EXECsref}/.
gmake clean

# sref_biasestimate.fd
cd ${SORCsref}/sref_biasestimate.fd
gmake clean
gmake
mv sref_estimate_bias ${EXECsref}/.
gmake clean

# sref_bufr.fd
cd ${SORCsref}/sref_bufr.fd
gmake clean
gmake
mv sref_bufr ${EXECsref}/.
gmake clean

# ENSADD
cd ${SORCsref}/global_ensadd.fd
gmake clobber
gmake
mv global_ensadd ${EXECsref}/.
gmake clean

# POSTGP
cd ${SORCsref}/global_postgp.fd
gmake clobber
./make.sh
mv global_postgs ${EXECsref}/.
gmake clobber

###
### WRFBUCKET
###
cd ${SORCsref}/sref_wrfbucket.fd
gmake -f Makefile clean
gmake -f Makefile
mv sref_pcpbucket_g212 ${EXECsref}/.
mv sref_pcpbucket_g216 ${EXECsref}/.
mv sref_pcpbucket_g221 ${EXECsref}/.
mv sref_pcpbucket_g243 ${EXECsref}/.
mv sref_pcpbucket_g132 ${EXECsref}/.
gmake -f Makefile clean

# sref_ens_gen.grib2.fd
cd ${SORCsref}/sref_ens_gen.grib2.fd
gmake clean
gmake
mv sref_ensprod ${EXECsref}/.
gmake clean

# sref_ens_gen.fd
cd ${SORCsref}/sref_ens_gen.fd
gmake -f makefile_REG clean
gmake -f makefile_REG
mv sref_ens_gen ${EXECsref}/.
gmake -f makefile_REG clean
gmake -f makefile_DS clean
gmake -f makefile_DS
mv sref_ens_gen_DS ${EXECsref}/.
gmake -f makefile_DS clean

# sref_dwnsfcst.fd
cd ${SORCsref}/sref_dwnsfcst.fd
gmake clean
gmake
mv sref_dvrtma_debias ${EXECsref}/.
gmake clean

# sref_dwnsvect.fd
cd ${SORCsref}/sref_dwnsvect.fd
gmake clean
gmake
mv sref_dvrtma_bias ${EXECsref}/.
gmake clean

# global_postgp.fd
cd ${SORCsref}/global_postgp.fd
gmake clean
./make.sh
mv global_postgs ${EXECsref}/.
gmake clean

# sref_fastcopygb.fd
cd ${SORCsref}/sref_fastcopygb.fd
gmake clean
gmake
mv fastcopygb ${EXECsref}/.
gmake clean

# sref_meansndp.fd
cd ${SORCsref}/sref_meansndp.fd
gmake clobber
gmake
mv sref_meansndp ${EXECsref}/.
gmake clean

# sref_sndp.fd
cd ${SORCsref}/sref_sndp.fd
rm -f *.o
./make_all
mv sref_*sndp ${EXECsref}/.
rm -f *.o

cd ${SORCsref}/sref_stnmlist.fd
gmake clean
gmake
mv sref_stnmlist ${EXECsref}/.
gmake clean

# sref_memberranking.fd
cd ${SORCsref}/sref_memberranking.fd
gmake clean
gmake
mv sref_ranking ${EXECsref}/.
gmake clean

# sref_qpfbiasestimate.fd
cd ${SORCsref}/sref_qpfbiasestimate.fd
gmake -f makefile_meanqpf clean
gmake -f makefile_meanqpf
mv sref_estimate_meanqpfbias ${EXECsref}/.
gmake -f makefile_meanqpf clean
gmake -f makefile_qpf clean
gmake -f makefile_qpf
mv sref_estimate_qpfbias ${EXECsref}/.
gmake -f makefile_qpf clean

# sref_qpfcalfcsterr.fd
cd ${SORCsref}/sref_qpfcalfcsterr.fd
gmake -f makefile_meanqpf clean
gmake -f makefile_meanqpf
mv sref_cal_meanqpffcsterr ${EXECsref}/.
gmake -f makefile_meanqpf clean
gmake -f makefile_qpf clean
gmake -f makefile_qpf
mv sref_cal_qpffcsterr ${EXECsref}/.
gmake -f makefile_qpf clean

cd ${SORCsref}/sref_wgtmkr_nmmb.fd
gmake clean
gmake
mv wgtmkr_nmmb.x ${EXECsref}/.
gmake clean

cd ${SORCsref}/sref_wgtmkr.fd
gmake -f Makefile_dell clean
gmake -f Makefile_dell
mv wgtmkr.x ${EXECsref}/.
gmake -f Makefile_dell clean

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
./compile nps

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
./compile wps

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

if ${BUILD_COLDSTART}; then
###
### COLDSTART
###
  cd ${SORCsref}/sref_coldstart.fd
  gmake clean
  ./make.sh
  mv sref_coldstart_wrf     ${EXECsref}/.
fi

if ${BUILD_POST}; then
###
### POST
###
cd ${SORCsref}/sref_post.fd
make -f makefile_dell clean
make -f makefile_dell
mv sref_post     ${EXECsref}/.
make -f makefile_dell clean

fi

if ${BUILD_PRDGEN}; then
###
### PRDGEN
###
  cd ${SORCsref}/sref_prdgen.fd
  gmake -f Makefile_dell clean
  gmake -f Makefile_dell
  mv sref_prdgen     ${EXECsref}/.
  gmake -f Makefile_dell clean
fi


date
