#!/bin/sh

set -eux

base=`pwd`/..

#module load /nwpara2/modulefiles/SREF/v7.1.0_wcoss
module load $base/modulefiles/SREF/v7.1.0_wcoss
module list

date

#######################################

SORCsref=`pwd`
EXECsref=`pwd`/../exec


if [ 1 == 1 ]; then
###
### ARW WRF
###

cd ${SORCsref}/sref_wrf_v3.5.1.fd
./clean -a
cd ..

fi

if [ 1 == 1 ]; then
###
### WPS
###

cd ${SORCsref}/sref_wps_v3.5.1.fd
./clean -a
rm -f fort_netcdf.f f_test.f90 cf_test fort_netcdf c_test.c netcdf.inc
cd ..

fi


if [ 1 == 1 ]; then
###
### NPS
###
cd ${SORCsref}/sref_nps.fd
./clean -a
cd ..

fi

if [ 1 == 1 ]; then
###
### NMMB
###
cd ${SORCsref}/sref_nmmb.fd/src
cp conf/configure.nems.Wcoss.intel conf/configure.nems
gmake clean
cd ..

fi


if [ 1 == 1 ]; then
###
### PERTURB
###
cd ${SORCsref}/sref_perturb.fd
cd dio
rm -f *.o *.mod
cd ..
make -f Makefile clean
cd ..

fi


if [ 1 == 1 ]; then
###
### POST
###
cd ${SORCsref}/sref_post.fd
make -f makefile clean
cd ..

fi

if [ 1 == 1 ]; then
###
### COLDSTART
###
cd ${SORCsref}/sref_coldstart.fd
make -f makefile clean
cd ..

fi

if [ 1 == 1 ]; then
###
### WGTMKR
###
cd ${SORCsref}/sref_wgtmkr.fd
make -f Makefile_wcoss clean
cd ..

fi


if [ 1 == 1 ]; then
###
### PRDGEN
###
cd ${SORCsref}/sref_prdgen.fd
make -f Makefile_wcoss clean
cd ..

fi


if [ 1 == 1 ]; then
###
### WRFBUCKET
###
cd ${SORCsref}/sref_wrfbucket.fd
make -f Makefile_wcoss clean
cd ..

fi

# ENSADD
if [ 1 == 1 ]; then
cd ${SORCsref}/global_ensadd.fd
make -f makefile_wcoss clean

cd ..
fi

# global_postgs
if [ 1 == 1 ]; then
cd ${SORCsref}/global_postgp.fd
make -f makefile_wcoss clean

cd ..
fi

# sref_biasestimate.fd

cd ${SORCsref}/sref_biasestimate.fd
make -f makefile_wcoss clean
cd ..

# sref_bufr.fd
cd ${SORCsref}/sref_bufr.fd
make -f makefile_wcoss clean
cd ..

# sref_calfcsterr.fd
cd ${SORCsref}/sref_calfcsterr.fd
make -f makefile_wcoss clean
cd ..

# sref_cluster_NCEP.fd
cd ${SORCsref}/sref_cluster_NCEP.fd
rm -f *.o
make -f makefile_cluster_wcoss clean
cd ..

# sref_cluster_OU.fd
cd ${SORCsref}/sref_cluster_OU.fd
rm -f *.o
make -f makefile_wcoss clean
cd ..

# sref_dwnsfcst.fd
cd ${SORCsref}/sref_dwnsfcst.fd
make -f makefile_wcoss clean
cd ..

# sref_dwnsvect.fd
cd ${SORCsref}/sref_dwnsvect.fd
make -f makefile_wcoss clean
cd ..

# sref_ens_gen.fd
cd ${SORCsref}/sref_ens_gen.fd
rm -f *.o *.mod
make -f makefile_REG_wcoss clean
make -f makefile_DS_wcoss clean
cd ../

# sref_fastcopygb.fd
cd ${SORCsref}/sref_fastcopygb.fd
make -f makefile_wcoss clean
cd ..

# sref_meansndp.fd
cd ${SORCsref}/sref_meansndp.fd
make -f makefile_wcoss clean
cd ..

# sref_memberranking.fd
cd ${SORCsref}/sref_memberranking.fd
make -f makefile_wcoss clean
cd ..

# sref_qpfbiasestimate.fd
cd ${SORCsref}/sref_qpfbiasestimate.fd
make -f makefile_meanqpf_wcoss clean
make -f makefile_qpf_wcoss clean
cd ..

# sref_qpfcalfcsterr.fd
cd ${SORCsref}/sref_qpfcalfcsterr.fd
rm -f *.o
make -f makefile_meanqpf_wcoss clean
make -f makefile_qpf_wcoss clean
cd ..

# sref_sndp.fd
cd ${SORCsref}/sref_sndp.fd
rm -f *.o
make -f makefile_em_wcoss clean
make -f makefile_nmmb_wcoss clean
rm -f *.o
cd ..


rm -f ${EXECsref}/*

date
