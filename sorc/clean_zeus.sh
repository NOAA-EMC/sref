#!/bin/sh

set -aeux

date

SORCsref=`pwd`
EXECsref=`pwd`/../exec


if [ 1 == 1 ]; then
###
### ARW and NMM WRF
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
cd ${SORCsref}/sref_nems.fd/src
cp conf/configure.nems.Zeus.intel_12 conf/configure.nems
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
make -f Makefile_zeus clean
cd ..

fi


if [ 1 == 1 ]; then
###
### PRDGEN
###
cd ${SORCsref}/sref_prdgen.fd
make -f Makefile_zeus clean
cd ..

fi


if [ 1 == 1 ]; then
###
### WRFBUCKET
###
cd ${SORCsref}/sref_wrfbucket.fd
make -f Makefile clean
cd ..

fi


date
