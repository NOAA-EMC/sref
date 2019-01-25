#!/bin/sh

set -eux

date

SORCsref=`pwd`
EXECsref=`pwd`/../exec


if [ 1 == 1 ]; then
###
### ARW and NMM WRF
###

cd ${SORCsref}/sref_wrf_v3.5.1.fd
./clean -a
./clean -aa
rm -f compile_arw.sc.log
cd ..

fi

if [ 1 == 1 ]; then
###
### WPS
###

cd ${SORCsref}/sref_wps_v3.5.1.fd
./clean -a
rm -f configure.wps.backup
cd ..

fi


if [ 1 == 1 ]; then
###
### NPS
###
cd ${SORCsref}/sref_nps.fd
./clean -a
rm -f configure.nps.backup
cd ..

fi

if [ 1 == 1 ]; then
###
### NMMB
###
cd ${SORCsref}/sref_nmmb.fd/src
cp conf/configure.nems.dell.intel conf/configure.nems
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
make -f Makefile_dell clean
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
make -f Makefile_dell clean
cd ..

fi


if [ 1 == 1 ]; then
###
### PRDGEN
###
cd ${SORCsref}/sref_prdgen.fd
make -f Makefile_dell clean
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

if [ 1 == 1 ]; then
###
### TRANSUTIL
###
cd ${SORCsref}/transutil.fd/v1.0.0/
rm -f libtransutil_4.a  libtransutil_8.a  libtransutil_d.a
cd ../..

fi


date
