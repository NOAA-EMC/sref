#! /bin/ksh

export WRF_NMM_CORE=1
export WRF_EM_CORE=0
#export NETCDF=/usrx/local/NetCDF/3.6.3
export WRFIO_NCD_LARGE_FILE_SUPPORT=1

TARGDIR=../../exec

./clean -a

( echo 32 ; echo 1 ) | ./configure

./compile nmm_real > compile_nmm.sc.log 2>&1

cp ./main/real_nmm.exe $TARGDIR/sref_real_nmm
cp ./main/wrf.exe  $TARGDIR/sref_wrf_nmm
