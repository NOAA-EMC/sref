#! /bin/ksh

export WRF_NMM_CORE=0
export WRF_EM_CORE=1
export WRFIO_NCD_LARGE_FILE_SUPPORT=1

TARGDIR=../../exec

./clean -a

( echo 32 ; echo 1 ) | ./configure

./compile em_real > compile_arw.log 2>&1

cp ./main/real.exe $TARGDIR/sref_real_arw
cp ./main/wrf.exe  $TARGDIR/sref_wrf_arw
