#!/bin/sh

set -aeux

SREFDIR=/gpfs/dell2/emc/modeling/noscrub/${LOGNAME}/sref.v7_dell
cd $SREFDIR/run

cyc=$1
FLENGTH=87
INCR=3
INCRBUFR=1
ymdh=`cat /gpfs/hps/nco/ops/com/date/t${cyc}z | cut -c7-16`

ymdh=20200130${1}

export SMSBIN=${HOME}/sms
export MACHINE=wcoss
export RES=16km
export IOFORM=2
export runflag=3hrly
#export runflag=hrly
export OUTGRD=255
#export OUTGRD=132

PDY=`echo $ymdh | cut -c1-8`
CYC=`echo $ymdh | cut -c9-10`

#rm -rf /gpfs/dell2/ptmp/$LOGNAME/sref
#mkdir -p /gpfs/dell2/ptmp/$LOGNAME/tmpnwprd
#rm -rf /gpfs/dell2/ptmp/$LOGNAME/tmpnwprd/jlogfile_sref

#rm -f *.bsub *.log

for MODEL in SREF ; do

for MEMBER in mean; do

cat SREF_MEANBUFR.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_FLENGTH_:$FLENGTH:g | \
    sed s:_INCR_:$INCRBUFR:g | \
    sed s:_SREFDIR_:$SREFDIR:g | \
    sed s:_MODEL_:$MODEL:g | \
    sed s:_MEMBER_:$MEMBER:g | \
    sed s:_RES_:$RES:g | \
    sed s:_MACHINE_:$MACHINE:g > SREF_BUFR_${MODEL}_${MEMBER}.bsub

bsub < SREF_BUFR_${MODEL}_${MEMBER}.bsub
echo "SREF waiting ${MODEL} ${MEMBER} BUFR" >> /gpfs/dell2/ptmp/$LOGNAME/tmpnwprd/jlogfile_sref

done

done
