#!/bin/sh

set -aeux

SREFDIR=/gpfs/dell2/emc/modeling/noscrub/${LOGNAME}/sref.v7_dell
cd $SREFDIR/run

cyc=$1
ymdh=`cat /gpfs/hps/nco/ops/com/date/t${cyc}z | cut -c7-16`

ymdh=20190318${1}

export SMSBIN=${HOME}/sms

PDY=`echo $ymdh | cut -c1-8`
CYC=`echo $ymdh | cut -c9-10`

#rm -rf /gpfs/dell2/ptmp/$LOGNAME/sref
#mkdir -p /gpfs/dell2/ptmp/$LOGNAME/tmpnwprd
#rm -rf /gpfs/dell2/ptmp/$LOGNAME/tmpnwprd/jlogfile_sref

#rm -f *.bsub *.log

cat SREF_DOWNSCALING.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_SREFDIR_:$SREFDIR:g > SREF_DOWNSCALING.bsub

bsub < SREF_DOWNSCALING.bsub
echo "SREF waiting DOWNSCALING" >> /gpfs/dell2/ptmp/$LOGNAME/tmpnwprd/jlogfile_sref

