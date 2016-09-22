#!/bin/sh

set -aeux

SREFDIR=/meso/save/${LOGNAME}/sref.v7.0.0
cd $SREFDIR/run

cyc=$1
ymdh=`cat /com/date/t${cyc}z | cut -c7-16`

ymdh=20150121${1}

export SMSBIN=${HOME}/sms

PDY=`echo $ymdh | cut -c1-8`
CYC=`echo $ymdh | cut -c9-10`

#rm -rf /ptmp/$LOGNAME/sref
#mkdir -p /ptmp/$LOGNAME/tmpnwprd
#rm -rf /ptmp/$LOGNAME/tmpnwprd/jlogfile_sref

#rm -f *.bsub *.log

cat SREF_GEFS2SREF.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_SREFDIR_:$SREFDIR:g > SREF_GEFS2SREF.bsub

bsub < SREF_GEFS2SREF.bsub
echo "SREF waiting GEFS2SREF" >> /ptmpp1/$LOGNAME/tmpnwprd/jlogfile_sref

