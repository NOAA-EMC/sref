#!/bin/sh

set -aeux

SREFDIR=/meso/save/${LOGNAME}/sref.v7.0.0
cd $SREFDIR/run

cyc=$1
ymdh=`cat /com/date/t${cyc}z | cut -c7-16`

ymdh=20150331${1}

export SMSBIN=${HOME}/sms

PDY=`echo $ymdh | cut -c1-8`
CYC=`echo $ymdh | cut -c9-10`

#rm -rf /ptmpp1/$LOGNAME/sref
#mkdir -p /ptmpp1/$LOGNAME/tmpnwprd
#rm -rf /ptmpp1/$LOGNAME/tmpnwprd/jlogfile_sref

#rm -f *.bsub *.log

cat SREF_CLUSTER.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_SREFDIR_:$SREFDIR:g > SREF_CLUSTER.bsub

bsub < SREF_CLUSTER.bsub
echo "SREF waiting CLUSTER" >> /ptmpp1/$LOGNAME/tmpnwprd/jlogfile_sref

