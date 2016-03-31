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

#rm -rf /ptmp/$LOGNAME/sref
#mkdir -p /ptmp/$LOGNAME/tmpnwprd
#rm -rf /ptmp/$LOGNAME/tmpnwprd/jlogfile_sref

#rm -f *.bsub *.log

cat SREF_CALFCSTERR.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_SREFDIR_:$SREFDIR:g > SREF_CALFCSTERR.bsub

cat SREF_BIAS.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_SREFDIR_:$SREFDIR:g > SREF_BIAS.bsub

cat SREF_QPFBIAS.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_SREFDIR_:$SREFDIR:g > SREF_QPFBIAS.bsub
bsub < SREF_CALFCSTERR.bsub
echo "SREF waiting CALFCSTERR" >> /ptmpp1/$LOGNAME/tmpnwprd/jlogfile_sref

#bsub < SREF_BIAS.bsub
#echo "SREF waiting BIAS" >> /ptmp/$LOGNAME/tmpnwprd/jlogfile_sref

#bsub < SREF_QPFBIAS.bsub
#echo "SREF waiting QPFBIAS" >> /ptmp/$LOGNAME/tmpnwprd/jlogfile_sref

