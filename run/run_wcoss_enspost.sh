#!/bin/sh

set -aeux

SREFDIR=/meso/save/${LOGNAME}/sref.v7.0.0
cd $SREFDIR/run

#type=1hrly
 type=3hrly
#type=16km
#type=all

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

if [ $type = 3hrly -o $type = all ]; then
cat SREF_ENSPOST.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_SREFDIR_:$SREFDIR:g > SREF_ENSPOST.bsub
bsub < SREF_ENSPOST.bsub
echo "SREF waiting ENSPOST" >> /ptmpp1/$LOGNAME/tmpnwprd/jlogfile_sref
fi

if [ $type = 1hrly -o $type = all ]; then
cat SREF_ENSPOST_1hrly.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_SREFDIR_:$SREFDIR:g > SREF_ENSPOST_1hrly.bsub
bsub < SREF_ENSPOST_1hrly.bsub
echo "SREF waiting ENSPOST_1hrly" >> /ptmpp1/$LOGNAME/tmpnwprd/jlogfile_sref
fi

if [ $type = 16km -o $type = all ]; then
 for ORDER in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30; do
#for ORDER in 1 2 3 4 5 6 7 8 9 10; do
#for ORDER in 1 2 3 4 5 6 7 8; do
#for ORDER in 1; do
cat SREF_ENSPOST132.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_ORDER_:$ORDER:g | \
    sed s:_SREFDIR_:$SREFDIR:g > SREF_ENSPOST132_$ORDER.bsub
bsub < SREF_ENSPOST132_$ORDER.bsub
echo "SREF waiting ENSPOST132_$ORDER" >> /ptmpp1/$LOGNAME/tmpnwprd/jlogfile_sref
done
fi

