#!/bin/sh

set -aeux

SREFDIR=/meso/save/${LOGNAME}/sref.v7.0.0
cd $SREFDIR/run

cyc=$1
#FLENGTH=87
#INCR=3
ymdh=`cat /com/date/t${cyc}z | cut -c7-16`

ymdh=20150324${1}

export SMSBIN=${HOME}/sms
export MACHINE=wcoss

PDY=`echo $ymdh | cut -c1-8`
CYC=`echo $ymdh | cut -c9-10`

#rm -rf /ptmpp1/$LOGNAME/sref
#mkdir -p /ptmpp1/$LOGNAME/tmpnwprd
#rm -rf /ptmpp1/$LOGNAME/tmpnwprd/jlogfile_sref

#rm -f *.bsub *.log

# In operation, ctl and others are run separately: ctl after GFS 90hr finishes; others after GEFS 90hr finishes.
for MEMBER in ctl c00 p01 p02 p03 p04 p05 p06 p07 p08 p09 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20; do

#cat SREF_GFSPREP.bsub.in | \
#    sed s:_PDY_:$PDY:g | \
#    sed s:_CYC_:$CYC:g | \
#    sed s:_FLENGTH_:$FLENGTH:g | \
#    sed s:_INCR_:$INCR:g | \
#    sed s:_SREFDIR_:$SREFDIR:g | \
#    sed s:_MEMBER_:$MEMBER:g | \
#    sed s:_MACHINE_:$MACHINE:g > SREF_GFSPREP_${MEMBER}.bsub
#bsub < SREF_GFSPREP_${MEMBER}.bsub

cat SREF_GFSPREP.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_SREFDIR_:$SREFDIR:g | \
    sed s:_MEMBER_:$MEMBER:g | \
    sed s:_MACHINE_:$MACHINE:g > SREF_GFSPREP_${MEMBER}.bsub
bsub < SREF_GFSPREP_${MEMBER}.bsub

#cat SREF_GFSPREP1.bsub.in | \
#    sed s:_PDY_:$PDY:g | \
#    sed s:_CYC_:$CYC:g | \
#    sed s:_SREFDIR_:$SREFDIR:g | \
#    sed s:_MEMBER_:$MEMBER:g | \
#    sed s:_MACHINE_:$MACHINE:g > SREF_GFSPREP1_${MEMBER}.bsub
#bsub < SREF_GFSPREP1_${MEMBER}.bsub

#cat SREF_GFSPREP2.bsub.in | \
#    sed s:_PDY_:$PDY:g | \
#    sed s:_CYC_:$CYC:g | \
#    sed s:_SREFDIR_:$SREFDIR:g | \
#    sed s:_MEMBER_:$MEMBER:g | \
#    sed s:_MACHINE_:$MACHINE:g > SREF_GFSPREP2_${MEMBER}.bsub
#bsub < SREF_GFSPREP2_${MEMBER}.bsub

done

