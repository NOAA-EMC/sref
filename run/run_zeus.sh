#!/bin/sh

set -aeux

SREFDIR=/scratch2/portfolios/NCEPDEV/meso/save/Dusan.Jovic/sref/zeus_sref.v6.0.0
cd $SREFDIR/run

rm -rf /ptmp/$LOGNAME/sref

cyc=$1
FLENGTH=24
INCR=3
ymdh=`cat /com/date/t${cyc}z | cut -c7-16`

ymdh=20120606${1}

export SMSBIN=${HOME}/sms
export MACHINE=zeus
export RES=25km
export IOFORM=2

PDY=`echo $ymdh | cut -c1-8`
CYC=`echo $ymdh | cut -c9-10`

rm -rf /ptmp/$LOGNAME/sref
mkdir -p /ptmp/Dusan.Jovic/tmpnwprd
rm -rf /ptmp/Dusan.Jovic/tmpnwprd/jlogfile_sref

rm -f *.qsub *.log

for MODEL in NMM ARW NMB ; do
#for MODEL in NMB; do

for MEMBER in ctl n01 p01 n02 p02 n03 p03; do
#for MEMBER in ctl n01 p01 n02 p02; do
#for MEMBER in ctl n02 p02 n03 p03; do
#for MEMBER in ctl n03 p03; do
#for MEMBER in ctl n01 p01 ; do
#for MEMBER in ctl ; do

cat SREF_PREP.qsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_FLENGTH_:$FLENGTH:g | \
    sed s:_INCR_:$INCR:g | \
    sed s:_SREFDIR_:$SREFDIR:g | \
    sed s:_MODEL_:$MODEL:g | \
    sed s:_MEMBER_:$MEMBER:g | \
    sed s:_IOFORM_:$IOFORM:g | \
    sed s:_RES_:$RES:g | \
    sed s:_MACHINE_:$MACHINE:g > SREF_PREP_${MODEL}_${MEMBER}.qsub

if [ $MEMBER != ctl ]; then
echo "qsub $SREFDIR/run/SREF_PREP_${MODEL}_${MEMBER}.qsub" >> SREF_PREP_${MODEL}_ctl.qsub
echo "echo \"SREF waiting ${MODEL} ${MEMBER} PREP\" >> /ptmp/Dusan.Jovic/tmpnwprd/jlogfile_sref" >> SREF_PREP_${MODEL}_ctl.qsub
fi


cat SREF_REAL.qsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_FLENGTH_:$FLENGTH:g | \
    sed s:_INCR_:$INCR:g | \
    sed s:_SREFDIR_:$SREFDIR:g | \
    sed s:_MODEL_:$MODEL:g | \
    sed s:_MEMBER_:$MEMBER:g | \
    sed s:_IOFORM_:$IOFORM:g | \
    sed s:_RES_:$RES:g | \
    sed s:_MACHINE_:$MACHINE:g > SREF_REAL_${MODEL}_${MEMBER}.qsub

cat SREF_FCST.qsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_FLENGTH_:$FLENGTH:g | \
    sed s:_INCR_:$INCR:g | \
    sed s:_SREFDIR_:$SREFDIR:g | \
    sed s:_MODEL_:$MODEL:g | \
    sed s:_MEMBER_:$MEMBER:g | \
    sed s:_IOFORM_:$IOFORM:g | \
    sed s:_RES_:$RES:g | \
    sed s:_MACHINE_:$MACHINE:g > SREF_FCST_${MODEL}_${MEMBER}.qsub

cat SREF_POST.qsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_FLENGTH_:$FLENGTH:g | \
    sed s:_INCR_:$INCR:g | \
    sed s:_SREFDIR_:$SREFDIR:g | \
    sed s:_MODEL_:$MODEL:g | \
    sed s:_MEMBER_:$MEMBER:g | \
    sed s:_IOFORM_:$IOFORM:g | \
    sed s:_RES_:$RES:g | \
    sed s:_MACHINE_:$MACHINE:g > SREF_POST_${MODEL}_${MEMBER}.qsub


#if [ ${MEMBER} == ctl ]; then
#CTLJOB=`qsub SREF_PREP_${MODEL}_${MEMBER}.qsub`
#echo "SREF waiting ${MODEL} ${MEMBER} PREP" >> /ptmp/Dusan.Jovic/tmpnwprd/jlogfile_sref
#else
#qsub -W depend=afterok:$CTLJOB SREF_PREP_${MODEL}_${MEMBER}.qsub
#echo "SREF waiting ${MODEL} ${MEMBER} PREP" >> /ptmp/Dusan.Jovic/tmpnwprd/jlogfile_sref
#fi

done

qsub SREF_PREP_${MODEL}_ctl.qsub
echo "SREF waiting ${MODEL} ctl PREP" >> /ptmp/Dusan.Jovic/tmpnwprd/jlogfile_sref

done
