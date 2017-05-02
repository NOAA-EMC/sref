#!/bin/sh

set -aeux

SREFDIR=/gpfs/hps/emc/meso/save/${LOGNAME}/sref.v7.1.0
cd $SREFDIR/run

cyc=$1
FLENGTH=87
INCR=3
##ymdh=`cat /gpfs/gp2/nco/ops/com/date/t${cyc}z | cut -c7-16`

ymdh=20170502${cyc}

export SMSBIN=${HOME}/sms
export MACHINE=wcoss
export RES=16km
export IOFORM=2

PDY=`echo $ymdh | cut -c1-8`
CYC=`echo $ymdh | cut -c9-10`

rm -rf /gpfs/hps/ptmp/$LOGNAME/sref
mkdir -p /gpfs/hps/ptmp/$LOGNAME/tmpnwprd
rm -rf /gpfs/hps/ptmp/$LOGNAME/tmpnwprd/jlogfile_sref
rm -f *.bsub *.log

for MODEL in ARW NMB; do

for MEMBER in ctl n01 p01 n02 p02 n03 p03 n04 p04 n05 p05 n06 p06; do
#for MEMBER in ctl n01 p01 n02 p02 n03 p03 n04 p04; do

cat SREF_PREP.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_FLENGTH_:$FLENGTH:g | \
    sed s:_INCR_:$INCR:g | \
    sed s:_SREFDIR_:$SREFDIR:g | \
    sed s:_MODEL_:$MODEL:g | \
    sed s:_MEMBER_:$MEMBER:g | \
    sed s:_IOFORM_:$IOFORM:g | \
    sed s:_RES_:$RES:g | \
    sed s:_MACHINE_:$MACHINE:g > SREF_PREP_${MODEL}_${MEMBER}.bsub

if [ $MEMBER != ctl ]; then
echo "bsub < $SREFDIR/run/SREF_PREP_${MODEL}_${MEMBER}.bsub" >> SREF_PREP_${MODEL}_ctl.bsub
echo "echo \"SREF waiting ${MODEL} ${MEMBER} PREP\" >> /gpfs/hps/ptmp/$LOGNAME/tmpnwprd/jlogfile_sref" >> SREF_PREP_${MODEL}_ctl.bsub
fi

cat SREF_REAL.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_FLENGTH_:$FLENGTH:g | \
    sed s:_INCR_:$INCR:g | \
    sed s:_SREFDIR_:$SREFDIR:g | \
    sed s:_MODEL_:$MODEL:g | \
    sed s:_MEMBER_:$MEMBER:g | \
    sed s:_IOFORM_:$IOFORM:g | \
    sed s:_RES_:$RES:g | \
    sed s:_MACHINE_:$MACHINE:g > SREF_REAL_${MODEL}_${MEMBER}.bsub
cat SREF_FCST.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_FLENGTH_:$FLENGTH:g | \
    sed s:_INCR_:$INCR:g | \
    sed s:_SREFDIR_:$SREFDIR:g | \
    sed s:_MODEL_:$MODEL:g | \
    sed s:_MEMBER_:$MEMBER:g | \
    sed s:_IOFORM_:$IOFORM:g | \
    sed s:_RES_:$RES:g | \
    sed s:_MACHINE_:$MACHINE:g > SREF_FCST_${MODEL}_${MEMBER}.bsub

# 3 hourly 212, 216, 234 grids
export runflag=3hrly
export OUTGRD=255
cat SREF_POST.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_FLENGTH_:$FLENGTH:g | \
    sed s:_INCR_:$INCR:g | \
    sed s:_SREFDIR_:$SREFDIR:g | \
    sed s:_MODEL_:$MODEL:g | \
    sed s:_MEMBER_:$MEMBER:g | \
    sed s:_IOFORM_:$IOFORM:g | \
    sed s:_RES_:$RES:g | \
    sed s:_MACHINE_:$MACHINE:g | \
    sed s:_runflag_:$runflag:g | \
    sed s:_OUTGRD_:$OUTGRD:g > SREF_POST_${MODEL}_${MEMBER}.bsub

# 1 hourly 212, 216, 243 grids
export runflag=hrly
export OUTGRD=255
cat SREF_POST2.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_FLENGTH_:$FLENGTH:g | \
    sed s:_INCR_:$INCR:g | \
    sed s:_SREFDIR_:$SREFDIR:g | \
    sed s:_MODEL_:$MODEL:g | \
    sed s:_MEMBER_:$MEMBER:g | \
    sed s:_IOFORM_:$IOFORM:g | \
    sed s:_RES_:$RES:g | \
    sed s:_MACHINE_:$MACHINE:g | \
    sed s:_runflag_:$runflag:g | \
    sed s:_OUTGRD_:$OUTGRD:g > SREF_POST2_${MODEL}_${MEMBER}.bsub

# 3 hourly 132 grid
export runflag=3hrly
export OUTGRD=132
cat SREF_POST3.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_FLENGTH_:$FLENGTH:g | \
    sed s:_INCR_:$INCR:g | \
    sed s:_SREFDIR_:$SREFDIR:g | \
    sed s:_MODEL_:$MODEL:g | \
    sed s:_MEMBER_:$MEMBER:g | \
    sed s:_IOFORM_:$IOFORM:g | \
    sed s:_RES_:$RES:g | \
    sed s:_MACHINE_:$MACHINE:g | \
    sed s:_runflag_:$runflag:g | \
    sed s:_OUTGRD_:$OUTGRD:g > SREF_POST3_${MODEL}_${MEMBER}.bsub

# 1 hourly 132 grid
export runflag=hrly
export OUTGRD=132
cat SREF_POST4.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_FLENGTH_:$FLENGTH:g | \
    sed s:_INCR_:$INCR:g | \
    sed s:_SREFDIR_:$SREFDIR:g | \
    sed s:_MODEL_:$MODEL:g | \
    sed s:_MEMBER_:$MEMBER:g | \
    sed s:_IOFORM_:$IOFORM:g | \
    sed s:_RES_:$RES:g | \
    sed s:_MACHINE_:$MACHINE:g | \
    sed s:_runflag_:$runflag:g | \
    sed s:_OUTGRD_:$OUTGRD:g > SREF_POST4_${MODEL}_${MEMBER}.bsub


#if [ ${MEMBER} == ctl ]; then
#CTLJOB=`bsub SREF_PREP_${MODEL}_${MEMBER}.bsub`
#echo "SREF waiting ${MODEL} ${MEMBER} PREP" >> /gpfs/hps/ptmp/Dusan.Jovic/tmpnwprd/jlogfile_sref
#else
#bsub -W depend=afterok:$CTLJOB SREF_PREP_${MODEL}_${MEMBER}.bsub
#echo "SREF waiting ${MODEL} ${MEMBER} PREP" >> /gpfs/hps/ptmp/Dusan.Jovic/tmpnwprd/jlogfile_sref
#fi

done
bsub < SREF_PREP_${MODEL}_ctl.bsub
echo "SREF waiting ${MODEL} ctl PREP" >> /gpfs/hps/ptmp/$LOGNAME/tmpnwprd/jlogfile_sref

done
