#!/bin/bash

set -eux

FV3HREFDIR=/gpfs/hps/emc/meso/save/${LOGNAME}/fv3href/fv3href_trunk
cd $FV3HREFDIR/run

cyc=$1
FLENGTH=48
##ymdh=`cat /gpfs/gp2/nco/ops/com/date/t${cyc}z | cut -c7-16`

ymdh=20170704${cyc}

PDY=`echo $ymdh | cut -c1-8`
CYC=`echo $ymdh | cut -c9-10`

rm -rf /gpfs/hps/ptmp/$LOGNAME/fv3href
mkdir -p /gpfs/hps/ptmp/$LOGNAME/fv3href
rm -f *.bsub *.log

export MACHINE=wcoss_cray

#for MEMBER in ctl mem1 mem2 mem3 mem4 mem5 mem6 mem7 mem8 mem9 ; do
for MEMBER in ctl mem1 mem2 ; do

cat FV3HREF_PREP.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_FV3HREFDIR_:$FV3HREFDIR:g | \
    sed s:_MEMBER_:$MEMBER:g | \
    sed s:_MACHINE_:$MACHINE:g > FV3HREF_PREP_${MEMBER}.bsub

cat FV3HREF_FCST.bsub.in | \
    sed s:_PDY_:$PDY:g | \
    sed s:_CYC_:$CYC:g | \
    sed s:_FLENGTH_:$FLENGTH:g | \
    sed s:_FV3HREFDIR_:$FV3HREFDIR:g | \
    sed s:_MEMBER_:$MEMBER:g | \
    sed s:_MACHINE_:$MACHINE:g > FV3HREF_FCST_${MEMBER}.bsub

bsub                                       < FV3HREF_PREP_${MEMBER}.bsub
bsub -w "done(\"FV3HREF_PREP_${MEMBER}\")" < FV3HREF_FCST_${MEMBER}.bsub

echo "FV3HREF waiting ${MEMBER} PREP" >> /gpfs/hps/ptmp/$LOGNAME/fv3href/jlogfile_sref
echo "FV3HREF waiting ${MEMBER} FCST" >> /gpfs/hps/ptmp/$LOGNAME/fv3href/jlogfile_sref

done
