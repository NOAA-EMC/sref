#! /bin/ksh
################################################################################
####  UNIX Script Documentation Block
#                      .                                             .
# Script name:         sref_wrfprdgen.sh
# Script description:  Runs WRF-NMM or WRF-EM product generator
#
# Author:        Eric Rogers       Org: NP22         Date: 2004-07-07
#
# Script history log:
# 2003-11-01  Matt Pyle - Original script for parallel
# 2004-07-07  Eric Rogers - Preliminary modifications for production.
# 2005-08-02  Jun Du - Modified to fit into SREF system
# 2005-01-25  Jun Du - changed to use a unified version of prdgen
# 2008-07-18  Jun Du - upgraded to match WRFv2.2
# 2008-08-13  Jun Du - (1) change output frequency from 3hrly to 1hrly for 
#                      grid212 in 00-39hr period;
#                      (2) corrected a bug: mistakenly using grid243 files as for 
#                      grid212 files in converting to grib2 format
# 2010-05-07  Jun Du - generalized the script to run NEMS ensemble (MODEL=NMB)
# 2010-07-09  Jun Du - output grid221 data in grib2 and alert
# 2012-04-16  Jun Du - Added a new 16km NA grid 132: double resolution of grid 221
# 2015-02-20  Jun Du - Added individual member ID to pgrb1 file
####################################################################
set -aux

export EXECGLOBAL=${NWROOTp1}/exec

fhr=$1
CYC=$2
MODEL=$3
MEMBER=$4
PDY=$5
OUTGRD=$6

#if [ $OUTGRD -ne 132 ]; then
# OUTGRD=
#fi

PAIR=$MEMBER
if [ $MODEL = NMM ];then CORE=nmm;fi
#if [ $MODEL = ARW ];then CORE=em;fi
if [ $MODEL = ARW ];then CORE=arw;fi
if [ $MODEL = NMB ];then CORE=nmb;fi

if [ $MEMBER = ctl ];then 
 PAIR=ctl;
(( e1 = 1 ))
(( e2 = 2 ))
fi
if [ $MEMBER = n01 ];then 
 PAIR=n1;
(( e1 = 2 ))
(( e2 = 1 ))
fi
if [ $MEMBER = n02 ];then 
 PAIR=n2;
(( e1 = 2 ))
(( e2 = 2 ))
fi
if [ $MEMBER = n03 ];then 
 PAIR=n3;
(( e1 = 2 ))
(( e2 = 3 ))
fi
if [ $MEMBER = n04 ];then 
 PAIR=n4;
(( e1 = 2 ))
(( e2 = 4 ))
fi
if [ $MEMBER = n05 ];then 
 PAIR=n5;
(( e1 = 2 ))
(( e2 = 5 ))
fi
if [ $MEMBER = n06 ];then 
 PAIR=n6;
(( e1 = 2 ))
(( e2 = 6 ))
fi
if [ $MEMBER = p01 ];then 
 PAIR=p1;
(( e1 = 3 ))
(( e2 = 1 ))
fi
if [ $MEMBER = p02 ];then 
 PAIR=p2;
(( e1 = 3 ))
(( e2 = 2 ))
fi
if [ $MEMBER = p03 ];then 
 PAIR=p3;
(( e1 = 3 ))
(( e2 = 3 ))
fi
if [ $MEMBER = p04 ];then 
 PAIR=p4;
(( e1 = 3 ))
(( e2 = 4 ))
fi
if [ $MEMBER = p05 ];then 
 PAIR=p5;
(( e1 = 3 ))
(( e2 = 5 ))
fi
if [ $MEMBER = p06 ];then 
 PAIR=p6;
(( e1 = 3 ))
(( e2 = 6 ))
fi

YYYY=`echo $PDY | cut -c1-4`
MM=`echo $PDY | cut -c5-6`
DD=`echo $PDY | cut -c7-8`
CYCLE=$YYYY$MM$DD$CYC

export XLFRTEOPTS="unit_vars=yes"

msg="JOB $job FOR WRF-${MODEL} NEST=$MODEL HAS BEGUN"
postmsg "$jlogfile" "$msg"

filedir=$WORK_DIR
#bindir=$utilexec

startd=$YYYY$MM$DD
startdate=$CYCLE

mkdir -p $COMOUT

date=`$NDATE $fhr $CYCLE`

mkdir -p $WORK_DIR/wrf_prdgen
cd $WORK_DIR/wrf_prdgen

export fhr
export tmmark=tm00

#
# WRF-EM post output for Hawaii/Puerto Rico are on mercator grid, while WRF-NMM is on lat/lon grid.
# WRF-EM post uses the same grid ID number as the WRF-NMM output. To make lat/lon grids from 
# WRF-EM runs do the following:
#
# 1) Use overgridnum_grib to make a temporary GRIB file with a new grid ID number
# (252 for Hawaii, 253 for Puerto Rico)
#
# 2) Run wrf_prdgen to interpolate the Hawaii/Puerto Rico mercator grid to the
# WRF-NMM Hawaii (grid #250) and WRF-NMM (grid #248) lat/lon grids
#
# 3) These WRF-EM lat/lon grids will be used to compute the ensemble mean
#

# make GRIB file with pressure data every 25 mb
cp $PARMsref/sref_wrf_master.ctl${OUTGRD}_$MODEL  master.ctl

#if [ $MODEL = NMB ];then
#ln -s $WORK_DIR/NMMPRS.GrbF${fhr} .
#cat >input${fhr}.prd <<EOF5
#NMMPRS.GrbF${fhr}
#EOF5
#else
ln -s $WORK_DIR/WRFPRS.GrbF${fhr} .
cat >input${fhr}.prd <<EOF5
WRFPRS.GrbF${fhr}
EOF5
#fi

rm -f fort.*
ln -sf master.ctl                     fort.10
if [ $OUTGRD -eq 132 ]; then
 ln -sf $FIXsref/sref_wgt_${CORE}_132_$RES  fort.21
else
 ln -sf $FIXsref/sref_wgt_${CORE}_212_$RES  fort.21
 ln -sf $FIXsref/sref_wgt_${CORE}_216_$RES  fort.22
 ln -sf $FIXsref/sref_wgt_${CORE}_221_$RES  fort.23
 ln -sf $FIXsref/sref_wgt_${CORE}_243_$RES  fort.24
fi
$EXECsref/sref_prdgen < input${fhr}.prd > prdgen.out${fhr}

echo "done executing prdgen" > $WORK_DIR/prdgendone${fhr}

# if need to save native grid in grib
#cp $WORK_DIR/wrfpost/WRFPRS${fhr}.tm00 $COMOUT/$CORE.t${CYC}z.wrfprs${fhr}.tm00

# compute 3-h precip buckets for ARW core (NMM calculates within the model)
###fhrprev=`expr $fhr - 3`
fhrprev=$(( $fhr - 3 ))

if [ $fhrprev -lt 10 ]
then
fhrprev=0$fhrprev
fi

if [ $OUTGRD -eq 132 ]; then

if [ $fhr -eq 00 ] ; then
# save grids
# cp meso.AWIP16${fhr}.tm00 $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}
  rm -f temp.pgrb temp.ipgrb epgrb
  cp meso.AWIP16${fhr}.tm00 temp.pgrb
  ${GRBINDEX} temp.pgrb temp.ipgrb
  $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb
  cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}
  sync

else

 if [ $fhr -eq 01 -o $fhr -eq 02 -o $fhr -eq 04 -o $fhr -eq 05 -o $fhr -eq 07 -o $fhr -eq 08 -o $fhr -eq 10 -o $fhr -eq 11 -o \
      $fhr -eq 13 -o $fhr -eq 14 -o $fhr -eq 16 -o $fhr -eq 17 -o $fhr -eq 19 -o $fhr -eq 20 -o $fhr -eq 22 -o $fhr -eq 23 -o \
      $fhr -eq 25 -o $fhr -eq 26 -o $fhr -eq 28 -o $fhr -eq 29 -o $fhr -eq 31 -o $fhr -eq 32 -o $fhr -eq 34 -o $fhr -eq 35 -o \
      $fhr -eq 37 -o $fhr -eq 38 ];then
# 132 grid
  if [ $MODEL = em -o $MODEL = ARW ]; then
    rm -f input.card
    echo "$WORK_DIR/wrf_prdgen" > input.card
    echo meso.AWIP16 >> input.card
    echo $fhrprev >> input.card
    echo $fhr >> input.card
    $EXECsref/sref_pcpbucket_g132 < input.card
#   cat meso.AWIP16${fhr}.tm00 WRFPCP${fhr}.tm00 > sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}
#   cp sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr} $COMOUT/.
    rm -f temp.pgrb temp.ipgrb epgrb
    cat meso.AWIP16${fhr}.tm00 WRFPCP${fhr}.tm00 > temp.pgrb
    rm WRFPCP${fhr}.tm00
    ${GRBINDEX} temp.pgrb temp.ipgrb
    $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb
    cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}
    sync
  else
#   cp meso.AWIP16${fhr}.tm00 $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}
    rm -f temp.pgrb temp.ipgrb epgrb
    cp meso.AWIP16${fhr}.tm00 temp.pgrb
    ${GRBINDEX} temp.pgrb temp.ipgrb
    $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb
    cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}
    sync
  fi

 else

# 132 grid
  if [ $MODEL = em -o $MODEL = ARW ]; then
    rm -f input.card
    echo "$WORK_DIR/wrf_prdgen" > input.card
    echo meso.AWIP16 >> input.card
    echo $fhrprev >> input.card
    echo $fhr >> input.card
    $EXECsref/sref_pcpbucket_g132 < input.card
#   cat meso.AWIP16${fhr}.tm00 WRFPCP${fhr}.tm00 > sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}
    rm -f temp.pgrb temp.ipgrb epgrb
    cat meso.AWIP16${fhr}.tm00 WRFPCP${fhr}.tm00 > temp.pgrb
    rm WRFPCP${fhr}.tm00
    ${GRBINDEX} temp.pgrb temp.ipgrb
    $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb
    cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}
    sync
  else
#   cp meso.AWIP16${fhr}.tm00 $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}
    rm -f temp.pgrb temp.ipgrb epgrb
    cp meso.AWIP16${fhr}.tm00 temp.pgrb
    ${GRBINDEX} temp.pgrb temp.ipgrb
    $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb 
    cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}
    sync
  fi

 fi
fi

# convert to GRIB2
if [ $SENDCOM = YES ]
then
  if [ `expr $fhr % 3` -ne 0 ]; then
   ${CNVGRIB} -g12 -p40 $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr} $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}.grib2
   ${WGRIB2} $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}.grib2 -s > $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}.grib2.idx
  else
   ${CNVGRIB} -g12 -p40 $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr} $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}.grib2
   ${WGRIB2} $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}.grib2 -s > $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}.grib2.idx
  fi

   if [ "$SENDDBN" = "YES" ]; then
    if [ `expr $fhr % 3` -ne 0 ]; then
      $DBNROOT/bin/dbn_alert MODEL SREF_${MODEL}_132_PGB_GB2 $job $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}.grib2
      $DBNROOT/bin/dbn_alert MODEL SREF_${MODEL}_132_PGB_GB2_WIDX $job $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}.grib2.idx
    else
      $DBNROOT/bin/dbn_alert MODEL SREF_${MODEL}_132_PGB_GB2 $job $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}.grib2
      $DBNROOT/bin/dbn_alert MODEL SREF_${MODEL}_132_PGB_GB2_WIDX $job $COMOUT/sref_$CORE.t${CYC}z.pgrb132.$PAIR.f${fhr}.grib2.idx
    fi
   fi
fi

####### above is 132 grid #########################
else
####### below is non-132 grid #########################

if [ $fhr -eq 00 ] ; then
# save grids
# cp meso.AWIP3D${fhr}.tm00 $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}
  rm -f temp.pgrb temp.ipgrb epgrb
  cp  meso.AWIP3D${fhr}.tm00 temp.pgrb
  ${GRBINDEX} temp.pgrb temp.ipgrb
  $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb
  cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}
# cp nam.AWIPAK${fhr}.tm00 $COMOUT/sref_$CORE.t${CYC}z.pgrb216.$PAIR.f${fhr}
  rm -f temp.pgrb temp.ipgrb epgrb
  cp  nam.AWIPAK${fhr}.tm00 temp.pgrb
  ${GRBINDEX} temp.pgrb temp.ipgrb
  $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb
  cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb216.$PAIR.f${fhr}
# cp nam.AWIP32${fhr}.tm00 $COMOUT/sref_$CORE.t${CYC}z.pgrb221.$PAIR.f${fhr}
  rm -f temp.pgrb temp.ipgrb epgrb
  cp  nam.AWIP32${fhr}.tm00 temp.pgrb
  ${GRBINDEX} temp.pgrb temp.ipgrb
  $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb
  cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb221.$PAIR.f${fhr}
# cp nam.AWIPHI${fhr}.tm00 $COMOUT/sref_$CORE.t${CYC}z.pgrb243.$PAIR.f${fhr}
  rm -f temp.pgrb temp.ipgrb epgrb
  cp  nam.AWIPHI${fhr}.tm00 temp.pgrb
  ${GRBINDEX} temp.pgrb temp.ipgrb
  $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb
  cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb243.$PAIR.f${fhr}
  sync

else

 if [ $fhr -eq 01 -o $fhr -eq 02 -o $fhr -eq 04 -o $fhr -eq 05 -o $fhr -eq 07 -o $fhr -eq 08 -o $fhr -eq 10 -o $fhr -eq 11 -o \
      $fhr -eq 13 -o $fhr -eq 14 -o $fhr -eq 16 -o $fhr -eq 17 -o $fhr -eq 19 -o $fhr -eq 20 -o $fhr -eq 22 -o $fhr -eq 23 -o \
      $fhr -eq 25 -o $fhr -eq 26 -o $fhr -eq 28 -o $fhr -eq 29 -o $fhr -eq 31 -o $fhr -eq 32 -o $fhr -eq 34 -o $fhr -eq 35 -o \
      $fhr -eq 37 -o $fhr -eq 38 ];then
# 212 grid
  if [ $MODEL = em -o $MODEL = ARW ]; then
    rm -f input.card
    echo "$WORK_DIR/wrf_prdgen" > input.card
    echo meso.AWIP3D >> input.card
    echo $fhrprev >> input.card
    echo $fhr >> input.card
    $EXECsref/sref_pcpbucket_g212 < input.card
#   cat meso.AWIP3D${fhr}.tm00 WRFPCP${fhr}.tm00 > sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}
#   cp sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr} $COMOUT/.
    rm -f temp.pgrb temp.ipgrb epgrb
    cat meso.AWIP3D${fhr}.tm00 WRFPCP${fhr}.tm00 > temp.pgrb
    rm WRFPCP${fhr}.tm00
    ${GRBINDEX} temp.pgrb temp.ipgrb
    $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb
    cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}
    sync
  else
#   cp meso.AWIP3D${fhr}.tm00 $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}
    rm -f temp.pgrb temp.ipgrb epgrb
    cp meso.AWIP3D${fhr}.tm00 temp.pgrb
    ${GRBINDEX} temp.pgrb temp.ipgrb
    $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb
    cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}
    sync
  fi

 else

# 212 grid
  if [ $MODEL = em -o $MODEL = ARW ]; then
    rm -f input.card
    echo "$WORK_DIR/wrf_prdgen" > input.card
    echo meso.AWIP3D >> input.card
    echo $fhrprev >> input.card
    echo $fhr >> input.card
    $EXECsref/sref_pcpbucket_g212 < input.card
#   cat meso.AWIP3D${fhr}.tm00 WRFPCP${fhr}.tm00 > sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}
#   cp sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr} $COMOUT/.
    rm -f temp.pgrb temp.ipgrb epgrb
    cat meso.AWIP3D${fhr}.tm00 WRFPCP${fhr}.tm00 > temp.pgrb
    rm WRFPCP${fhr}.tm00
    ${GRBINDEX} temp.pgrb temp.ipgrb
    $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb
    cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}
    sync
  else
#   cp meso.AWIP3D${fhr}.tm00 $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}
    rm -f temp.pgrb temp.ipgrb epgrb
    cp meso.AWIP3D${fhr}.tm00 temp.pgrb
    ${GRBINDEX} temp.pgrb temp.ipgrb
    $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb 
    cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}
    sync
  fi
# 216 grid
  if [ $MODEL = em -o $MODEL = ARW ]; then
    rm -f input.card
    echo "$WORK_DIR/wrf_prdgen" > input.card
    echo nam.AWIPAK >> input.card
    echo $fhrprev >> input.card
    echo $fhr >> input.card
    $EXECsref/sref_pcpbucket_g216 < input.card
#   cat nam.AWIPAK${fhr}.tm00 WRFPCP${fhr}.tm00 > sref_$CORE.t${CYC}z.pgrb216.$PAIR.f${fhr}
#   cp sref_$CORE.t${CYC}z.pgrb216.$PAIR.f${fhr} $COMOUT/.
    rm -f temp.pgrb temp.ipgrb epgrb
    cat nam.AWIPAK${fhr}.tm00 WRFPCP${fhr}.tm00 > temp.pgrb
    rm WRFPCP${fhr}.tm00
    ${GRBINDEX} temp.pgrb temp.ipgrb
    $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb
    cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb216.$PAIR.f${fhr}
    sync
  else
#   cp nam.AWIPAK${fhr}.tm00 $COMOUT/sref_$CORE.t${CYC}z.pgrb216.$PAIR.f${fhr}
    rm -f temp.pgrb temp.ipgrb epgrb
    cp nam.AWIPAK${fhr}.tm00 temp.pgrb
    ${GRBINDEX} temp.pgrb temp.ipgrb
    $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb 
    cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb216.$PAIR.f${fhr}
    sync
  fi
# 221 grid
  if [ $MODEL = em -o $MODEL = ARW ]; then
    rm -f input.card
    echo "$WORK_DIR/wrf_prdgen" > input.card
    echo nam.AWIP32 >> input.card
    echo $fhrprev >> input.card
    echo $fhr >> input.card
    $EXECsref/sref_pcpbucket_g221 < input.card
#   cat nam.AWIP32${fhr}.tm00 WRFPCP${fhr}.tm00 > sref_$CORE.t${CYC}z.pgrb221.$PAIR.f${fhr}
#   cp sref_$CORE.t${CYC}z.pgrb221.$PAIR.f${fhr} $COMOUT/.
    rm -f temp.pgrb temp.ipgrb epgrb
    cat nam.AWIP32${fhr}.tm00 WRFPCP${fhr}.tm00 > temp.pgrb
    rm WRFPCP${fhr}.tm00
    ${GRBINDEX} temp.pgrb temp.ipgrb
    $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb
    cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb221.$PAIR.f${fhr}
    sync
  else
#   cp nam.AWIP32${fhr}.tm00 $COMOUT/sref_$CORE.t${CYC}z.pgrb221.$PAIR.f${fhr}
    rm -f temp.pgrb temp.ipgrb epgrb
    cp nam.AWIP32${fhr}.tm00 temp.pgrb
    ${GRBINDEX} temp.pgrb temp.ipgrb
    $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb
    cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb221.$PAIR.f${fhr}
    sync
  fi
# 243 grid
  if [ $MODEL = em -o $MODEL = ARW ]; then
    rm -f input.card
    echo "$WORK_DIR/wrf_prdgen" > input.card
    echo nam.AWIPHI >> input.card
    echo $fhrprev >> input.card
    echo $fhr >> input.card
    $EXECsref/sref_pcpbucket_g243 < input.card
#   cat nam.AWIPHI${fhr}.tm00 WRFPCP${fhr}.tm00 > sref_$CORE.t${CYC}z.pgrb243.$PAIR.f${fhr}
#   cp sref_$CORE.t${CYC}z.pgrb243.$PAIR.f${fhr} $COMOUT/.
    rm -f temp.pgrb temp.ipgrb epgrb
    cat nam.AWIPHI${fhr}.tm00 WRFPCP${fhr}.tm00 > temp.pgrb
    rm WRFPCP${fhr}.tm00
    ${GRBINDEX} temp.pgrb temp.ipgrb
    $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb
    cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb243.$PAIR.f${fhr}
    sync
  else
    cp nam.AWIPHI${fhr}.tm00 $COMOUT/sref_$CORE.t${CYC}z.pgrb243.$PAIR.f${fhr}
    rm -f temp.pgrb temp.ipgrb epgrb
    cp nam.AWIPHI${fhr}.tm00 temp.pgrb
    ${GRBINDEX} temp.pgrb temp.ipgrb
    $USHsref/global_ensadd.sh $e1 $e2 temp.pgrb temp.ipgrb epgrb
    cp epgrb $COMOUT/sref_$CORE.t${CYC}z.pgrb243.$PAIR.f${fhr}
    sync
  fi
 fi
fi

# convert to GRIB2
if [ $SENDCOM = YES ]
then
  if [ `expr $fhr % 3` -ne 0 ]; then
   ${CNVGRIB} -g12 -p40 $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr} $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}.grib2
   ${WGRIB2} $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}.grib2 -s > $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}.grib2.idx
  else
   ${CNVGRIB} -g12 -p40 $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr} $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}.grib2
   ${WGRIB2} $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}.grib2 -s > $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}.grib2.idx
   ${CNVGRIB} -g12 -p40 $COMOUT/sref_$CORE.t${CYC}z.pgrb221.$PAIR.f${fhr} $COMOUT/sref_$CORE.t${CYC}z.pgrb221.$PAIR.f${fhr}.grib2
   ${WGRIB2} $COMOUT/sref_$CORE.t${CYC}z.pgrb221.$PAIR.f${fhr}.grib2 -s > $COMOUT/sref_$CORE.t${CYC}z.pgrb221.$PAIR.f${fhr}.grib2.idx
   ${CNVGRIB} -g12 -p40 $COMOUT/sref_$CORE.t${CYC}z.pgrb243.$PAIR.f${fhr} $COMOUT/sref_$CORE.t${CYC}z.pgrb243.$PAIR.f${fhr}.grib2
   ${WGRIB2} $COMOUT/sref_$CORE.t${CYC}z.pgrb243.$PAIR.f${fhr}.grib2 -s > $COMOUT/sref_$CORE.t${CYC}z.pgrb243.$PAIR.f${fhr}.grib2.idx
   ${CNVGRIB} -g12 -p40 $COMOUT/sref_$CORE.t${CYC}z.pgrb216.$PAIR.f${fhr} $COMOUT/sref_$CORE.t${CYC}z.pgrb216.$PAIR.f${fhr}.grib2
   ${WGRIB2} $COMOUT/sref_$CORE.t${CYC}z.pgrb216.$PAIR.f${fhr}.grib2 -s > $COMOUT/sref_$CORE.t${CYC}z.pgrb216.$PAIR.f${fhr}.grib2.idx
  fi

   if [ "$SENDDBN" = "YES" ]; then
    if [ `expr $fhr % 3` -ne 0 ]; then
      $DBNROOT/bin/dbn_alert MODEL SREF_${MODEL}_212_PGB_GB2 $job $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}.grib2
      $DBNROOT/bin/dbn_alert MODEL SREF_${MODEL}_212_PGB_GB2_WIDX $job $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}.grib2.idx
    else
      $DBNROOT/bin/dbn_alert MODEL SREF_${MODEL}_212_PGB_GB2 $job $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}.grib2
      $DBNROOT/bin/dbn_alert MODEL SREF_${MODEL}_212_PGB_GB2_WIDX $job $COMOUT/sref_$CORE.t${CYC}z.pgrb212.$PAIR.f${fhr}.grib2.idx
      $DBNROOT/bin/dbn_alert MODEL SREF_${MODEL}_221_PGB_GB2 $job $COMOUT/sref_$CORE.t${CYC}z.pgrb221.$PAIR.f${fhr}.grib2
      $DBNROOT/bin/dbn_alert MODEL SREF_${MODEL}_221_PGB_GB2_WIDX $job $COMOUT/sref_$CORE.t${CYC}z.pgrb221.$PAIR.f${fhr}.grib2.idx
    fi
   fi
fi

fi

exit
