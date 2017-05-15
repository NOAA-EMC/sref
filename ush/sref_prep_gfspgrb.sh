#!/bin/ksh
#
######################################################################################
# UNIX Script Documentation Block                                                    #
#                                                                                    #
# Script name:         sref_prep_gfspgrb.sh                                          #
# Script description:  grab and prepare global gprb data for SREF                    #
#                      ensemble                                                      #
#                                                                                    #
# Author:        Jun Du       Org: NP22         Date: 2017-01-10                     #
#                                                                                    #
# Abstract: The script gets global model pgrb data and prepares it for SREF to use   #
#                                                                                    #
# Script history log:                                                                #
# 2017-01-10  Jun Du  - Initial script                                               #
#                                                                                    #
######################################################################################

#set -aux
set -x

DATE=$1
CYC=$2
GLBPAIR=$3
SREFMEM=$4
export fhr=$5
resolution=$6
model=$7

cd $DATA/${GLBPAIR}_$fhr

#
## GFS Files
#
if [ $model = GFS ]; then

 str=`ls -s $COMINgfs/gfs.$DATE/gfs.t${CYC}z.pgrbf90`  #get file size in block (1024 bytes)
     set -A fsize $str
     if [ -s ${COMINgfs}/gfs.$DATE/gfs.t${CYC}z.master.grb2f90 ] && [ ${fsize[0]} -gt 800 ] ; then # fcst file is finished?
       echo  ${COMINgfs}/gfs.$DATE/gfs.t${CYC}z.master.grb2f90 " exist"
       export COMIN_SIG=$COMINgfs/gfs.$DATE
       newCYC=$CYC
       newfhr=$fhr
       echo 'using the most current GFS data'
     else
       if [ $CYC -eq 00 ]; then
        newCYC=18
        export COMIN_SIG=$COMINgfs/gfs.$PDYm1
       fi
       if [ $CYC -eq 06 ]; then
        newCYC=00
        export COMIN_SIG=$COMINgfs/gfs.$DATE
       fi
       if [ $CYC -eq 12 ]; then
        newCYC=06
        export COMIN_SIG=$COMINgfs/gfs.$DATE
       fi
       if [ $CYC -eq 18 ]; then
        newCYC=12
        export COMIN_SIG=$COMINgfs/gfs.$DATE
       fi
       newfhr=`expr $fhr + 6`
       if [ $newfhr -lt 10 ];then newfhr=0$newfhr;fi
       echo 'using an older cycle GFS data'
     fi

#convert from grib2 to grib1 here
       if [ $newfhr -lt 100 ];then
cat ${COMIN_SIG}/gfs.t${newCYC}z.pgrb2.0p25.f0${newfhr} ${COMIN_SIG}/gfs.t${newCYC}z.pgrb2b.0p25.f0$newfhr > temp_grib2
       else
cat ${COMIN_SIG}/gfs.t${newCYC}z.pgrb2.0p25.f${newfhr} ${COMIN_SIG}/gfs.t${newCYC}z.pgrb2b.0p25.f$newfhr > temp_grib2
       fi
$CNVGRIB -g21 temp_grib2 $GFSOUT/gfs.t${newCYC}z.master.grbf$newfhr
#$CNVGRIB -g21 ${COMIN_SIG}/gfs.t${newCYC}z.master.grb2f${newfhr} $GFSOUT/gfs.t${newCYC}z.master.grbf$newfhr
#cp ${COMIN_SIG}/gfs.t${newCYC}z.pgrbf${newfhr} $GFSOUT/gfs.t${newCYC}z.master.grbf$newfhr

fi

if [ $model = gens ]; then
# Note: Currently not used because the GEFS has less vertical pressure levels in pgrb files (25mb interval is needed)

 str=`ls -s $COMINgens/gefs.$DATE/${CYC}/pgrb2bp5/gep20.t${CYC}z.pgrb2b.0p50.f090`  #get file size in block (1024 bytes)
     set -A fsize $str
     if [ -s ${COMINgens}/gefs.$DATE/${CYC}/pgrb2bp5/gep20.t${CYC}z.pgrb2b.0p50.f090 ] && [ ${fsize[0]} -gt 800 ] ; then # fcst file is finished?
       echo  ${COMINgens}/gefs.$DATE/${CYC}/pgrb2bp5/gep20.t${CYC}z.pgrb2b.0p50.f090 " exist"

       export COMIN_SIG=$COMINgens/gefs.$DATE
       newCYC=$CYC
       newfhr=$fhr
       echo 'using the most current GEFS data'
     else
       if [ $CYC -eq 00 ]; then
        newCYC=18
        export COMIN_SIG=$COMINgens/gefs.$PDYm1
       fi
       if [ $CYC -eq 06 ]; then
        newCYC=00
        export COMIN_SIG=$COMINgens/gefs.$DATE
       fi
       if [ $CYC -eq 12 ]; then
        newCYC=06
        export COMIN_SIG=$COMINgens/gefs.$DATE
       fi
       if [ $CYC -eq 18 ]; then
        newCYC=12
        export COMIN_SIG=$COMINgens/gefs.$DATE
       fi
       newfhr=`expr $fhr + 6`
       if [ $newfhr -lt 10 ];then newfhr=0$newfhr;fi
       echo 'using an older cycle GEFS data'
     fi

       if [ $newfhr -lt 100 ];then
  cat ${COMIN_SIG}/${newCYC}/pgrb2ap5/ge${GLBPAIR}.t${newCYC}z.pgrb2a.0p50.f0$newfhr ${COMIN_SIG}/${newCYC}/pgrb2bp5/ge${GLBPAIR}.t${newCYC}z.pgrb2b.0p50.f0$newfhr > tempfile_grib2
       else
  cat ${COMIN_SIG}/${newCYC}/pgrb2ap5/ge${GLBPAIR}.t${newCYC}z.pgrb2a.0p50.f$newfhr ${COMIN_SIG}/${newCYC}/pgrb2bp5/ge${GLBPAIR}.t${newCYC}z.pgrb2b.0p50.f$newfhr > tempfile_grib2
       fi

#convert from grib2 to grib1 here
$CNVGRIB -g21 tempfile_grib2 $GFSOUT/${model}_${SREFMEM}.t${newCYC}z.pgrbf$newfhr

fi

