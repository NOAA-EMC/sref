#!/bin/ksh
#
######################################################################################
# UNIX Script Documentation Block                                                    #
#                                                                                    #
# Script name:         sref_prep_sig2pgb.sh                                          #
# Script description:  convert global ensemble sigma data into pgrb format for WRF   #
#                      ensemble                                                      #
#                                                                                    #
# Author:        Jun Du       Org: NP22         Date: 2005-07-15                     #
#                                                                                    #
# Abstract: The script gets global ensemble sigma data and flux grib data to convert #
#           them into pgrib data needed for LBCs used by the WRF-NMM and WRF-EM      #
#           ensemble run's.                                                          #
#                                                                                    #
# Script history log:                                                                #
# 2005-07-?  David Michaud  - cut from a global post script                          #
# 2005-07-?  Jun Du - write it into a script using for SREF WRF members              #
# 2006-04-05 Jun Du - modified directory name and file names to match the new global #
#                     ensemble output                                                #
# 2007-08-13 Jun Du - add capability of selecting either 1.0 degree or 0.5 degree to
#                     be converted to
# 2008-05-20 Jun Du - (1) $GLBPAIR's name is modified to include gefs lowres ctl run (gec00)
#                     (2) gens output's vertical levels set to be the same as gfs'
# 2011-07-18 Jun Du - Output GEFS in 0.5Deg due to the increased resolution of GEFS to T254
# 2013-02-27 Jun Du - Added capability to process GFS's EnKF data
# 2015-03-24 Jun Du - Added a capability to search for older GEFS or GFS data
######################################################################################

set -aux

DATE=$1
CYC=$2
GLBPAIR=$3
SREFMEM=$4
export fhr=$5
resolution=$6
model=$7

cd $DATA/${GLBPAIR}_$fhr
#
## GFS 0.5 Degree Files
#
if [ $resolution -eq 5 -a $model = GFS ]; then

 str=`ls -s $COMINgfs/gfs.$DATE/gfs.t${CYC}z.sf90`  #get file size in block (1024 bytes)
     set -A fsize $str
     if [ -s ${COMINgfs}/gfs.$DATE/gfs.t${CYC}z.sf90 ] && [ ${fsize[0]} -gt 800 ] ; then # fcst file is finished?
       echo  ${COMINgfs}/gfs.$DATE/gfs.t${CYC}z.sf90 " exist"
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

   export POSTGPVARS="IO=720,JO=361,KPO=47,PO=1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,725.,700.,675.,650.,625.,600.,575.,550.,525.,500.,475.,450.,425.,400.,375.,350.,325.,300.,275.,250.,225.,200.,175.,150.,125.,100.,70.,50.,30.,20.,10.,7.,5.,3.,2.,1.,POB(51)=1000.,KZZ=8,ZZ=305.,457.,610.,914.,1829.,2743.,3658.,4572.,KPTO=6,"
   export IGEN_ANL=81
   export IGEN_FCST=96
   export IO=720
   export JO=361
fi

if [ $resolution -eq 5 -a $model = gens ]; then

 str=`ls -s $COMINgens/gefs.$DATE/${CYC}/sfcsig/gep20.t${CYC}z.sf90`  #get file size in block (1024 bytes)
     set -A fsize $str
     if [ -s ${COMINgens}/gefs.$DATE/${CYC}/sfcsig/gep20.t${CYC}z.sf90 ] && [ ${fsize[0]} -gt 800 ] ; then # fcst file is finished?
       echo  ${COMINgens}/gefs.$DATE/${CYC}/sfcsig/gep20.t${CYC}z.sf90 " exist"
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

   export POSTGPVARS="IO=720,JO=361,KPO=47,PO=1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,725.,700.,675.,650.,625.,600.,575.,550.,525.,500.,475.,450.,425.,400.,375.,350.,325.,300.,275.,250.,225.,200.,175.,150.,125.,100.,70.,50.,30.,20.,10.,7.,5.,3.,2.,1.,POB(51)=1000.,KZZ=8,ZZ=305.,457.,610.,914.,1829.,2743.,3658.,4572.,KPTO=6,"
   export IGEN_ANL=81
   export IGEN_FCST=96
   export IO=720
   export JO=361
fi

if [ $resolution -eq 5 -a $model = EnKF ]; then
   export POSTGPVARS="IO=720,JO=361,KPO=47,PO=1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,725.,700.,675.,650.,625.,600.,575.,550.,525.,500.,475.,450.,425.,400.,375.,350.,325.,300.,275.,250.,225.,200.,175.,150.,125.,100.,70.,50.,30.,20.,10.,7.,5.,3.,2.,1.,POB(51)=1000.,KZZ=8,ZZ=305.,457.,610.,914.,1829.,2743.,3658.,4572.,KPTO=6,"
   export IGEN_ANL=81
   export IGEN_FCST=96
   export COMIN_SIG=$COMINgfs/enkf.$DATE
   export IO=720
   export JO=361
fi

#
## GEFS 1.0 Degree Files
#
if [ $resolution -eq 1 -a $model = gens ]; then
   export POSTGPVARS="IO=360,JO=181,KPO=47,PO=1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,725.,700.,675.,650.,625.,600.,575.,550.,525.,500.,475.,450.,425.,400.,375.,350.,325.,300.,275.,250.,225.,200.,175.,150.,125.,100.,70.,50.,30.,20.,10.,7.,5.,3.,2.,1.,POB(51)=1000.,KZZ=8,ZZ=305.,457.,610.,914.,1829.,2743.,3658.,4572.,KPTO=6,"
   export IGEN_ANL=81
   export IGEN_FCST=96
   export COMIN_SIG=$COMINgens/gefs.$DATE
   export IO=360
   export JO=181
fi

if [ $model = GFS ]; then
   until test -s ${COMIN_SIG}/gfs.t${newCYC}z.sf${newfhr}
   do
      sleep 20 
   done
   until test -s ${COMIN_SIG}/gfs.t${newCYC}z.sfluxgrbf${newfhr}
   do
      sleep 20 
   done
   cp ${COMIN_SIG}/gfs.t${newCYC}z.sf$newfhr sigfile
   cp ${COMIN_SIG}/gfs.t${newCYC}z.sfluxgrbf$newfhr flxfile
fi

if [ $model = gens ]; then
   until test -s ${COMIN_SIG}/${newCYC}/sfcsig/ge${GLBPAIR}.t${newCYC}z.sf${newfhr}
   do
      sleep 20 
   done
   until test -s ${COMIN_SIG}/${newCYC}/sflux/ge${GLBPAIR}.t${newCYC}z.sfluxgrbf${newfhr}
   do
      sleep 20 
   done
   cp ${COMIN_SIG}/${newCYC}/sfcsig/ge${GLBPAIR}.t${newCYC}z.sf$newfhr sigfile
   cp ${COMIN_SIG}/${newCYC}/sflux/ge${GLBPAIR}.t${newCYC}z.sfluxgrbf$newfhr flxfile
fi

if [ $model = EnKF ]; then
#  until test -s ${COMIN_SIG}/${CYC}/siganl_${DATE}${CYC}_mem080
   until test -s ${COMIN_SIG}/${CYC}/sfg_${DATE}${CYC}_fhr03_mem080
   do
      sleep 20 
   done
#  until test -s $COMINgfs/gfs.$DATE/gfs.t${CYC}z.sfluxgrbf00
   until test -s $COMINgfs/gfs.$DATE/gfs.t${CYC}z.sfluxgrbf03
   do
      sleep 20 
   done
#  cp ${COMIN_SIG}/${CYC}/siganl_${DATE}${CYC}_mem0$fhr sigfile
#  cp $COMINgfs/gfs.$DATE/gfs.t${CYC}z.sfluxgrbf00 flxfile
   cp ${COMIN_SIG}/${CYC}/sfg_${DATE}${CYC}_fhr03_mem0$fhr sigfile
   cp $COMINgfs/gfs.$DATE/gfs.t${CYC}z.sfluxgrbf03 flxfile
fi

export SIGINP=sigfile
export FLXINP=flxfile
#export FLXINP=/dev/null
export FLXIOUT=/dev/null
export PGBOUT=pgbfile
export PGIOUT=pgifile
export IGEN=${IGEN_FCST}
export ANOMCATSH=${USHglobal}/global_anomcat.sh
export POSTGPSH=${POSTGPSH:-${USHglobal}/global_postgp.sh}

$POSTGPSH

if [ $model = GFS ]; then
  mv pgbfile $GFSOUT/gfs.t${CYC}z.master.grbf$fhr
fi
if [ $model = gens ]; then
  mv pgbfile $GFSOUT/${model}_${SREFMEM}.t${CYC}z.pgrbf$fhr
fi
if [ $model = EnKF ]; then
# mv pgbfile $GFSOUT/${model}_mem0${fhr}.t${CYC}z.pgrbf00
  mv pgbfile $GFSOUT/${model}_mem0${fhr}.t${CYC}z.pgrbf03
fi

