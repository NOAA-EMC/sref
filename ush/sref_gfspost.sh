#!/bin/ksh
#
######################################################################################
# UNIX Script Documentation Block                                                    #
#                                                                                    #
# Script name:         sref_gfspost.sh                                               #
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
# 2008-05-20 Jun Du - (1) $GLBMEM's name is modified to include gefs lowres ctl run (gec00)
#                     (2) gens output's vertical levels set to be the same as gfs'
# 2011-07-18 Jun Du - Output GEFS in 0.5Deg due to the increased resolution of GEFS to T254
######################################################################################

set -eux

DATE=$1
CYC=$2
GLBMEM=$3
FHR=$4
GLB_MODEL=$5

mkdir -p $DATA/${GLB_MODEL}_${GLBMEM}_$FHR
cd $DATA/${GLB_MODEL}_${GLBMEM}_$FHR
#
## GFS 0.5 Degree Files
#
export POSTGPVARS="IO=720,JO=361,KPO=47,PO=1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,725.,700.,675.,650.,625.,600.,575.,550.,525.,500.,475.,450.,425.,400.,375.,350.,325.,300.,275.,250.,225.,200.,175.,150.,125.,100.,70.,50.,30.,20.,10.,7.,5.,3.,2.,1.,POB(51)=1000.,KZZ=8,ZZ=305.,457.,610.,914.,1829.,2743.,3658.,4572.,KPTO=6,"
export IGEN_ANL=81
export IGEN_FCST=96
export IO=720
export JO=361

if [ $GLB_MODEL = gfs ]; then
   COMIN_SIG=$COMINgfs/gfs.$DATE
   until test -s ${COMIN_SIG}/gfs.t${CYC}z.sf${FHR}
   do
      sleep 20 
   done
   until test -s ${COMIN_SIG}/gfs.t${CYC}z.sfluxgrbf${FHR}
   do
      sleep 20 
   done
   ln -sf ${COMIN_SIG}/gfs.t${CYC}z.sf$FHR sigfile
   ln -sf ${COMIN_SIG}/gfs.t${CYC}z.sfluxgrbf$FHR flxfile

elif [ $GLB_MODEL = gens ]; then
   COMIN_SIG=$COMINgens/gefs.$DATE
   until test -s ${COMIN_SIG}/${CYC}/sfcsig/ge${GLBMEM}.t${CYC}z.sf${FHR}
   do
      sleep 20 
   done
   until test -s ${COMIN_SIG}/${CYC}/sflux/ge${GLBMEM}.t${CYC}z.sfluxgrbf${FHR}
   do
      sleep 20 
   done
   ln -sf ${COMIN_SIG}/${CYC}/sfcsig/ge${GLBMEM}.t${CYC}z.sf$FHR sigfile
   ln -sf ${COMIN_SIG}/${CYC}/sflux/ge${GLBMEM}.t${CYC}z.sfluxgrbf$FHR flxfile

else
   echo "Unknown GLB_MODEL: $GLB_MODEL"
   exit 1
fi

export SIGINP=sigfile
export FLXINP=flxfile
export FLXIOUT=/dev/null
export PGBOUT=pgbfile
export PGIOUT=pgifile
export IGEN=${IGEN_FCST}
export ANOMCATSH=
export POSTGPSH=${POSTGPSH:-${USHglobal}/global_postgp.sh}
export XC=
export MP=s

$POSTGPSH

if [ $GLB_MODEL = gfs ]; then
  mv pgbfile $GFSOUT/gfs.t${CYC}z.master.grbf$FHR
fi
if [ $GLB_MODEL = gens ]; then
  mv pgbfile $GFSOUT/${GLB_MODEL}_${GLBMEM}.t${CYC}z.pgrbf$FHR
fi
