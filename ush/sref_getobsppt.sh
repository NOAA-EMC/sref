#!/bin/ksh -x
#
###################################################################
# This script gets 3-hourly Stage II precipitation analysis 
# i.e  00/03/06/09/12/15/18/21Z for a day (21z-21z)
# 
# Change Log:
# 10/01/10: Jun Du, initial script for precip bias correction job
# 04/22/11: Jun Du, added an alternative option of using CCPA data
#
###################################################################
set -x
export XLFRTEOPTS="unit_vars=yes"

option=${option:-CCPA3}

vday=`echo $1 | cut -c 1-8`
vdaym1=`$NDATE -24 ${vday}12 | cut -c 1-8`
vdayp1=`$NDATE +24 ${vday}12 | cut -c 1-8`

YYYY=`echo $vday | cut -c 1-4`
MM=`echo $vday   | cut -c 5-6`

#export CPGB=${COPYGB:-${NWROOTp1}/util/exec/copygb}
export CPGB=${CPGB:-${COPYGB:?}}
#export EXECnam=${EXECnam:-${NWROOTp1}/exec}
export EXECnam=${EXECnam:-${NWROOTp3}/nam.${nam.ver}/exec/dell.exec}

if [ $option = ST2 ]; then
 get Stage2 hourly data
#export COMINccpa=${COMROOTp1}/hourly/prod/nam_pcpn_anal.
 export COMINccpa=${COMROOT}/pcpanl/prod/pcpanl.
else
 get 3hrly CCPA data
 export COMINccpa=$(compath.py ccpa/prod)
fi

# Get Stage II data:
if [ $option = ST2 ]; then
# get Stage2 hourly accumulated precip data from /com directory 
# covering 21Z $vdaym1 - 21Z $vday
cp ${COMINccpa}${vdaym1}/ST2ml${vdaym1}22.Grb.Z .
cp ${COMINccpa}${vdaym1}/ST2ml${vdaym1}23.Grb.Z .
cp ${COMINccpa}${vday}/ST2ml${vday}0*.Grb.Z .
cp ${COMINccpa}${vday}/ST2ml${vday}1*.Grb.Z .
cp ${COMINccpa}${vday}/ST2ml${vday}20.Grb.Z .
cp ${COMINccpa}${vday}/ST2ml${vday}21.Grb.Z .

uncompress ST2ml*.Grb.Z

# process hourly data into 3hrly accumulated precip
hr=0
accdate=${vday}00
while [ $accdate -le ${vday}21 ]; do
  if [ $hr -lt 10 ];then hr=0$hr;fi
  time1=`$NDATE -2 $accdate`  
  time2=`$NDATE -1 $accdate`  
  time3=$accdate
  cat > input_acc3h << EOF
obs
ST2.
ST2ml$time1.Grb
ST2ml$time2.Grb
ST2ml$time3.Grb
EOF

$EXECnam/nam_stage4_acc < input_acc3h
export err=$?; err_chk

# convert to grid 212
#  for grid in 212 218 
 for grid in 212 
 do
   $CPGB -g${grid} -i3 -x ST2.$accdate.03h ppt$hr
 done

 accdate=`$NDATE +3 $accdate`  
 hr=`expr $hr + 3`
done 
fi    # finished accumulating 3-hourly Stage II data.

# If use CCPA data:
grid=212
if [ $option = CCPA3 ]; then
 hr=0
 while [ $hr -le 21 ]; do
  if [ $hr -lt 10 ];then hr=0$hr;fi
# $CPGB -g${grid} -i3 -x ${COMINccpa}/gefs.$vday/$hr/ccpa.t${hr}z.03h.hrap.conus ppt$hr

  if [ $hr -eq 00 ];then
   $CPGB -g${grid} -i3 -x ${COMINccpa}/ccpa.$vday/00/ccpa.t${hr}z.03h.hrap.conus ppt$hr
  fi
  if [ $hr -eq 03 -o $hr -eq 06 ];then
   $CPGB -g${grid} -i3 -x ${COMINccpa}/ccpa.$vday/06/ccpa.t${hr}z.03h.hrap.conus ppt$hr
  fi
  if [ $hr -eq 09 -o $hr -eq 12 ];then
   $CPGB -g${grid} -i3 -x ${COMINccpa}/ccpa.$vday/12/ccpa.t${hr}z.03h.hrap.conus ppt$hr
  fi
  if [ $hr -eq 15 -o $hr -eq 18 ];then
   $CPGB -g${grid} -i3 -x ${COMINccpa}/ccpa.$vday/18/ccpa.t${hr}z.03h.hrap.conus ppt$hr
  fi
  if [ $hr -eq 21 ];then
   $CPGB -g${grid} -i3 -x ${COMINccpa}/ccpa.$vdayp1/00/ccpa.t${hr}z.03h.hrap.conus ppt$hr
  fi

  hr=`expr $hr + 3`
 done
fi #finished getting CCPA data
 
ls -l ppt*
rm -f ST* input*

exit

