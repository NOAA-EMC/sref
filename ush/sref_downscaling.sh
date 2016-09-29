#!/bin/ksh
################################################################
# UNIX Script Documentation Block
#
# Script name:         sref_downscaling.sh
# Script description:  This program is to downscale a low-res forecast 
#                      into a hi-res forecast with a pre-calculated 
#                      dowbscaling vector (derived from exsref_dwnsvect.sh.sms)
#
# Author:      Jun Du       Org: NP21         Date: 2011-12
#
# Change log:
#  08/xx/2007: Bo Cui, original code works for GEFS
#  11/28/2011: Brian Etherton/DET, tested/evaluated the GEFS codes for SREF
#  12/09/2011: Jun Du, rewrote the scripts, excluded sfc PRES and added
#                      2m q, output raw SREF in NDGD grid, and tested and 
#                      implemented for the operational SREF
#
################################################################
set -x

#. ./prep_step

hourlist=$FHRLIST
modellist=$MDLLIST
memberlist=$MEMLIST
#hourlist=$1
#modellist=$2
#memberlist=$3

for nfhrs in $hourlist
do
 for model in $modellist
 do
   for member in $memberlist
   do

   echo $nfhrs $model $member $cyc

#
#  Interpolate raw forecasts into 5km ndgd grid 
#
    grid='255 3 1073 689 20191 238445 8 265000 5079 5079 0 64 25000 25000 0 0 0'
#   rawsref=sref_${model}.t${cyc}z.pgrb212.${member}.f${nfhrs} 
#   cp $COMIN/$rawsref $rawsref
#   infile=$rawsref
    infile=$COMIN/sref_${model}.t${cyc}z.pgrb212.${member}.f${nfhrs} 
    outfile=sref_${model}.t${cyc}z.pgrb212.${member}.f${nfhrs}_temp
    cfile=sref_raw_${model}.t${cyc}z.pgrb197.${member}.f${nfhrs}

    >$outfile
#   ${WGRIB} $infile | grep ":PRES:" | grep "sfc"  | $WGRIB -i $infile  -grib -append -o $outfile
    ${WGRIB} $infile | grep ":TMP:"  | grep "2 m"  | $WGRIB -i $infile  -grib -append -o $outfile
    ${WGRIB} $infile | grep ":SPFH:" | grep "2 m"  | $WGRIB -i $infile  -grib -append -o $outfile
    ${WGRIB} $infile | grep ":UGRD:" | grep "10 m" | $WGRIB -i $infile  -grib -append -o $outfile
    ${WGRIB} $infile | grep ":VGRD:" | grep "10 m" | $WGRIB -i $infile  -grib -append -o $outfile
#   ${WGRIB} $infile | grep ":TMP:kpds5=11:kpds6=105:kpds7=2:"   | $WGRIB -i -grib $infile -o $outfile
#   ${WGRIB} $infile | grep ":SPFH:kpds5=51:kpds6=105:kpds7=2:"  | $WGRIB -i -grib $infile -o $outfile
#   ${WGRIB} $infile | grep ":UGRD:kpds5=33:kpds6=105:kpds7=10:" | $WGRIB -i -grib $infile -o $outfile
#   ${WGRIB} $infile | grep ":VGRD:kpds5=34:kpds6=105:kpds7=10:" | $WGRIB -i -grib $infile -o $outfile
    ls -l $outfile
    ${COPYGB} -g"$grid" -i1,1 -x $outfile $cfile
#   $EXECsref/fastcopygb -g"$grid" -i1,1 -x $outfile $cfile
    ls -l $file
    mv $cfile ${COM_NDGD}/.
    rm -f sref*temp
#   mv sref*temp  ${COM_NDGD}/.

#
#  Interpolate bias corrected forecasts into 5km ndgd grid 
#
    grid='255 3 1073 689 20191 238445 8 265000 5079 5079 0 64 25000 25000 0 0 0'
    infile=$COMOUT/sref_${model}.t${cyc}z.pgrb212.${member}.f${nfhrs} 
    outfile=sref_${model}.t${cyc}z.pgrb212.${member}.f${nfhrs}_temp
    cfile=sref_bc_${model}.t${cyc}z.pgrb197.${member}.f${nfhrs}

    >$outfile
#   ${WGRIB} $infile | grep ":PRES:" | grep "sfc"  | $WGRIB -i $infile  -grib -append -o $outfile
    ${WGRIB} $infile | grep ":TMP:"  | grep "2 m"  | $WGRIB -i $infile  -grib -append -o $outfile
    ${WGRIB} $infile | grep ":SPFH:" | grep "2 m"  | $WGRIB -i $infile  -grib -append -o $outfile
    ${WGRIB} $infile | grep ":UGRD:" | grep "10 m" | $WGRIB -i $infile  -grib -append -o $outfile
    ${WGRIB} $infile | grep ":VGRD:" | grep "10 m" | $WGRIB -i $infile  -grib -append -o $outfile
    ${COPYGB} -g"$grid" -i1,1 -x $outfile $cfile
#   $EXECsref/fastcopygb -g"$grid" -i1,1 -x $outfile $cfile
    cp $cfile ${COM_NDGD}/.
    rm -f sref*temp

#
# Set the index to default 0 (downscaling vector exists)
#
    cstart=0

# Given a fcst hour to figure out a corresponding anl cycle or dwnsvect cyc
    if [ $cyc -eq 03 ]; then
     if [ $nfhrs -eq 00 -o $nfhrs -eq 24 -o $nfhrs -eq 48 -o $nfhrs -eq 72 ];then
      anlcyc=03
     fi
     if [ $nfhrs -eq 03 -o $nfhrs -eq 27 -o $nfhrs -eq 51 -o $nfhrs -eq 75 ];then
      anlcyc=06
     fi
     if [ $nfhrs -eq 06 -o $nfhrs -eq 30 -o $nfhrs -eq 54 -o $nfhrs -eq 78 ];then
      anlcyc=09
     fi
     if [ $nfhrs -eq 09 -o $nfhrs -eq 33 -o $nfhrs -eq 57 -o $nfhrs -eq 81 ];then
      anlcyc=12
     fi
     if [ $nfhrs -eq 12 -o $nfhrs -eq 36 -o $nfhrs -eq 60 -o $nfhrs -eq 84 ];then
      anlcyc=15
     fi
     if [ $nfhrs -eq 15 -o $nfhrs -eq 39 -o $nfhrs -eq 63 -o $nfhrs -eq 87 ];then
      anlcyc=18
     fi
     if [ $nfhrs -eq 18 -o $nfhrs -eq 42 -o $nfhrs -eq 66 ];then
      anlcyc=21
     fi
     if [ $nfhrs -eq 21 -o $nfhrs -eq 45 -o $nfhrs -eq 69 ];then
      anlcyc=00
     fi
    fi

    if [ $cyc -eq 09 ]; then
     if [ $nfhrs -eq 00 -o $nfhrs -eq 24 -o $nfhrs -eq 48 -o $nfhrs -eq 72 ];then
      anlcyc=09
     fi
     if [ $nfhrs -eq 03 -o $nfhrs -eq 27 -o $nfhrs -eq 51 -o $nfhrs -eq 75 ];then
      anlcyc=12
     fi
     if [ $nfhrs -eq 06 -o $nfhrs -eq 30 -o $nfhrs -eq 54 -o $nfhrs -eq 78 ];then
      anlcyc=15
     fi
     if [ $nfhrs -eq 09 -o $nfhrs -eq 33 -o $nfhrs -eq 57 -o $nfhrs -eq 81 ];then
      anlcyc=18
     fi
     if [ $nfhrs -eq 12 -o $nfhrs -eq 36 -o $nfhrs -eq 60 -o $nfhrs -eq 84 ];then
      anlcyc=21
     fi
     if [ $nfhrs -eq 15 -o $nfhrs -eq 39 -o $nfhrs -eq 63 -o $nfhrs -eq 87 ];then
      anlcyc=00
     fi
     if [ $nfhrs -eq 18 -o $nfhrs -eq 42 -o $nfhrs -eq 66 ];then
      anlcyc=03
     fi
     if [ $nfhrs -eq 21 -o $nfhrs -eq 45 -o $nfhrs -eq 69 ];then
      anlcyc=06
     fi
    fi

    if [ $cyc -eq 15 ]; then
     if [ $nfhrs -eq 00 -o $nfhrs -eq 24 -o $nfhrs -eq 48 -o $nfhrs -eq 72 ];then
      anlcyc=15
     fi
     if [ $nfhrs -eq 03 -o $nfhrs -eq 27 -o $nfhrs -eq 51 -o $nfhrs -eq 75 ];then
      anlcyc=18
     fi
     if [ $nfhrs -eq 06 -o $nfhrs -eq 30 -o $nfhrs -eq 54 -o $nfhrs -eq 78 ];then
      anlcyc=21
     fi
     if [ $nfhrs -eq 09 -o $nfhrs -eq 33 -o $nfhrs -eq 57 -o $nfhrs -eq 81 ];then
      anlcyc=00
     fi
     if [ $nfhrs -eq 12 -o $nfhrs -eq 36 -o $nfhrs -eq 60 -o $nfhrs -eq 84 ];then
      anlcyc=03
     fi
     if [ $nfhrs -eq 15 -o $nfhrs -eq 39 -o $nfhrs -eq 63 -o $nfhrs -eq 87 ];then
      anlcyc=06
     fi
     if [ $nfhrs -eq 18 -o $nfhrs -eq 42 -o $nfhrs -eq 66 ];then
      anlcyc=09
     fi
     if [ $nfhrs -eq 21 -o $nfhrs -eq 45 -o $nfhrs -eq 69 ];then
      anlcyc=12
     fi
    fi

    if [ $cyc -eq 21 ]; then
     if [ $nfhrs -eq 00 -o $nfhrs -eq 24 -o $nfhrs -eq 48 -o $nfhrs -eq 72 ];then
      anlcyc=21
     fi
     if [ $nfhrs -eq 03 -o $nfhrs -eq 27 -o $nfhrs -eq 51 -o $nfhrs -eq 75 ];then
      anlcyc=00
     fi
     if [ $nfhrs -eq 06 -o $nfhrs -eq 30 -o $nfhrs -eq 54 -o $nfhrs -eq 78 ];then
      anlcyc=03
     fi
     if [ $nfhrs -eq 09 -o $nfhrs -eq 33 -o $nfhrs -eq 57 -o $nfhrs -eq 81 ];then
      anlcyc=06
     fi
     if [ $nfhrs -eq 12 -o $nfhrs -eq 36 -o $nfhrs -eq 60 -o $nfhrs -eq 84 ];then
      anlcyc=09
     fi
     if [ $nfhrs -eq 15 -o $nfhrs -eq 39 -o $nfhrs -eq 63 -o $nfhrs -eq 87 ];then
      anlcyc=12
     fi
     if [ $nfhrs -eq 18 -o $nfhrs -eq 42 -o $nfhrs -eq 66 ];then
      anlcyc=15
     fi
     if [ $nfhrs -eq 21 -o $nfhrs -eq 45 -o $nfhrs -eq 69 ];then
      anlcyc=18
     fi
    fi

    ifile=dvrtma.t${anlcyc}z.ndgd_conus 
    if [ -s $COM_IN/sref.$PDY/${cyc}/misc/$ifile ]; then
      cp $COM_IN/sref.$PDY/${cyc}/misc/$ifile $ifile
    elif [ -s $COM_IN/sref.$PDYm1/${cyc}/misc/$ifile ]; then
      cp $COM_IN/sref.$PDYm1/${cyc}/misc/$ifile $ifile
    elif [ -s $COM_IN/sref.$PDYm2/${cyc}/misc/$ifile ]; then
      cp $COM_IN/sref.$PDYm2/${cyc}/misc/$ifile $ifile
    elif [ -s $COM_IN/sref.$PDYm3/${cyc}/misc/$ifile ]; then
      cp $COM_IN/sref.$PDYm3/${cyc}/misc/$ifile $ifile
    elif [ -s $COM_IN/sref.$PDYm4/${cyc}/misc/$ifile ]; then
      cp $COM_IN/sref.$PDYm4/${cyc}/misc/$ifile $ifile
    elif [ -s $COM_IN/sref.$PDYm5/${cyc}/misc/$ifile ]; then
      cp $COM_IN/sref.$PDYm5/${cyc}/misc/$ifile $ifile
    elif [ -s $COM_IN/sref.$PDYm6/${cyc}/misc/$ifile ]; then
      cp $COM_IN/sref.$PDYm6/${cyc}/misc/$ifile $ifile
    elif [ -s $COM_IN/sref.$PDYm7/${cyc}/misc/$ifile ]; then
      cp $COM_IN/sref.$PDYm7/${cyc}/misc/$ifile $ifile
    else
      echo " Note: There is no downscaling vector available at " t${anlcyc}z 
      cstart=1
    fi

    echo "&message"  >input.$nfhrs.$model.$member
    echo " icstart=${cstart}," >> input.$nfhrs.$model.$member
    echo "/" >>input.$nfhrs.$model.$member

    ofile=sref_ds_${model}.t${cyc}z.pgrb197.${member}.f${nfhrs}

#
#   Downscaling a forecast
#
      rm fort.*
      ln -sf $ifile fort.11
      ln -sf $cfile fort.12
      ln -sf $ofile fort.51

#     startmsg
      $EXECsref/sref_dvrtma_debias   <input.$nfhrs.$model.$member  > $pgmout.$nfhrs.$model.$member  2> errfile
      cp $ofile ${COM_NDGD}/.
#     export err=$?;$DATA/err_chk

  done
 done
done

set +x
echo " "
echo "Leaving ush script sref_downscaling.sh"
echo " "
set -x
