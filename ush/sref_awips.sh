#!/bin/sh
######################################################################
#  UTILITY SCRIPT NAME :  sref_awips.sh
#         DATE WRITTEN :  06/07/2006
#
#  Abstract:  This utility script produces the  (mean, spread and
#             probability) of selected variables for grids 212, 216
#             and 243 in GRIB2 format for AWIPS from the SREF model.
#
#
#####################################################################
echo "-----------------------------------------------------"
echo "JSREF_AWIPS ( 03Z,09Z,15Z AND 21Z) SREF postprocessing"
echo "-----------------------------------------------------"
echo "History: JUNE 2006 - First implementation of this new script."
echo "         AUG  2013 - Add process to produce MMEFS for AWIPS  "
echo "                     from the SREF model  "
echo "03/31/2015 - Jun Du, Expanded for SREF.v7.0.0"
echo " "
#####################################################################

########################################
  msg="HAS BEGUN!"
  postmsg "$jlogfile" "$msg"
########################################

set -x 

echo " ------------------------------------------------------------"
echo "  BEGIN MAKING SREF (for grids 212, 216, 243) GRIB2 PRODUCTS "
echo " ------------------------------------------------------------"

set +x
echo " "
echo "###############################################################"
echo "                                                               "
echo "   Process SREF GRIB2 (Mean, Spread and Probability) PRODUCTS  "
echo "                                                               "
echo "###############################################################"
echo " "
set -x

for grid in 212 216 243; do
  for type in mean spread prob ; do
    export GRID=$grid
    export type
    mkawpgrb.sh "00"
  done
done

echo " ------------------------------------------"
echo " BEGIN MAKING  SREF (MMEFS)  GRIB2 PRODUCTS"
echo " ------------------------------------------"

set +x
echo " "
echo "###########################################"
echo "      Process SREF GRIB2 (MMEFS) PRODUCTS  "
echo " FORECAST HOURS 00 - 87.                   "
echo "###########################################"
echo " "
set -x

#for type in em nmb nmm
for type in arw nmb
do
   for id in ctl n1 n2 n3 n4 n5 n6 p1 p2 p3 p4 p5 p6
   do
      for hr in 00 03 06 09 12 15 18 21 24 27 30 33 36 39 42 45 48 \
                51 54 57 60 63 66 69 72 75 78 81 84 87
      do
         cp ${COMROOT}/${NET}/${envir}/${RUN}.${PDY}/${cyc}/pgrb/sref_${type}.t${cyc}z.pgrb221.${id}.f${hr}.grib2  .
         ${WGRIB2} sref_${type}.t${cyc}z.pgrb221.${id}.f${hr}.grib2 | grep "APCP"  |  \
         ${WGRIB2} -i sref_${type}.t${cyc}z.pgrb221.${id}.f${hr}.grib2 -grib sref_${type}.t${cyc}z.pgrb221.${id}.f${hr}.grib2_apcp
         ${WGRIB2}  sref_${type}.t${cyc}z.pgrb221.${id}.f${hr}.grib2 | grep "TMP:2 m" | \
         ${WGRIB2} -i sref_${type}.t${cyc}z.pgrb221.${id}.f${hr}.grib2 -grib sref_${type}.t${cyc}z.pgrb221.${id}.f${hr}.grib2_tmp_2m
         cat sref_${type}.t${cyc}z.pgrb221.${id}.f${hr}.grib2_apcp sref_${type}.t${cyc}z.pgrb221.${id}.f${hr}.grib2_tmp_2m >> sref_${type}.t${cyc}z.pgrb221.${id}_all
     done
   done
done

############################################
# Processing GRIB2 SREF grid 221 for MMEFS
############################################

  pgm=tocgrib2
  export pgm;. prep_step
  startmsg

#for type in em nmm nmb
for type in arw nmb
do
  for id in ctl n1 n2 n3 n4 n5 n6 p1 p2 p3 p4 p5 p6
  do

   export FORT11=sref_${type}.t${cyc}z.pgrb221.${id}_all
   export FORT31=" "
   export FORT51=grib2_sref_${type}.t${cyc}z.${id}_pgrb221
   ${TOCGRIB2} <  $PARMsref/grib2_sref_${type}_pgrb221

   err=$?;export err ;err_chk
   echo " error from tocgrib2=",$err

  if [ $SENDCOM = "YES" ] ; then

   ##############################
   # Post Files to PCOM
   ##############################

     cp  grib2_sref_${type}.t${cyc}z.${id}_pgrb221  $pcom/grib2_sref_${type}.t${cyc}z.${id}_pgrb221

  fi

  ##########################
  # Distribute Data to NCF
  #########################

  if [ $SENDDBN = "YES" ] ; then
#
#    Distribute Data to TOC (AWIPS)
#
     $DBNROOT/bin/dbn_alert NTC_LOW $NET $job   $pcom/grib2_sref_${type}.t${cyc}z.${id}_pgrb221
     echo " "
  fi

  done
done

################################################################################
# GOOD RUN
set +x
echo "**************JOB JSREF_AWIPS  COMPLETED NORMALLY ON THE WCOSS"
set -x
################################################################################

msg="HAS COMPLETED NORMALLY!"
postmsg "$jlogfile" "$msg"

############## END OF SCRIPT #######################
