#!/bin/ksh

###############################
# Note:  This job is for bias corrected SREF only (another job called ens_post.sh 
#        is for raw SREF). 
#
# change log:
# 06/08/05, Jun Du: add multi-grids capability for SREF outputs as well as
#                   the precipitable water (PWAT) field; grib product definition
#                   is changed for the lifted index in RSM forecasts due to the use
#                   WRF post. Forecast length is extended from 63hr to 87hr.
# 12/08/05, Jun Du: add WRF_nmm and WRF_em members
# 10/17/07, Xiaoxue Wang: Changed to accomadate sref bias corrected product.
# 11/19/10, Jun Du: Added precipitation and "-nv" option in g12 conversion
# 09/28/11, Jun Du: Added 221 grid (hold back now)
# 06/01/15, SPA SH: comment out alert of grib1 products
###############################

echo "------------------------------------------------"
echo "Ensemble Postprocessing"
echo "------------------------------------------------"

set -x

msg="Starting exsemble.sh for memeber $memberlist"
postmsg "$jlogfile" "$msg"

#####################################
# Define Script/Exec Variables
#####################################

#####################################
# Set up hours to dump                        
#####################################

model=$1
member=$2
FSTART_HIGH=03
FEND_HIGH=87
FINC_HIGH=3
FILE_HIGH1="\$COMOUT/\${RUN}_\$model.\${cycle}.pgrb212.\${member}.f\${HOUR}"
#FILE_HIGH2="\$COMOUT/\${RUN}_\$model.\${cycle}.pgrb216.\${member}.f\${HOUR}"
#FILE_HIGH3="\$COMOUT/\${RUN}_\$model.\${cycle}.pgrb243.\${member}.f\${HOUR}"
#FILE_HIGH4="\$COMOUT/\${RUN}_\$model.\${cycle}.pgrb221.\${member}.f\${HOUR}"


WGRIB_LIST=":HGT:1000.mb|:HGT:850.mb|:HGT:700.mb|:HGT:500.mb|:HGT:250.mb|UGRD:10.m.a|:UGRD:850.mb|:UGRD:500.mb|:UGRD:250.mb|:VGRD:10.m.a|:VGRD:850.mb|:VGRD:500.mb|:VGRD:250.mb|:TMP:2.m.a|:TMP:850.mb|:TMP:700.mb|:RH:850.mb|:RH:700.mb|:PRMSL:MSL|:APCP:sfc"

OUTPUT_FILE1="\$COMOUT/\${RUN}_\${model}.\${cycle}.pgrb212.\${member}"
#OUTPUT_FILE2="\$COMOUT/\${RUN}_\${model}.\${cycle}.pgrb216.\${member}"
#OUTPUT_FILE3="\$COMOUT/\${RUN}_\${model}.\${cycle}.pgrb243.\${member}"
#OUTPUT_FILE4="\$COMOUT/\${RUN}_\${model}.\${cycle}.pgrb221.\${member}"

if [ $model = em -o $model = arw ]; then
  MODEL_TYPE=ARW
else
  MODEL_TYPE=`echo $model| tr [a-z] [A-Z]`
fi
alert_type="${RUN}_${MODEL_TYPE}_PGB"
ALERT_TYPE=`echo $alert_type | tr [a-z] [A-Z]`

#####################################
# Concatenate the grib files for all     
#  hours into a single file.  Begin 
#  by deleting old files, if present.
#####################################
for member in $memberlist
do 
   if [ -s pgrbf${member} ]
   then
       rm  pgrbf${member}
   fi
   if [ -s pgrb1f${member} ]
   then
       rm  pgrb1f${member}
   fi
   if [ -s pgrb2f${member} ]
   then
       rm  pgrb2f${member}
   fi
   if [ -s pgrb3f${member} ]
   then
       rm  pgrb3f${member}
   fi
   if [ -s pgrb4f${member} ]
   then
       rm  pgrb4f${member}
   fi

   typeset -Z2 HOUR

   HOUR=$FSTART_HIGH

   while [ $HOUR -le $FEND_HIGH ]
   do
    filename1=`eval echo $FILE_HIGH1`
    cat $filename1  >>  pgrb1f${member}  

   # filename2=`eval echo $FILE_HIGH2`
   # cat $filename2  >>  pgrb2f${member}  

   # filename3=`eval echo $FILE_HIGH3`
   # cat $filename3  >>  pgrb3f${member}  

   # filename4=`eval echo $FILE_HIGH4`
   # cat $filename4  >>  pgrb4f${member}  

    let HOUR=$HOUR+$FINC_HIGH
   done 

#   HOUR=$FSTART_LOW
#   #
#   # IF FSTART_LOW IS SET THEN...
#   #
#   if [ "$FSTART_LOW" != ""  ]
#   then
#      while [ $HOUR -le $FEND_LOW ]
#      do
#         filename=`eval echo $FILE_LOW`
#         cat $filename  >>  pgrbf${file}    
#
#         let HOUR=$HOUR+$FINC_LOW
#      done
#   fi

#####################################
# Run wgrib to get the inventory of     
#  this new concatenated file.       
#####################################
    
  $WGRIB -s pgrb1f${member} > pgrb1f${member}.inv
#  $WGRIB -s pgrb2f${member} > pgrb2f${member}.inv
#  $WGRIB -s pgrb3f${member} > pgrb3f${member}.inv
#  $WGRIB -s pgrb4f${member} > pgrb4f${member}.inv
#####################################
# Use egrep to filter the fields         
#####################################

  egrep $WGRIB_LIST pgrb1f${member}.inv | $WGRIB pgrb1f${member} -s -grib -i -o ensemble1.${member}
#  egrep $WGRIB_LIST pgrb2f${member}.inv | $WGRIB pgrb2f${member} -s -grib -i -o ensemble2.${member}
#  egrep $WGRIB_LIST pgrb3f${member}.inv | $WGRIB pgrb3f${member} -s -grib -i -o ensemble3.${member}
#  egrep $WGRIB_LIST pgrb4f${member}.inv | $WGRIB pgrb4f${member} -s -grib -i -o ensemble4.${member}

#####################################
# Write the files out to permanent         
#  diskspace
#####################################
  filename1=`eval echo $OUTPUT_FILE1`
  cp ensemble1.${member} $filename1

#  filename2=`eval echo $OUTPUT_FILE2`
#  cp ensemble2.${member} $filename2

#  filename3=`eval echo $OUTPUT_FILE3`
#  cp ensemble3.${member} $filename3

#  filename4=`eval echo $OUTPUT_FILE4`
#  cp ensemble4.${member} $filename4

   #################################
   # Convert to grib2 format
   #################################
   $CNVGRIB -g12 -p40 -nv $filename1 ${filename1}.grib2
#   $CNVGRIB -g12 -p40 -nv $filename2 ${filename2}.grib2
#   $CNVGRIB -g12 -p40 -nv $filename3 ${filename3}.grib2
#   $CNVGRIB -g12 -p40 -nv $filename4 ${filename4}.grib2
   $WGRIB2 ${filename1}.grib2 -s >${filename1}.grib2.idx
#   $WGRIB2 ${filename2}.grib2 -s >${filename2}.grib2.idx
#   $WGRIB2 ${filename3}.grib2 -s >${filename3}.grib2.idx
#   $WGRIB2 ${filename4}.grib2 -s >${filename4}.grib2.idx
   
#   cp ${filename1}.grib2 $COMOUT

   if [ "$SENDDBN" = "YES" ]
   then
#     $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE} $job $filename1
#     $DBNROOT/bin/dbn_alert MODEL $ALERT_TYPE $job $filename2
#     $DBNROOT/bin/dbn_alert MODEL $ALERT_TYPE $job $filename3
#     $DBNROOT/bin/dbn_alert MODEL $ALERT_TYPE $job $filename4

#     if [ $SENDDBN_GB2 = YES ]
#     then

       $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2 $job ${filename1}.grib2
#       $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2 $job ${filename2}.grib2
#       $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2 $job ${filename3}.grib2
#       $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2 $job ${filename4}.grib2

       $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2_WIDX $job ${filename1}.grib2.idx
#       $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2_WIDX $job ${filename2}.grib2.idx
#       $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2_WIDX $job ${filename3}.grib2.idx
#       $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2_WIDX $job ${filename4}.grib2.idx
          
#     fi
  fi
done        

msg="ENDING exsemble.sh for memeber $memberlist"
postmsg "$jlogfile" "$msg"
