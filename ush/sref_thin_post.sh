#!/bin/ksh

###############################
# Note:  This job is for raw SREF. There is another job called sref_ens_post.sh is 
#        particularly for bias corrected SREF for the same purpose.
#
# change log:
# 06/08/05, Jun Du: add multi-grids capability for SREF outputs as well as
#                   the precipitable water (PWAT) field; grib product definition
#                   is changed for the lifted index in RSM forecasts due to the use
#                   WRF post. Forecast length is extended from 63hr to 87hr.
# 12/08/05, Jun Du: add WRF_nmm and WRF_em members
# 09/28/11, Jun Du: added model NMMB, 221 grid and "-nv" flag in grib conversion 
# 05/03/12, Jun Du: added 16km NA grid 132 and changed the old job script name
#                   from ens_post.sh to sref_thin_post.sh
# 03/24/15, Jun Du: Added PRATE and visbility as requested by WPC
# 05/18/15, SPA SH: Comment out the dbn_alert of grib1 products

echo "------------------------------------------------"
echo "Ensemble Postprocessing"
echo "------------------------------------------------"

set -eux

msg="Starting sref_thin_post.sh for memeber $memberlist"
postmsg "$jlogfile" "$msg"

#####################################
# Define Script/Exec Variables
#####################################

# JY WGRIB=$WGRIB

#####################################
# Set up hours to dump                        
#####################################

case $RUN in
   sref) FSTART_HIGH=00
         FEND_HIGH=87
         FINC_HIGH=3
       if [ $OUTGRD -eq 132 ]; then
         FILE_HIGH1="\$COMOUT/\${RUN}_\$model.\${cycle}.pgrb132.\${member}.f\${HOUR}"
       else
         FILE_HIGH1="\$COMOUT/\${RUN}_\$model.\${cycle}.pgrb212.\${member}.f\${HOUR}"
         FILE_HIGH2="\$COMOUT/\${RUN}_\$model.\${cycle}.pgrb216.\${member}.f\${HOUR}"
         FILE_HIGH3="\$COMOUT/\${RUN}_\$model.\${cycle}.pgrb243.\${member}.f\${HOUR}"
         FILE_HIGH4="\$COMOUT/\${RUN}_\$model.\${cycle}.pgrb221.\${member}.f\${HOUR}"
       fi

         if [ "$model" = "eta" -o "$model" = "kfeta" ]
         then
            WGRIB_LIST=":PWAT:atmos.col|:HGT:1000.mb|:HGT:850.mb|:HGT:700.mb|:HGT:500.mb|:HGT:250.mb|UGRD:10.m.a|:UGRD:850.mb|:UGRD:700.mb|:UGRD:500.mb|:UGRD:250.mb|:VGRD:10.m.a|:VGRD:850.mb|:VGRD:700.mb|:VGRD:500.mb|:VGRD:250.mb|:TMP:2.m.a|:TMP:850.mb|:TMP:700.mb|:RH:850.mb|:RH:700.mb|:APCP|:PRMSL:MSL|:ABSV:500.mb|:ABSV:250.mb|:CAPE:sfc|:CIN:sfc|:PLI.30.mb.a|:TMAX:2.m.a|:TMIN:2.m.a|:CSNOW:sfc|:CICEP:sfc|:CFRZR:sfc|:CRAIN:sfc|:VIS|:PRATE"
         fi

         if [ "$model" = "rsm" -o "$model" = "nmm" -o "$model" = "em" -o "$model" = "arw" -o "$model" = "nmb" ]
         then
            WGRIB_LIST=":PWAT:atmos.col|:HGT:1000.mb|:HGT:850.mb|:HGT:700.mb|:HGT:500.mb|:HGT:250.mb|UGRD:10.m.a|:UGRD:850.mb|:UGRD:700.mb|:UGRD:500.mb|:UGRD:250.mb|:VGRD:10.m.a|:VGRD:850.mb|:VGRD:700.mb|:VGRD:500.mb|:VGRD:250.mb|:TMP:2.m.a|:TMP:850.mb|:TMP:700.mb|:RH:850.mb|:RH:700.mb|:APCP|:PRMSL:MSL|:ABSV:500.mb|:ABSV:250.mb|:CAPE:sfc|:CIN:sfc|:LFTX:500.1000.mb|:TMAX:2.m.a|:TMIN:2.m.a|:CSNOW:sfc|:CICEP:sfc|:CFRZR:sfc|:CRAIN:sfc|:VIS|:PRATE"
         fi

       if [ $OUTGRD -eq 132 ]; then
         OUTPUT_FILE1="\$COMOUT/\${RUN}_\${model}.\${cycle}.pgrb132.\${member}"
       else
         OUTPUT_FILE1="\$COMOUT/\${RUN}_\${model}.\${cycle}.pgrb212.\${member}"
         OUTPUT_FILE2="\$COMOUT/\${RUN}_\${model}.\${cycle}.pgrb216.\${member}"
         OUTPUT_FILE3="\$COMOUT/\${RUN}_\${model}.\${cycle}.pgrb243.\${member}"
         OUTPUT_FILE4="\$COMOUT/\${RUN}_\${model}.\${cycle}.pgrb221.\${member}"
       fi
         alert_type="${RUN}_${MODEL}_PGB"
         ALERT_TYPE=`echo $alert_type | tr [a-z] [A-Z]`
         ;;

esac

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
   if [ -s ensemble.${member} ]
   then
       rm  ensemble.${member}
   fi
   if [ -s ensemble1.${member} ]
   then
       rm  ensemble1.${member}
   fi
   if [ -s ensemble2.${member} ]
   then
       rm  ensemble2.${member}
   fi
   if [ -s ensemble3.${member} ]
   then
       rm  ensemble3.${member}
   fi
   if [ -s ensemble4.${member} ]
   then
       rm  ensemble4.${member}
   fi

   typeset -Z2 HOUR

   HOUR=$FSTART_HIGH

   file=''

   while [ $HOUR -le $FEND_HIGH ]
   do
     if [ "$RUN" = "sref" ]
     then
      if [ $OUTGRD -eq 132 ]; then
       filename1=`eval echo $FILE_HIGH1`
       cat $filename1  >>  pgrb1f${file}    
      else
       filename1=`eval echo $FILE_HIGH1`
       cat $filename1  >>  pgrb1f${file}    

       filename2=`eval echo $FILE_HIGH2`
       cat $filename2  >>  pgrb2f${file}    

       filename3=`eval echo $FILE_HIGH3`
       cat $filename3  >>  pgrb3f${file}    

       filename4=`eval echo $FILE_HIGH4`
       cat $filename4  >>  pgrb4f${file}    
      fi
     else
      filename=`eval echo $FILE_HIGH`
      cat $filename  >>  pgrbf${file}    
     fi

      let HOUR=$HOUR+$FINC_HIGH
   done 

   #
   # IF FSTART_LOW IS SET THEN...
   #
   FSTART_LOW=''
   if [ "$FSTART_LOW" != ""  ]
   then
      HOUR=$FSTART_LOW
      while [ $HOUR -le $FEND_LOW ]
      do
         filename=`eval echo $FILE_LOW`
         cat $filename  >>  pgrbf${file}    

         let HOUR=$HOUR+$FINC_LOW
      done
   fi

#####################################
# Run wgrib to get the inventory of     
#  this new concatenated file.       
#####################################
    
   if [ "$RUN" = "sref" ]
   then
    if [ $OUTGRD -eq 132 ]; then
     $WGRIB -s pgrb1f${file} > pgrb1f${file}.inv
    else
     $WGRIB -s pgrb1f${file} > pgrb1f${file}.inv
     $WGRIB -s pgrb2f${file} > pgrb2f${file}.inv
     $WGRIB -s pgrb3f${file} > pgrb3f${file}.inv
     $WGRIB -s pgrb4f${file} > pgrb4f${file}.inv
    fi
   else
    $WGRIB -s pgrbf${file} > pgrbf${file}.inv
   fi
#####################################
# Use egrep to filter the fields         
#####################################

   if [ "$RUN" = "sref" ]
   then
    if [ $OUTGRD -eq 132 ]; then
     egrep $WGRIB_LIST pgrb1f${file}.inv | $WGRIB pgrb1f${file} -s -grib -i -o ensemble1.${member}
    else
     egrep $WGRIB_LIST pgrb1f${file}.inv | $WGRIB pgrb1f${file} -s -grib -i -o ensemble1.${member}
     egrep $WGRIB_LIST pgrb2f${file}.inv | $WGRIB pgrb2f${file} -s -grib -i -o ensemble2.${member}
     egrep $WGRIB_LIST pgrb3f${file}.inv | $WGRIB pgrb3f${file} -s -grib -i -o ensemble3.${member}
     egrep $WGRIB_LIST pgrb4f${file}.inv | $WGRIB pgrb4f${file} -s -grib -i -o ensemble4.${member}
    fi
   else
    egrep $WGRIB_LIST pgrbf${file}.inv | $WGRIB pgrbf${file} -s -grib -i -o ensemble.${member}
   fi

   #####################################
   # Write the files out to permanent         
   #  diskspace
   #####################################
    if  [ ${member} = "c0" ]
    then
       member="cnt"
    fi

   if [ "$RUN" = "sref" ]
   then
    if [ $OUTGRD -eq 132 ]; then
     filename1=`eval echo $OUTPUT_FILE1`
     cp ensemble1.${member} $filename1
    else 
     filename1=`eval echo $OUTPUT_FILE1`
     cp ensemble1.${member} $filename1

     filename2=`eval echo $OUTPUT_FILE2`
     cp ensemble2.${member} $filename2

     filename3=`eval echo $OUTPUT_FILE3`
     cp ensemble3.${member} $filename3

     filename4=`eval echo $OUTPUT_FILE4`
     cp ensemble4.${member} $filename4
    fi

   else
    filename=`eval echo $OUTPUT_FILE`
    cp ensemble.${member} $filename
   fi

   #################################
   # Convert to grib2 format
   #################################
#  CNVGRIB=$utilexec/cnvgrib
#  WGRIB2=$WGRIB2
   if [ $OUTGRD -eq 132 ]; then
    $CNVGRIB -g12 -p40 -nv $filename1 ${filename1}.grib2
    $WGRIB2 ${filename1}.grib2 -s >${filename1}.grib2.idx
   else
    $CNVGRIB -g12 -p40 -nv $filename1 ${filename1}.grib2
    $CNVGRIB -g12 -p40 -nv $filename2 ${filename2}.grib2
    $CNVGRIB -g12 -p40 -nv $filename3 ${filename3}.grib2
    $CNVGRIB -g12 -p40 -nv $filename4 ${filename4}.grib2
    $WGRIB2 ${filename1}.grib2 -s >${filename1}.grib2.idx
    $WGRIB2 ${filename2}.grib2 -s >${filename2}.grib2.idx
    $WGRIB2 ${filename3}.grib2 -s >${filename3}.grib2.idx
    $WGRIB2 ${filename4}.grib2 -s >${filename4}.grib2.idx
   fi
   
   #cp ${filename1}.grib2 $COMOUT

    if [ "$SENDDBN" = "YES" ]
    then
       if [ "$RUN" = "sref" ]
       then
#        if [ $OUTGRD -eq 132 ]; then
#         $DBNROOT/bin/dbn_alert MODEL $ALERT_TYPE $job $filename1
#        else
#         $DBNROOT/bin/dbn_alert MODEL $ALERT_TYPE $job $filename1
#         $DBNROOT/bin/dbn_alert MODEL $ALERT_TYPE $job $filename2
#         $DBNROOT/bin/dbn_alert MODEL $ALERT_TYPE $job $filename3
#         $DBNROOT/bin/dbn_alert MODEL $ALERT_TYPE $job $filename4
#        fi

          #if [ $SENDDBN_GB2 = YES ]
          #then

           if [ $OUTGRD -eq 132 ]; then
            $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2 $job ${filename1}.grib2
            $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2_WIDX $job ${filename1}.grib2.idx
           else
            $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2 $job ${filename1}.grib2
            $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2 $job ${filename2}.grib2
            $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2 $job ${filename3}.grib2
            $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2 $job ${filename4}.grib2

            $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2_WIDX $job ${filename1}.grib2.idx
            $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2_WIDX $job ${filename2}.grib2.idx
            $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2_WIDX $job ${filename3}.grib2.idx
            $DBNROOT/bin/dbn_alert MODEL ${ALERT_TYPE}_GB2_WIDX $job ${filename4}.grib2.idx
           fi
          
          #fi
       fi
    fi
done        

if [ "$model" = "eta" ]
then
  echo done > $FCSTDIR/${model}ens.done
fi

msg="ENDING sref_thin_post.sh for memeber $memberlist"
postmsg "$jlogfile" "$msg"
