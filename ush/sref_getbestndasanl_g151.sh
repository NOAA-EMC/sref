set -eux

DATE=$1

echo $DATE > datem00

DATEB03=`$NDATE -3 $DATE`
echo $DATEB03 > dateb03
DATEM03=`$NDATE +3 $DATE`
echo $DATEM03 > datem03
DATEM06=`$NDATE +6 $DATE`
echo $DATEM06 > datem06
DATEM09=`$NDATE +9 $DATE`
echo $DATEM09 > datem09
DATEM12=`$NDATE +12 $DATE`
echo $DATEM12 > datem12

PDYb3=`cut -c 1-8 dateb03`
HHb3=`cut -c 9-10 dateb03`
PDY00=`cut -c 1-8 datem00`
HH00=`cut -c 9-10 datem00`
PDY03=`cut -c 1-8 datem03`
HH03=`cut -c 9-10 datem03`
PDY06=`cut -c 1-8 datem06`
HH06=`cut -c 9-10 datem06`
PDY09=`cut -c 1-8 datem09`
HH09=`cut -c 9-10 datem09`
PDY12=`cut -c 1-8 datem12`
HH12=`cut -c 9-10 datem12`

case $HH00 in
 00) ndas1=$COMINnam/ndas.$PDY12/ndas.t${HH12}z.awp15100.tm12
     ndas2=$COMINnam/ndas.$PDY06/ndas.t${HH06}z.awp15100.tm06
     ndas3=$COMINnam/nam.$PDY00/nam.t${HH00}z.awp15100.tm00;;
 03) ndas1=$COMINnam/ndas.$PDY09/ndas.t${HH09}z.awp15100.tm09
     ndas2=$COMINnam/ndas.$PDY03/ndas.t${HH03}z.awp15100.tm03
     ndas3=$COMINnam/nam.$PDYb3/nam.t${HHb3}z.awp15103.tm00;;
 06) ndas1=$COMINnam/ndas.$PDY12/ndas.t${HH12}z.awp15100.tm12
     ndas2=$COMINnam/ndas.$PDY06/ndas.t${HH06}z.awp15100.tm06
     ndas3=$COMINnam/nam.$PDY00/nam.t${HH00}z.awp15100.tm00;;
 09) ndas1=$COMINnam/ndas.$PDY09/ndas.t${HH09}z.awp15100.tm09
     ndas2=$COMINnam/ndas.$PDY03/ndas.t${HH03}z.awp15100.tm03
     ndas3=$COMINnam/nam.$PDYb3/nam.t${HHb3}z.awp15103.tm00;;
 12) ndas1=$COMINnam/ndas.$PDY12/ndas.t${HH12}z.awp15100.tm12
     ndas2=$COMINnam/ndas.$PDY06/ndas.t${HH06}z.awp15100.tm06
     ndas3=$COMINnam/nam.$PDY00/nam.t${HH00}z.awp15100.tm00;;
 15) ndas1=$COMINnam/ndas.$PDY09/ndas.t${HH09}z.awp15100.tm09
     ndas2=$COMINnam/ndas.$PDY03/ndas.t${HH03}z.awp15100.tm03
     ndas3=$COMINnam/nam.$PDYb3/nam.t${HHb3}z.awp15103.tm00;;
 18) ndas1=$COMINnam/ndas.$PDY12/ndas.t${HH12}z.awp15100.tm12
     ndas2=$COMINnam/ndas.$PDY06/ndas.t${HH06}z.awp15100.tm06
     ndas3=$COMINnam/nam.$PDY00/nam.t${HH00}z.awp15100.tm00;;
 21) ndas1=$COMINnam/ndas.$PDY09/ndas.t${HH09}z.awp15100.tm09
     ndas2=$COMINnam/ndas.$PDY03/ndas.t${HH03}z.awp15100.tm03
     ndas3=$COMINnam/nam.$PDYb3/nam.t${HHb3}z.awp15103.tm00;;
esac

if [ -s $ndas1 ]; then
  cp $ndas1 awp15100.${DATE}
elif [ -s $ndas2 ]; then
  cp $ndas2 awp15100.${DATE}
elif [ -s $ndas3 ]; then
  cp $ndas3 awp15100.${DATE}
else
  set +x
  echo "**********************************************************"
  echo "NO NDAS DATA AVAILABLE, NEED AT LEAST ONE OF FOLLOWING:"
  echo $ndas1
  echo $ndas2
  echo $ndas3
  echo "PLEASE CHECK IF RAP FINISHES FINE!"
  echo "**********************************************************"
  set -x
  export err=911
  err_chk
fi

exit

