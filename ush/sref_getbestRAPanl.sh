set -eux

DATE=$1

echo $DATE > datem00

DATEM01=`$NDATE -1 $DATE`
echo $DATEM01 > datem01
DATEM02=`$NDATE -2 $DATE`
echo $DATEM02 > datem02
DATEM03=`$NDATE -3 $DATE`
echo $DATEM03 > datem03

PDY00=`cut -c 1-8 datem00`
HH00=`cut -c 9-10 datem00`
PDY01=`cut -c 1-8 datem01`
HH01=`cut -c 9-10 datem01`
PDY02=`cut -c 1-8 datem02`
HH02=`cut -c 9-10 datem02`
PDY03=`cut -c 1-8 datem03`
HH03=`cut -c 9-10 datem03`

case $HH00 in
 00) rap1=$COMINrap/rap.$PDY00/rap.t${HH00}z.awp130pgrbf00.grib2
     rap2=$COMINrap/rap.$PDY01/rap.t${HH01}z.awp130pgrbf01.grib2
     rap3=$COMINrap/rap.$PDY02/rap.t${HH02}z.awp130pgrbf02.grib2
     rap4=$COMINrap/rap.$PDY03/rap.t${HH03}z.awp130pgrbf03.grib2;;
 03) rap1=$COMINrap/rap.$PDY00/rap.t${HH00}z.awp130pgrbf00.grib2
     rap2=$COMINrap/rap.$PDY01/rap.t${HH01}z.awp130pgrbf01.grib2
     rap3=$COMINrap/rap.$PDY02/rap.t${HH02}z.awp130pgrbf02.grib2
     rap4=$COMINrap/rap.$PDY03/rap.t${HH03}z.awp130pgrbf03.grib2;;
 06) rap1=$COMINrap/rap.$PDY00/rap.t${HH00}z.awp130pgrbf00.grib2
     rap2=$COMINrap/rap.$PDY01/rap.t${HH01}z.awp130pgrbf01.grib2
     rap3=$COMINrap/rap.$PDY02/rap.t${HH02}z.awp130pgrbf02.grib2
     rap4=$COMINrap/rap.$PDY03/rap.t${HH03}z.awp130pgrbf03.grib2;;
 09) rap1=$COMINrap/rap.$PDY00/rap.t${HH00}z.awp130pgrbf00.grib2
     rap2=$COMINrap/rap.$PDY01/rap.t${HH01}z.awp130pgrbf01.grib2
     rap3=$COMINrap/rap.$PDY02/rap.t${HH02}z.awp130pgrbf02.grib2
     rap4=$COMINrap/rap.$PDY03/rap.t${HH03}z.awp130pgrbf03.grib2;;
 12) rap1=$COMINrap/rap.$PDY00/rap.t${HH00}z.awp130pgrbf00.grib2
     rap2=$COMINrap/rap.$PDY01/rap.t${HH01}z.awp130pgrbf01.grib2
     rap3=$COMINrap/rap.$PDY02/rap.t${HH02}z.awp130pgrbf02.grib2
     rap4=$COMINrap/rap.$PDY03/rap.t${HH03}z.awp130pgrbf03.grib2;;
 15) rap1=$COMINrap/rap.$PDY00/rap.t${HH00}z.awp130pgrbf00.grib2
     rap2=$COMINrap/rap.$PDY01/rap.t${HH01}z.awp130pgrbf01.grib2
     rap3=$COMINrap/rap.$PDY02/rap.t${HH02}z.awp130pgrbf02.grib2
     rap4=$COMINrap/rap.$PDY03/rap.t${HH03}z.awp130pgrbf03.grib2;;
 18) rap1=$COMINrap/rap.$PDY00/rap.t${HH00}z.awp130pgrbf00.grib2
     rap2=$COMINrap/rap.$PDY01/rap.t${HH01}z.awp130pgrbf01.grib2
     rap3=$COMINrap/rap.$PDY02/rap.t${HH02}z.awp130pgrbf02.grib2
     rap4=$COMINrap/rap.$PDY03/rap.t${HH03}z.awp130pgrbf03.grib2;;
 21) rap1=$COMINrap/rap.$PDY00/rap.t${HH00}z.awp130pgrbf00.grib2
     rap2=$COMINrap/rap.$PDY01/rap.t${HH01}z.awp130pgrbf01.grib2
     rap3=$COMINrap/rap.$PDY02/rap.t${HH02}z.awp130pgrbf02.grib2
     rap4=$COMINrap/rap.$PDY03/rap.t${HH03}z.awp130pgrbf03.grib2;;

# Note this new larger domain data has several extra pressure levels including 75mb
#00) rap1=$COMINrap/rap.$PDY00/rap.t${HH00}z.awip32f00.grib2
#    rap2=$COMINrap/rap.$PDY01/rap.t${HH01}z.awip32f01.grib2
#    rap3=$COMINrap/rap.$PDY02/rap.t${HH02}z.awip32f02.grib2
#    rap4=$COMINrap/rap.$PDY03/rap.t${HH03}z.awip32f03.grib2;;
#03) rap1=$COMINrap/rap.$PDY00/rap.t${HH00}z.awip32f00.grib2
#    rap2=$COMINrap/rap.$PDY01/rap.t${HH01}z.awip32f01.grib2
#    rap3=$COMINrap/rap.$PDY02/rap.t${HH02}z.awip32f02.grib2
#    rap4=$COMINrap/rap.$PDY03/rap.t${HH03}z.awip32f03.grib2;;
#06) rap1=$COMINrap/rap.$PDY00/rap.t${HH00}z.awip32f00.grib2
#    rap2=$COMINrap/rap.$PDY01/rap.t${HH01}z.awip32f01.grib2
#    rap3=$COMINrap/rap.$PDY02/rap.t${HH02}z.awip32f02.grib2
#    rap4=$COMINrap/rap.$PDY03/rap.t${HH03}z.awip32f03.grib2;;
#09) rap1=$COMINrap/rap.$PDY00/rap.t${HH00}z.awip32f00.grib2
#    rap2=$COMINrap/rap.$PDY01/rap.t${HH01}z.awip32f01.grib2
#    rap3=$COMINrap/rap.$PDY02/rap.t${HH02}z.awip32f02.grib2
#    rap4=$COMINrap/rap.$PDY03/rap.t${HH03}z.awip32f03.grib2;;
#12) rap1=$COMINrap/rap.$PDY00/rap.t${HH00}z.awip32f00.grib2
#    rap2=$COMINrap/rap.$PDY01/rap.t${HH01}z.awip32f01.grib2
#    rap3=$COMINrap/rap.$PDY02/rap.t${HH02}z.awip32f02.grib2
#    rap4=$COMINrap/rap.$PDY03/rap.t${HH03}z.awip32f03.grib2;;
#15) rap1=$COMINrap/rap.$PDY00/rap.t${HH00}z.awip32f00.grib2
#    rap2=$COMINrap/rap.$PDY01/rap.t${HH01}z.awip32f01.grib2
#    rap3=$COMINrap/rap.$PDY02/rap.t${HH02}z.awip32f02.grib2
#    rap4=$COMINrap/rap.$PDY03/rap.t${HH03}z.awip32f03.grib2;;
#18) rap1=$COMINrap/rap.$PDY00/rap.t${HH00}z.awip32f00.grib2
#    rap2=$COMINrap/rap.$PDY01/rap.t${HH01}z.awip32f01.grib2
#    rap3=$COMINrap/rap.$PDY02/rap.t${HH02}z.awip32f02.grib2
#    rap4=$COMINrap/rap.$PDY03/rap.t${HH03}z.awip32f03.grib2;;
#21) rap1=$COMINrap/rap.$PDY00/rap.t${HH00}z.awip32f00.grib2
#    rap2=$COMINrap/rap.$PDY01/rap.t${HH01}z.awip32f01.grib2
#    rap3=$COMINrap/rap.$PDY02/rap.t${HH02}z.awip32f02.grib2
#    rap4=$COMINrap/rap.$PDY03/rap.t${HH03}z.awip32f03.grib2;;
esac

if [ -s $rap1 ]; then
  $CNVGRIB -g21 $rap1 rap.t${HH00}z.awp130f00
elif [ -s $rap2 ]; then
  $CNVGRIB -g21 $rap2 rap.t${HH00}z.awp130f00
elif [ -s $rap3 ]; then
  $CNVGRIB -g21 $rap3 rap.t${HH00}z.awp130f00
elif [ -s $rap4 ]; then
  $CNVGRIB -g21 $rap4 rap.t${HH00}z.awp130f00
else
  set +x
  echo "**********************************************************"
  echo "NO RAP DATA AVAILABLE, NEED AT LEAST ONE OF FOLLOWING:"
  echo $rap1
  echo $rap2
  echo $rap3
  echo $rap4
  echo "PLEASE CHECK IF RAP FINISHES FINE!"
  echo "**********************************************************"
  set -x
  export err=911
  err_chk
fi

exit

