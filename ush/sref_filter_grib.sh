#!/bin/sh

set -eux

if [ $MACHINE == wcoss ]; then

##  . /usrx/local/Modules/default/init/sh
##  module load cfp
###. /opt/modules/default/init/bash
###module load PrgEnv-intel

  CMDFILE=cmdfile
  rm -f $CMDFILE

  n=0
  for GRIBFILE in $*; do

cat <<- END_OF_JOB > job.$n
#!/bin/sh
set -eux
echo "$0: filtering ${GRIBFILE}"
$WGRIB -s ${GRIBFILE} | grep -v ":1 mb:" \\
                               | grep -v ":2 mb:" \\
                               | grep -v ":3 mb:" \\
                               | grep -v ":5 mb:" \\
                               | grep -v ":7 mb:" \\
                               | grep -v ":75 mb:" \\
                               | $WGRIB -i -grib ${GRIBFILE} -o ${GRIBFILE}_new
mv ${GRIBFILE}_new ${GRIBFILE}
END_OF_JOB

  echo "sh `pwd`/job.$n" >> $CMDFILE
  (( n=n+1 ))

  done

#export NODES=2
#module load cfp-intel-sandybridge

##module load PrgEnv-intel                      #not needed for static exe
##module list
chmod +x $CMDFILE
##  mpirun.lsf cfp $CMDFILE
# aprun -n 16 -N 8 -d 3  cfp $CMDFILE
$CMDFILE

else

  for GRIBFILE in $*; do

  echo "$0: filtering ${GRIBFILE}"
   
  ${WGRIB} -s ${GRIBFILE} | grep -v ":1 mb:" \
                                 | grep -v ":2 mb:" \
                                 | grep -v ":3 mb:" \
                                 | grep -v ":5 mb:" \
                                 | grep -v ":7 mb:" \
                                 | grep -v ":75 mb:" \
                                 | ${WGRIB} -i -grib ${GRIBFILE} -o new_grib_file
  
  mv new_grib_file ${GRIBFILE}
  
  done
fi
