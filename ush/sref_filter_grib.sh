#!/bin/sh

set -eux

if [ $MACHINE == dell ]; then

#  . /usrx/local/Modules/default/init/sh
#  module load cfp

  CMDFILE=cmdfile
  rm -f $CMDFILE

  n=0
  for GRIBFILE in $*; do

cat <<- END_OF_JOB > job.$n
#!/bin/sh
set -eux
echo "$0: filtering ${GRIBFILE}"
${WGRIB:?} -s ${GRIBFILE} | grep -v ":0 mb:" \\
                               | grep -v ":1 mb:" \\
                               | grep -v ":2 mb:" \\
                               | grep -v ":3 mb:" \\
                               | grep -v ":5 mb:" \\
                               | grep -v ":7 mb:" \\
                               | grep -v ":15 mb:" \\
                               | grep -v ":40 mb:" \\
                               | grep -v ":75 mb:" \\
                               | grep -v ":125 mb:" \\
                               | grep -v ":175 mb:" \\
                               | grep -v ":225 mb:" \\
                               | grep -v ":275 mb:" \\
                               | grep -v ":325 mb:" \\
                               | grep -v ":375 mb:" \\
                               | grep -v ":425 mb:" \\
                               | grep -v ":475 mb:" \\
                               | grep -v ":525 mb:" \\
                               | grep -v ":575 mb:" \\
                               | grep -v ":625 mb:" \\
                               | grep -v ":675 mb:" \\
                               | grep -v ":725 mb:" \\
                               | grep -v ":775 mb:" \\
                               | grep -v ":825 mb:" \\
                               | grep -v ":875 mb:" \\
                               | ${WGRIB:?} -i -grib ${GRIBFILE} -o ${GRIBFILE}_new
mv ${GRIBFILE}_new ${GRIBFILE}
END_OF_JOB

  echo "sh `pwd`/job.$n" >> $CMDFILE
  (( n=n+1 ))

  done

# mpirun.lsf cfp $CMDFILE
  mpirun cfp $CMDFILE

else

  for GRIBFILE in $*; do

  echo "$0: filtering ${GRIBFILE}"
   
  ${WGRIB:?} -s ${GRIBFILE} | grep -v ":0 mb:" \
                                 | grep -v ":1 mb:" \
                                 | grep -v ":2 mb:" \
                                 | grep -v ":3 mb:" \
                                 | grep -v ":5 mb:" \
                                 | grep -v ":7 mb:" \
                                 | grep -v ":15 mb:" \
                                 | grep -v ":40 mb:" \
                                 | grep -v ":75 mb:" \
                                 | grep -v ":125 mb:" \
                                 | grep -v ":175 mb:" \
                                 | grep -v ":225 mb:" \
                                 | grep -v ":275 mb:" \
                                 | grep -v ":325 mb:" \
                                 | grep -v ":375 mb:" \
                                 | grep -v ":425 mb:" \
                                 | grep -v ":475 mb:" \
                                 | grep -v ":525 mb:" \
                                 | grep -v ":575 mb:" \
                                 | grep -v ":625 mb:" \
                                 | grep -v ":675 mb:" \
                                 | grep -v ":725 mb:" \
                                 | grep -v ":775 mb:" \
                                 | grep -v ":825 mb:" \
                                 | grep -v ":875 mb:" \
                                 | ${WGRIB:?} -i -grib ${GRIBFILE} -o new_grib_file
  
  mv new_grib_file ${GRIBFILE}
  
  done
fi
