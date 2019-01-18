#!/bin/sh

#---------------------------------------------------------------------------------
#  The driver script for compiling the emcsfc_coldstart program.  Loads
#  module files and exports environment variables required by the makefile
#  Then, invokes the makefile.  Only works on WCOSS Dell.
#
#  To invoke: type 'make.sh' from the command line.  If successfully built, 
#  the executable will be installed the ../../exec subdirectory.
#
#  See the README.build file for more details.
#---------------------------------------------------------------------------------

#set -x

mac=$(hostname | cut -c1-1)

if [ $mac = v -o $mac = m ] ; then  # WCOSS Dell

  echo
  echo "BUILD EMCSFC_COLDSTART PROGRAM ON WCOSS DELL."
  echo

  module purge
  module load EnvVars/1.0.2
  module load ips/18.0.1.163
  module load impi/18.0.1
  module load NetCDF/3.6.3
  module load sfcio/1.0.0
  module load ip/3.0.1
  module load sp/2.0.2
  module load w3nco/2.0.6
  module load bacio/2.0.2
  module load jasper/1.900.1
  module load zlib/1.2.11
  module load libpng/1.2.44
  module load g2/3.1.0
  module load nemsio/2.2.3
  module load landsfcutil/2.1.0
  module list

  export FCOMP=mpif90
  export FFLAGS="-O0 -r8 -i4 -FR -convert big_endian"
  export FPPFLAGS="-fpp -DGFS=0 -save-temps"
  export LFLAGS="-qopenmp"

else

  echo "MACHINE OPTION NOT FOUND. EXIT."
  exit

fi

# invoke the makefile

make clean
make
rc=$?

if ((rc != 0));then
  echo "BUILD FAILED. EXIT."
  exit
else
  make install
fi

exit
