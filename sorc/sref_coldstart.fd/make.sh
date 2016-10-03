#!/bin/sh

#---------------------------------------------------------------------------------
#  The driver script for compiling the emcsfc_coldstart program.  Loads
#  module files and exports environment variables required by the makefile
#  Then, invokes the makefile.  
#
#  Only tested on Zeus and the NCEP WCOSS machines.
#
#  To invoke: type 'make.sh' from the command line.  If successfully built, 
#  the executable will be installed the ../../exec subdirectory.
#
#  See the README.build file for more details.
#---------------------------------------------------------------------------------

#set -x

mac=$(hostname | cut -c1-1)

if [ $mac = t -o $mac = g ] ; then  #For WCOSS

  echo
  echo "BUILD EMCSFC_COLDSTART PROGRAM ON WCOSS"
  echo

  module purge

# load intel compiler

  module load ics/12.1
  export FCOMP=mpfort
  export FFLAGS="-O0 -r8 -i4 -FR -convert big_endian -compiler intel"
  export FPPFLAGS="-fpp -DGFS=0 -save-temps"
  export LFLAGS="-openmp"

# load ncep library modules

  module load NetCDF/3.6.3
  module load sfcio/v1.0.0
  module load ip/v2.0.0
  module load sp/v2.0.2
  module load w3nco/v2.0.6
  module load bacio/v2.0.1
  module load jasper/v1.900.1
  module load z/v1.2.6
  module load png/v1.2.44
  module load g2/v2.5.0
  module load nemsio/v2.2.1
  module load landsfcutil/v2.0.0

# ensure libraries and include directories exist.
  export CHECK_LIBS=/nwprod/spa_util/check_libs.bash
  make check_prereqs
  rc=$?
  if ((rc != 0));then
    echo "MISSING LIBRARY AND/OR INCLUDE DIRECTORY. EXIT."
    exit
  fi

elif [ $mac = s -o $mac = l ] ; then  #For WCOSS_CRAY

  echo
  echo "BUILD EMCSFC_COLDSTART PROGRAM ON WCOSS_CRAY"
  echo


# load intel compiler
  module load PrgEnv-intel
  module load cray-mpich

  export FCOMP=ftn
  export FFLAGS="-O0 -r8 -i4 -FR -convert big_endian"
  export FPPFLAGS="-fpp -DGFS=0 -save-temps"
  export LFLAGS="-openmp"

# load ncep library modules

  module load NetCDF-intel-sandybridge/3.6.3
  module load sfcio-intel/1.0.0
  module load ip-intel/2.0.0
  module load sp-intel/2.0.2
  module load w3nco-intel/2.0.6
  module load bacio-intel/2.0.1
  module load jasper-gnu-sandybridge/1.900.1
  module load zlib-intel-sandybridge/1.2.7
  module load png-intel-sandybridge/1.2.44
  module load g2-intel/2.5.0
  module load nemsio-intel/2.2.2
  module load landsfcutil-intel/2.0.0

# ensure libraries and include directories exist.
  export CHECK_LIBS=/gpfs/hps/nco/ops/nwprod/spa_util/check_libs.bash
  make check_prereqs
  rc=$?
  if ((rc != 0));then
    echo "MISSING LIBRARY AND/OR INCLUDE DIRECTORY. EXIT."
    exit
  fi

elif [ $mac = f ] ; then  #For Zeus

  echo
  echo "BUILD EMCSFC_COLDSTART PROGRAM ON ZEUS"
  echo

  module purge
  module load mpt/2.06
  module load intel/12-12.0.4.191
  export FCOMP=ifort
  export FFLAGS="-O0 -convert big_endian -r8 -i4 -FR"
  export FPPFLAGS="-fpp -DGFS=0 -save-temps"
  export LFLAGS="-openmp -lmpi"

# load ncep library modules

  module load netcdf/3.6.3-intel
  export NETCDF_INCLUDE=-I${NETCDF}/include
  export NETCDF_LDFLAGS_F="-L${NETCDF}/lib -lnetcdf"
  module load sfcio/v1.1.0
  module load ip/v2.0.0
  module load sp/v2.0.2
  module load w3nco/v2.0.6
  module load bacio/v2.0.1
  module load jasper/1.900.1
  module load z/1.2.6
  module load png/12
  module load g2/v2.5.0
  module load nemsio/v2.2.1
  module load landsfcutil/v2.0.0

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
