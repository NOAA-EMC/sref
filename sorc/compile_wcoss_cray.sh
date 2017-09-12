#!/bin/bash
set -eux

base=`pwd`/..


date
hostname

#######################################

SORCfv3href=`pwd`
EXECfv3href=`pwd`/../exec
mkdir -m 775 -p $EXECfv3href

#
# build preprocessing tools
#
(
    set +x
    . $MODULESHOME/init/bash
    #module load /nwpara2/modulefiles/FV3HREF/v1.0.0_wcoss_cray
    module load $base/modulefiles/FV3HREF/v1.0.0_wcoss_cray
    module list
    set -x

    if true; then
    (
        ###
        ### fre-nctools
        ###
        cd ${SORCfv3href}/fre-nctools.fd
        ./BUILD_TOOLS.csh cray
    )
    fi


    if true; then
    (
        ###
        ### orog
        ###
        cd ${SORCfv3href}/orog.fd
        ./makefile.sh_cray
    )
    fi


    if true; then
    (
        ###
        ### global_chgres
        ###
        cd ${SORCfv3href}/global_chgres.fd
        ./makefile.sh
        if [ -d nst_mask_namchg.fd ]; then
          cd nst_mask_namchg.fd
          ./makefile.sh
        fi
    )
    fi

)

if true; then
(
  ###
  ### nemsfv3gfs
  ###
  # apply atmos_model.F90 patch until FV3 trunk is updated
#  cp ${SORCfv3href}/patches/atmos_model.F90_fixed_for_nest ${SORCfv3href}/nemsfv3gfs.fd/FV3/atmos_model.F90
  cd ${SORCfv3href}/nemsfv3gfs.fd
  cd tests
  ./compile.sh ${SORCfv3href}/nemsfv3gfs.fd/FV3 wcoss_cray "32BIT=Y" 32bit >& make.out.32bit
  cp fv3_32bit.exe ${EXECfv3href}
)
fi

date
