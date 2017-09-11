#!/bin/sh
set -eux

if [ $MACHINE == wcoss_cray ]; then
  . /opt/modules/default/init/bash
  module load PrgEnv-intel
  module unload intel
  module load intel/16.3.210
  module use /gpfs/hps/nco/ops/nwprod/modulefiles
  module load prod_util
  module list
fi

       JJOB=$1
export PDY=$2
export cyc=$3
export MEMBER=$4

export FLENGTH=48

export KMP_AFFINITY=disabled
export OMP_STACKSIZE=1024m
export PMI_LABEL_ERROUT=1
export MKL_CBWR=AVX2

export MP_EUIDEVICE=sn_all
export MP_EUILIB=us
export MP_MPILIB=mpich2

export RUN_ENVIR=para
export envir=prod
export job=$$

export CONFIG_FILE=${FV3HREFDIR}/parm/fv3href_para_config_cray

export RUN_ENVIR=para
export envir=prod
export job=${JJOB}_${MEMBER}_${cyc}

echo "Will now execute ${FV3HREFDIR}/jobs/${JJOB}"
${FV3HREFDIR}/jobs/${JJOB}
