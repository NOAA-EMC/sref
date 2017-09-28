#!/bin/bash

#BSUB -oo _FV3HREFDIR_/run/fv3href_fcst__MEMBER_.log
#BSUB -eo _FV3HREFDIR_/run/fv3href_fcst__MEMBER_.log
#BSUB -J FV3HREF_FCST__MEMBER_
#BSUB -W 02:00
#BSUB -q dev
#BSUB -M 256
#BSUB -P SREF-T2O
#BSUB -extsched 'CRAYLINUX[]'

set -eux

MY_PTMP=/gpfs/hps3/ptmp/$LOGNAME
trap 'echo "FV3HREF error _MEMBER_ FCST" >> $MY_PTMP/fv3href/jlogfile_sref' INT TERM ERR
echo "FV3HREF started _MEMBER_ FCST" >> $MY_PTMP/fv3href/jlogfile_sref

. /opt/modules/default/init/bash
module load PrgEnv-intel
module unload intel
module load intel/16.3.210
module load cray-mpich
module load craype-hugepages4M
module load NetCDF-intel-sandybridge/4.2
module swap pmi/5.0.6-1.0000.10439.140.2.ari pmi/5.0.11
module use /gpfs/hps/nco/ops/nwprod/modulefiles
module load prod_util
module list

ulimit -s unlimited
ulimit -a

export NODES=144
export KMP_AFFINITY=disabled
export OMP_NUM_THREADS=2
export OMP_STACKSIZE=1024m
export PMI_LABEL_ERROUT=1
export MKL_CBWR=AVX2

export RUN_ENVIR=para
export envir=prod
export job=$$

export CONFIG_FILE=_FV3HREFDIR_/parm/fv3href_para_config_cray

export PDY=_PDY_
export cyc=_CYC_
export FLENGTH=_FLENGTH_
export MEMBER=_MEMBER_
export MACHINE=_MACHINE_

export APRUN="aprun -n 1728 -N 12 -j 1 -d 2 -cc depth"

_FV3HREFDIR_/jobs/JFV3HREF_FCST

echo "FV3HREF finished _MEMBER_ FCST" >> $MY_PTMP/fv3href/jlogfile_sref
