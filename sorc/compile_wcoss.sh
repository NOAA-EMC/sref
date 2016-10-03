#!/bin/sh
set -eux

base=`pwd`/..

set +x
module purge
#module load /nwpara2/modulefiles/SREF/v7.1.0_wcoss
module load $base/modulefiles/SREF/v7.1.0_wcoss
module list
set -x

date
hostname

#######################################

SORCsref=`pwd`
EXECsref=`pwd`/../exec
mkdir -m 775 -p $EXECsref

if [ 1 == 1 ]; then
###
### ARW WRF
###

cd ${SORCsref}/sref_wrf_v3.5.1.fd
cp ${SORCsref}/sref_wrf_v3.5.1.fd/phys/module_sf_noahlsm.F_driersoil ${SORCsref}/sref_wrf_v3.5.1.fd/phys/module_sf_noahlsm.F
./compile_arw_wcoss.sh
ls -l ${EXECsref}/sref_wrf_arw
mv ${EXECsref}/sref_wrf_arw ${EXECsref}/sref_wrf_arw_driersoil

cp ${SORCsref}/sref_wrf_v3.5.1.fd/phys/module_sf_noahlsm.F_normalsoil ${SORCsref}/sref_wrf_v3.5.1.fd/phys/module_sf_noahlsm.F
./compile_arw_wcoss.sh
ls -l ${EXECsref}/sref_wrf_arw
mv ${EXECsref}/sref_wrf_arw ${EXECsref}/sref_wrf_arw_normalsoil

cd ..

fi


if [ 1 == 1 ]; then
###
### WPS
###

cd ${SORCsref}/sref_wps_v3.5.1.fd
./clean -a
./configure <<EOF
23
EOF
./compile wps

cp ungrib.exe  ${EXECsref}/sref_ungrib
cp geogrid.exe ${EXECsref}/sref_geogrid
cp metgrid.exe ${EXECsref}/sref_metgrid

./compile mod_levs
cp util/src/mod_levs.exe  ${EXECsref}/sref_ARW_levs

cd ..

fi


if [ 1 == 1 ]; then
###
### NPS
###

cd ${SORCsref}/sref_nps.fd
./conf wcoss
./clean -a
./configure <<EOF
2
EOF
./compile nps

cp ungrib.exe     ${EXECsref}/sref_ungrib_nems
cp ungribp.exe    ${EXECsref}/sref_ungribp_nems
cp geogrid.exe    ${EXECsref}/sref_geogrid_nems
cp metgrid.exe    ${EXECsref}/sref_metgrid_nems
cp nemsinterp.exe ${EXECsref}/sref_real_nmb

./compile mod_levs
cp util/src/mod_levs.exe  ${EXECsref}/sref_NMB_levs

cd ..

fi

if [ 1 == 1 ]; then
###
### NMMB
###
cd ${SORCsref}/sref_nmmb.fd/src
cp ${SORCsref}/sref_nmmb.fd/src/atmos/phys/module_LS_NOAHLSM.F90_driersoil ${SORCsref}/sref_nmmb.fd/src/atmos/phys/module_LS_NOAHLSM.F90
./configure 6.3r_nmm_wcoss
source conf/modules.nems
gmake clean
gmake nmm
mv ../exe/NEMS.x        ${EXECsref}/sref_wrf_nmb_driersoil
gmake clean

cd ${SORCsref}/sref_nmmb.fd/src
cp ${SORCsref}/sref_nmmb.fd/src/atmos/phys/module_LS_NOAHLSM.F90_normalsoil ${SORCsref}/sref_nmmb.fd/src/atmos/phys/module_LS_NOAHLSM.F90
./configure 6.3r_nmm_wcoss
source conf/modules.nems
gmake clean
gmake nmm
mv ../exe/NEMS.x        ${EXECsref}/sref_wrf_nmb_normalsoil
gmake clean
fi

set +x
module purge
#module load /nwpara2/modulefiles/SREF/v7.1.0_wcoss
module load $base/modulefiles/SREF/v7.1.0_wcoss
module list
set -x

if [ 1 == 1 ]; then
###
### PERTURB
###

cd ${SORCsref}/sref_perturb.fd
cd dio
sh compile_wcoss.sh
cd ..
make -f Makefile clean
make -f Makefile
cp breeding_arw.exe    ${EXECsref}/sref_breeding_arw
#cp breeding_nmm.exe    ${EXECsref}/sref_breeding_nmm
cp breeding_nmb.exe    ${EXECsref}/sref_breeding_nmb
cp lbc_perturb_wrf.exe ${EXECsref}/sref_lbc_perturb_wrf
cp lbc_perturb_nmb.exe ${EXECsref}/sref_lbc_perturb_nmb

cd ..

fi

if [ 1 == 1 ]; then
###
### POST
###
cd ${SORCsref}/sref_post.fd
make -f makefile clean
make -f makefile

mv sref_post     ${EXECsref}/.
rm *.o *.mod

cd ..

fi

if [ 1 == 1 ]; then
###
### COLDSTART
###
cd ${SORCsref}/sref_coldstart.fd
./make.sh

mv sref_coldstart_wrf     ${EXECsref}/.

cd ..

fi

if [ 1 == 1 ]; then
###
### PRDGEN
###
cd ${SORCsref}/sref_prdgen.fd
make -f Makefile_wcoss

mv sref_prdgen     ${EXECsref}/.

cd ..
fi


if [ 1 == 1 ]; then
###
### WRFBUCKET
###
cd ${SORCsref}/sref_wrfbucket.fd
make -f Makefile_wcoss
mv sref_pcpbucket_g212 ${EXECsref}/.
mv sref_pcpbucket_g216 ${EXECsref}/.
mv sref_pcpbucket_g221 ${EXECsref}/.
mv sref_pcpbucket_g243 ${EXECsref}/.
mv sref_pcpbucket_g132 ${EXECsref}/.

cd ..

fi

# ENSADD
if [ 1 == 1 ]; then
cd ${SORCsref}/global_ensadd.fd
make -f makefile_wcoss
mv global_ensadd ${EXECsref}/.
cd ..
fi

# global_postgs
if [ 1 == 1 ]; then
cd ${SORCsref}/global_postgp.fd
make -f makefile_wcoss clean
make -f makefile_wcoss global_postgs
mv global_postgs ${EXECsref}/.
make -f makefile_wcoss clean
cd ..
fi

# sref_biasestimate.fd

cd ${SORCsref}/sref_biasestimate.fd
make -f makefile_wcoss clean
make -f makefile_wcoss
mv sref_estimate_bias ${EXECsref}/.
make -f makefile_wcoss clean
cd ..

# sref_bufr.fd
cd ${SORCsref}/sref_bufr.fd
make -f makefile_wcoss clean
make -f makefile_wcoss
mv sref_bufr ${EXECsref}/.
make -f makefile_wcoss clean
cd ..

# sref_calfcsterr.fd
cd ${SORCsref}/sref_calfcsterr.fd
make -f makefile_wcoss clean
make -f makefile_wcoss
mv sref_calfcsterr ${EXECsref}/.
make -f makefile_wcoss clean
cd ..

# sref_cluster_NCEP.fd
cd ${SORCsref}/sref_cluster_NCEP.fd
rm -f *.o
make -f makefile_cluster_wcoss
mv sref_cluster_NCEP ${EXECsref}/.
rm -f *.o
make -f makefile_weight_wcoss
mv sref_clusterweight ${EXECsref}/.
rm -f *.o
cd ..

# sref_cluster_OU.fd
cd ${SORCsref}/sref_cluster_OU.fd
rm -f *.o
make -f makefile_wcoss
mv sref_cluster_OU ${EXECsref}/.
make -f makefile_wcoss clean
cd ..

# sref_dwnsfcst.fd
cd ${SORCsref}/sref_dwnsfcst.fd
make -f makefile_wcoss clean
make -f makefile_wcoss
mv sref_dvrtma_debias ${EXECsref}/.
make -f makefile_wcoss clean
cd ..

# sref_dwnsvect.fd
cd ${SORCsref}/sref_dwnsvect.fd
make -f makefile_wcoss clean
make -f makefile_wcoss
mv sref_dvrtma_bias ${EXECsref}/.
make -f makefile_wcoss clean
cd ..

# sref_ens_gen.fd
cd ${SORCsref}/sref_ens_gen.fd
rm -f *.o *.mod
make -f makefile_REG_wcoss
mv sref_ens_gen ${EXECsref}/.
rm -f *.o *.mod
make -f makefile_DS_wcoss
mv sref_ens_gen_DS ${EXECsref}/.
rm -f *.o *.mod
cd ../

# sref_fastcopygb.fd
cd ${SORCsref}/sref_fastcopygb.fd
make -f makefile_wcoss
mv fastcopygb ${EXECsref}/.
make -f makefile_wcoss clean
cd ..

# sref_meansndp.fd
cd ${SORCsref}/sref_meansndp.fd
make -f makefile_wcoss clean
make -f makefile_wcoss
mv sref_meansndp ${EXECsref}/.
make -f makefile_wcoss clean
cd ..

# sref_memberranking.fd
cd ${SORCsref}/sref_memberranking.fd
make -f makefile_wcoss clean
make -f makefile_wcoss
mv sref_ranking ${EXECsref}/.
make -f makefile_wcoss clean
cd ..

# sref_qpfbiasestimate.fd
cd ${SORCsref}/sref_qpfbiasestimate.fd
make -f makefile_meanqpf_wcoss clean
make -f makefile_meanqpf_wcoss
mv sref_estimate_meanqpfbias ${EXECsref}/.
make -f makefile_meanqpf_wcoss clean
make -f makefile_qpf_wcoss
mv sref_estimate_qpfbias ${EXECsref}/.
make -f makefile_qpf_wcoss clean
cd ..

# sref_qpfcalfcsterr.fd
cd ${SORCsref}/sref_qpfcalfcsterr.fd
rm -f *.o
make -f makefile_meanqpf_wcoss
mv sref_cal_meanqpffcsterr ${EXECsref}/.
rm -f *.o
make -f makefile_qpf_wcoss
mv sref_cal_qpffcsterr ${EXECsref}/.
rm -f *.o
cd ..

# sref_sndp.fd
cd ${SORCsref}/sref_sndp.fd
rm -f *.o
make -f makefile_em_wcoss
make -f makefile_nmmb_wcoss
mv sref_*sndp ${EXECsref}/.
rm -f *.o
cd ..



date
