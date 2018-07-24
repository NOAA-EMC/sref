#!/bin/sh
set -aux

base=$(dirname $PWD)

module purge
module use ${base}/modulefiles
module load SREF/v7.0
module list

date

#######################################

SORCsref=${base}/sorc
EXECsref=${base}/exec
mkdir -m 775 -p $EXECsref

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
rm *.o *.mod
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
#ln -sf makefile.wcoss.conf makefile.conf
#make -f makefile
make.sh

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
make -f Makefile
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
make -f makefile
mv global_ensadd ${EXECsref}/.
cd ..
fi

# global_postgs
if [ 1 == 1 ]; then
cd ${SORCsref}/global_postgp.fd
make.sh
mv global_postgs ${EXECsref}/.
cd ..
fi

# sref_biasestimate.fd

cd ${SORCsref}/sref_biasestimate.fd
make clean
make
mv sref_estimate_bias ${EXECsref}/.
make clean

cd ..

# sref_bufr.fd
cd ${SORCsref}/sref_bufr.fd
make clean
make
mv sref_bufr ${EXECsref}/.
make clean

cd ..

# sref_calfcsterr.fd
cd ${SORCsref}/sref_calfcsterr.fd
make clean
make
mv sref_calfcsterr ${EXECsref}/.
make clean
cd ..

# sref_cluster_NCEP.fd
cd ${SORCsref}/sref_cluster_NCEP.fd
rm *.o
make -f makefile_cluster
mv sref_cluster_NCEP ${EXECsref}/.
rm *.o
make -f makefile_weight
mv sref_clusterweight ${EXECsref}/.
rm *.o
cd ..

# sref_cluster_OU.fd
cd ${SORCsref}/sref_cluster_OU.fd
rm *.o
make -f makefile
mv sref_cluster_OU ${EXECsref}/.
rm *.o
cd ..

# sref_dwnsfcst.fd
cd ${SORCsref}/sref_dwnsfcst.fd
make clean
make
mv sref_dvrtma_debias ${EXECsref}/.
make clean
cd ..

# sref_dwnsvect.fd
cd ${SORCsref}/sref_dwnsvect.fd
make clean
make
mv sref_dvrtma_bias ${EXECsref}/.
make clean
cd ..

# sref_ens_gen.fd
cd ${SORCsref}/sref_ens_gen.fd
rm -f *.o *.mod
make -f makefile_REG
mv sref_ens_gen ${EXECsref}/.
rm *.o *.mod
make -f makefile_DS
mv sref_ens_gen_DS ${EXECsref}/.
rm *.o *.mode
cd ../

# sref_fastcopygb.fd
cd ${SORCsref}/sref_fastcopygb.fd
make
mv fastcopygb ${EXECsref}/.
rm *.o 
cd ..

# sref_meansndp.fd
cd ${SORCsref}/sref_meansndp.fd
make clean
make
mv sref_meansndp ${EXECsref}/.
make clean
cd ..

# sref_memberranking.fd
cd ${SORCsref}/sref_memberranking.fd
make clean
make
mv sref_ranking ${EXECsref}/.
make clean
cd ..

# sref_qpfbiasestimate.fd
cd ${SORCsref}/sref_qpfbiasestimate.fd
rm *.o
make -f makefile_meanqpf
mv sref_estimate_meanqpfbias ${EXECsref}/.
rm *.o
make -f makefile_qpf
mv sref_estimate_qpfbias ${EXECsref}/.
rm *.o
cd ..

# sref_qpfcalfcsterr.fd
cd ${SORCsref}/sref_qpfcalfcsterr.fd
rm *.o
make -f makefile_meanqpf
mv sref_cal_meanqpffcsterr ${EXECsref}/.
rm *.o
make -f makefile_qpf
mv sref_cal_qpffcsterr ${EXECsref}/.
rm *.o
cd ..

# sref_sndp.fd
cd ${SORCsref}/sref_sndp.fd
rm *.o
chmod +x make_all
./make_all
mv sref_*sndp ${EXECsref}/.
rm *.o
cd ..

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
### old NMMB
###
#cd ${SORCsref}/sref_nems.fd/src
#cp ${SORCsref}/sref_nems.fd/src/atmos/phys/module_LS_NOAHLSM.F90_driersoil ${SORCsref}/sref_nems.fd/src/atmos/phys/module_LS_NOAHLSM.F90
#./esmf_version 3_wcoss
#source conf/modules.nems.wcoss
#gmake clean
#gmake nmm
#cd ../exe
#cp NEMS.x        ${EXECsref}/sref_wrf_nmb_driersoil
#
#cd ${SORCsref}/sref_nems.fd/src
#cp ${SORCsref}/sref_nems.fd/src/atmos/phys/module_LS_NOAHLSM.F90_normalsoil ${SORCsref}/sref_nems.fd/src/atmos/phys/module_LS_NOAHLSM.F90
#./esmf_version 3_wcoss
#source conf/modules.nems.wcoss
#gmake clean
#gmake nmm
#cd ../exe
#cp NEMS.x        ${EXECsref}/sref_wrf_nmb_normalsoil

###
### new NMMB
###
cd ${SORCsref}/sref_nmmb.fd/src
cp ${SORCsref}/sref_nmmb.fd/src/atmos/phys/module_LS_NOAHLSM.F90_driersoil ${SORCsref}/sref_nmmb.fd/src/atmos/phys/module_LS_NOAHLSM.F90
./configure 6.3r_nmm_wcoss
source conf/modules.nems
gmake clean
gmake nmm
cd ../exe
cp NEMS.x        ${EXECsref}/sref_wrf_nmb_driersoil

cd ${SORCsref}/sref_nmmb.fd/src
cp ${SORCsref}/sref_nmmb.fd/src/atmos/phys/module_LS_NOAHLSM.F90_normalsoil ${SORCsref}/sref_nmmb.fd/src/atmos/phys/module_LS_NOAHLSM.F90
./configure 6.3r_nmm_wcoss
source conf/modules.nems
gmake clean
gmake nmm
cd ../exe
cp NEMS.x        ${EXECsref}/sref_wrf_nmb_normalsoil

cd ..
gmake clean
fi

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

date
