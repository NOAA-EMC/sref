#######################################
# Script: ensadd.sh
# ABSTRACT:  This scripts produces a gribindex of 
#            files.  Theses files are then feed 
#            into the executable ensaddx 
#######################################
#cd ${WORK_DIR}
#EXECGLOBAL=${NWROOTp1}/exec
#EXECsref=/meso/save/Jun.Du/sref.v7.0.0/exec

set +x
echo " "
echo "Entering sub script  ensadd.sh"
echo "File 1 =$1"
echo "File 2 =$2"
echo "File 3 =$3"
echo "File 4 =$4"
echo "File 5 =$5"
echo " "
set -x

if [[ $# != 5 ]];then
  echo "Usage: $0 ienst iensi gribin indexin gribout"
  echo "       inserts ensemble PDS extensions in GRIB file"
  exit 1
fi

#eval $EXECGLOBAL/global_ensadd <<EOF
eval $EXECsref/global_ensadd <<EOF
 &namin
 ienst=$1,iensi=$2,cpgb='$3',cpgi='$4',cpge='$5' /
EOF

set +x
echo " "
echo "Leaving sub script  ensadd.sh"
echo " "
set -x
