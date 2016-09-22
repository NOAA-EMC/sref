#!/bin/sh

set -aeux

res=16km

for MODEL in em raw nmm nmb
do


for OUTGRID in 212 216 221 243 132
do

rm -rf fort.*

ln -sf griddef.out_${res}${MODEL} fort.15
./wgtmkr.x <<EOF
255
${OUTGRID}
EOF

cp fort.51 sref_wgt_${MODEL}_${OUTGRID}_${res}
#cp sref_wgt_${MODEL}_${OUTGRID}_${res} ../../fix

done
done
